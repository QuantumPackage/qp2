#if defined(__linux__)
#include <sys/syscall.h>   /* SYS_getcpu — defined for every Linux arch */
#endif
#include <unistd.h>        /* syscall()                                  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <limits.h>

#define MAX_CPUS  8192

static int        ngpus = 0;
static int        cpu_to_gpu[MAX_CPUS];   /* zero-initialised by C runtime */
static volatile int affinity_ready = 0;

#if defined(__linux__)

/* ── PCI vendor IDs ────────────────────────────────────────────────── */
#define VENDOR_NVIDIA  0x10de
#define VENDOR_AMD     0x1002
#define VENDOR_INTEL   0x8086   /* Xe / Arc / Ponte Vecchio             */

/* ── PCI class codes (upper 2 bytes of the 3-byte class register) ─── */
#define CLASS_VGA      0x0300   /* VGA-compatible controller            */
#define CLASS_3D       0x0302   /* 3D controller (headless compute GPU) */
#define CLASS_DISPLAY  0x0380   /* Display controller                   */

#define MAX_GPUS  64

typedef struct {
    unsigned int domain, bus, slot, func;
    int          device_index;   /* 0-based, sorted by BDF               */
    int          numa_node;      /* -1 when NUMA is not exposed          */
} gpu_info_t;

static gpu_info_t gpus[MAX_GPUS];

/* ── tiny helpers ──────────────────────────────────────────────────── */

static int read_hex(const char *path, unsigned long *val)
{
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    int ok = (fscanf(f, "%lx", val) == 1);
    fclose(f);
    return ok ? 0 : -1;
}

/* Returns the NUMA node that owns logical CPU core 'cpu',
   by parsing /sys/devices/system/node/nodeN/cpulist.
   Returns -1 when NUMA sysfs is absent (UMA system).            */
static int cpu_numa_node(int cpu)
{
    char path[PATH_MAX], buf[8192];
    for (int n = 0; n < 256; n++) {
        snprintf(path, sizeof(path),
                 "/sys/devices/system/node/node%d/cpulist", n);
        FILE *f = fopen(path, "r");
        if (!f) break;
        if (fgets(buf, sizeof(buf), f) == NULL) {
            fclose(f);
            break;
        }
        fclose(f);

        /* cpulist is a range string: "0-23,48-71"  */
        char *tok = strtok(buf, ",\n");
        while (tok) {
            int lo, hi;
            if (sscanf(tok, "%d-%d", &lo, &hi) == 2) {
                if (cpu >= lo && cpu <= hi) return n;
            } else if (sscanf(tok, "%d", &lo) == 1) {
                if (cpu == lo) return n;
            }
            tok = strtok(NULL, ",\n");
        }
    }
    return -1;
}

static int bdf_cmp(const void *a, const void *b)
{
    const gpu_info_t *ga = a, *gb = b;
    if (ga->domain != gb->domain) return (int)(ga->domain - gb->domain);
    if (ga->bus    != gb->bus   ) return (int)(ga->bus    - gb->bus);
    if (ga->slot   != gb->slot  ) return (int)(ga->slot   - gb->slot);
    return (int)(ga->func - gb->func);
}

/* ── public API ────────────────────────────────────────────────────── */

/*
 * Scan /sys/bus/pci/devices/, collect GPU entries, sort by BDF,
 * then build the cpu -> gpu map.
 *
 * Fallback strategy (in order):
 *  1. Prefer NUMA-local GPU.
 *  2. If no NUMA info is available on either side (UMA systems, VMs,
 *     containers), spread CPUs round-robin across GPUs by CPU index.
 *  3. If no GPU is detected at all, warn and map everything to 0 so
 *     the caller can still function.
 *
 * REQUIREMENT (CUDA): export CUDA_DEVICE_ORDER=PCI_BUS_ID
 * HIP already defaults to PCI bus order.
 *
 * NOT thread-safe by itself — call from a single-threaded context, or
 * let get_closest_gpu_() handle the guarded lazy init below.
 *
 * Fortran visible as: build_gpu_affinity_map()
 */
void build_gpu_affinity_map_(void)
{
    DIR *dir = opendir("/sys/bus/pci/devices");
    if (!dir) {
        fprintf(stderr,
                "[gpu_affinity] WARNING: cannot open /sys/bus/pci/devices"
                " — mapping all CPUs to device 0\n");
        affinity_ready = 1;
        return;   /* cpu_to_gpu already zero-initialised */
    }

    ngpus = 0;
    struct dirent *e;
    while ((e = readdir(dir)) && ngpus < MAX_GPUS) {
        if (e->d_name[0] == '.') continue;

        char base[PATH_MAX], path[PATH_MAX];
        snprintf(base, sizeof(base), "/sys/bus/pci/devices/%s", e->d_name);

        /* Filter by vendor ─────────────────────────────────────────── */
        unsigned long vendor;
        path[0] = '\0';
        strncat(path, base, sizeof(path) - 1);
        strncat(path, "/vendor", sizeof(path) - strlen(path) - 1);
        if (read_hex(path, &vendor) < 0) continue;
        if (vendor != VENDOR_NVIDIA &&
            vendor != VENDOR_AMD    &&
            vendor != VENDOR_INTEL) continue;

        /* Filter by PCI class (drop the programming-interface byte) ─ */
        unsigned long pci_class;
        path[0] = '\0';
        strncat(path, base, sizeof(path) - 1);
        strncat(path, "/class", sizeof(path) - strlen(path) - 1);
        if (read_hex(path, &pci_class) < 0) continue;
        unsigned int cls = (unsigned int)(pci_class >> 8);
        if (cls != CLASS_VGA && cls != CLASS_3D && cls != CLASS_DISPLAY)
            continue;
        /* Intel CLASS_VGA is the integrated GPU — not a compute device.
           Intel discrete compute GPUs (Arc, Ponte Vecchio) use CLASS_3D
           or CLASS_DISPLAY, so we can safely drop Intel CLASS_VGA here. */
        if (vendor == VENDOR_INTEL && cls == CLASS_VGA)
            continue;

        /* NUMA node (/sys value is -1 when not a NUMA system) ──────── */
        long numa = -1;
        path[0] = '\0';
        strncat(path, base, sizeof(path) - 1);
        strncat(path, "/numa_node", sizeof(path) - strlen(path) - 1);
        FILE *f = fopen(path, "r");
        if (f) {
            if (fscanf(f, "%ld", &numa) != 1) numa = -1;
            fclose(f);
        }

        /* BDF from directory name:  dddd:bb:ss.f ───────────────────── */
        unsigned int dom, bus, slot, func;
        if (sscanf(e->d_name, "%x:%x:%x.%x",
                   &dom, &bus, &slot, &func) != 4) continue;

        gpus[ngpus++] = (gpu_info_t){
            .domain = dom, .bus = bus, .slot = slot, .func = func,
            .numa_node = (int)numa
        };
    }
    closedir(dir);

    /* ── No GPU detected ─────────────────────────────────────────── */
    if (ngpus == 0) {
        fprintf(stderr,
                "[gpu_affinity] WARNING: no GPU found in /sys/bus/pci/devices"
                " — mapping all CPUs to device 0\n");
        affinity_ready = 1;
        return;   /* cpu_to_gpu already zero-initialised → device 0 */
    }

    /* Sort by BDF so device_index matches the runtime's PCI-order ─── */
    qsort(gpus, ngpus, sizeof(gpu_info_t), bdf_cmp);
    for (int g = 0; g < ngpus; g++) gpus[g].device_index = g;

    /* ── Build cpu -> closest gpu map ──────────────────────────────── *
     *                                                                   *
     * Fill every entry with a valid round-robin default first, then    *
     * overwrite with the NUMA-local GPU where topology is known.       *
     * This avoids any sentinel/-1 window that a racing reader could    *
     * observe, and handles UMA / no-NUMA-sysfs systems automatically.  *
     * ─────────────────────────────────────────────────────────────── */

    /* Pass 1: round-robin default — always a valid device index -------- */
    for (int c = 0; c < MAX_CPUS; c++)
        cpu_to_gpu[c] = c % ngpus;

    /* Pass 2: NUMA-aware overwrite where topology is available --------- */
    int numa_mapped = 0;
    for (int c = 0; c < MAX_CPUS; c++) {
        int node = cpu_numa_node(c);
        if (node < 0) continue;
        for (int g = 0; g < ngpus; g++) {
            if (gpus[g].numa_node == node) {
                cpu_to_gpu[c] = gpus[g].device_index;
                numa_mapped++;
                break;
            }
        }
    }

    if (numa_mapped == 0) {
        fprintf(stderr,
                "[gpu_affinity] WARNING: no NUMA topology found"
                " — assigning CPUs round-robin across %d GPU(s)\n", ngpus);
    }

    affinity_ready = 1;

    printf("[gpu_affinity] %d GPU(s) found (sorted by PCI BDF):\n", ngpus);
    for (int g = 0; g < ngpus; g++)
        printf("  device %d  %04x:%02x:%02x.%x  NUMA node %d\n",
               gpus[g].device_index,
               gpus[g].domain, gpus[g].bus, gpus[g].slot, gpus[g].func,
               gpus[g].numa_node);
}

#else

/* Non-Linux (macOS): no /sys PCI or NUMA topology.  Map every CPU to
   device 0 — on such platforms gpu_num is 0 (simple_gpu CPU emulation)
   so the affinity map is inert. */
void build_gpu_affinity_map_(void)
{
    affinity_ready = 1;
}

#endif /* defined(__linux__) */


static int get_current_cpu(void)
{
#if defined(__linux__)
    unsigned int cpu, node;
    if (syscall(SYS_getcpu, &cpu, &node, NULL) == 0)
        return (int)cpu;
    return 0;   /* fallback: should never happen on Linux >= 2.6.19 */
#else
    return 0;   /* no getcpu() on macOS */
#endif
}

/*
 * Returns the GPU device index closest to the calling thread.
 * Intended to be called from within an OpenMP parallel region.
 *
 * Thread safety: uses an OMP critical section for the one-time
 * lazy initialisation so that concurrent first callers don't race
 * on cpu_to_gpu while it is being written by build_gpu_affinity_map_().
 * After the first call sets affinity_ready = 1 the critical section
 * is bypassed entirely on every subsequent call.
 */
int get_closest_gpu_(void)
{
    /* Double-checked locking pattern with explicit flushes — required  *
     * for correctness in OpenMP.  The flush before the outer test      *
     * forces the thread to read affinity_ready from memory rather than *
     * a cached register.  The implicit flush at the end of the         *
     * critical section ensures cpu_to_gpu/ngpus are visible to all     *
     * threads before affinity_ready is observed as 1.                  */
#pragma omp flush(affinity_ready)
    if (!affinity_ready) {
#pragma omp critical (gpu_affinity_init)
        {
#pragma omp flush(affinity_ready)
            if (!affinity_ready)
                build_gpu_affinity_map_();
        } /* implicit flush here — cpu_to_gpu and ngpus now visible    */
    }

    int cpu = get_current_cpu();
    if (cpu < 0 || cpu >= MAX_CPUS) {
        /* cpu index out of our table — clamp to round-robin fallback  */
        return (ngpus > 0) ? (cpu % ngpus) : 0;
    }

    int gpu = cpu_to_gpu[cpu];

    /* Final safety clamp: should never trigger, but defends against   *
     * future changes that could leave a stale -1 or an out-of-range   *
     * value in the table.                                             */
    if (gpu < 0 || (ngpus > 0 && gpu >= ngpus))
        gpu = (ngpus > 0) ? (cpu % ngpus) : 0;

    return gpu;
}
