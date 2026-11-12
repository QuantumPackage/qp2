#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __linux__
/* Best-effort pin of the BLAS pool, called before the first BLAS
 * call.  Verified 2026-11-11 (Ampere Xeon 6730P): the rPT2 drift
 * across thread counts is caused by the OMP default team
 * (OMP_NUM_THREADS at launch), which sizes the pool OpenBLAS uses
 * to split the K reduction of the small-K GEMMs of the one-electron
 * AO->MO transform -- not by OPENBLAS_NUM_THREADS.  The launch
 * contract OMP_NUM_THREADS=1 + QP_NTHREADS_PT2=N (see the pt2.irp.f
 * header) fixes the root cause.  This pin is hardening for BLAS
 * builds where OPENBLAS_NUM_THREADS is the lever; on the Ampere
 * nodes it is inert (the env vars are read at dlopen, before this
 * constructor). */
/* Not exported by reference BLAS, MKL or Accelerate.  Declared weak so
 * the link succeeds against any BLAS; the call is skipped when absent. */
extern void openblas_set_num_threads(int) __attribute__((weak));
static void qp_blas_pin_pool(void)
{
  if (openblas_set_num_threads != NULL)
    openblas_set_num_threads(1);
}
/* Bit-reproducibility of the pt2 stochastic walk across thread
 * counts: every integral-table build (one_e AO->MO transform, jj
 * and jj_exchange tables, the MO-integrals cache, the Cholesky
 * transform, H|u>, Hii) must run serially, i.e. with the OMP
 * default at 1; the walk itself is parallelized by its explicit
 * NUM_THREADS(nproc_target) clause and needs no OMP default.
 *
 * The constructor below runs at process start (before the Fortran
 * main, hence before the provider calls that IRPF90 hoists above it
 * build the tables): it pins OPENBLAS_NUM_THREADS and
 * VECLIB_MAXIMUM_THREADS to 1, and warns when the OMP default is
 * not 1, which would silently reintroduce the drift.
 *
 * Only the pt2 program gets this (the walk is the only part of QP2
 * whose bit-reproducibility is a contract); other programs keep
 * their BLAS threading.  Skippable with QP2_NO_BLAS_FORCE=1 for
 * diagnosis.
 */
__attribute__((constructor))
static void qp_blas_env_init(void)
{
  char comm[64];
  FILE *f;
  const char *omp;

  if (getenv("QP2_NO_BLAS_FORCE") != NULL) return;
  /* /proc/self/comm is a regular file, not a symlink: readlink(2)
   * fails on it (EINVAL) and the whole gate would never run. */
  f = fopen("/proc/self/comm", "r");
  if (f == NULL) return;
  if (fgets(comm, 64, f) == NULL) { fclose(f); return; }
  fclose(f);
  comm[strcspn(comm, "\n")] = '\0';
  if (strcmp(comm, "pt2") != 0) return;
  setenv("OPENBLAS_NUM_THREADS", "1", 1);
  setenv("VECLIB_MAXIMUM_THREADS", "1", 1);
  qp_blas_pin_pool();
  omp = getenv("OMP_NUM_THREADS");
  if (omp == NULL || atoi(omp) > 1)
    fprintf(stderr,
      "pt2: OMP_NUM_THREADS is not 1; rPT2 may not be bit-reproducible\n"
      "across thread counts. Launch with OMP_NUM_THREADS=1 and set the\n"
      "walk team size with QP_NTHREADS_PT2 instead; the walk stays\n"
      "parallel through its explicit OpenMP clause.\n");
}
#endif

void usleep_c(int s)
{
  usleep((useconds_t) s);
}

void sscanf_ssds_c(const char* str, char* s1, char* s2, int* i, char* s3)
{
  sscanf(str, "%s %s %d %s", s1, s2, i, s3);
  s1[strlen(s1)] = ' ';
  s2[strlen(s2)] = ' ';
  s3[strlen(s3)] = ' ';
}

void sscanf_dd_c(const char* str, int* i1, int* i2)
{
  sscanf(str, "%d %d", i1, i2);
}

void sscanf_ddd_c(const char* str, int* i1, int* i2, int* i3)
{
  sscanf(str, "%d %d %d", i1, i2, i3);
}

void sscanf_ss_c(const char* str, char* s1, char* s2)
{
  sscanf(str, "%s %s", s1, s2);
  s1[strlen(s1)] = ' ';
  s2[strlen(s2)] = ' ';
}

void sscanf_sd_c(const char* str, char* s1, int* i)
{
  sscanf(str, "%s %d", s1, i);
  s1[strlen(s1)] = ' ';
}

#ifdef __APPLE__
#include <sys/sysctl.h>
#include <mach/mach.h>

/* Total physical memory of the host, in bytes (macOS). */
long long qp_phys_mem_bytes(void)
{
  long long size = 0;
  size_t len = sizeof(size);
  if (sysctlbyname("hw.memsize", &size, &len, NULL, 0) != 0)
    return 0;
  return size;
}

static mach_task_basic_info_data_t qp_task_basic(void)
{
  mach_task_basic_info_data_t info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                (task_info_t)&info, &count) != KERN_SUCCESS) {
    info.resident_size = 0;
    info.virtual_size = 0;
  }
  return info;
}

/* Resident set size of the current process, in bytes (macOS). */
long long qp_rss_bytes(void)
{
  return (long long) qp_task_basic().resident_size;
}

/* Virtual size of the current process, in bytes (macOS). */
long long qp_vsize_bytes(void)
{
  return (long long) qp_task_basic().virtual_size;
}
#endif

