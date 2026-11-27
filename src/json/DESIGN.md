# JSON result logging in QP3

Status: design, ready to implement. Nothing in this document is implemented
yet.

## 1. Goal

After every calculation, QP writes one JSON document per run holding

* the EZFIO input parameters the calculation inspected,
* the numbers the calculation produced (energies, PT2, diagnostics).

The document is consumed by other codes (TREX2 and friends) instead of
grepping QP stdout. It is valid JSON at every instant, including after a
`kill -9`, because validity is a property of how the file is written and
never depends on the program reaching an orderly exit.

Adding a number to the document must be one line of Fortran:

```fortran
call json_add('/results/scf/energy', SCF_energy)
```

Non-goals: replacing EZFIO, logging wave-function vectors, keeping a
per-iteration append-only trace (see §14), Windows support.

## 2. Why the current mechanism is replaced

`src/json/json_formats.irp.f` provides 20 Fortran format strings and the call
sites write JSON text directly to a unit with `write(json_unit, fmt)`. The
scheme is CIPSI-only, and what it produces in practice is not a document:

* every file on this machine is a two-byte stub. The four
  `~/qp2/*.ezfio/json/00001.json` files are each exactly `{`. The file is
  opened with `'w'` at the first provider touch and only terminated by
  `call json_close`, so any exit before the last executable statement of the
  program — `stop 1`, `qp_stop`, a crash, a killed job — leaves it open.
* opening and closing are matched by hand at distant call sites, and the
  "trailing comma" / "no trailing comma" variants (`_fmtx`) are paired
  incorrectly. `scf_utils/roothaan_hall_scf.irp.f` writes `saved` (with a
  trailing comma) and the MO save outside the `lock_io` region that guards the
  rest of the iteration report, and closes the loop's dictionary with the
  `_fmtx` variant after opening it with the `uopen` one;
  `cipsi/write_cipsi_json.irp.f` re-emits a `"states"` key into a dictionary
  that already has one on every iteration; `kohn_sham/ks_scf.irp.f` and
  `kohn_sham_rs/rs_ks_scf.irp.f` write the bare lines `'"scf" : ['` and `']'`
  around the SCF.
* the formats themselves are unsound for the job: reals go out with `E22.15`,
  i.e. 15 significant digits where 17 are needed to round-trip a `double`, so
  what a consumer reads is not the number QP computed; and `json_str_fmt`
  emits its string through a bare `A` descriptor, so the first `"` or `\` in a
  `mo_label` or a functional name will break the file (no call site uses it
  yet).

Supervising that text layout from the outside has no chance of staying
correct: there is no object that knows what is open. Reconstructing anything
from the result is guesswork — `scripts/utility/qp_json.py` ships a `fix_json`
function that tries to re-balance the braces of the broken files it reads.
This design throws the text-writer away, builds the document in memory, and
serialises it.

## 3. Layers

```
Fortran call sites        call json_add('/results/scf/energy', e)
                          call json_init('scf')
        |
src/json/json.irp.f       IRPF90 glue: json_filename, json_init
src/json/qp_json.f90      plain module: generic json_add over 6 specifics
        |                 (bind(C) interfaces, no state)
src/json/qp_json_shim.c   C shim: owns the yyjson document, path resolution,
        |                 envelope, timestamps, throttle, writer thread
external/yyjson/yyjson.c/.h  vendored yyjson 0.13.0 (MIT)
```

No JSON state is visible to Fortran. Fortran passes a value and an absolute
path; everything else is the shim's business. The idiom is the same one
`src/utils/c_functions.f90` + `src/utils/c_funcs.c` already use: plain `.f90`
and `.c` files next to the `.irp.f` files of a module directory are compiled
and linked by the existing build, with no build-system change.

**The shim must not be named `qp_json.c`.** `get_l_file_for_module()` in
`scripts/compilation/qp_create_ninja` names the object after the file stem
(`obj = '{0}.o'.format(os.path.splitext(f)[0])`, line 478), and the object
directory is shared by every source file of the module
(`src/json/IRPF90_temp/`). `qp_json.c` and `qp_json.f90` would both produce
`qp_json.o` and one would silently overwrite the other. Hence
`qp_json_shim.c`.

## 4. Vendoring yyjson

* Version 0.13.0 (tag `0.13.0`, `YYJSON_VERSION_STRING "0.13.0"`, MIT).
* Two files, the flat amalgamation, in `external/yyjson/`: `yyjson.c` and
  `yyjson.h`, with the upstream `LICENSE` and a `README` pointing at the
  release. Not a submodule.
* An external library, treated like ezfio: the build system compiles
  `external/yyjson/yyjson.c` with the global `CFLAGS` and archives it once
  per tree as `external/yyjson/libyyjson.a`, exposed as
  `$QP_ROOT/lib/libyyjson.a`. `scripts/compilation/qp_create_ninja` emits
  the `build_yyjson` edge (ninja mode) and the `$(YYJSON_LIB)` recipe
  (make mode); every `build_binaries` edge and every `module_*` target
  depends on it.
* `src/json` consumes it like any module consumes an external library:
  `src/json/LIB` lists `$QP_ROOT/lib/libyyjson.a`, and `src/json/CFLAGS`
  adds `-I$QP_ROOT/external/yyjson`, which the shim's
  `#include "yyjson.h"` now needs (the header no longer sits in the module
  directory). No module compiles yyjson itself.

**Trap.** `config/build.ninja` sets `CFLAGS = -Ofast`, which implies
`-ffast-math`. yyjson's `dtoa` must print round-trip-exact doubles; a
`-ffast-math` miscompile there would silently corrupt every energy in the
file. Phase 0 (§13) gates this with a round-trip test before anything else.
Fallback if it fails: set `CFLAGS = -O3` in the build configuration (only two
small `.c` files in the whole tree) and record the measurement in the README.

## 5. Document

### 5.1 Root shape

Seven root keys, all created by the shim when `qp_json_begin` runs:

```json
{
  "schema": "qp3/1",
  "creator": { "name": "Quantum Package", "version": "unknown" },
  "provenance": {
    "ezfio_file": "/home/x/water.ezfio",
    "host": "lcpqsv225",
    "pid": 4242,
    "started": "2026-11-24T14:03:11Z"
  },
  "units": { "energy": "hartree", "length": "bohr", "time": "second" },
  "run": { "kind": "fci", "timestamp": "2026-11-24T14:07:52Z", "elapsed": 281.4 },
  "input": {
    "determinants": { "n_states": 2, "s2_eig": true },
    "hartree_fock": { "thresh_scf": 1e-8 }
  },
  "results": {
    "scf": { "iteration": 12, "energy": -75.4123, "converged": true },
    "cipsi": {
      "n_det": 120345, "n_cfg": 8123,
      "states": {
        "energy": [-75.9123, -75.8011],
        "s2": [0.0, 0.0],
        "pt2": [-0.4123, -0.3512],
        "pt2_err": [0.0001, 0.0002],
        "rpt2": [-0.4110, -0.3500],
        "rpt2_err": [0.0002, 0.0003],
        "variance": [0.02, 0.03],
        "variance_err": [0.001, 0.001],
        "extrapolated_energy": [-76.2011, -76.1002]
      },
      "history": {
        "n_iter": 8, "n_iter_max": 8,
        "shape": [2, 8],
        "energy": [ ... 16 numbers ... ],
        "rpt2":   [ ... 16 numbers ... ]
      }
    }
  }
}
```

Rules:

* `schema` is the compatibility handle. A consumer that does not know
  `"qp3/1"` must refuse the file.
* `creator.version` is a placeholder (`"unknown"`) until a build-time
  `git describe` string exists (§14).
* `units` is declared once, for the whole document. Every energy is hartree,
  every length bohr, every time second, in `input` and in `results`. Call
  sites never carry a unit in a key.
* `run.kind` is the program name (`'scf'`, `'ks_scf'`, `'rs_ks_scf'`,
  `'fci'`, `'pt2'`, `'mrci'`, `'casscf'`, `'cipsi_orb_opt'`), passed once by
  `json_init`. `run.timestamp` and `run.elapsed` are maintained by the shim at
  each actual write, never by Fortran.
* `input` is nested like the EZFIO groups: `/input/<ezfio_group>/<ezfio_name>`.
* `results` is grouped by method, not by program: `/results/scf/...`,
  `/results/ks/...`, `/results/cipsi/...`, `/results/casscf/...`,
  `/results/mo_opt/...`. `program fci` and `program mrci` both call
  `run_cipsi`, so both fill `/results/cipsi`; `run.kind` says which program
  ran. `ks_scf` and `rs_ks_scf` fill `/results/scf` from the SCF driver, plus
  `/results/ks` for the RS decomposition.

### 5.2 Per-state data

QP holds per-state quantities as rank-1 arrays over `1:N_states`
(`psi_energy_with_nucl_rep`, `pt2_data % pt2`, `psi_s2`, ...). They are
logged as arrays under a `states` object, one key per quantity, not as an
array of per-state objects. This is what makes the Fortran call one line: no
index has to be formatted into a path.

Rank-2 arrays are not supported by the API (§7.4). When a 2-D block must be
logged, the call site flattens it in Fortran order and records the shape in a
sibling array:

```fortran
call json_add('/results/cipsi/history/shape', (/ n, N_iter /))
```

with the JSON array `[n1, n2]` meaning the data has `n1` elements in the first
Fortran dimension. A reader reshapes with
`numpy.reshape(data, shape, order='F')`. This is the only shape convention in
the document.

## 6. File model

One file per run, rewritten in place, atomically:

* Path: `<EZFIO>/json/NNNNN.json`, `NNNNN` being the first unused number
  (existing provider logic, §8.1). A restart gets the next number; the
  previous run's file is never touched again.
* Every write replaces the whole file. The document in memory is cumulative,
  so the newest file is a complete, self-contained record of the run: a
  `program fci` run contains the SCF results computed in-process as well as
  the CIPSI results.
* A temp file `<path>.tmp.<pid>` is written in the same directory,
  `fflush`ed, `fsync`ed, `fclose`d and then `rename(2)`d onto `<path>`.
  `rename` within a directory is atomic, so a reader (or a post-mortem) sees
  either the whole previous version or the whole new one, never a truncated
  document.
* Nothing else is done for crash safety: no signal handler, no locking
  protocol between reader and writer. Validity is a property of the layout. A
  `.tmp.<pid>` file left by a `SIGKILL` is garbage; readers only accept
  `^[0-9]{5}\.json$`.
* Power-loss durability (directory `fsync`) is out of scope.

## 7. C shim `qp_json_shim.c`

### 7.1 State

```c
#define QP_JSON_PATH_MAX 4096
#define QP_JSON_TICK_NS  100000000LL      /* 100 ms */

static pthread_mutex_t mu    = PTHREAD_MUTEX_INITIALIZER; /* the document */
static pthread_mutex_t io_mu = PTHREAD_MUTEX_INITIALIZER; /* the file      */

static yyjson_mut_doc *doc;                  /* NULL until first add       */
static char     out_path[QP_JSON_PATH_MAX];  /* empty until qp_json_begin  */
static char     run_kind[64];
static char     ezfio_file[QP_JSON_PATH_MAX];
static char     host[256];
static int      enabled = 1;                 /* 0 on non-master ranks      */
static int      dirty;                       /* document newer than file   */
static int      stop_thread;
static int      thread_started;
static int      atexit_done;
static double   min_interval = -1.0;         /* <0: not read from env yet  */
static double   t_last_write = -1e18;        /* CLOCK_MONOTONIC seconds    */
static struct timespec t_start;
static pthread_t writer_tid;
```

Lock order, everywhere: `io_mu` then `mu`. The document is only touched under
`mu`; the file is only touched under `io_mu`. `mu` is never held across a
`read`/`write`/`fsync`/`rename`.

### 7.2 Path language

```c
static void put_scalar(const char *path, int type,
                       long long iv, double dv, const char *sv);
static void put_array(const char *path, int type, const void *data, int n);
```

with `type` in `{T_INT, T_REAL, T_BOOL, T_STR, T_ARR_INT, T_ARR_REAL}`.
Both take `mu`, create the document lazily, walk the path and insert. All
four public scalar setters and the two array setters are one-liners over
them, so the whole mutation logic lives in one place.

Parsing rules, in this order; a rule violation prints one `stderr` line
(`qp_json: ...`) and returns without mutating:

* the path must start with `/` and be shorter than `QP_JSON_PATH_MAX`;
* it must not contain a space (this catches a stray `write` format);
* components are separated by `/`; empty components are ignored, so `/a//b`
  and `/a/b/` are accepted;
* intermediate components create objects. If an existing intermediate node
  is not an object, it is replaced by one and a warning is printed;
* the last component is the key. `yyjson_mut_obj_put` replaces an existing
  value for that key, so re-logging the same path is idempotent.

Consequences for call sites: absolute paths from the root, lowercase
`snake_case`, no spaces, no leading `results/` shortcut. There is no
per-call handle, no open/close, no unit number.

### 7.3 Value rules

* `T_INT` → `yyjson_mut_int`; `T_REAL` → `yyjson_mut_real` (shortest
  round-trip representation, so an energy reads back bit-exactly);
  `T_BOOL` → `yyjson_mut_bool`; `T_STR` → `yyjson_mut_strcpy` (the string is
  copied, so Fortran temporaries are safe).
* Non-finite reals become `null`. This is done by the write flag
  `YYJSON_WRITE_INF_AND_NAN_AS_NULL`, not by an `isfinite` test — `isfinite`
  is folded to `true` under `-ffast-math`. Without the flag one NaN makes
  `yyjson_mut_write` return `NULL` and the whole file stops being updated,
  which is exactly the failure mode we must avoid.
* `null` therefore means "not finite / not available". Sentinel values such
  as `pt2_data % pt2 = -huge(1.d0)` are finite and are logged literally; call
  sites write a path only when the value is meaningful.
* Empty arrays produce `[]`.

### 7.4 Arrays

Rank-1 only, `integer` and `real(c_double)`. No rank-2: a rank-2 actual
argument must fail to compile rather than be silently flattened.

### 7.5 Initialisation

```c
void qp_json_begin(const char *kind, const char *filename,
                   const char *ezfio_file, int master);
```

1. If `master == 0`: set `enabled = 0`, free the document if one was built by
   an earlier `qp_json_add`, and return. Every other entry point returns
   immediately when `enabled == 0`, so non-master ranks do no work and no I/O.
2. Read `QP_JSON_MIN_INTERVAL` (§7.6) once.
3. Ensure the document exists (create it if `json_add` was never called),
   then set `schema`, `creator`, `provenance` (`ezfio_file`, `gethostname`,
   `getpid`, `started` from `CLOCK_REALTIME` → `gmtime_r` →
   `%Y-%m-%dT%H:%M:%SZ`), `units`, an empty `input` and `results` object if
   they are absent, and `run.kind`. These puts are unconditional: `qp_json_add`
   makes them idempotent, and the document may already exist because a
   `json_add` ran first (§7.2).
4. Store `out_path`. If it is already set to the same string, return
   silently; if it is set to a different string, warn and keep the first one.
5. Start the writer thread if `min_interval > 0` (§7.6), register the `atexit`
   handler once (§7.7), and force one synchronous write, so the file exists
   and is valid from the first millisecond with provenance and schema even if
   the program dies before the first throttled write.

`json_add` before `qp_json_begin` is legal: the document is built lazily and
nothing is written until `out_path` is known. A program that never calls
`qp_json_begin` loses its observations silently except for one warning
printed by the exit handler (§7.7).

### 7.6 Writer thread and throttle

`QP_JSON_MIN_INTERVAL` (seconds, floating point, default `1.0`):

| value | behaviour |
|---|---|
| `> 0` | writer thread; at most one write per interval, and a dirty document is written once the interval has elapsed since the previous write |
| `0` | no thread; every mutation is followed by a synchronous write — deterministic mode for tests |
| invalid / NaN / negative | default `1.0` |
| `pthread_create` fails | falls back to `0` and warns |

Thread loop (100 ms tick, ten wake-ups per second, negligible for a QP run):

```c
for (;;) {
    struct timespec tick = {0, QP_JSON_TICK_NS};
    nanosleep(&tick, NULL);
    pthread_mutex_lock(&mu);
    int stop = stop_thread;
    pthread_mutex_unlock(&mu);
    if (stop) return NULL;
    do_write(0);
}
```

`do_write(force)` is the only writer:

```c
lock(io_mu);
  lock(mu);
    if (!enabled || !doc || !out_path[0])             { unlock; unlock; return 0; }
    if (!dirty && !force)                             { unlock; unlock; return 0; }
    if (!force && since(t_last_write) < min_interval) { unlock; unlock; return 0; }
    dirty = 0;              /* clear *before* the snapshot: any later
                               mutation re-dirties and is picked up next tick */
    stamp_run();            /* run.timestamp, run.elapsed */
    text = yyjson_mut_write(doc, QP_JSON_WRITE_FLAGS, NULL);
  unlock(mu);
  if (!text) { lock(mu); dirty = 1; unlock(mu); warn; unlock(io_mu); return -1; }
  /* fopen "<out_path>.tmp.<pid>", fwrite, fflush, fsync, fclose, rename */
  lock(mu); t_last_write = now_monotonic(); if (!ok) dirty = 1; unlock(mu);
unlock(io_mu);
```

`QP_JSON_WRITE_FLAGS` = `YYJSON_WRITE_PRETTY_TWO_SPACES |
YYJSON_WRITE_NEWLINE_AT_END | YYJSON_WRITE_INF_AND_NAN_AS_NULL`. Pretty
printing costs microseconds once per second and makes the file readable with
`jq`, `less` and `diff`.

Clear-then-snapshot under `mu` is what makes the coalescing correct without a
generation counter: a mutation that arrives while the text is being built (or
while the file is being written) sets `dirty` again and is picked up by the
next tick.

Because `io_mu` is taken before the snapshot, a write that starts later also
serialises a later document, so the last write always wins with the freshest
content. No ordering flag is needed.

In interval-0 mode the mutation path calls `do_write(0)` **after** releasing
`mu`, never while holding it: `do_write` takes `io_mu` then `mu`, `mu` is not
recursive, so a nested call would hang the process.

The 1 s delay is a *maximum staleness*: results are on disk at most one
second after Fortran produced them. `run.elapsed` and `run.timestamp` advance
once per actual write, not per `json_add`.

### 7.7 Exit path

`qp_json_begin` registers, once:

```c
static void on_exit(void) {
    lock(mu);
    if (!enabled || !doc) { unlock; return; }
    stop_thread = 1;
    unlock(mu);
    if (thread_started) {
        struct timespec tick = {0, QP_JSON_TICK_NS};
        pthread_join(writer_tid, NULL);   /* at most one 100 ms tick */
    }
    if (!out_path[0]) {
        fprintf(stderr, "qp_json: qp_json_begin was never called; "
                        "the collected values are discarded\n");
        return;
    }
    do_write(1);                          /* force, ignore the throttle */
}
```

so a normal `end program` or `stop 1` leaves a file containing every value
logged so far, with no Fortran-side discipline needed. The thread is created
joinable (not detached) precisely so the exit handler can join it: after the
join, no other thread touches the document or the file, and the process pays
at most 100 ms. The document is not freed — the process is exiting.

A `SIGKILL`/`SIGTERM`/`SIGSEGV` skips the handler. The file on disk is then
the last throttled write, at most one second old, and valid.

No signal handlers are installed. Writing JSON from a signal handler would
run `malloc` and `stdio` in an async-signal context, and a "dirty" flag set
from the handler is observationally identical to the one-second throttle we
already have.

### 7.8 C API (complete)

```c
void qp_json_begin(const char *kind, const char *filename,
                   const char *ezfio_file, int master);
void qp_json_add_int(const char *path, long long v);
void qp_json_add_real(const char *path, double v);
void qp_json_add_bool(const char *path, int v);
void qp_json_add_str(const char *path, const char *v);
void qp_json_add_arr_int(const char *path, const long long *v, int n);
void qp_json_add_arr_real(const char *path, const double *v, int n);
```

Seven entries, plus the internal `on_exit`/`do_write`/`put_*`/`stamp_run`.
No header file is shipped: the prototypes live at the top of
`qp_json_shim.c` and are mirrored one-to-one by the `bind(C)` interfaces in
`qp_json.f90`.

### 7.9 Two consequences of the yyjson mutable allocator

A `yyjson_mut_doc` allocates every value and every key in per-document memory
pools; `yyjson_mut_doc_free` is the only way to release anything. So:

* **Overwriting a path abandons the previous value in the pool.** Nothing
  leaks outside the document, but memory grows monotonically with the number
  of `json_add` calls, not with the number of distinct keys. Order of
  magnitude for the hottest site (§10.4): ~80 values and ~10 keys per CIPSI
  iteration, ~3 kB, so ~30 MB over 10^4 iterations. Accepted: a CIPSI run of
  that length is hours of compute. If a consumer ever needs a 10^5-iteration
  run, the fix is a 220-entry path→value cache in the shim, not a redesign;
  it is deliberately not built now. P2 measures RSS at the end of a long run
  and records the number in the README.
* **`yyjson_mut_obj_put` moves a replaced key to the end of its object.** Key
  order therefore changes between two writes of the same file. This is
  invisible to a JSON consumer (`jq` reorders nothing) and only shows up if
  someone diffs two versions of the same file. Do not "fix" it with a custom
  in-place setter.

## 8. Fortran side

### 8.1 `src/json/json.irp.f` (rewritten)

`json_filename` keeps today's numbering loop (which already implements "one
file per run") but takes over the `ezfio_set_json_empty(.False.)` call that
today lives in the deleted `json_unit` provider, and grows from 128 to 2048
characters: `ezfio_filename` is `character*(1024)`, so a 128-character buffer
silently truncates long EZFIO paths.

```fortran
BEGIN_PROVIDER [ character*(2048), json_filename ]
 implicit none
 BEGIN_DOC
 ! Path of the JSON result file of this run: <EZFIO>/json/NNNNN.json
 END_DOC
 integer                      :: counter
 character*(2048)             :: prefix
 logical                      :: exists

 prefix = trim(ezfio_filename)//'/json/'

 call lock_io
 call ezfio_set_json_empty(.False.)
 exists = .True.
 counter = 0
 do while (exists)
   counter += 1
   write(json_filename, '(A,I5.5,A)') trim(prefix), counter, '.json'
   INQUIRE(FILE=trim(json_filename), EXIST=exists)
 enddo
 call unlock_io
END_PROVIDER
```

`json_unit` and `json_close` are deleted at the end of P2, not in P1: until
the last old call site is ported, every program must keep building, so
`json.irp.f` carries them unchanged alongside the new `json_init` (§13). The C
`atexit` handler replaces `json_close`, so Fortran has no finalisation call.

The `ezfio_set_json_empty(.False.)` call sits in this provider rather than in
`json_init` because it does two things at once: it marks the EZFIO file as
carrying JSON output, and it creates `<EZFIO>/json/` (the EZFIO setter runs
`create_dir_all` on the parent). The directory must exist before the C side
`fopen`s into it, and the `INQUIRE` loop just below runs in it already, so this
is the only ordering that is guaranteed.

`subroutine json_init(kind)` is the single entry point:

```fortran
subroutine json_init(kind)
  use iso_c_binding, only: c_int
  use qp_json, only: qp_json_begin, json_add
  implicit none
  BEGIN_DOC
  ! Start the JSON logging of this run and write the first version of the
  ! file. Called once per program, before any json_add.
  END_DOC
  character(len=*), intent(in) :: kind

  PROVIDE json_filename
  call qp_json_begin(trim(kind), trim(json_filename), trim(ezfio_filename), &
                     merge(1_c_int, 0_c_int, mpi_master))
  if (mpi_master) then
    ! ezfio_files is excluded from the automatic input capture (§9.2): a NEED
    ! entry on json would close a dependency cycle. Its one input that the
    ! calculation really inspects (ao_basis checks it against a format
    ! version) is logged here by hand.
    call json_add('/input/ezfio_files/ezfio_convention', ezfio_convention)
  endif
end
```

`mpi_master` is passed to the shim, which drops every later call on non-master
ranks. No call site needs an `if (mpi_master)` guard, which is the whole point
of "easy in Fortran".

`json_init` is called once, as the first executable statement of the
calculation routine of each program that produces results (§10.2). Programs
that never call it never create a `json/` directory and never write:
`qp_run save_natorb` and the other utilities stay silent.

### 8.2 `src/json/qp_json.f90` (new, ~150 lines)

Plain module, house style copied from `src/utils/c_functions.f90` (a `use
iso_c_binding` inside each interface body, `character(kind=c_char) :: x(*)`,
`trim(s)//c_null_char` at the call).

```fortran
module qp_json
  use iso_c_binding
  implicit none
  private
  public :: json_add, qp_json_begin

  interface
    subroutine qp_json_begin(kind, filename, ezfio_file, master) &
        bind(C, name="qp_json_begin")
      use iso_c_binding
      character(kind=c_char), intent(in) :: kind(*), filename(*), ezfio_file(*)
      integer(c_int), value, intent(in)  :: master
    end subroutine
    ! ... one interface per C entry of §7.8 ...
  end interface

  interface json_add
     module procedure json_add_int, json_add_real, json_add_bool, &
                      json_add_str, json_add_arr_int, json_add_arr_real
  end interface
contains
  ...
end module
```

Wrapper semantics:

| Fortran actual | specific | conversion |
|---|---|---|
| `integer` | `json_add_int` | `int(value, c_long_long)` |
| `real(c_double)` (= `double precision`) | `json_add_real` | none |
| `logical` | `json_add_bool` | `merge(1_c_int, 0_c_int, value)` |
| `character(len=*)` | `json_add_str` | `trim(value)//c_null_char` |
| `integer, dimension(:)` | `json_add_arr_int` | `buf = int(values, c_long_long)` |
| `real(c_double), dimension(:)` | `json_add_arr_real` | `buf = values` |

Two details that matter:

* The public array specifics take **assumed-shape, non-contiguous** dummies
  (`dimension(:)`, no `contiguous`), because call sites pass strided sections
  such as `extrapolated_energy(N_iter, 1:n)`. They copy into a local
  contiguous buffer, which is also what `bind(C)` requires.
* Use `integer(c_long_long)` and `real(c_double)` in the interface bodies,
  not `c_int64`: `src/utils/c_functions.f90` records that gfortran ≥ 16 fails
  to resolve `c_int64` inside a module's interface body.

An actual argument of an unsupported type (a single-precision `real`, an
`integer(8)`, a rank-2 array) is a compile error by design: no silent
narrowing, no silent flattening.

No runtime error of these wrappers ever stops a run. A bad path or a failed
allocation prints one `stderr` line and returns.

## 9. Input capture

### 9.1 Generator hook

`scripts/ezfio_interface/ezfio_generate_provider.py` is the only durable
place: the generated `*/ezfio_interface.irp.f` files are gitignored
(`src/.gitignore:9`).

Three edits to `EZFIO_Provider`:

1. `data`: add a `%(json_use)s` line between the `BEGIN_PROVIDER` header and
   `implicit none`, and a `%(json_log)s` line right after the
   `call ezfio_get_...` line, inside `if (has)`, inside `if (mpi_master)`:

```
BEGIN_PROVIDER [ %(type)s, %(name)s %(size)s ]
%(json_use)s
  implicit none
  ...
    call ezfio_has_%(ezfio_dir)s_%(ezfio_name)s(has)
    if (has) then
      call ezfio_get_%(ezfio_dir)s_%(ezfio_name)s(%(name)s)
%(json_log)s
    else
```

`use` must precede `implicit none` in Fortran, hence that position.

2. a new setter, called from `__repr__` right after `self.set_test_null_size()`
   (which is what guarantees `size` exists):

```python
    json_type_ok = ("logical", "integer", "double precision")

    def set_json_hook(self):
        if "size" not in self.__dict__:
            self.__dict__["size"] = ""
        scalar = (self.size == "")
        typed = (self.type in self.json_type_ok
                 or self.type.startswith("character"))
        if scalar and typed:
            self.json_use = "  use qp_json, only: json_add"
            self.json_log = ("      call json_add('/input/%(ezfio_dir)s/%(ezfio_name)s',"
                             " %(name)s)" % self.__dict__)
        else:
            self.json_use = ""
            self.json_log = ""
```

3. no change to `self.values`: `__repr__` only uses it for a check that tests
   the loop variable, so adding fields there is unnecessary; `self.data %
   self.__dict__` is satisfied as long as `set_json_hook` ran.

Why these tests:

* `size == ""` is the scalar test — `set_size` maps `"1"` to `""`, and arrays
  carry `", <dims>"`. 43 of the 246 providers are arrays; a multi-megabyte
  integral array has no place in the document.
* the type whitelist keeps the generated code compilable if a provider with a
  non-interoperable type (`integer(bit_kind)`, a single-precision `real`) ever
  appears. The 203 scalar providers of today are exactly `logical` 65,
  `integer` 51, `double precision` 43, `character*(32)` 42,
  `character*(256)` 2 — all supported, and nothing is left out.
* `Disk_access` providers never reach this code: they are generated by
  `gen_ezfio_provider_disk_access`, a different template.

`ezfio_dir`/`ezfio_name` are defaulted to the lowercased module name and the
provider name (`ei_handler.py:255-256`) and lowercased again by the setters,
so `/input/hartree_fock/thresh_scf` is the shape of every key.

`use` inside a `BEGIN_PROVIDER` is legal IRPF90 (111 occurrences in the tree).
Per-provider emission means `aux_quantities` (providers, zero scalars) gets no
`use` and needs no `NEED` edit.

### 9.2 NEED entries

The rule is one line: **a module whose generated provider file, or whose own
source, executes `use qp_json` needs `json` in its `NEED`**, so that
`json/qp_json.mod` is compiled before that module.

Applying the §9.1 filter to the whole tree through
`ei_handler.get_dict_config_file()` — the real code path, not a copy of its
logic — gives 32 modules holding at least one scalar typed provider, and 203
such providers in total (246 providers carry the `provider` token; the other
43 are array-valued). Two of the 32 are excluded below for dependency cycles
(`ezfio_files`, `utils`), leaving 30 modules and 201 hooked providers. Two of
those 30 already list `json`:

```
cipsi         scf_utils
```

The remaining 28 gain the line:

```
ao_basis             ao_extra_basis        ao_one_e_ints    ao_two_e_ints
basis                becke_numerical_grid  bitmask          casscf_cipsi
davidson_keywords    density_for_dft       determinants     dft_keywords
electrons            extra_nuclei          hamiltonian      mo_one_e_ints
mo_optimization_utils mo_two_e_ints        mol_properties   mrci
mu_of_r              nuclei                perturbation     pseudo
trexio               two_body_rdm          utils_cc         utils_trust_region
```

and one more module gains it for a hand-written call: `kohn_sham_rs`, whose
`rs_ks_scf.irp.f` does `use qp_json, only: json_add` for the `/results/ks`
block (§10.3). That is 29 `NEED` edits.

Exceptions and near-misses, so the build agent does not add them by reflex:

* `src/ezfio_files` and `src/utils` are excluded: both are in `json`'s
  dependency closure (`json → ezfio_files → utils → mpi`), so a `NEED` entry
  would close a cycle. Each has exactly one scalar provider, and only one of
  them matters: `ezfio_files/ezfio_convention` is inspected by `ao_basis`
  (format-version check) and is logged by hand in `json_init` (§8.1);
  `utils/restore_symm` is declared but never read by any source file — see
  §14 — so it is not logged.
* `hartree_fock` already lists `json`, but its `EZFIO.cfg` has no entry with
  the `provider` interface token, so nothing is hooked there. The existing
  entry is harmless after the port (it is what ordered `json` before
  `hartree_fock` in the first place); leave it.
* `kohn_sham`, `fci`, `cas_based_on_top`, `csf` and the other cfg-bearing
  directories with no `provider` keyword need no edit: no generated provider →
  no `use`. `kohn_sham` only calls `json_init` (a subroutine), and
  `cipsi_orb_opt` does not touch the JSON layer at all.
* `src/json/NEED` gains `mpi` (for `mpi_master` in `json_init`); it becomes
  `ezfio_files mpi`. This is redundant with the transitive dependency but
  makes the requirement explicit; `mpi/NEED` is empty, so no cycle.

Cycles were checked for the whole graph with these entries added: none.

`call json_init(...)` needs no `NEED` entry wherever it appears: it is an
IRPF90 subroutine, resolved globally like the current `call json_close` in
`src/mrci`. `call json_add(...)` is different — it is a module generic, so it
requires the `use` and therefore the `NEED` entry. Every `json_init` caller is
already ordered after `json` transitively (`fci` and `mrci` through `cipsi`,
`kohn_sham` through `scf_utils`, `casscf_cipsi` and `mo_optimization_utils`
through their own new entry).

## 10. Migration

### 10.1 Call sites that write JSON today

| file | today | after |
|---|---|---|
| `src/hartree_fock/scf.irp.f` (83, 92-93) | `write(json_unit,json_array_open_fmt) 'scf'` … `json_close` | `call json_init('scf')` at the top of `run`; the three `write`s and `json_close` deleted |
| `src/scf_utils/roothaan_hall_scf.irp.f` (188-220) | 12 `write(json_unit, …)` calls per iteration, with `lock_io` regions that do not cover the `saved` write | 7 `json_add` calls at the end of the loop body (§10.3) |
| `src/kohn_sham/ks_scf.irp.f` (104-108) | `write(json_unit,*) '"scf" : ['`, `']'`, `json_close` | `call json_init('ks_scf')`; the three lines deleted. The SCF driver fills `/results/scf` |
| `src/kohn_sham_rs/rs_ks_scf.irp.f` (104-107) | same raw writes | `call json_init('rs_ks_scf')`; the three lines deleted; `/results/ks/{one_e_energy,two_e_energy,exchange_dft,correlation_dft}` after the SCF |
| `src/fci/fci.irp.f` (42, 58-61) | `json_array_open_fmt 'fci'` … `json_close`, all empty | `call json_init('fci')` after `END_DOC`; 5 lines deleted |
| `src/mrci/mrci.irp.f` (42, 52-55) | same, empty | `call json_init('mrci')`; 5 lines deleted |
| `src/fci/pt2.irp.f` | nothing (results only in stdout) | `json_init('pt2')` in `run`; `call write_cipsi_json(pt2_data, pt2_data_err)` after `print_summary` |
| `src/cipsi/cipsi.irp.f` (113, 162) | 2 × `call write_cipsi_json(...)` | unchanged call sites; `json_add('/results/cipsi/algorithm', 'deterministic')` next to the first one |
| `src/cipsi/stochastic_cipsi.irp.f` (104, 143) | 2 × `call write_cipsi_json(...)` | unchanged call sites; `/results/cipsi/algorithm` = `'stochastic'` |
| `src/cipsi/write_cipsi_json.irp.f` | 20 text writes | rewritten body (§10.4) |
| `src/casscf_cipsi/casscf.irp.f` | nothing | `json_init('casscf')` in `run`; `/results/casscf/...` per iteration (§10.5) |
| `src/mo_optimization_utils/routine_opt_mos.irp.f` | `result_opt` text file | `json_init('cipsi_orb_opt')` at the top of `run_optimization_mos_CIPSI`; `/results/mo_opt/...` per cycle (§10.5) |

Deleted outright: `src/json/json_formats.irp.f` (all 20 format providers),
the `json_unit` provider, `json_close`, and their `PROVIDE` traces.

`program cis`, `cisd`, `ccsd`, `mp2` log nothing today and are left alone;
adding `/results/<method>/energy` to them is a follow-up (§14).

### 10.2 `json_init` sites

One call each, as the first executable statement of the calculation routine,
after the `END_DOC` block and before the first `PROVIDE` of a calculation
provider:

| program | where |
|---|---|
| `scf` | top of `subroutine run` (`hartree_fock/scf.irp.f:70`), replacing the `'scf'` array open |
| `ks_scf` | top of `subroutine run` (`kohn_sham/ks_scf.irp.f:87`), replacing the raw open |
| `rs_ks_scf` | top of `subroutine run` (`kohn_sham_rs/rs_ks_scf.irp.f`), replacing the raw open |
| `fci` | `program fci` body, line 38 (`fci/fci.irp.f`) |
| `pt2` | top of `subroutine run` (`fci/pt2.irp.f:71`) |
| `mrci` | `program mrci` body, line 38 |
| `casscf` | top of `subroutine run` (`casscf_cipsi/casscf.irp.f:21`) |
| `cipsi_orb_opt` | top of `run_optimization_mos_CIPSI`, which is the only caller (`mo_optimization/cipsi_orb_opt.irp.f` calls it) |

`run_cipsi` / `run_stochastic_cipsi` are also reached from `casscf` and
`routine_opt_mos`; those programs therefore get `/results/cipsi` blocks in
addition to their own, which is intended.

### 10.3 SCF

In `Roothaan_Hall_SCF`, replace the whole JSON block of the loop body
(lines 188-220: the two `lock_io` regions, the
`json_dict_uopen`/`json_dict_close_uopen` dance, the `saved` write and the
final `write(json_unit, json_dict_close_fmtx)` after the loop) with:

```fortran
    call json_add('/results/scf/iteration', iteration_SCF)
    call json_add('/results/scf/energy', energy_SCF)
    call json_add('/results/scf/delta_energy', Delta_energy_SCF)
    call json_add('/results/scf/max_error_diis', max_error_DIIS)
    call json_add('/results/scf/level_shift', level_shift)
    call json_add('/results/scf/dim_diis', dim_DIIS)
    call json_add('/results/scf/converged', converged)
```

placed after the `write(6,'(I4, 1X, ...)')` line, i.e. once `converged` has
been computed for this iteration. Each iteration overwrites the previous
values, so these scalars hold the last completed iteration, which is the
useful thing to have after a crash; the whole per-iteration series is
carried in `/results/scf/history/`, one element per iteration:
`energy`, `delta_energy`, `max_error_diis`, `level_shift`, `dim_diis`,
`converged`, `epoch`, `wall_time`, `memory` (§10.7). `saved` is dropped:
the MO state is EZFIO's business, not a result. `scf_utils/NEED` already contains `json`;
`scf_utils/roothaan_hall_scf.irp.f` needs `use qp_json, only: json_add` above
`implicit none`.

For `ks_scf` and `rs_ks_scf` this table is the converged energy as well: the
KS SCF loop is the same `Roothaan_Hall_SCF`. `run.kind` says which Hamiltonian
produced it. `rs_ks_scf` adds its decomposition after the call:

```fortran
  call json_add('/results/ks/one_e_energy', one_e_energy)
  call json_add('/results/ks/two_e_energy', two_e_energy)
  call json_add('/results/ks/exchange_dft', e_exchange_dft)
  call json_add('/results/ks/correlation_dft', e_correlation_dft)
```

### 10.4 CIPSI

`write_cipsi_json` keeps its signature and its call sites, and becomes:

```fortran
subroutine write_cipsi_json(pt2_data, pt2_data_err)
  use selection_types
  use qp_json, only: json_add
  implicit none
  BEGIN_DOC
!  Logs the CIPSI results of the current iteration.
  END_DOC
  type(pt2_type), intent(in) :: pt2_data, pt2_data_err
  integer :: n

  n = min(N_states,N_det)
  if (n < 1) return

  call json_add('/results/cipsi/n_det', N_det)
  if (s2_eig) call json_add('/results/cipsi/n_cfg', N_Configuration)
  call json_add('/results/cipsi/states/energy', psi_energy_with_nucl_rep(1:n))
  call json_add('/results/cipsi/states/s2', psi_s2(1:n))
  call json_add('/results/cipsi/states/pt2', pt2_data % pt2(1:n))
  call json_add('/results/cipsi/states/pt2_err', pt2_data_err % pt2(1:n))
  call json_add('/results/cipsi/states/rpt2', pt2_data % rpt2(1:n))
  call json_add('/results/cipsi/states/rpt2_err', pt2_data_err % rpt2(1:n))
  call json_add('/results/cipsi/states/variance', pt2_data % variance(1:n))
  call json_add('/results/cipsi/states/variance_err', pt2_data_err % variance(1:n))
  if (N_iter > 0) then
    call json_add('/results/cipsi/history/n_iter', N_iter)
    call json_add('/results/cipsi/history/n_iter_max', N_iter_max)
    call json_add('/results/cipsi/states/extrapolated_energy', &
                  extrapolated_energy(N_iter,1:n))
  endif
  ! history/shape and the per-iteration series keys: see §10.7
end
```

`/results/cipsi/states/*` is the final state of the calculation (every
call overwrites it with the current iteration's values), while the
per-iteration series is carried in `/results/cipsi/history/*`, with
`shape` = `[n, nh]` and the flat arrays indexed state-fastest, i.e.
element `(s, i)` at flat index `(i-1)*n + (s-1)` (§10.7).

Notes that a build agent must honour:

* `n = min(N_states, N_det)` follows the current guard (a tiny wave function
  can hold fewer determinants than states).
* `pt2_type` (`src/cipsi_utils/selection_types.f90`) has exactly four
  components: `pt2`, `rpt2`, `variance`, `overlap`. The `_err` values are the
  same components of the `pt2_data_err` object — there are no `rpt2_err` or
  `variance_err` components.
* `history/rpt2` stores `pt2_data % rpt2`, not `pt2` (the quantity QP
  extrapolates); it is the same key the old sliding `N_iter_max = 8`
  window used, with the corrected content of §10.7.
  `extrapolated_energy` is declared `(N_iter_max, N_states)`, so
  `extrapolated_energy(N_iter,1:n)` is the strided section holding the
  current value for each state.
* `history/n_iter` is `N_iter` as maintained by `increment_n_iter`: the number
  of iterations the extrapolation window holds, **not** the number of
  iterations the run performed. It saturates at `N_iter_max = 8` and starts
  from 0 in a new process (`N_iter` is not an EZFIO variable, `src/iterations`
  has no `EZFIO.cfg`). `history/n_iter_max` is logged beside it so a reader can
  tell a saturated window from a short run. QP keeps no lifetime iteration
  counter; if a consumer needs one, it has to be added upstream (§14).
* `n_csf` is **dropped** on purpose. `N_csf` is a provider that diagonalises
  S² per configuration and builds the determinant→CSF transformation
  (`src/csf/csf.irp.f:10`). IRPF90 emits `PROVIDE` for every referenced
  provider at subroutine entry, unconditionally, so keeping `N_csf` here
  would make `program pt2` — which starts calling this routine in §10.1 —
  build the whole CSF transformation for a number nobody needs.
* `N_det` is a plain 4-byte `integer` (`src/determinants/determinants.irp.f:3`),
  `N_Configuration` likewise, so the generic resolves.

### 10.5 CASSCF and orbital optimisation

`casscf_cipsi/casscf.irp.f`: inside the `if (.not. converged)` branch, right
after the `criterion_casscf` block computes `converged` (line 124, before
`if(.not.converged)then ... iteration += 1`), where `e_av`, `Ev`, `E_PT2`,
`PT2` and `norm_grad_vec2` are all live:

```fortran
     call json_add('/results/casscf/iteration', iteration)
     call json_add('/results/casscf/energy', e_av)
     call json_add('/results/casscf/states/energy', Ev(1:N_states))
     call json_add('/results/casscf/states/energy_pt2', E_PT2(1:N_states))
     call json_add('/results/casscf/states/pt2', PT2(1:N_states))
     call json_add('/results/casscf/norm_grad', norm_grad_vec2)
     call json_add('/results/casscf/converged', converged)
```

`mo_optimization_utils/routine_opt_mos.irp.f`, right after `delta_e` is
computed (before the `result_opt` append):

```fortran
      call json_add('/results/mo_opt/iteration', nb_iter)
      call json_add('/results/mo_opt/n_det', n_det)
      call json_add('/results/mo_opt/energy_cipsi', e_cipsi)
      call json_add('/results/mo_opt/energy_opt', e_opt)
      call json_add('/results/mo_opt/delta_energy', delta_e)
```

Both files need `use qp_json, only: json_add` above their `implicit none`
(`casscf.irp.f:30`, `routine_opt_mos.irp.f:12`), and both modules gain the
`json` NEED entry of §9.2. The `json_init` calls need nothing else.

The `mo_optimization/result_opt` text file stays as it is (it is an EZFIO
artifact other scripts read). `iteration` in `casscf.irp.f` is 1-based
(initialised to 1); `nb_iter` is incremented at the *end* of the loop body, so
the value logged here is 0-based and matches the `Cipsi step: <nb_iter>` banner
of that cycle. Say so in the README: consumers must not assume a numbering
convention for `iteration` across methods, and the SCF `iteration` (1-based) is
a different counter from a different loop.

### 10.6 In-tree Python consumers

Three scripts in `scripts/` read these files and break on the new schema. They
are not part of the build, but they must not be left silently returning wrong
answers (their `try/except: pass` blocks swallow every mismatch).

`scripts/utility/qp_json.py` — rewrite as a thin, strict reader:

```python
import os, json, re

_NAME = re.compile(r'^\d{5}\.json$')

def files(ezfio_filename):
    prefix = os.path.join(ezfio_filename, 'json')
    if not os.path.isdir(prefix):
        return []
    return sorted(n for n in os.listdir(prefix) if _NAME.match(n))

def load(filename):
    with open(filename) as f:
        return json.load(f)

def load_all(ezfio_filename):
    prefix = os.path.join(ezfio_filename, 'json')
    out = {}
    for n in files(ezfio_filename):
        try:
            d = load(os.path.join(prefix, n))
        except (OSError, ValueError):
            continue
        if d.get('schema') == 'qp3/1':
            out[n] = d
    return out

def load_last(ezfio_filename):
    names = files(ezfio_filename)
    if not names:
        raise FileNotFoundError(
            os.path.join(ezfio_filename, 'json', 'NNNNN.json'))
    return load(os.path.join(ezfio_filename, 'json', names[-1]))
```

`fix_json` and `fix` are deleted: files are valid by construction, and
`fix_json`'s `s.replace(' ', '')` corrupts the strings of a pretty-printed
document. Zero-padded names sort lexicographically like numerically, and the
`.tmp.<pid>` files are excluded by `_NAME`.

`scripts/qp_extract_cipsi_data.py` — change the per-file access:

```python
for f in sorted(d):
    try:
        j = d[f]['results']['cipsi']
    except KeyError:
        continue
    st = j['states']
    n = len(st['energy'])
    out = f" {j['n_det']:8d}"
    ex = st.get('extrapolated_energy', [0.0]*n)
    for i in range(n):
        out += (f" {st['energy'][i]:16.8f}  {st['pt2'][i]:e}  {st['pt2_err'][i]:e}"
                f"  {st['rpt2'][i]:e}  {st['rpt2_err'][i]:e}  {ex[i]:16.8f}")
    print(out)
```

`scripts/qp_exc_energy.py` — the only change is the record builder at the top
of `read_data`; the fitting code below it is untouched. Rebuild the
per-iteration records from `history`:

```python
data = qp_json.load_last(ezfio_filename)['results']['cipsi']
h = data.get('history')
if h is None:
    st = data['states']
    n = len(st['energy'])
    lines = [{'n_det': data['n_det'],
              'states': [{'energy': st['energy'][i], 'pt2': st['rpt2'][i]}
                         for i in range(n)]}]
else:
    n, m = h['shape']
    lines = [{'n_det': data['n_det'],
              'states': [{'energy': h['energy'][i + j*n],
                          'pt2': h['rpt2'][i + j*n]} for i in range(n)]}
             for j in range(m)]
```

One semantic note for the commit message: the trajectory stored in the
document is `rpt2` (the quantity QP extrapolates), so the script's `pt2` slot
is fed with `rpt2`. That is the same choice `extrapolate_data` makes.

### 10.7 Per-iteration history and symmetry keys

Convention (stated once): one array per quantity, index-aligned; the
iteration number is index + 1 (no `iteration` array, except Davidson, where
the printed step number is not an index). `epoch` = Unix seconds,
`wall_time` = seconds since `human_t0`, `memory` = resident GiB.

```
/results/cipsi/history/   shape[n,nh] energy pt2 rpt2 variance pt2_err   (n,nh)
                          n_det n_cfg pt2_windows pt2_tasks samples
                          pt2_time pt2_epoch_start pt2_epoch_end
                          epoch wall_time memory                      (nh)
                          n_iter n_iter_max   (rolling extrapolation window)
/results/scf/history/     energy delta_energy max_error_diis level_shift
                          dim_diis converged epoch wall_time memory    (nit)
/results/davidson/calls/NNNN/  algorithm n_state epoch itertot converged
                          wall_time shape[n,nstep] iteration energy
                          residual s2[n,nstep]
/results/symmetry/        group n_irrep  states/<s>/irrep  states/<s>/purity
```

* CIPSI. `write_cipsi_json` appends one row per call to SAVE arrays and
  rewrites the history keys every call; the document is the state of the
  run. The series is continuous across the orbital-optimisation cycles of
  `cipsi_orb_opt` (the SAVE arrays persist in one process; there is no
  reset and no cycle index). If the CI space changes under a running
  program (`n = min(N_states,N_det)` differs), the series restarts, by
  design. `samples` is one value per iteration (all states share the
  total). `pt2_err` is `+huge` on the first stochastic iteration by the
  existing seeding in `pt2_stochastic`; that is the record. `history/rpt2`
  holds `rpt2` (the same key the rolling window used, corrected content).
* PT2 counters. `pt2_type` carries `windows`, `tasks`, `samples`, `time`,
  `epoch_start`, `epoch_end`; `pt2_alloc` zeroes them. They are per-call
  scalars, so `pt2_type_size` (which counts the allocatable doubles only)
  is unchanged; `pt2_data_I` holds `N_det_generators` instances, so the
  cost is ~40 B per generator, acceptable. The counters are per-process
  (no MPI broadcast exists). In `pt2_stochastic` the walk is run once per
  state, so `windows`, `tasks` and `samples` are summed over the states:
  they count the work of the whole call, not of the last state. The walk
  also stops as soon as the statistical criterion is met (`do_exit` and
  `conv_crit <= relative_error`, or the last position), so `windows` and
  `tasks` measure the work actually done and stay below the full sweep
  once the determinant space is large enough for the comb; below the comb
  threshold the sweep is exhaustive and one (generator, subset) task is
  run per generator. On the
  deterministic path (`omp_selection`) there is one window pass for all
  states, `samples = 0` (no estimator).
* SCF. The history is appended in the loop and written once after it,
  guarded by `nhist > 0`; scalar series are 1-D, no `shape` key.
* Davidson. One block per solve; `calls/NNNN` is monotone for the whole
  process and shared by the four variants (`davidson_diag_hjj_sjj`,
  `davidson_diag_nonsym_hjj`, `davidson_diag_csf_hjj`,
  `davidson_general_ext_rout_diag_dressed`). A skipped step is neither
  printed nor recorded. The state-following retry -- the routine returns
  when not every state can be attributed, so that the caller may raise
  `N_st_diag` -- closes its block before returning, with `converged =
  false`, its `itertot` and the steps recorded so far: an aborted attempt
  is a call, and if it were dropped the counter would show a hole. On the
  `stop -1` divergence of
  `dav_diag_dressed_ext_rout` a half-written block (algorithm, epoch,
  itertot, converged, wall_time) is acceptable. The `energy` of the
  generic variant excludes the nuclear repulsion: the record says what is
  printed, it does not "fix" it.
* Symmetry. `results/symmetry/group` and `n_irrep` come from
  `pg_report_mo_symmetry_work` (SCF and symmetry runs) and
  `pg_print_state_irreps_work` (CI runs); `states/<s>/irrep` reuses the
  exact variable the table prints, so JSON == stdout.

## 11. Reconstruction from several files

The rule consumers must implement:

1. Enumerate `<EZFIO>/json/[0-9][0-9][0-9][0-9][0-9].json`, sorted by number.
   Ignore `.tmp.<pid>` and anything without `"schema": "qp3/1"`.
2. Take the **highest-numbered** file as the state of the calculation. It is
   complete: `results` carries every number the run computed, including
   in-process SCF for an `fci` run, and `input` carries what the restart
   actually read.
3. Never merge fields from two files. Older files are history only: the
   sequence of `results.cipsi.states` across files, ordered by number and
   labelled by `run.timestamp`, is the CIPSI trajectory across restarts.
4. `run.kind` tells which program produced a file, `provenance.ezfio_file`
   which calculation it belongs to (relevant when several EZFIO files are
   compared), `run.elapsed` how long that run took.

Because a restart re-reads the wave function from EZFIO, the new file's
`input` block records what the restart actually used, and its first
`results.cipsi` block is the state at restart. No cross-file bookkeeping is
needed in QP itself.

## 12. Failure modes and guarantees

| event | file on disk |
|---|---|
| normal exit, `stop 1` | complete up to the last `json_add` (atexit forces a write) |
| `SIGKILL`, `SIGTERM`, `SIGSEGV`, OOM killer | last throttled write: valid, at most 1 s old |
| disk full / `rename` fails | previous version intact, one warning on stderr, run continues |
| NaN or Inf in a logged value | `null` at that key, document still written |
| bad path in a call site | that value is absent, one warning, run continues |
| non-master MPI rank | no thread, no document, no I/O |
| `json_init` never called | nothing is written, one warning at exit |
| two QP processes on one EZFIO file | both pick the same number and race on it; each write is still atomic and valid, so the file always decodes, but the last process wins. Concurrent runs on one EZFIO file are already unsupported (they corrupt the EZFIO data itself). |

The run is never stopped by the logging layer. A consumer must treat missing
keys as "not available" and never assume a key exists.

## 13. Phasing

Each phase leaves the tree building and a working program.

**P0 — vendored library, shim, Fortran module. No call sites.**

Files: `yyjson.c`, `yyjson.h`, `qp_json_shim.c`, `qp_json.f90`. Acceptance in
a throwaway directory under `$HOME/OpenCode/tmp/` (not committed), compiled by
hand with the same flags as the build (`gcc -Ofast`, `gfortran -std=legacy
-ffree-line-length-none`):

1. round-trip test: log `0.1d0`, `1.d0/3.d0`, `1.d-300`, `huge(1.d0)`,
   `-0.d0`, `1.d-308`, read back with `jq` and compare to `%.17g` — this is
   the `-Ofast` gate for yyjson's `dtoa` (§4);
2. NaN/Inf land as `null`;
3. a 3-second loop of `json_add` (thousands of calls) rewrites the file about
   three times — at most once per second — not once per call;
4. `kill -9` during that loop leaves `jq . <file>` succeeding;
5. `QP_JSON_MIN_INTERVAL=0` writes on every call and starts no thread;
6. `qp_json_begin(..., 0)` produces no file at all.

**P1 — glue and SCF.** Rewrite `json.irp.f` (add `json_init`, grow
`json_filename` to 2048 and move the `ezfio_set_json_empty` call into it;
leave the old `json_unit`/`json_close` providers in place), set `json/NEED` to
`ezfio_files mpi`, port `scf.irp.f` and `roothaan_hall_scf.irp.f`.
`json_formats.irp.f` stays until P2, so `program fci`, `mrci`, `ks_scf` and
`rs_ks_scf` still build and still behave as before. Acceptance: `qp_run scf
water.ezfio` then `jq -e . water.ezfio/json/00001.json` succeeds and
`.results.scf.energy` equals the energy printed by QP.

**P2 — CIPSI, CASSCF, orbital optimisation, KS.** `write_cipsi_json`,
`cipsi.irp.f`, `stochastic_cipsi.irp.f`, `fci.irp.f`, `pt2.irp.f`,
`mrci.irp.f`, `casscf.irp.f`, `routine_opt_mos.irp.f`, `ks_scf.irp.f`,
`rs_ks_scf.irp.f`, plus the `json_init` sites and the `NEED` entries those
modules require. Once nothing references the old writers any more
(`grep -rn 'json_unit\|json_formats' src --include='*.irp.f'` returns
nothing), delete `json_formats.irp.f`, `json_unit` and `json_close`.
Acceptance: an `fci` run on a small system, `jq
'.results.cipsi.states.energy'` matches the printed energies for every state,
the file survives `kill -9` mid-run, `qp_run pt2` produces
`/results/cipsi/states/pt2` without computing `N_csf` (compare the wall time
with the previous build on the same input), and the run's peak RSS is
recorded next to the number of CIPSI iterations so §7.9's estimate can be
checked against reality.

**P3 — input capture.** Hook in `ezfio_generate_provider.py`, then the 29
`NEED` edits of §9.2 (plus `mpi` in `src/json/NEED`). Acceptance: after a build,
`grep -r 'call json_add' src/*/ezfio_interface.irp.f | wc -l` counts 201
calls (`git grep` sees none of them: the generated files are gitignored); a
real run shows `/input/<group>/<key>` for each input printed in the QP header,
plus the hand-logged `/input/ezfio_files/ezfio_convention`.

**P4 — tools and documentation.** The three Python files of §10.6 and
`src/json/README.rst`: schema, units, the atomic-replace rule,
`QP_JSON_MIN_INTERVAL`, the reconstruction rule of §11, yyjson version and
MIT licence, the `-Ofast` round-trip measurement.

## 14. Deferred

* Build-time `creator.version` (a generated `git describe` string) — the
  field exists and reads `"unknown"` until then.
* Append-only per-iteration trace (`jsonl`), for people who want to plot a
  convergence that spans several restarts. The per-run document now carries
  the whole per-iteration series of the run in its history (§10.7).
* PT2 counters are per-process (no MPI reduction exists); in
  `pt2_stochastic` the walk runs once per state, so `windows`, `tasks` and
  `samples` are summed over the states (§10.7).
* A Davidson block left half-written by an abnormal exit stays
  half-written (§10.7).
* `n_csf` (§10.4) and array-valued EZFIO inputs (43 of 246 providers), if a
  consumer ever asks for them.
* `/results/<method>/energy` for `cis`, `cisd`, `ccsd`, `mp2`, which log
  nothing today.
* Directory `fsync` for power-loss durability.
* `utils/restore_symm` is declared in `utils/EZFIO.cfg` (default `.False.`) but
  no source file reads it: the `restore_symmetry` subroutines take their
  arguments directly. Either wire it into `restore_symmetry` or drop it from
  the config; until then it is not an input anyone inspected, so it stays out
  of the document (§9.2).
* A lifetime CIPSI iteration counter: `history/n_iter` saturates at
  `N_iter_max = 8` (§10.4). Adding a real one means a new counter in
  `increment_n_iter`, not a change to the logging layer.
* Unit tests for the Fortran wrappers: QP3 has no harness for individual
  Fortran procedures, so P0 stays a manual recipe. The bats suite
  (`tests/bats/common.bats.sh`, one `NN.name.bats` per program) runs whole
  programs and is the natural place for a `jq -e` regression check on
  `results.scf.energy` (`10.hf.bats`), the CIPSI state arrays (`40.fci.bats`)
  and the RS decomposition (`61.rsks.bats`), once the schema stops moving.

## 15. File inventory

| action | file |
|---|---|
| add | `src/json/yyjson.c`, `src/json/yyjson.h` (yyjson 0.13.0, MIT) |
| add | `src/json/qp_json_shim.c` (~300 lines) |
| add | `src/json/qp_json.f90` (~150 lines) |
| add | `src/json/DESIGN.md` (this document) |
| delete | `src/json/json_formats.irp.f` |
| rewrite | `src/json/json.irp.f`, `src/json/README.rst` |
| edit | `src/json/NEED` (`+ mpi`) |
| edit | `scripts/ezfio_interface/ezfio_generate_provider.py` |
| edit | `NEED` files: the 28 modules of §9.2, plus `kohn_sham_rs` (§10.3) |
| edit | 12 Fortran call-site files (§10.1) |
| edit | `scripts/utility/qp_json.py`, `scripts/qp_extract_cipsi_data.py`, `scripts/qp_exc_energy.py` (§10.6) |
