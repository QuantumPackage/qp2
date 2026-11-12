# qp-cli — notes for the next agent (P4, wave 3)

Status of this crate, 2026-11-12. It exists because `qp_run` needs the four
pieces of the OCaml tools that `qp_run` shared with no other crate:

| module | state | tests |
|---|---|---|
| `root` | **done**, `Qpackage.root` + #224 fallback | 3 |
| `executables` | **done**, `data/executables` parser | 3 |
| `launch` | **done**, `Sys.command` equivalent (job line, exit codes) | 5 |
| `ffi` | **done**, the only two `unsafe` blocks of the workspace | 2 |
| `git` | **done**, banner stamp from the environment | 2 |
| `signal` | **done**, SIGUSR1 poll + `qp_stop` orchestration | 0 (e2e) |
| `stop` | **done**, `ignore @@ Sys.command ("qp_stop "^ezfio)` | 0 (e2e) |
| `lib` | **done**, `Error`/`Result` | 1 |

16 unit tests. `signal` and `stop` have none on purpose: they need a
private `PATH` and a process-global signal flag, so `bins/qp_run/tests/`
covers them end to end (`signal.rs`, 4 tests: a real child, a `qp_stop`
shim, the real `bin/qp_stop` where it can run, and a job killed by a
signal).

The `qp_run` binary itself is **not here**: it lives in
`rust/bins/qp_run` (main flow, 15 integration tests). This crate holds
only what a second tool (`qp_edit`, `qp_create_ezfio`, …) will also need.
Nothing in `qp-cli` knows the job-launch details or the command line.

Build and test what exists:

```
cd rust
cargo test -j4 -p qp-cli
cargo clippy -j4 -p qp-cli --all-targets -- -D warnings
cargo fmt -p qp-cli -- --check
```

2026-11-12 result: 16 tests pass, clippy clean, fmt clean.

## The two `unsafe` blocks (PI decision, wave 3)

The wave gate says "zero own `unsafe`". `qp_run` cannot be ported faithfully
without these two calls: `std` has no local time and no signal API, and the
plan allows no extra crate for them (`chrono` would only solve the first).
The PI approved exactly two blocks, both in `ffi.rs`, each with its own
`// SAFETY:` comment and no `unsafe` anywhere else in the workspace (tests
included). Do not add a third; if a new call seems to need one, stop and ask.

| function | reference site | `std` gap |
|---|---|---|
| `local_time` | `Unix.localtime`, `qp_run.ml:102` | no local-time conversion |
| `install_sigusr1` | `Sys.set_signal Sys.sigusr1`, `qp_run.ml:88` | no signal API |

The former `hostname` (`Unix.gethostname`, `TaskServer.ml:89`) and
`ipv4_for_interface` (`QP_NIC`) blocks lived with the multi-node layer
(task-server address) and were removed with it, 2026-11.

Invariants each block relies on:

1. **`localtime_r`** writes only through `&mut broken`, which is an owned,
   zeroed local; the return pointer is null-checked and `broken` is discarded
   on failure. `time` is an owned `time_t`. `localtime_r` is the reentrant
   form, so no static buffer is read concurrently.
2. **`sigaction`** takes a zeroed owned `sigaction` (padding), an
   `extern "C" fn(c_int)` matching `sighandler_t`, an `sa_mask` initialised by
   `sigemptyset` on the same object, and `SA_RESTART` (`waitpid` must not fail
   with `EINTR`). The old action is not needed (null third argument).

**Signal safety.** The handler does exactly one relaxed `AtomicBool::store`
and nothing else (async-signal-safe). The reference runs
`Sys.command ("qp_stop "^ezfio)` *inside* its handler; `system(3)`/`fork(3)`
are not async-signal-safe, so a handler that allocates and forks can
deadlock on the allocator or on an internal lock. The main thread polls
`take_sigusr1()` every 200 ms and runs `qp_stop` from normal code. `take_sigusr1` is a `swap`, so a signal arriving between "read" and
"clear" is not lost.

## Deviations from the reference (all deliberate)

1. **SIGUSR1 polling.** `qp_stop` starts up to 200 ms after the signal and
   runs on the main thread (delaying job-exit detection by its duration). Two
   signals inside one interval collapse into one `qp_stop` (idempotent
   `touch`). A signal after the job's exit is dropped. See `signal.rs`.
2. **Git stamp from the environment** (`QP_GIT_SHA1`, `QP_GIT_DATE`,
   `QP_GIT_MESSAGE`, default `unknown`). The reference bakes the values into
   `Git.ml` at build time; baking them into a Rust binary would print the
   build-time commit of a *copied* binary and make the build
   non-reproducible. `qp_run`'s banner therefore shows `unknown` unless the
   wrapper exports the three variables.
3. **`Sys.command` is `/bin/sh -c`** with inherited standard streams, in
   both `launch` and `stop`. Exit status is `status.code()`, or 255 when the
   *shell* was signalled (`caml_sys_system_command`). A grandchild killed by a
   signal is reported by `/bin/sh` as `128+n`, where `system()` would report
   255; `qp_run` only tests `= 0`, so the difference is invisible. The thin
   launcher passes no extra environment to the job (the task-server address
   variables are gone).
4. **`$QP_ROOT` fallback (issue #224).** The reference reads the variable and
   fails when unset. `root()` reads it first and, when unset or empty, takes
   the first ancestor of `current_exe` (4 levels, `SEARCH_DEPTH`) containing
   `data/executables`. No path is baked at compile time. An empty variable is
   treated as unset (the reference would use `""` and fail opening
   `data/executables` with a different message).
5. **Error text.** Rust tools print `<tool>: <message>` on stderr with the
   reference's exit status (2 for `failwith`, 1 for the two explicit `exit 1`
   sites) instead of OCaml's `Fatal error: exception Failure("…")`. Every
   *message* is the reference's string, including the second line of the
   `QP_ROOT` failure.
6. **`Random.self_init ()` is not ported.** The reference calls it at module
   initialisation; nothing in `qp_run` draws a random number, so the call has
   no observable effect and `std` needs no seed.

## The test seam (in `bins/qp_run/src/main.rs`)

| variable | effect | why |
|---|---|---|
| `QP_RUN_QP_EDIT` | program substituted in `"<editor> -c <ezfio>"`; default `qp_edit` | the `-c` pre-check must be an injectable subprocess: the reference `qp_edit` is not guaranteed on `PATH`, and a test needs to observe the call and force a non-zero status |

It is read once, by the binary, and documented on `main.rs`'s module
header. `qp-cli` does not read it. Nothing else in the environment is a
seam: `QP_ROOT` is the tree root, and the reference's address variables
(`QP_NIC`, `QP_TASK_DEBUG`, `QP_RUN_ADDRESS*`) went with the multi-node
layer, 2026-11. (The former `QP_RUN_MAX_SCANS` seam capped the port-window
scan; the scan is gone.)

The integration tests build a throwaway tree: a temp `$QP_ROOT` with
`data/executables`, a `bin/` of shim scripts, and an EZFIO directory. The
child's output goes to files, not pipes: a hung job would hold a pipe open
and defeat the hard deadline that kills and reaps the child with its output
attached.

## Reference facts worth not re-deriving

- **History (multi-node layer removed 2026-11):** `qp_run.ml:39-77` created
  the port window **before** the EZFIO existence check (`qp_run.ml:83`), so
  a bad directory still failed after a bind scan. The port was kept by the
  pre-check sockets themselves (they became the REP and PUB sockets), which
  was the #393 fix; the probes `P+2..P+9` were bound and released.
- `qp_run.ml:88-89` installs SIGUSR1 before the executable lookup, i.e. even
  a run that dies on `Executable … not found` has the handler installed.
- `qp_run.ml:102` prints `tm_hour + isdst` (a manual DST adjustment);
  `ffi::Tm::hour` is the raw `tm_hour`, and `main.rs` adds `isdst`.
- `qp_run.ml:121` runs `qp_edit -c` on the input; `qp_run.ml:195` runs the
  executable through `Sys.command`, so a program with a space in its path
  works on both sides only because the `data/executables` path is
  concatenated unquoted.
- `qp_stop` is `$QP_ROOT/bin/qp_stop`, a bash script touching
  `<ezfio>/work/qpkill`; it is not ported (the file protocol belongs to the
  Fortran half), so `stop::qp_stop` runs the bare name through `PATH` and
  ignores the status, like `ignore @@ Sys.command`.
- `Sys.command` on a missing program is `/bin/sh`'s 127; the reference prints
  `sh: …: command not found` on the inherited stderr and continues for
  `qp_stop`, and fails the run for the executable (non-zero).

## h2 SCF interop, measured 2026-11-05 (do not re-derive the numbers)

Two `cp -R` copies of `rust/tests/golden/fixtures/h2`, `n_it_scf_max` 50,
`source ~/qp2/quantum_package.rc`, `OMP_NUM_THREADS=1 QP_NIC=lo0`; reference
`~/qp2/ocaml/qp_run.native scf`, Rust `rust/target/debug/qp_run scf` with
`QP_ROOT=~/qp2` (same `scf`, same pre-fix `qp_edit`). Both exit 0, energy
`-1.005924963287844` (target `-1.005924963288527`, 6.8e-13), 284 decompressed
EZFIO files identical, 289-file raw list identical, 115-line stdout identical
after the harness `NOISE` and path normalisation, `work/qp_run_address` gone
in both. Wall 3.42 s (reference) vs 3.51 s (Rust), user/sys 0.07/0.11 s
both; `scf` itself reports 1.04 s wall in both streams.

`qp_edit -c` on h2 (pre-fix native) creates only `work/.version` and
`work/empty` (`F\n\n`), no data file: `contracts/check-mutated.json` lists
`[]` for h2, confirmed empirically.

**Environment problem, not worked around:** `~/qp2/bin/qp_reset` calls GNU
`getopt -l`; macOS BSD `getopt` turns the argv into garbage and the script
exits 2 with `EZFIO directory not found`. `qp reset --mos` therefore cannot
run on this macmini. The h2 fixture is already in the `--mos` state (no
`mo_basis`, `io_*` = `None`), so the reset is a no-op; the gate skipped it.
Do not patch `~/qp2`.

## Where the plan/brief is wrong

- The P4 brief says `qp_run` uses `$QP_PORT` (default 43210). There is no
  `$QP_PORT`; the window is `10 * (pid mod 2823) + 32769`, `+100` per retry.
  Verified in `qp_run.ml`. (The port window itself is gone with the
  multi-node layer, 2026-11-12; this section stays as plan history.)
- The brief said the modulus is 2821; the OCaml is 2823.
- The #393 test wording ("pre-bind the window, assert a *fail-fast* with an
  unbounded scan") is self-contradictory: an unbounded scanner retries
  forever by design. The test uses `QP_RUN_MAX_SCANS=1` for the failure arm
  and no cap for the "moves to P+100" arm.
- The wave gate "zero own `unsafe`" cannot hold for a faithful `qp_run`; the
  PI approved the two blocks above. Report, do not work around.

## Formatting gate

`cargo fmt -p qp-cli -- --check` and `cargo fmt -p qp_run -- --check` are
clean. `cargo fmt --all -- --check` fails on
`rust/crates/qp-types/src/lib.rs:86` (`pub mod gto;` ordering), inherited
from wave 3's predecessors; it was left alone on purpose.
