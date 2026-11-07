//! SIGUSR1 → `qp_stop` while the job runs (`ocaml/qp_run.ml:87-89`).
//!
//! `qp_run` installs a SIGUSR1 handler before launching the job; the handler
//! runs `qp_stop <ezfio>`, which creates `<ezfio>/work/qpkill` and makes the
//! Fortran executable stop at its next polling point. A user therefore
//! interrupts a run with `kill -USR1 <pid of qp_run>`.
//!
//! # Why the handler only sets a flag (documented deviation)
//!
//! The reference runs `Sys.command ("qp_stop "^ezfio)` *inside* the handler.
//! `system(3)`/`fork(3)` are not async-signal-safe, and `qp_run` keeps several
//! threads alive (the task server, its PUB thread, the 0MQ I/O threads), so a
//! handler that allocates and forks can deadlock on the allocator or on the
//! 0MQ internal locks. Here the handler does one relaxed atomic store
//! ([`crate::ffi`]) and the main thread polls it every [`POLL_INTERVAL`] while
//! it waits for the job, then runs `qp_stop` from normal code.
//!
//! Consequences, all acceptable for an interrupt:
//!
//! * `qp_stop` starts up to [`POLL_INTERVAL`] after the signal, not
//!   immediately.
//! * It runs on the main thread, so job-exit detection is delayed by its
//!   duration (a `touch` behind one `/bin/sh`).
//! * Two signals in the same interval collapse into one `qp_stop`; the
//!   reference would run one per signal. Idempotent by construction (`touch`).
//! * A signal that arrives between the job's exit and the process's exit is
//!   dropped instead of running `qp_stop` after the job is gone.
//!
//! # Test coverage
//!
//! This module has no unit test: the flag is a process-global static shared
//! with [`crate::ffi`]'s SIGUSR1 test, and `cargo test` runs the two in
//! threads of one process, so each could consume the other's signal. The
//! behaviour is covered end to end by `bins/qp_run/tests/signal.rs`, which
//! signals a real child process (its own flag) and observes a `qp_stop` shim
//! on its `PATH`.

use std::process::ExitStatus;
use std::time::Duration;

use crate::launch;
use crate::{Error, ffi, stop};

/// How often the main thread checks for SIGUSR1 while the job runs.
///
/// 200 ms bounds the interrupt latency and costs five `try_wait` calls per
/// second — the reference pays nothing here because its handler runs at signal
/// delivery, but it cannot be made async-signal-safe (see the module docs).
pub const POLL_INTERVAL: Duration = Duration::from_millis(200);

/// Install the SIGUSR1 handler (`Sys.set_signal Sys.sigusr1`).
///
/// Call once, before the job is spawned. The flag is not cleared: a signal
/// that arrives between this call and the first poll is still seen.
///
/// # Errors
///
/// [`Error::Os`] when `sigaction` fails.
pub fn install() -> crate::Result<()> {
    ffi::install_sigusr1()?;
    Ok(())
}

/// Wait for `child`, running `qp_stop <ezfio>` each time SIGUSR1 was received.
///
/// The reference's handler replaces the blocking `Sys.command` of the job, so
/// this is the port of "wait for the job while being interruptible": it polls
/// every [`POLL_INTERVAL`] instead of blocking, which is the deviation
/// documented in the module header. `ezfio` is the EZFIO directory passed to
/// `qp_stop` verbatim.
///
/// # Errors
///
/// [`Error::Io`] when `try_wait` fails. A `qp_stop` failure is printed as
/// `qp_run: <message>` and does not stop the job: the reference's handler
/// ignores the command's status and cannot fail either.
pub fn wait_for_job(child: &mut std::process::Child, ezfio: &str) -> Result<ExitStatus, Error> {
    launch::wait_with_ticks(child, POLL_INTERVAL, &mut || {
        if ffi::take_sigusr1()
            && let Err(error) = stop::qp_stop(ezfio)
        {
            eprintln!("qp_run: {error}");
        }
    })
}
