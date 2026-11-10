//! The `qp_stop` command that SIGUSR1 runs (`ocaml/qp_run.ml:89`).
//!
//! `qp_stop` is the reference bash script `$QP_ROOT/bin/qp_stop`, which
//! touches `<ezfio>/work/qpkill`; the running Fortran executable polls that
//! file and stops. The script is not ported (its `qpkill`/`qpstop` file
//! protocol belongs to the Fortran half), so [`qp_stop`] runs the shell
//! command [`qp_zmq::port_window::qp_stop_command`] builds — a bare `qp_stop`
//! resolved through `PATH`, with the EZFIO directory concatenated without
//! quoting, exactly like the reference.
//!
//! # Exit status
//!
//! The reference calls `ignore @@ Sys.command (...)`: the status is
//! discarded, so a missing `qp_stop` (`/bin/sh` exits `127`, printing
//! `qp_stop: command not found` on the inherited stderr) is not an error
//! either. Only a failure to *start* the shell is reported, which is the
//! reference's `Sys_error`.
//!
//! There is no unit test here: the command needs a `qp_stop` on `PATH`, and
//! only a child process can have a private `PATH` (edition 2024 forbids
//! `std::env::set_var`). `bins/qp_run/tests/signal.rs` covers it end to end
//! with a shim script.

use crate::Result;
use crate::launch;

/// Ask the running job to stop (`ignore @@ Sys.command ("qp_stop "^ezfio)`).
///
/// Runs synchronously: the shell and the script are waited for, so this takes
/// as long as the script does (a few milliseconds for `touch`). The SIGUSR1
/// poller blocks on it, which delays the next job-exit check by the same
/// amount ([`crate::signal`]).
///
/// # Errors
///
/// [`crate::Error::Io`] when `/bin/sh` cannot be spawned; the script's exit
/// status is deliberately ignored, as in the reference.
pub fn qp_stop(ezfio: &str) -> Result<()> {
    let command = qp_zmq::port_window::qp_stop_command(ezfio);
    launch::run_command(&command, &[]).map(|_status| ())
}
