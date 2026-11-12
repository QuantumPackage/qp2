//! The `Sys.command` equivalent used to run jobs and helper tools.
//!
//! `Unix.system` runs `/bin/sh -c "<command>"`, inherits the standard streams,
//! and returns `WEXITSTATUS`; a child killed by a signal returns `255`
//! (`caml_sys_system_command`). This module reproduces that: the command
//! string is passed verbatim to `/bin/sh -c`, so shell quoting, `PATH`
//! lookup, redirections and a missing program (exit `127`) behave exactly as
//! in the reference, including the `gdb --args ` prefix of `qp_run -p gdb`.
//!
//! # Why a shell and not `Command::new(program)`
//!
//! `qp_run` builds `<prefix><program path> <ezfio>` as one string
//! (`ocaml/qp_run.ml:185-195`); the prefix may itself contain arguments
//! (`gdb --args `), so splitting it here would re-implement the shell. The
//! price is one `/bin/sh` per job, which is what the reference pays too.

use std::process::{Child, Command, ExitStatus};
use std::time::Duration;

use crate::Result;

/// The shell `Unix.system` uses.
pub const SHELL: &str = "/bin/sh";

/// Exit status reported for a child that was killed by a signal
/// (`WIFSIGNALED` in `caml_sys_system_command`).
pub const SIGNALED_EXIT_CODE: i32 = 255;

/// The command line `qp_run` runs for one job (`ocaml/qp_run.ml:185-195`).
///
/// `prefix` is the `-p` value: `Some("gdb")` becomes `gdb --args `, any other
/// value becomes `<value> `, and `None` becomes the empty string. `program` is
/// the path from `data/executables` and `ezfio` the directory, both
/// concatenated without quoting, exactly like the reference.
#[must_use]
pub fn job_command(prefix: Option<&str>, program: &str, ezfio: &str) -> String {
    let prefix = match prefix {
        Some("gdb") => "gdb --args ",
        Some(other) => return format!("{other} {program} {ezfio}"),
        None => "",
    };
    format!("{prefix}{program} {ezfio}")
}

/// Spawn `command` through `/bin/sh -c`, inheriting the standard streams.
///
/// `env` pairs are added to the child's environment only; the parent's
/// environment is inherited and untouched. The reference used the same
/// mechanism for `QP_RUN_ADDRESS`, with `Unix.putenv` on the parent,
/// which would need `unsafe` in edition 2024 and would race the
/// SIGUSR1 handler; the thin launcher passes no env pairs.
///
/// # Errors
///
/// [`crate::Error::Io`] when the shell cannot be spawned.
pub fn spawn_command(command: &str, env: &[(&str, &str)]) -> Result<Child> {
    let child = Command::new(SHELL)
        .arg("-c")
        .arg(command)
        .envs(env.iter().copied())
        .spawn()?;
    Ok(child)
}

/// Run `command` to completion and return its exit status code
/// (`Sys.command` without polling).
///
/// # Errors
///
/// [`crate::Error::Io`] when the shell cannot be spawned or waited for.
pub fn run_command(command: &str, env: &[(&str, &str)]) -> Result<i32> {
    let mut child = spawn_command(command, env)?;
    Ok(exit_code(child.wait()?))
}

/// The reference's status mapping: the exit code, or
/// [`SIGNALED_EXIT_CODE`] when the shell itself was killed.
///
/// A *grandchild* killed by a signal is reported by `/bin/sh` as `128 + n`,
/// which is what the shell exits with; the reference's `system()` reports the
/// same case as `255` because it observes the shell directly. Both are
/// non-zero, and `qp_run` only tests `= 0`.
#[must_use]
pub fn exit_code(status: ExitStatus) -> i32 {
    status.code().unwrap_or(SIGNALED_EXIT_CODE)
}

/// Wait for `child`, calling `on_tick` every `interval` while it runs.
///
/// `qp_run` needs this instead of a blocking `wait` because SIGUSR1 must run
/// `qp_stop` while the job runs (`ocaml/qp_run.ml:88`). The cost is one
/// `try_wait` plus `on_tick` every `interval`, i.e. 5 wake-ups per second at
/// the default 200 ms, and up to 200 ms of extra latency before the job's
/// exit is observed.
///
/// # Errors
///
/// [`crate::Error::Io`] when `try_wait` fails.
pub fn wait_with_ticks(
    child: &mut Child,
    interval: Duration,
    on_tick: &mut dyn FnMut(),
) -> Result<ExitStatus> {
    loop {
        if let Some(status) = child.try_wait()? {
            return Ok(status);
        }
        on_tick();
        std::thread::sleep(interval);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The prefix mapping of `qp_run.ml:185-189`.
    #[test]
    fn job_command_prefix() {
        assert_eq!(
            job_command(None, "/qp/src/hartree_fock/scf", "/tmp/h2"),
            "/qp/src/hartree_fock/scf /tmp/h2"
        );
        assert_eq!(
            job_command(Some("gdb"), "/qp/scf", "/tmp/h2"),
            "gdb --args /qp/scf /tmp/h2"
        );
        assert_eq!(
            job_command(Some("valgrind"), "/qp/scf", "/tmp/h2"),
            "valgrind /qp/scf /tmp/h2"
        );
    }

    /// Exit codes: `0`, an arbitrary code, and a shell killed by a signal.
    #[test]
    fn exit_codes_are_propagated() {
        assert_eq!(run_command("exit 0", &[]).expect("exit 0"), 0);
        assert_eq!(run_command("exit 3", &[]).expect("exit 3"), 3);
        assert_eq!(
            run_command("kill -TERM $$", &[]).expect("signalled"),
            SIGNALED_EXIT_CODE
        );
    }

    /// A program that does not exist is the shell's 127, as `Sys.command`
    /// reports it.
    #[test]
    fn missing_program_is_127() {
        assert_eq!(
            run_command("qp_no_such_program_0a1b2c", &[]).expect("missing"),
            127
        );
    }

    /// The extra environment pairs reach the child and only the child.
    #[test]
    fn environment_pairs_are_added() {
        let dir = tempfile::tempdir().expect("tempdir");
        let out = dir.path().join("env.txt");
        let command = format!("printenv QP_CLI_TEST_MARKER > {}", out.display());
        assert_eq!(
            run_command(&command, &[("QP_CLI_TEST_MARKER", "tcp://1.2.3.4:5")]).expect("run"),
            0
        );
        assert_eq!(
            std::fs::read_to_string(&out).expect("read").trim(),
            "tcp://1.2.3.4:5"
        );
        assert!(
            std::env::var("QP_CLI_TEST_MARKER").is_err(),
            "the parent environment is untouched"
        );
    }

    /// The polling wait returns as soon as the child is gone and ticks at
    /// least once on a child that outlives one interval.
    #[test]
    fn wait_with_ticks_polls() {
        let mut child = spawn_command("sleep 0.3", &[]).expect("spawn");
        let mut ticks = 0_u32;
        let status = wait_with_ticks(&mut child, Duration::from_millis(50), &mut || {
            ticks += 1;
        })
        .expect("wait");
        assert_eq!(exit_code(status), 0);
        assert!(ticks >= 1, "ticks = {ticks}");
        // Once the status has been reaped, `std` caches it and the polling
        // wait returns without a tick. Reaping here first is deterministic,
        // unlike spawning `true` and hoping it has already exited.
        let mut child = spawn_command("true", &[]).expect("spawn");
        while child.try_wait().expect("try_wait").is_none() {
            std::thread::sleep(Duration::from_millis(10));
        }
        let mut again = 0_u32;
        let status = wait_with_ticks(&mut child, Duration::from_millis(50), &mut || {
            again += 1;
        })
        .expect("wait");
        assert_eq!(exit_code(status), 0);
        assert_eq!(again, 0, "the cached status needs no tick");
    }
}
