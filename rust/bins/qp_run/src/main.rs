//! `qp_run` — the Quantum Package launcher.
//!
//! `qp_run` is the launcher every user types: it checks the input, prints
//! the banner, and runs one Quantum Package executable.
//!
//! The multi-node layer (task server, port window, `work/qp_run_address`)
//! was removed in 2026-11: parallelism is now in-process OpenMP inside
//! each executable, so the launcher only launches.
//!
//! # What it computes
//!
//! Nothing physical: it is a launcher. The numbers it produces are the wall
//! time of the run (whole seconds, printed as `h:mm:ss`) and the exit status
//! of the executable it ran.
//!
//! # Flow (order is observable; do not reorder)
//!
//! 1. `QP_ROOT` (environment, else derived from the executable path —
//!    [`qp_cli::root`]) → `data/executables`. Resolved before the command
//!    line, as the reference's `Qpackage.root` module initialiser runs.
//! 2. Command line ([`qp_types::command_line::RunArgs`]); `-h`/`--help` is
//!    `clap`'s built-in (help once, exit 0), a wrong number of positionals
//!    is a `clap` usage error (exit 2).
//! 3. EZFIO directory exists (`the EZFIO directory does not exist: <d>`).
//! 4. SIGUSR1 handler (→ `qp_stop`, see [`qp_cli::signal`]).
//! 5. Program in `data/executables`, else `Possible choices:` + failure.
//! 6. Banner: `Date:` (local time), Quantum Package, Git stamp, EZFIO dir.
//! 7. `qp_edit -c <dir>`; non-zero → `the input is inconsistent
//!    (qp_edit -c failed)`.
//! 8. The executable runs with its inherited environment (no extra
//!    variables).
//! 9. `Program exited with code %d.` (non-zero only), `Wall time:`,
//!    one-second pause, exit with the executable's status.
//!
//! # Environment variables
//!
//! | variable | meaning |
//! |---|---|
//! | `QP_ROOT` | QP2 tree (else derived from the executable) |
//! | `QP_GIT_SHA1`/`QP_GIT_DATE`/`QP_GIT_MESSAGE` | banner stamp, `unknown` when unset |
//! | `QP_RUN_QP_EDIT` | program used for the `-c` pre-check (default `qp_edit`) |
//!
//! The last one is a documented test seam; nothing else in the tool reads
//! it. It exists because the reference's `qp_edit` path cannot otherwise be
//! exercised from a test.
//!
//! # Deviations from the reference
//!
//! * `Sys.command` is `/bin/sh -c` with inherited standard streams
//!   ([`qp_cli::launch`]), so the `Program exited with code %d.` line and
//!   the exit code match, including the `255` of a signalled shell.
//! * SIGUSR1 sets a flag; `qp_stop` runs from the main thread at most 200 ms
//!   later instead of inside the handler (see `qp_cli::signal`).
//! * The Git stamp comes from the environment instead of being compiled in.
//! * Error text: `<tool>: <message>` on stderr with the reference's exit
//!   status (2 for a `failwith`) instead of
//!   `Fatal error: exception Failure("…")`.
//! * A broken stdout pipe panics Rust's `println!` (exit 101) where the
//!   reference raises `Sys_error` (exit 2).
//!
//! # Running it
//!
//! ```text
//! cargo run -p qp_run -- scf /path/to/h2.ezfio
//! cargo test -p qp_run -j4 -- --test-threads=4         # see tests/
//! ```

#![deny(missing_docs)]

use std::io::Write;
use std::path::Path;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use clap::Parser;
use qp_cli::{GitStamp, Programs, launch, signal};
use qp_types::command_line::{RunArgs, RunConfig};

/// Exit status of a fatal error (an uncaught `Failure`/`Sys_error` in the
/// reference, `CONTRACTS.md` §2.4).
const EXIT_FATAL: i32 = 2;

/// Program used for the `qp_edit -c` pre-check. Test seam; the reference
/// hard-codes `qp_edit`.
const QP_EDIT_VAR: &str = "QP_RUN_QP_EDIT";

/// The reference's `qp_edit` path in `Sys.command ("qp_edit -c "^dir)`.
const DEFAULT_QP_EDIT: &str = "qp_edit";

/// The failed `qp_edit -c` pre-check.
const INCONSISTENT_INPUT: &str = "the input is inconsistent (qp_edit -c failed)";

/// Pause after the `Wall time:` line (`Unix.sleep 1`), so a terminal that
/// launched `qp_run` shows the final output before the prompt.
const FINAL_SLEEP: Duration = Duration::from_secs(1);

/// A fatal error and the status the process exits with.
#[derive(Debug)]
struct Fatal {
    /// Message printed as `qp_run: <message>` on stderr.
    message: String,
    /// Exit status.
    code: i32,
}

impl Fatal {
    /// The reference's `failwith`: uncaught exception, exit 2.
    fn failure(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            code: EXIT_FATAL,
        }
    }
}

impl From<qp_cli::Error> for Fatal {
    fn from(error: qp_cli::Error) -> Self {
        Self::failure(error.to_string())
    }
}

impl From<qp_cli::ffi::Error> for Fatal {
    fn from(error: qp_cli::ffi::Error) -> Self {
        Self::from(qp_cli::Error::Os(error))
    }
}

impl From<std::io::Error> for Fatal {
    fn from(error: std::io::Error) -> Self {
        Self::from(qp_cli::Error::from(error))
    }
}

fn main() {
    let code = match run() {
        Ok(code) => code,
        Err(fatal) => {
            eprintln!("qp_run: {}", fatal.message);
            fatal.code
        }
    };
    std::process::exit(code);
}

/// Parse the command line and run the tool.
///
/// Returns the exit status of the whole process: the executable's status,
/// or `EXIT_FATAL` when a check fails.
fn run() -> Result<i32, Fatal> {
    // `QP_ROOT` is resolved before the command line is parsed: the reference's
    // `Qpackage.root` runs from a module initialiser, so a missing root
    // fails even when the command line is also bad (both are exit 2, the
    // root message wins).
    let root = qp_cli::root()?;
    let programs = Programs::load(&root)?;

    // `clap` owns the command line: `-h`/`--help` prints its help once and
    // exits 0; a usage error (unknown option, missing or extra positional)
    // prints its own plain diagnostic and exits 2.
    let args = RunArgs::parse();
    let config = args.validate();
    drive(&config, &programs)
}

/// The launcher body.
fn drive(config: &RunConfig, programs: &Programs) -> Result<i32, Fatal> {
    let time_start = unix_time();

    // 1. The EZFIO directory must exist.
    let ezfio = config.ezfio_dir.as_str();
    if !Path::new(ezfio).exists() {
        return Err(Fatal::failure(format!(
            "the EZFIO directory does not exist: {ezfio}"
        )));
    }

    // 2. SIGUSR1 → qp_stop, before the job can be interrupted.
    signal::install()?;

    // 3. The program must be one of `data/executables`.
    let Some(program) = programs.find(&config.program) else {
        println!("\nPossible choices:");
        for name in programs.visible() {
            println!("* {name}");
        }
        // Keep stdout ahead of the error on a terminal.
        let _ = std::io::stdout().flush();
        return Err(Fatal::failure(format!(
            "program not found in data/executables: {}",
            config.program
        )));
    };

    // 4. Banner.
    print_banner(ezfio, time_start)?;

    // 5. Input consistency.
    let editor = std::env::var(QP_EDIT_VAR).unwrap_or_else(|_| DEFAULT_QP_EDIT.to_string());
    let status = launch::run_command(&format!("{editor} -c {ezfio}"), &[])?;
    if status != 0 {
        return Err(Fatal::failure(INCONSISTENT_INPUT));
    }

    // 6. Run the executable. The job inherits the environment: the launcher
    //    sets no extra variables (the task-server address is gone).
    let command = launch::job_command(config.prefix.as_deref(), program, ezfio);
    let mut child = launch::spawn_command(&command, &[])?;
    let status = signal::wait_for_job(&mut child, ezfio)?;
    let job_code = launch::exit_code(status);
    if job_code != 0 {
        println!("Program exited with code {job_code}.");
    }

    // 7. Timing and the final status.
    print_wall_time(time_start);
    std::thread::sleep(FINAL_SLEEP);
    Ok(job_code)
}

/// Print the startup block.
///
/// The hour is `tm_hour + (if tm_isdst then 1 else 0)`: the reference adds an
/// hour when daylight saving is in effect, although `localtime` already
/// applied it. Reproduced as-is — the line is informational, and the QP2
/// logs it is compared against contain the extra hour.
///
/// # Errors
///
/// [`Fatal`] when `localtime_r` fails (an epoch outside `time_t`'s range).
fn print_banner(ezfio: &str, start: i64) -> Result<(), Fatal> {
    let tm = qp_cli::ffi::local_time(start)?;
    let stamp = GitStamp::from_env();
    let hour = tm.hour + u32::from(tm.isdst);
    println!(
        "Date: {:02}/{:02}/{:04} {hour:02}:{:02}:{:02}",
        tm.mday, tm.mon, tm.year, tm.min, tm.sec
    );
    println!("===============");
    println!("Quantum Package");
    println!("===============");
    println!();
    println!("Git Commit: {}", stamp.message);
    println!("Git Date  : {}", stamp.date);
    println!("Git SHA1  : {}", stamp.sha1);
    println!("EZFIO Dir : {ezfio}");
    println!();
    println!();
    Ok(())
}

/// Print `Wall time: <h>:<mm>:<ss>` and the blank line.
///
/// The duration is in whole seconds (`Unix.time` resolution), truncated
/// toward zero like `Unix.gmtime`.
fn print_wall_time(start: i64) {
    let elapsed = (unix_time() - start).max(0);
    println!(
        "Wall time: {}:{:02}:{:02}",
        elapsed / 3600,
        (elapsed % 3600) / 60,
        elapsed % 60
    );
    println!();
}

/// Seconds since the Unix epoch, truncated to whole seconds (`Unix.time`).
fn unix_time() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |elapsed| {
            i64::try_from(elapsed.as_secs()).unwrap_or(i64::MAX)
        })
}
