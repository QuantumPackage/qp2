//! `qp_run` — Rust port of `ocaml/qp_run.ml` (P4, wave 3).
//!
//! `qp_run` is the launcher every QP2 user types: it reserves a ZMQ port
//! window, starts the task server (the 18-message plain-text protocol, byte
//! compatible with the f77zmq Fortran workers and with `qp_run -s` slaves),
//! runs one Quantum Package executable, and stops the server when the
//! executable is done.
//!
//! # What it computes
//!
//! Nothing physical: it is a launcher. The numbers it produces are the port
//! window (`10 * (pid mod 2823) + 32769`, `+100` per retry, see
//! [`ports`]), the wall time of the run (whole seconds, printed as
//! `h:mm:ss`), and the exit status of the executable it ran.
//!
//! # Flow (order is observable; do not reorder)
//!
//! 1. `QP_ROOT` (environment, else derived from the executable path —
//!    [`qp_cli::root`]; issue #224 forbids baked paths) → `data/executables`.
//!    Resolved before the command line, as the reference's `Qpackage.root`
//!    module initialiser runs.
//! 2. Command line ([`qp_types::command_line::RunArgs`]); `-h`/`--help` is
//!    `clap`'s built-in (help once, exit 0), a wrong number of positionals
//!    is a `clap` usage error (exit 2).
//! 3. Port window (before any check, as in `qp_run.ml:79`).
//! 4. EZFIO directory exists (`the EZFIO directory does not exist: <d>`).
//! 5. SIGUSR1 handler (→ `qp_stop`, see [`qp_cli::signal`]).
//! 6. Program in `data/executables`, else `Possible choices:` + failure.
//! 7. Banner: `Date:` (local time), Quantum Package, Git stamp, EZFIO dir.
//! 8. `qp_edit -c <dir>` unless `-s`; non-zero → `the input is inconsistent
//!    (qp_edit -c failed)`.
//! 9. Task server address resolved **once** in the main thread (issue #393),
//!    then the server thread starts.
//! 10. `<work>/qp_run_address` written (master) or read (slave).
//! 11. The executable runs with `QP_RUN_ADDRESS` (and
//!     `QP_RUN_ADDRESS_MASTER` for a slave) in *its* environment.
//! 12. `TaskServer.stop`, server joined, address file removed, `Wall time:`,
//!     one-second pause, exit with the executable's status.
//!
//! # Environment variables
//!
//! | variable | meaning | reference |
//! |---|---|---|
//! | `QP_ROOT` | QP2 tree (else derived from the executable) | read by `Qpackage.root` |
//! | `QP_NIC` | interface whose IPv4 address the server advertises | same |
//! | `QP_TASK_DEBUG` | print the task-server wire form | same |
//! | `QP_GIT_SHA1`/`QP_GIT_DATE`/`QP_GIT_MESSAGE` | banner stamp, `unknown` when unset | baked by `create_git_sha1.sh` |
//! | `QP_RUN_ADDRESS` | set **for the executable only** | `Unix.putenv` in the parent |
//! | `QP_RUN_ADDRESS_MASTER` | same, slave mode | same |
//! | `QP_RUN_QP_EDIT` | program used for the `-c` pre-check (default `qp_edit`) | hard-coded |
//! | `QP_RUN_MAX_SCANS` | cap on port-window attempts; unset = unlimited (reference) | none |
//!
//! The last two are documented test seams; nothing else in the tool reads
//! them. They exist because the reference's `qp_edit` path and its unbounded
//! port scan cannot both be exercised from a test otherwise.
//!
//! # Deviations from the reference
//!
//! * `Sys.command` is `/bin/sh -c` with inherited standard streams
//!   ([`qp_cli::launch`]), so the `Program exited with code %d.` line and the
//!   exit code match, including the `255` of a signalled shell.
//! * `set_io_threads 16` is not portable (see `qp-zmq/NOTES.md`).
//! * SIGUSR1 sets a flag; `qp_stop` runs from the main thread at most 200 ms
//!   later instead of inside the handler (see `qp_cli::signal`).
//! * The Git stamp comes from the environment instead of being compiled in.
//! * The address file is written after the server thread starts, as in the
//!   reference; a slave that reads it between `write` and the first bind sees
//!   the same race the reference has.
//! * Error text: `<tool>: <message>` on stderr with the reference's exit status
//!   (2 for a `failwith`, 1 for the two explicit `exit 1` sites) instead of
//!   `Fatal error: exception Failure("…")`.
//! * A broken stdout pipe panics Rust's `println!` (exit 101) where the
//!   reference raises `Sys_error` (exit 2). Shared with `qp-zmq`.
//!
//! # Running it
//!
//! ```text
//! cargo run -p qp_run -- scf /path/to/h2.ezfio
//! cargo run -p qp_run -- -s scf /path/to/h2.ezfio      # slave
//! cargo test -p qp_run -j4 -- --test-threads=4         # see tests/
//! ```

#![deny(missing_docs)]

mod ports;

use std::io::Write;
use std::path::Path;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use clap::Parser;
use qp_cli::host::SystemHost;
use qp_cli::{GitStamp, Programs, launch, signal};
use qp_types::command_line::{EXIT_USAGE, RunArgs, RunConfig};
use qp_zmq::task_server;

/// Exit status of a fatal error (an uncaught `Failure`/`Sys_error` in the
/// reference, `CONTRACTS.md` §2.4).
const EXIT_FATAL: i32 = EXIT_USAGE;

/// Exit status of the two explicit `exit 1` sites of `qp_run.ml`: a failed
/// `TaskServer.stop`, and a dead state channel after a successful run.
const EXIT_ONE: i32 = 1;

/// Cap on the port-window scan (`ports` module docs). Test seam.
const MAX_SCANS_VAR: &str = "QP_RUN_MAX_SCANS";

/// Program used for the `qp_edit -c` pre-check. Test seam; the reference
/// hard-codes `qp_edit`.
const QP_EDIT_VAR: &str = "QP_RUN_QP_EDIT";

/// The reference's `qp_edit` path in `Sys.command ("qp_edit -c "^dir)`.
const DEFAULT_QP_EDIT: &str = "qp_edit";

/// Address of the task server, read by the slaves (`qp_run.ml:129`).
const ADDRESS_FILE: &str = "qp_run_address";

/// Address of *this* task server; set for the executable (`qp_run.ml:177`).
const ADDRESS_VAR: &str = "QP_RUN_ADDRESS";

/// Address of the master's task server; set for the executable in slave mode
/// (`qp_run.ml:138`).
const MASTER_VAR: &str = "QP_RUN_ADDRESS_MASTER";

/// The slave's missing/empty master address (`qp_run.ml:140`, reference
/// `failwith "No master is not running"`; wave S2.5 wording).
const NO_MASTER: &str = "no master is running (work/qp_run_address missing or empty)";

/// The master's failed `qp_edit -c` pre-check (`qp_run.ml:123`, reference
/// `failwith "Error: Input inconsistent\n"`; wave S2.5 wording, the
/// trailing-newline quirk goes).
const INCONSISTENT_INPUT: &str = "the input is inconsistent (qp_edit -c failed)";

/// Pause after the `Wall time:` line (`Unix.sleep 1`, `qp_run.ml:216`), so a
/// terminal that launched `qp_run` shows the final output before the prompt.
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

    /// One of the reference's explicit `exit 1` sites.
    fn exit_one(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            code: EXIT_ONE,
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
/// Returns the exit status of the whole process: the executable's status, or
/// `0`/`1` from the task-server outcome.
fn run() -> Result<i32, Fatal> {
    // `QP_ROOT` is resolved before the command line is parsed: the reference's
    // `Qpackage.root` runs from a module initialiser, so a missing root
    // fails even when the command line is also bad (both are exit 2, the
    // root message wins — wave S2.1).
    let root = qp_cli::root()?;
    let programs = Programs::load(&root)?;

    // `clap` owns the command line (wave S2.1): `-h`/`--help` prints its
    // help once and exits 0; a usage error (unknown option, missing or
    // extra positional) prints its own plain diagnostic and exits 2. The
    // reference's double-print of the help on a bad argument count is gone.
    let args = RunArgs::parse();
    let config = args.validate();
    drive(&config, &programs)
}

/// The body of the reference's `run` (`qp_run.ml:19-227`).
fn drive(config: &RunConfig, programs: &Programs) -> Result<i32, Fatal> {
    // 1. Reserve the port window before any other check.
    let context = zmq::Context::new();
    let max_scans = max_scans_from_env()?;
    let window = ports::bind_window(&context, std::process::id(), max_scans)
        .map_err(|error| Fatal::failure(error.to_string()))?;
    let port = window.port;

    let time_start = unix_time();

    // 2. The EZFIO directory must exist.
    let ezfio = config.ezfio_dir.as_str();
    if !Path::new(ezfio).exists() {
        return Err(Fatal::failure(format!(
            "the EZFIO directory does not exist: {ezfio}"
        )));
    }

    // 3. SIGUSR1 → qp_stop, before the job can be interrupted.
    signal::install()?;

    // 4. The program must be one of `data/executables`.
    let Some(program) = programs.find(&config.program) else {
        println!("\nPossible choices:");
        for name in programs.visible() {
            println!("* {name}");
        }
        // Keep stdout ahead of the error on a terminal (`%!` in the reference).
        let _ = std::io::stdout().flush();
        return Err(Fatal::failure(format!(
            "program not found in data/executables: {}",
            config.program
        )));
    };

    // 5. Banner.
    print_banner(ezfio, time_start)?;

    // 6. Input consistency, master only.
    if !config.slave {
        let editor = std::env::var(QP_EDIT_VAR).unwrap_or_else(|_| DEFAULT_QP_EDIT.to_string());
        let status = launch::run_command(&format!("{editor} -c {ezfio}"), &[])?;
        if status != 0 {
            return Err(Fatal::failure(INCONSISTENT_INPUT));
        }
    }

    // 7. `work/`, and the slave's master address.
    let work = qp_cli::ezfio_work(Path::new(ezfio))?;
    let address_file = work.join(ADDRESS_FILE);
    let master = if config.slave {
        Some(read_master_address(&address_file)?)
    } else {
        None
    };

    // 8. Resolve the advertised address exactly once, in this thread, before
    //    the server thread exists (issue #393: a concurrent lazy force in the
    //    OCaml runtime killed the server thread).
    let ip = task_server::resolve_address_from_env(&SystemHost)
        .map_err(|error| Fatal::failure(error.to_string()))?;
    let address = format!("tcp://{ip}:{port}");

    // 9. Start the task server on the pre-bound sockets.
    let server_address = address.clone();
    let server = std::thread::Builder::new()
        .name("task_server".to_string())
        .spawn(move || {
            let result = task_server::run(window.rep, window.pub_socket, &server_address);
            if let Err(error) = &result {
                // The reference prints the exception and its backtrace, then
                // re-raises; `stop` reports the dead server (issue #393).
                // Wave S2.5: plain one-line form.
                eprintln!("qp_run: the task server stopped: {error}");
            }
            result
        })?;

    // 10. Publish the address for the slaves (master only).
    if !config.slave {
        std::fs::write(&address_file, format!("{address}\n"))?;
    }

    // 11. Run the executable.
    let mut env = vec![(ADDRESS_VAR, address.as_str())];
    if let Some(master) = master.as_deref() {
        env.push((MASTER_VAR, master));
    }
    let command = launch::job_command(config.prefix.as_deref(), program, ezfio);
    let mut child = launch::spawn_command(&command, &env)?;
    let status = signal::wait_for_job(&mut child, ezfio)?;
    let job_code = launch::exit_code(status);
    if job_code != 0 {
        println!("Program exited with code {job_code}.");
    }

    // 12. Stop the server and release its ports. A failure here is the
    //     reference's `exit 1`: no clean-up, no `Wall time:` line.
    task_server::stop(port).map_err(|error| Fatal::exit_one(error.to_string()))?;
    let outcome = join_server(server)?;
    // Dropping the context terminates it (`Zmq.Context.terminate`), which
    // releases the window before another `qp_run` can try to take it.
    drop(context);
    if !config.slave {
        std::fs::remove_file(&address_file)?;
    }

    // 13. Timing and the final status.
    print_wall_time(time_start);
    std::thread::sleep(FINAL_SLEEP);
    if job_code != 0 {
        return Ok(job_code);
    }
    match outcome.pub_error {
        None => Ok(0),
        Some(message) => Err(Fatal::exit_one(format!(
            "{message} -- the calculation finished, but the task server state channel was down"
        ))),
    }
}

/// Join the task-server thread and return its outcome.
///
/// A thread that returned an error (a fatal request, a bind failure) is a
/// `failwith` in the reference; a panicked thread is the reference's re-raise
/// in `Thread.join`. Both exit 2.
fn join_server(
    handle: std::thread::JoinHandle<Result<task_server::RunOutcome, task_server::ServerError>>,
) -> Result<task_server::RunOutcome, Fatal> {
    match handle.join() {
        Ok(Ok(outcome)) => Ok(outcome),
        Ok(Err(error)) => Err(Fatal::failure(error.to_string())),
        Err(panic) => {
            let payload = panic
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| {
                    panic
                        .downcast_ref::<&str>()
                        .map(|payload| payload.to_string())
                })
                .unwrap_or_else(|| "non-string payload".to_string());
            Err(Fatal::failure(format!(
                "the task server thread panicked: {payload}"
            )))
        }
    }
}

/// Read the master's address from `<work>/qp_run_address` (`qp_run.ml:131-141`).
///
/// The reference reads one line, closes the file and trims it; any read error
/// selects `No master is not running`. An empty file raises `End_of_file`
/// there (also exit 2, different text); this port reports [`NO_MASTER`] for it
/// as well, because it means the same thing: there is no usable master.
fn read_master_address(path: &Path) -> Result<String, Fatal> {
    let text = std::fs::read_to_string(path).map_err(|_| Fatal::failure(NO_MASTER))?;
    let address = text.lines().next().unwrap_or("").trim();
    if address.is_empty() {
        return Err(Fatal::failure(NO_MASTER));
    }
    Ok(address.to_string())
}

/// `QP_RUN_MAX_SCANS` as a scan cap (`ports` module docs).
///
/// # Errors
///
/// [`Fatal`] when the variable is set to something that is not a `u32`.
fn max_scans_from_env() -> Result<Option<u32>, Fatal> {
    match std::env::var(MAX_SCANS_VAR) {
        Err(_) => Ok(None),
        Ok(text) => text
            .parse::<u32>()
            .map(Some)
            .map_err(|_| Fatal::failure(format!("{MAX_SCANS_VAR} is not a port count: {text:?}"))),
    }
}

/// Print the startup block (`qp_run.ml:106-118`).
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

/// Print `Wall time: <h>:<mm>:<ss>` and the reference's blank line
/// (`qp_run.ml:212-217`).
///
/// `Unix.gmtime (end - start)` is broken into `tm_yday`, `tm_hour`, `tm_min`,
/// `tm_sec`, and the reference prints `tm_yday * 24 + tm_hour`, which is
/// `seconds / 3600` for a non-negative duration. The duration is in whole
/// seconds (`Unix.time` resolution), truncated toward zero like `Unix.gmtime`.
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
