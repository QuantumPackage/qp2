//! SIGUSR1 → `qp_stop` while the job runs (`ocaml/qp_run.ml:87-89`).
//!
//! The reference runs `Sys.command ("qp_stop "^ezfio)` *inside* its signal
//! handler. This port's handler only stores a flag and the main thread polls
//! it every 200 ms (`qp_cli::signal`), so the observable contract is: after
//! `kill -USR1 <pid of qp_run>`, `qp_stop <ezfio>` runs at least once while
//! the job is alive, and the run still ends normally.
//!
//! The test cannot use `std::env::set_var` (edition 2024 makes it `unsafe`),
//! so `qp_stop` is a shim in the fixture's `bin/`, which the fixture already
//! prepends to the child's `PATH`.

mod common;

use std::os::unix::fs::PermissionsExt;
use std::path::Path;
use std::process::Command;
use std::time::Duration;

use common::{Fixture, Running, TIMEOUT, wait_for_file};

/// Job body: announce that it runs, then wait for the test's `go` file (with
/// a bounded loop, so a broken test cannot leave it forever). The launcher
/// creates no `work/`; the job makes its own, as the Fortran programs do.
const SLOW_JOB: &str = "\
mkdir -p \"$1/work\"
printf 'ready\\n' > \"$1/work/ready\"
i=0
while [ ! -f \"$1/work/go\" ] && [ \"$i\" -lt 600 ]; do
  sleep 0.05
  i=$((i + 1))
done
exit 0";

/// A SIGUSR1 sent while the job runs makes `qp_stop <ezfio>` run, and the run
/// still completes with the job's status.
#[test]
fn sigusr1_runs_qp_stop() {
    let fixture = Fixture::new();
    fixture.program("slow", SLOW_JOB);
    fixture.shim("qp_stop", "printf '%s\\n' \"$1\" > \"$QP_TEST_MARKER\"");
    let marker = fixture.dir().join("stop_marker");

    let mut command = fixture.command(&["slow", fixture.ezfio_str()]);
    command.env("QP_TEST_MARKER", &marker);
    let running = Running::spawn(command, &fixture.logs());
    let pid = running.id();

    if !wait_for_file(&fixture.work().join("ready"), Duration::from_secs(20)) {
        let outcome = running.kill();
        panic!("the job never started\n{outcome}");
    }

    // The handler is installed before the job is spawned (`drive` step 3),
    // so the signal cannot arrive too early.
    let status = Command::new("kill")
        .args(["-USR1", &pid.to_string()])
        .status()
        .expect("run kill(1)");
    assert!(status.success(), "kill -USR1 {pid} failed");

    let saw_stop = wait_for_file(&marker, Duration::from_secs(10));
    // Unblock the job in every case, so the launcher can exit and be reaped.
    std::fs::write(fixture.work().join("go"), "").expect("release the job");
    let outcome = running.wait(TIMEOUT);

    assert!(saw_stop, "qp_stop was not run after SIGUSR1\n{outcome}");
    assert_eq!(
        std::fs::read_to_string(&marker).expect("the marker").trim(),
        fixture.ezfio_str(),
        "qp_stop receives the EZFIO directory verbatim"
    );
    assert_eq!(outcome.code(), 0, "{outcome}");
    assert!(outcome.stdout.contains("Wall time: "), "{outcome}");
}

/// Without a signal, `qp_stop` is never run.
#[test]
fn no_signal_means_no_qp_stop() {
    let fixture = Fixture::new();
    fixture.program("quick", "exit 0");
    fixture.shim("qp_stop", "printf 'called\\n' > \"$QP_TEST_MARKER\"");
    let marker = fixture.dir().join("stop_marker");
    let mut command = fixture.command(&["quick", fixture.ezfio_str()]);
    command.env("QP_TEST_MARKER", &marker);
    let outcome = Running::spawn(command, &fixture.logs()).wait(TIMEOUT);
    assert_eq!(outcome.code(), 0, "{outcome}");
    assert!(!marker.exists(), "qp_stop must not run without a signal");
}

/// Job body: wait for one of its children, which it kills with SIGTERM
/// after a short delay. The script's last command is the `wait`, so the
/// shell exits with `128 + 15`, the POSIX status of a signalled child.
const SIGNALED_JOB: &str = "\
mkdir -p \"$1/work\"
sleep 30 &
child=$!
sleep 0.3
kill -TERM \"$child\"
wait \"$child\"";

/// A job killed by a signal is reported by `/bin/sh` as `128 + n` (unlike a
/// signalled *shell*, which `qp_cli::launch` maps to `255`), and the
/// launcher prints it and exits with it.
#[test]
fn a_job_killed_by_a_signal_exits_128_plus_the_signal() {
    let fixture = Fixture::new();
    fixture.program("signalled", SIGNALED_JOB);
    let outcome = fixture.run(&["signalled", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 128 + 15, "{outcome}");
    assert!(
        outcome.stdout.contains("Program exited with code 143."),
        "{outcome}"
    );
    assert!(outcome.stdout.contains("Wall time: "), "{outcome}");
}

/// SIGUSR1 runs the *real* `bin/qp_stop` (not a shim): `<ezfio>/work/qpkill`
/// must appear. The real script needs GNU `getopt`; on macOS the system
/// `getopt` (BSD) mangles its arguments, so the test probes the script
/// first and skips itself when it cannot run there (the same environment
/// problem the `qp_reset` note in `crates/qp-cli/NOTES.md` records).
#[test]
fn sigusr1_runs_the_real_qp_stop() {
    // `rust/bins/qp_run` → the repository root.
    let repo = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("..");
    let real = std::fs::canonicalize(repo.join("bin").join("qp_stop"))
        .expect("the repository must carry bin/qp_stop");

    // Probe: does the real script work on this platform (GNU getopt)?
    // It only needs a directory with `work/`; the script checks `-d`.
    let probe = tempfile::tempdir().expect("tempdir");
    let probe_ezfio = probe.path().join("probe");
    std::fs::create_dir_all(probe_ezfio.join("work")).expect("probe work/");
    let output = Command::new(&real)
        .arg(probe_ezfio.to_str().expect("UTF-8 probe path"))
        .env("QP_ROOT", probe.path())
        .output()
        .expect("run the real qp_stop");
    if !output.status.success() || !probe_ezfio.join("work").join("qpkill").exists() {
        eprintln!(
            "skipping: the real bin/qp_stop cannot run on this platform \
             (needs GNU getopt); stderr = {}",
            String::from_utf8_lossy(&output.stderr)
        );
        return;
    }

    let fixture = Fixture::new();
    fixture.program("slow", SLOW_JOB);
    // The real script on the child's PATH: the fixture prepends bin/.
    let installed = fixture.root().join("bin").join("qp_stop");
    std::fs::copy(&real, &installed).expect("copy the real qp_stop");
    let mut permissions = std::fs::metadata(&installed)
        .expect("qp_stop metadata")
        .permissions();
    permissions.set_mode(0o755);
    std::fs::set_permissions(&installed, permissions).expect("chmod qp_stop");

    let running = Running::spawn(
        fixture.command(&["slow", fixture.ezfio_str()]),
        &fixture.logs(),
    );
    let pid = running.id();

    // The job creates `work/` before it waits (SLOW_JOB), as the Fortran
    // executables do: the real script only `touch`es, it never `mkdir -p`.
    if !wait_for_file(&fixture.work().join("ready"), Duration::from_secs(20)) {
        let outcome = running.kill();
        panic!("the job never started\n{outcome}");
    }

    let status = Command::new("kill")
        .args(["-USR1", &pid.to_string()])
        .status()
        .expect("run kill(1)");
    assert!(status.success(), "kill -USR1 {pid} failed");

    let saw_qpkill = wait_for_file(&fixture.work().join("qpkill"), Duration::from_secs(10));
    // Unblock the job in every case, so the launcher can exit and be reaped.
    std::fs::write(fixture.work().join("go"), "").expect("release the job");
    let outcome = running.wait(TIMEOUT);

    assert!(
        saw_qpkill,
        "the real qp_stop did not create work/qpkill\n{outcome}"
    );
    assert_eq!(outcome.code(), 0, "{outcome}");
}
