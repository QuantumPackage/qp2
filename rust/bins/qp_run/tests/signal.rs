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

use std::process::Command;
use std::time::Duration;

use common::{Fixture, Running, TIMEOUT, wait_for_file};

/// Job body: announce that it runs, then wait for the test's `go` file (with
/// a bounded loop, so a broken test cannot leave it forever).
const SLOW_JOB: &str = "\
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
