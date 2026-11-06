//! The issue #393 property: the port window is reserved *before* the server
//! starts, so a concurrent `qp_run` cannot steal a port after the pre-check.
//!
//! Both tests need the pid before the binary computes its window, so they
//! spawn `sh -c 'read x || exit 1; exec "$@"' sh <qp_run> …`: the shell blocks
//! on one line of standard input, the test binds a candidate port for that
//! pid, and only then releases the shell. `exec` keeps the pid, so the first
//! candidate is `first_port(pid)`.
//!
//! * The first test pre-binds `P` and checks that the run advertises `P+100`:
//!   the taken window was skipped, not retried later.
//! * The second pre-binds `P`, caps the scan at one attempt
//!   (`QP_RUN_MAX_SCANS=1`, the documented test seam) and checks that the tool
//!   fails *fast* with the bind error, under a hard deadline, instead of
//!   scanning forever.
//!
//! The pre-bound socket is a wildcard bind on `tcp://*:P`, which is exactly
//! what the rival `qp_run` does; a second wildcard bind on the same port fails
//! with `EADDRINUSE` (measured in `qp-zmq/NOTES.md`).

mod common;

use std::time::Duration;

use qp_zmq::port_window::first_port;
use zmq::{Context, Socket, SocketType};

use common::{Fixture, Running, TIMEOUT};

/// Bind `port` on all interfaces and keep the socket (and its context) alive.
fn reserve(context: &Context, port: u32) -> Socket {
    let socket = context.socket(SocketType::REP).expect("socket");
    socket.bind(&format!("tcp://*:{port}")).expect("bind");
    socket
}

/// The run's advertised address line, exactly (`QP_NIC` is the loopback
/// interface in every fixture command).
fn expected_line(port: u32) -> String {
    format!("Task server running : tcp://127.0.0.1:{port}")
}

/// A window taken by a rival makes the scan move by exactly 100.
#[test]
fn a_taken_window_is_skipped_by_100() {
    let fixture = Fixture::new();
    fixture.program("scf", "exit 0");
    let mut running = Running::spawn(
        fixture.gated_command(&["scf", fixture.ezfio_str()]),
        &fixture.logs(),
    );
    let taken = first_port(running.id());
    let context = Context::new();
    let _reservation = reserve(&context, taken);

    // The child cannot have exec'd yet (it is blocked on the gate), so the
    // rival bind above always precedes the binary's first bind attempt.
    running.release_gate();
    let outcome = running.wait(TIMEOUT);

    assert_eq!(outcome.code(), 0, "{outcome}");
    let expected = expected_line(taken + 100);
    assert!(
        outcome.stdout.lines().any(|line| line == expected),
        "expected {expected:?} in:\n{outcome}"
    );
}

/// A scan of one attempt against a taken window fails immediately with
/// `TaskServer.bind_socket`'s wording, and does not start anything.
#[test]
fn an_exhausted_scan_fails_fast() {
    let fixture = Fixture::new();
    fixture.program("scf", "exit 0");
    let mut command = fixture.gated_command(&["scf", fixture.ezfio_str()]);
    // Read before the scan; one attempt cannot reach a free window.
    command.env("QP_RUN_MAX_SCANS", "1");
    let mut running = Running::spawn(command, &fixture.logs());
    let taken = first_port(running.id());
    let context = Context::new();
    let _reservation = reserve(&context, taken);

    running.release_gate();
    // Hard deadline: the reference's unbounded scan would retry forever, so a
    // regression is a 10 s failure with the captured output, not a hang.
    let outcome = running.wait(Duration::from_secs(10));

    assert_eq!(outcome.code(), 2, "{outcome}");
    assert_eq!(
        outcome.stderr,
        format!("qp_run: Unable to bind the REP socket to port : {taken} \n"),
        "the wording (and its trailing space) is TaskServer.bind_socket's"
    );
    // The window is reserved first: no banner, no server, no job.
    assert!(!outcome.stdout.contains("Task server running"), "{outcome}");
    assert!(!outcome.stdout.contains("EZFIO Dir"), "{outcome}");
    assert!(!fixture.work().join("observed").exists());
}
