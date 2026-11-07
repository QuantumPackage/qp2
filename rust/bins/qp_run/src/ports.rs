//! Port-window reservation (`ocaml/qp_run.ml:21-77`).
//!
//! # What it does
//!
//! A window is ten consecutive ports `P..P+9`. Candidate first ports come from
//! [`qp_zmq::port_window::candidate_ports`]: `P0 = 10 * (pid mod 2823) + 32769`,
//! then `P + 100` at each retry, without end. For each candidate:
//!
//! 1. a **REP** socket is created and bound to `tcp://*:P` — and kept, because
//!    it becomes the task server's REP socket;
//! 2. a **PUB** socket is created and bound to `tcp://*:P+1` — and kept, for the
//!    same reason;
//! 3. throwaway REP sockets are bound to `tcp://*:P+2 .. tcp://*:P+9` and
//!    closed at once; the Fortran collector binds its PULL socket on `P+2`, so a
//!    window already in use by another `qp_run` is rejected here.
//!
//! Because the pre-check sockets *are* the real sockets, no bind can fail
//! later: a concurrent `qp_run` wanting the same window fails its own first
//! bind and rescan (issue #393). A failure at any step closes whatever this
//! attempt bound and moves to `P + 100`.
//!
//! `tcp://*` (all interfaces) is used, so [`qp_zmq::task_server::stop`] can
//! reach the REP socket through `tcp://localhost:<port>`.
//!
//! # Cost
//!
//! One 0MQ context, two long-lived sockets and eight throwaway sockets per
//! attempt; the scan runs once per process and normally succeeds on the first
//! candidate, so this is 10 binds and 10 small `String` endpoints
//! (`tcp://*:PORT`) plus one bound socket pair for the life of the run. The
//! sockets hold 0MQ I/O buffers in kernel and heap memory (hundreds of KiB),
//! which is why no socket is created per iteration of any other loop.
//!
//! # Deviation: `QP_RUN_MAX_SCANS`
//!
//! The reference retries forever. `QP_RUN_MAX_SCANS=N` (read by the binary)
//! caps the number of candidates at `N`, so a test can exhaust the scan in
//! milliseconds instead of hanging; unset (the default) keeps the reference's
//! unbounded behaviour. When the cap is reached the error names the socket and
//! port of the last attempt with the wording of
//! `TaskServer.bind_socket` — `Unable to bind the REP socket to port : <port> `
//! (trailing space included) — because the missing window is the same
//! condition. This seam exists only for the issue #393 test.

use std::fmt;

use qp_zmq::port_window::{PROBE_OFFSETS, candidate_ports};
use zmq::{Context, Socket, SocketType};

/// Linger of the REP socket, in milliseconds (`qp_run.ml:75`,
/// `Zmq.Socket.set_linger_period`).
pub const REP_LINGER_MS: i32 = 1_000_000;

/// Linger of the PUB socket, in milliseconds (`qp_run.ml:76`).
pub const PUB_LINGER_MS: i32 = 1_000;

/// Which socket of the window failed to bind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SocketKind {
    /// The task server's REP socket on `P`, or one of the throwaway probes on
    /// `P+2..P+9` (the reference creates the probes as REP sockets too).
    Rep,
    /// The task server's PUB socket on `P+1`.
    Pub,
}

impl SocketKind {
    /// The name `TaskServer.bind_socket` prints for this socket type.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Rep => "REP",
            Self::Pub => "PUB",
        }
    }
}

/// Why a window could not be reserved.
#[derive(Debug)]
pub enum WindowError {
    /// 0MQ refused to create a socket, which only happens when the context is
    /// being terminated.
    Create(zmq::Error),
    /// The scan is exhausted (`QP_RUN_MAX_SCANS`); this is the last attempt.
    Exhausted {
        /// Socket that failed to bind.
        socket: SocketKind,
        /// Port it failed to bind.
        port: u32,
    },
}

impl fmt::Display for WindowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Create(error) => write!(f, "Unable to create a socket: {error}"),
            // Same wording as `TaskServer.bind_socket`, trailing space included.
            Self::Exhausted { socket, port } => write!(
                f,
                "Unable to bind the {} socket to port : {port} ",
                socket.as_str()
            ),
        }
    }
}

impl std::error::Error for WindowError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Create(error) => Some(error),
            Self::Exhausted { .. } => None,
        }
    }
}

/// The reserved window: the first port and the two long-lived sockets.
///
/// Dropping the value closes the sockets, exactly like the reference's
/// process-exit abandonment (no message is ever queued on them at this point).
pub struct Window {
    /// First port of the window; the REP socket is bound to it.
    pub port: u32,
    /// REP socket bound to `tcp://*:<port>`; becomes the task server's.
    pub rep: Socket,
    /// PUB socket bound to `tcp://*:<port+1>`; becomes the task server's.
    pub pub_socket: Socket,
}

// `zmq::Socket` is not printable, so the ports are the only useful fields.
impl fmt::Debug for Window {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Window")
            .field("port", &self.port)
            .field("pub_port", &(self.port + 1))
            .finish_non_exhaustive()
    }
}

/// Reserve the first free window for process `pid`.
///
/// `max_scans` is the `QP_RUN_MAX_SCANS` cap: `None` scans until a window is
/// free (the reference), `Some(n)` tries at most `n` candidates and returns
/// [`WindowError::Exhausted`] for the last failure. `Some(0)` attempts
/// nothing and reports the first candidate.
///
/// The caller must keep `context` alive for the sockets' lifetime; every
/// socket created here holds a reference to it.
///
/// # Errors
///
/// [`WindowError::Exhausted`] when the scan cap is reached, or
/// [`WindowError::Create`] when 0MQ cannot create a socket.
pub fn bind_window(
    context: &Context,
    pid: u32,
    max_scans: Option<u32>,
) -> Result<Window, WindowError> {
    let mut last: Option<(SocketKind, u32)> = None;
    // `scan` is the number of candidates already tried, so the cap of 1 lets
    // exactly the first candidate through.
    for (scan, port) in candidate_ports(pid).enumerate() {
        if max_scans.is_some_and(|limit| scan >= limit as usize) {
            break;
        }
        match try_window(context, port) {
            Ok((rep, pub_socket)) => {
                rep.set_linger(REP_LINGER_MS).map_err(WindowError::Create)?;
                pub_socket
                    .set_linger(PUB_LINGER_MS)
                    .map_err(WindowError::Create)?;
                return Ok(Window {
                    port,
                    rep,
                    pub_socket,
                });
            }
            Err(WindowError::Exhausted { socket, port }) => last = Some((socket, port)),
            Err(error) => return Err(error),
        }
    }
    // `max_scans == Some(0)`: nothing was attempted, so report the candidate
    // the reference would have tried first.
    let (socket, port) = last.unwrap_or((SocketKind::Rep, qp_zmq::port_window::first_port(pid)));
    Err(WindowError::Exhausted { socket, port })
}

/// Bind the ten sockets of the window that starts at `port`.
///
/// On failure the sockets bound by this attempt are dropped (closed) and the
/// error names the socket that failed.
fn try_window(context: &Context, port: u32) -> Result<(Socket, Socket), WindowError> {
    let rep = context
        .socket(SocketType::REP)
        .map_err(WindowError::Create)?;
    let pub_socket = context
        .socket(SocketType::PUB)
        .map_err(WindowError::Create)?;
    let failed = |socket, port| WindowError::Exhausted { socket, port };
    if rep.bind(&endpoint(port)).is_err() {
        return Err(failed(SocketKind::Rep, port));
    }
    if pub_socket.bind(&endpoint(port + 1)).is_err() {
        return Err(failed(SocketKind::Pub, port + 1));
    }
    for offset in PROBE_OFFSETS {
        let probe_port = port + offset;
        let probe = context
            .socket(SocketType::REP)
            .map_err(WindowError::Create)?;
        let bound = probe.bind(&endpoint(probe_port));
        // Drop first: a failed probe must not keep the port bound while the
        // caller rescans (the reference closes it before raising).
        drop(probe);
        if bound.is_err() {
            return Err(failed(SocketKind::Rep, probe_port));
        }
    }
    Ok((rep, pub_socket))
}

/// `tcp://*:<port>`, the endpoint of every bind in the reference.
fn endpoint(port: u32) -> String {
    format!("tcp://*:{port}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use qp_zmq::port_window::first_port;

    /// Pids used by the tests. Each test gets its own value, so each scans its
    /// own window and the tests stay independent while `cargo test` runs them
    /// in parallel. Any pid is legal: `first_port` only maps it modulo 2823.
    const PID: u32 = 11;
    const PID_RESCAN: u32 = 22;
    const PID_PROBE: u32 = 33;
    const PID_CAP: u32 = 44;
    const PID_PUB: u32 = 55;
    const PID_ZERO: u32 = 66;

    /// Ports bound by the test, released when it returns.
    struct Reservation {
        /// Kept alive (and bound) for the duration of the test.
        _sockets: Vec<Socket>,
    }

    impl Reservation {
        /// Bind `ports` on all interfaces, the way another `qp_run` would.
        fn new(context: &Context, ports: &[u32]) -> Self {
            Self {
                _sockets: ports
                    .iter()
                    .map(|port| {
                        let socket = context.socket(SocketType::REP).expect("socket");
                        socket.bind(&endpoint(*port)).expect("bind");
                        socket
                    })
                    .collect(),
            }
        }
    }

    /// The window is `P..P+9`, and the two long-lived sockets are bound to
    /// `P` and `P+1`: a *second* wildcard bind on either port must fail, which
    /// is the issue #393 property.
    #[test]
    fn success_binds_the_first_candidate() {
        let context = Context::new();
        let window = bind_window(&context, PID, None).expect("window");
        assert_eq!(window.port, first_port(PID));
        assert_eq!(window.rep.get_linger().expect("linger"), REP_LINGER_MS);
        assert_eq!(
            window.pub_socket.get_linger().expect("linger"),
            PUB_LINGER_MS
        );
        // A concurrent qp_run would fail exactly here.
        let second = context.socket(SocketType::REP).expect("socket");
        assert!(
            second.bind(&endpoint(window.port)).is_err(),
            "a wildcard bind on {} must conflict",
            window.port
        );
    }

    /// A taken first window makes the scan jump by 100.
    #[test]
    fn taken_window_rescans_by_100() {
        let context = Context::new();
        let taken = first_port(PID_RESCAN);
        let _reservation = Reservation::new(&context, &[taken]);
        let window = bind_window(&context, PID_RESCAN, None).expect("window");
        assert_eq!(window.port, taken + 100);
    }

    /// A taken *probe* port (`P+5`) rejects the window too, not only a taken
    /// REP port, and the rejected attempt releases the REP and PUB sockets it
    /// had already bound.
    #[test]
    fn taken_probe_rescans() {
        let context = Context::new();
        let taken = first_port(PID_PROBE);
        let _reservation = Reservation::new(&context, &[taken + 5]);
        let window = bind_window(&context, PID_PROBE, None).expect("window");
        assert_eq!(window.port, taken + 100);
        let check = context.socket(SocketType::REP).expect("socket");
        check.bind(&endpoint(taken)).expect("P was released");
        check.bind(&endpoint(taken + 1)).expect("P+1 was released");
    }

    /// The scan cap makes the failure immediate and names the socket and port
    /// of the last attempt, with the reference wording and its trailing space.
    #[test]
    fn scan_cap_reports_the_last_failure() {
        let context = Context::new();
        let taken = first_port(PID_CAP);
        let _reservation = Reservation::new(&context, &[taken, taken + 100]);
        let error = bind_window(&context, PID_CAP, Some(1)).expect_err("cap reached");
        assert_eq!(
            error.to_string(),
            format!("Unable to bind the REP socket to port : {taken} ")
        );
        // With a cap of 2 the second candidate is the last one tried.
        let error = bind_window(&context, PID_CAP, Some(2)).expect_err("cap reached");
        assert_eq!(
            error.to_string(),
            format!("Unable to bind the REP socket to port : {} ", taken + 100)
        );
    }

    /// The PUB port (`P+1`) taken by someone else reports `PUB`, not `REP`.
    #[test]
    fn taken_pub_port_is_named() {
        let context = Context::new();
        let taken = first_port(PID_PUB);
        let _reservation = Reservation::new(&context, &[taken + 1]);
        let error = bind_window(&context, PID_PUB, Some(1)).expect_err("cap reached");
        assert_eq!(
            error.to_string(),
            format!("Unable to bind the PUB socket to port : {} ", taken + 1)
        );
    }

    /// `Some(0)` attempts nothing and names the first candidate.
    #[test]
    fn zero_scans_is_reported() {
        let context = Context::new();
        let error = bind_window(&context, PID_ZERO, Some(0)).expect_err("no scan");
        assert_eq!(
            error.to_string(),
            format!(
                "Unable to bind the REP socket to port : {} ",
                first_port(PID_ZERO)
            )
        );
    }
}
