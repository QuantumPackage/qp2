//! The port window that `qp_run` reserves before starting the task server
//! (`ocaml/qp_run.ml:21-77`).
//!
//! This module is the pure part of the successor's `qp_run` work: the
//! candidate sequence, the window layout and the `qp_stop` command string.
//! It binds nothing; the caller owns the sockets and the retry loop.
//!
//! # Algorithm
//!
//! The first candidate is
//!
//! ```text
//! P0 = 10 * (pid mod 2823) + 32769
//! ```
//!
//! and every rejected window advances by 100 (`P <- P + 100`, forever). A
//! window is ten consecutive ports `P..P+9`:
//!
//! * `P` is bound by the REP socket, which stays bound for the whole run;
//! * `P+1` is bound by the PUB socket, which stays bound as well;
//! * `P+2..P+9` are probed (bound, then released) to reject a window already
//!   used by another `qp_run`; the Fortran collector binds its PULL socket
//!   on `P+2`.
//!
//! The pre-check sockets *are* the real sockets, so a later bind cannot fail:
//! a concurrent `qp_run` that wants the same window fails its own pre-check
//! and re-scans. This is the issue #393 fix; `TaskServer` never scans or
//! binds the window, and `QP_NIC` is resolved before the server thread starts.
//!
//! Ports are `u32` like the reference's OCaml `int`. Candidates above 65535
//! cannot be bound; the caller's bind reports that as an ordinary bind
//! failure (`Unable to bind the ... socket to port : ...`), exactly as the
//! reference does. The largest first candidate, for `pid mod 2823 = 2822`,
//! is 60989, so a fresh process always starts inside the valid range.

/// Number of distinct pid values mapped onto the port space
/// (`2823` in `ocaml/qp_run.ml:72`).
pub const PORT_SLOTS: u32 = 2_823;

/// First port of the window for `pid = 0` (`32_769`).
pub const BASE_PORT: u32 = 32_769;

/// Number of ports in a reserved window, `P..P+9`.
pub const WINDOW_SIZE: u32 = 10;

/// Advance applied to the first port after a rejected window.
pub const RETRY_STRIDE: u32 = 100;

/// Offsets of the probed ports, relative to the first port of the window
/// (`[2; 3; ...; 9]` in `ocaml/qp_run.ml:63`).
pub const PROBE_OFFSETS: [u32; 8] = [2, 3, 4, 5, 6, 7, 8, 9];

/// First candidate port of the window for process id `pid`.
///
/// The reference computes `10 * (Unix.getpid () mod 2823) + 32_769`; the
/// modulo keeps the window clear of the ephemeral range for every pid.
#[must_use]
pub fn first_port(pid: u32) -> u32 {
    WINDOW_SIZE * (pid % PORT_SLOTS) + BASE_PORT
}

/// The ten ports of the window that starts at `port`, in ascending order.
///
/// The caller binds `[0]`, binds `[1]`, and probes `[2..]`; see the module
/// documentation for why.
#[must_use]
pub fn window(port: u32) -> [u32; WINDOW_SIZE as usize] {
    let mut ports = [0; WINDOW_SIZE as usize];
    for (offset, slot) in ports.iter_mut().enumerate() {
        *slot = port + u32::try_from(offset).unwrap_or(0);
    }
    ports
}

/// Iterator over the candidate first ports of the window.
///
/// The reference retries forever (`try_new_port (port + 100)`); this
/// iterator is infinite for the same reason. Construct it with
/// [`candidate_ports`].
#[derive(Clone, Debug)]
pub struct CandidatePorts {
    /// Next first port to yield.
    next: u32,
}

impl Iterator for CandidatePorts {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        let port = self.next;
        self.next = port.wrapping_add(RETRY_STRIDE);
        Some(port)
    }
}

/// Candidate first ports for process id `pid`: [`first_port`] then
/// `+100` at each step, without end.
///
/// The caller binds a window, and moves to the next candidate on failure.
#[must_use]
pub fn candidate_ports(pid: u32) -> CandidatePorts {
    CandidatePorts {
        next: first_port(pid),
    }
}

/// The shell command the SIGUSR1 handler runs to stop a `qp_run`
/// (`Sys.command ("qp_stop "^ezfio_file)`, `ocaml/qp_run.ml:89`).
///
/// `qp_stop` is the reference bash script (not ported): it touches
/// `<ezfio>/work/qpkill`, which the running executable polls. The argument
/// is concatenated verbatim, without quoting, exactly like the reference.
#[must_use]
pub fn qp_stop_command(ezfio: &str) -> String {
    format!("qp_stop {ezfio}")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Oracle: `10*(pid mod 2823)+32769` evaluated for the endpoints and
    /// for one full period.
    #[test]
    fn first_port_formula() {
        assert_eq!(first_port(0), 32_769);
        assert_eq!(first_port(1), 32_779);
        assert_eq!(first_port(2_822), 60_989);
        // One full period returns to the first window.
        assert_eq!(first_port(PORT_SLOTS), first_port(0));
        // The largest window stays inside the 16-bit port range.
        assert!(first_port(2_822) + WINDOW_SIZE - 1 <= u32::from(u16::MAX));
    }

    /// The candidate sequence is arithmetic with stride 100.
    #[test]
    fn candidate_sequence_stride() {
        let mut ports = candidate_ports(1);
        assert_eq!(ports.next(), Some(32_779));
        assert_eq!(ports.next(), Some(32_879));
        assert_eq!(ports.next(), Some(32_979));
    }

    /// The window is `P..P+9` and the probes are `P+2..P+9`.
    #[test]
    fn window_layout() {
        let ports = window(32_769);
        assert_eq!(ports[0], 32_769);
        assert_eq!(ports[9], 32_778);
        for (offset, port) in ports.iter().enumerate() {
            assert_eq!(*port, 32_769 + u32::try_from(offset).unwrap());
        }
        let probes: Vec<u32> = PROBE_OFFSETS.iter().map(|o| 32_769 + o).collect();
        assert_eq!(
            probes,
            vec![
                32_771, 32_772, 32_773, 32_774, 32_775, 32_776, 32_777, 32_778
            ]
        );
    }

    /// `Sys.command ("qp_stop "^ezfio)` concatenates without quoting.
    #[test]
    fn qp_stop_is_a_plain_concatenation() {
        assert_eq!(qp_stop_command("/tmp/a b/ezfio"), "qp_stop /tmp/a b/ezfio");
        assert_eq!(qp_stop_command(""), "qp_stop ");
    }
}
