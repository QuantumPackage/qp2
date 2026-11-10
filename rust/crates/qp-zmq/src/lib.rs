//! QP2 ZMQ message protocol and task server (Rust).
//!
//! Plain-text protocol, 18 message kinds; byte-compatible with the f77zmq
//! Fortran workers and remote `qp_run -s` slaves. Ported verbatim from
//! `ocaml/Message.ml`, `ocaml/Message_lexer.mll`, `ocaml/TaskServer.ml` and
//! `ocaml/Progress_bar.ml`.
//!
//! * [`message`] — byte-level codec for the wire protocol.
//! * [`queuing_system`] — the double-ended task queue.
//! * [`progress_bar`] — the `qp_run` progress display (TTY and log forms).
//! * [`port_window`] — pure port-window and `qp_stop` helpers for `qp_run`.
//! * [`task_server`] — the REP/PUB task server, its handlers and `stop`.
//!
//! The `qp_run` binary itself (SIGUSR1 handler installation, window binding,
//! slave mode) lives in `qp-cli`; this crate exposes the tested pieces and
//! the [`task_server::HostLookup`] seam that keeps OS calls out of it.
//!
//! # Ports
//!
//! The task server reserves a contiguous window of ten ports (P..P+9) that
//! starts at `10 * (pid mod 2823) + 32769` and moves forward in steps of 100
//! until every port of the window can be bound. There is no `$QP_PORT`
//! override; the reference (`ocaml/qp_run.ml:39-77`) derives the window from
//! the pid so that concurrent runs on one node land on different windows.
//! The window is reserved by binding the real REP and PUB sockets plus
//! throwaway REP probes on P+2..P+9 before the calculation starts, so a
//! second run that wants the same window fails its own bind and re-scans
//! instead of killing the first (issue #393). [`port_window`] holds the pure
//! candidate sequence; binding is the caller's job.

#![deny(missing_docs)]

pub mod port_window;
pub mod task_server;

pub(crate) mod message;
pub(crate) mod progress_bar;
pub(crate) mod queuing_system;
