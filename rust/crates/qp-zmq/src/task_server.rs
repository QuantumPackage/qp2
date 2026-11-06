//! The `qp_run` task server (`ocaml/TaskServer.ml`).
//!
//! The server is a REP socket that receives one request per client turn and
//! answers on the same socket, a PUB socket that broadcasts the job state
//! (`Waiting`, the job's state name, or `Stopped`), and an in-process PAIR
//! that feeds the PUB thread. Task payloads live in a [`QueuingSystem`]; the
//! handlers here are a direct port of the reference dispatch, including its
//! arm order (`ocaml/TaskServer.ml:920-941`), its exact `error` replies, and
//! its functional-update discipline (a failed request leaves the whole
//! program state untouched).
//!
//! # What the caller owns
//!
//! [`run`] does **not** choose or scan a port window: `qp_run` binds the REP
//! and PUB sockets on `P` and `P+1` (plus probes on `P+2..P+9`) and passes
//! them in, so a concurrent run cannot steal a port after the pre-check
//! (issue #393). See [`crate::port_window`] for the pure candidate sequence.
//! The TCP address is likewise resolved by the caller, once, before the
//! server thread starts: [`resolve_address`] is the pure part of
//! `TaskServer.ip_address`, including the `QP_NIC` read in
//! [`resolve_address_from_env`], but the single-force discipline is the
//! caller's.
//!
//! The `None` socket arguments are the reference's legacy path (bind with a
//! 60 × 1 s retry); `qp_run` always passes `Some`.
//!
//! # Timing
//!
//! Reference defaults are kept: 60 bind retries one second apart, a 1 s poll
//! timeout in the main loop and in the PUB thread, and a 30 s deadline for
//! [`stop`]. Unit tests override them through the private `Timing` value so
//! the suite stays fast; the error strings do not change.
//!
//! # Deviations (documented, not hidden)
//!
//! * `Zmq.Context.set_io_threads 16` is not portable with the declared `zmq`
//!   crate (no `ctx_set`); the default I/O thread count is used.
//! * Socket linger values are the reference's, in **milliseconds**: the OCaml
//!   binding passes the integer straight to `ZMQ_LINGER`. REP 1 000 000 ms,
//!   PUB 1 000 ms, `stop` REQ 1 000 000 ms then 1 000 ms. On the failure path
//!   of [`stop`] this port resets the linger to 0 before the socket is
//!   dropped, because Rust always runs `Drop` where the reference abandons
//!   the socket to process exit; otherwise a dead server would block the
//!   launcher for ~1000 s.
//! * A fatal request (parse or socket error) still sends `Stopped` to the PUB
//!   thread and joins it before returning. The reference lets the server
//!   thread die and leaks the PUB thread; returning an error from a Rust
//!   thread must not leak a thread waiting forever on its PAIR socket.
//! * The PUB thread records any socket error in [`RunOutcome::pub_error`] and
//!   marks the PAIR channel dead. The reference lets the thread die silently
//!   and keeps sending on the PAIR (which can then block forever); only
//!   external `set_stopped` traffic, which no QP client sends, can reach the
//!   difference.
//! * Data values are raw bytes (OCaml strings are byte strings); keys are
//!   protocol tokens and therefore UTF-8.
//! * `QP_TASK_DEBUG` prints the wire form ([`Message::encode_text`]) rather
//!   than the type-`t` debug printer. For every request kind the server can
//!   receive, the two strings are identical; they differ only for reply
//!   kinds, which fall into the `Invalid message` arm.
//! * An unparseable frame is fatal, exactly as `Message.of_string` is called
//!   outside the reference's error handler. The exception is returned as
//!   [`ServerError::Parse`] instead of killing the thread.
//! * `Assert_failure` sites are spelled `TaskServer.ml:<line>:<column>` (the
//!   same basename convention as `queuing_system.rs`). The reference embeds
//!   the path its compiler saw, which differs between build trees.
//! * The progress-bar `stderr` write errors and consecutive clock reads are
//!   documented in [`crate::progress_bar`].

use std::collections::HashMap;
use std::fmt;
use std::io::{self, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, OnceLock};
use std::thread;
use std::time::{Duration, Instant};

use zmq::{Context, Socket, SocketType};

use crate::message::{ConnectType, Id, Message, ParseError};
use crate::progress_bar::ProgressBar;
use crate::queuing_system::{QueueError, QueuingSystem};

/// inproc endpoint that carries the state updates from the main loop to the
/// PUB thread (`TaskServer.ml:754` and `:832`).
const PAIR_ENDPOINT: &str = "inproc://pair";

/// `end_job` state that skips the client check and resets the job
/// (`TaskServer.force_state`).
const FORCE_STATE: &str = "force";

/// Publication state of the task server (`TaskServer.pub_state`).
///
/// The state channel carries the bare text: `Waiting`, `Stopped`, or the
/// job's state name. `Running` is the reference's textual case: any name
/// other than the two keywords is a running job with that name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PubState {
    /// No task is queued; workers should ask for one.
    Waiting,
    /// A job is running in the named state.
    Running(String),
    /// The server is stopping.
    Stopped,
}

impl PubState {
    /// Parse a state frame (`TaskServer.pub_state_of_string`).
    ///
    /// `"Waiting"` and `"Stopped"` are keywords; every other string is a
    /// state name (including the empty string, as in the reference).
    #[must_use]
    pub fn parse(text: &str) -> Self {
        match text {
            "Waiting" => Self::Waiting,
            "Stopped" => Self::Stopped,
            name => Self::Running(name.to_string()),
        }
    }

    /// The text sent on the state channel (`string_of_pub_state`).
    #[must_use]
    pub fn as_str(&self) -> &str {
        match self {
            Self::Waiting => "Waiting",
            Self::Stopped => "Stopped",
            Self::Running(name) => name,
        }
    }
}

/// Operating-system facts needed to resolve the advertised TCP address.
///
/// `qp-cli` implements this trait with `gethostname`, `gethostbyname` and
/// `getifaddrs`; the task server only contains the (tested) decision logic,
/// so this crate stays free of `unsafe`.
pub trait HostLookup {
    /// The machine's host name (`Unix.gethostname`).
    ///
    /// The reference catches every exception and falls back to
    /// `"localhost"`; returning `Err` selects that fallback.
    fn hostname(&self) -> Result<String, String>;

    /// First IPv4 address of `host` (`Unix.gethostbyname` then
    /// `h_addr_list.(0)`), as a dotted quad.
    ///
    /// The reference raises `Not_found` when the name does not resolve and
    /// `Unix_error` on a resolver failure; both are returned as `Err` here.
    fn ipv4_for_host(&self, host: &str) -> Result<String, String>;

    /// First IPv4 address of a named interface
    /// (`get_ipv4_address_for_interface`).
    ///
    /// The C stub reports failure as a string starting with `"error"`; the
    /// reference tests `String.sub result 0 5 = "error"`. Any other string
    /// is used verbatim as the address.
    fn ipv4_for_interface(&self, interface: &str) -> String;
}

/// Why [`resolve_address`] could not produce an address.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ResolveError {
    /// The host name did not resolve to an IPv4 address.
    HostLookup,
    /// `QP_NIC` named an interface that has no IPv4 address.
    Interface(String),
}

impl fmt::Display for ResolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::HostLookup => f.write_str("Unable to find IP address from host name."),
            Self::Interface(interface) => {
                write!(f, "Unable to use network interface {interface}")
            }
        }
    }
}

impl std::error::Error for ResolveError {}

/// Resolve the TCP address advertised by the task server
/// (`TaskServer.ip_address`).
///
/// When `nic` is `None`, the host name is looked up and its first IPv4
/// address is returned; a host-name lookup failure selects the reference's
/// `"localhost"` fallback. When `nic` is `Some(interface)`, that interface's
/// address is returned, and a result starting with `"error"` becomes
/// [`ResolveError::Interface`].
///
/// The reference distinguishes `Unix_error` (turned into the
/// "Unable to find IP address..." failure) from `Not_found` (which escapes
/// uncaught); this port maps every host-lookup error to the failure message.
/// It also accepts an interface result shorter than five bytes, where the
/// reference's `String.sub` would raise `Invalid_argument`.
///
/// # Errors
///
/// [`ResolveError`] with the reference's exact message.
pub fn resolve_address<L: HostLookup + ?Sized>(
    lookup: &L,
    nic: Option<&str>,
) -> Result<String, ResolveError> {
    match nic {
        None => {
            let host = lookup
                .hostname()
                .unwrap_or_else(|_| "localhost".to_string());
            lookup
                .ipv4_for_host(&host)
                .map_err(|_| ResolveError::HostLookup)
        }
        Some(interface) => {
            let result = lookup.ipv4_for_interface(interface);
            if result.starts_with("error") {
                Err(ResolveError::Interface(interface.to_string()))
            } else {
                Ok(result)
            }
        }
    }
}

/// [`resolve_address`] with `QP_NIC` read from the environment.
///
/// This is the port of the reference site (`TaskServer.ml:106`): the
/// variable is read inside the address resolution, and an unset or
/// non-Unicode value is treated as unset (`Not_found`). The caller must still
/// invoke it at most once per process, from the main thread, before any
/// server thread starts (issue #393).
///
/// # Errors
///
/// The same errors as [`resolve_address`].
pub fn resolve_address_from_env<L: HostLookup + ?Sized>(
    lookup: &L,
) -> Result<String, ResolveError> {
    let nic = std::env::var("QP_NIC").ok();
    resolve_address(lookup, nic.as_deref())
}

/// A fatal task-server failure: the reference exception that would have
/// escaped `TaskServer.run`.
#[derive(Debug)]
pub enum ServerError {
    /// A `Failure msg` (`failwith`); the display form is the OCaml
    /// exception rendering `Failure("msg")`.
    Failure(String),
    /// A frame that the reference lexer / message constructor rejected.
    Parse(ParseError),
    /// A 0MQ call failed.
    Zmq(zmq::Error),
    /// A `stdout`/`stderr` write failed (the reference's `Sys_error`).
    Io(io::Error),
}

impl fmt::Display for ServerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Failure(message) => write!(f, "Failure({message:?})"),
            Self::Parse(error) => error.fmt(f),
            Self::Zmq(error) => write!(f, "ZMQ_exception({error:?}, {:?})", error.to_string()),
            Self::Io(error) => write!(f, "Sys_error({:?})", error.to_string()),
        }
    }
}

impl std::error::Error for ServerError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Parse(error) => Some(error),
            Self::Zmq(error) => Some(error),
            Self::Io(error) => Some(error),
            Self::Failure(_) => None,
        }
    }
}

impl From<ParseError> for ServerError {
    fn from(error: ParseError) -> Self {
        Self::Parse(error)
    }
}

impl From<zmq::Error> for ServerError {
    fn from(error: zmq::Error) -> Self {
        Self::Zmq(error)
    }
}

impl From<io::Error> for ServerError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

/// What [`run`] returns after a clean shutdown.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RunOutcome {
    /// The reference's `TaskServer.pub_error ()`: `Some(message)` when the
    /// state channel could not be bound (or died), `None` when it was
    /// healthy. `qp_run` exits with status 1 after a successful calculation
    /// when this is `Some`.
    pub pub_error: Option<String>,
}

/// Reference timing constants, overridable by unit tests.
#[derive(Clone, Copy, Debug)]
struct Timing {
    /// `bind_socket` attempts before giving up (`TaskServer.ml:84`).
    bind_retries: u32,
    /// Delay between bind attempts (`Unix.sleep 1`).
    bind_delay: Duration,
    /// Main-loop `Zmq.Poll.poll ~timeout` in milliseconds.
    loop_poll_ms: i64,
    /// PUB-thread `Zmq.Poll.poll ~timeout` in milliseconds.
    pub_poll_ms: i64,
    /// [`stop`] deadline (`TaskServer.ml:163`).
    stop_deadline: Duration,
    /// [`stop`] poll timeout in milliseconds (`TaskServer.ml:166`).
    stop_poll_ms: i64,
}

impl Default for Timing {
    fn default() -> Self {
        Self {
            bind_retries: 60,
            bind_delay: Duration::from_secs(1),
            loop_poll_ms: 1_000,
            pub_poll_ms: 1_000,
            stop_deadline: Duration::from_secs(30),
            stop_poll_ms: 1_000,
        }
    }
}

/// Bind `socket` to `tcp://*:<port>`, retrying `retries` times with
/// `delay` between attempts (`TaskServer.bind_socket`).
///
/// On exhaustion the message is the reference's, including its trailing
/// space: `Unable to bind the <type> socket to port : <port> `.
///
/// The reference retries every `Unix.Unix_error` and re-raises other
/// exceptions; `zmq::Error` does not expose that distinction, so every bind
/// error is retried here.
fn bind_socket(
    socket: &Socket,
    socket_type: &str,
    port: u32,
    retries: u32,
    delay: Duration,
) -> Result<(), String> {
    let endpoint = format!("tcp://*:{port}");
    let mut remaining = retries;
    loop {
        if remaining == 0 {
            return Err(format!(
                "Unable to bind the {socket_type} socket to port : {port} "
            ));
        }
        match socket.bind(&endpoint) {
            Ok(()) => return Ok(()),
            Err(_) => {
                thread::sleep(delay);
                remaining -= 1;
            }
        }
    }
}

/// Reset a socket's linger to 0 on drop unless disarmed.
///
/// [`stop_with`] needs it: the reference abandons its REQ socket to process
/// exit on the failure path, where a 1 000 000 ms linger would block Rust's
/// `Drop` for ~1000 s.
struct LingerReset<'a> {
    /// Socket whose linger is reset.
    socket: &'a Socket,
    /// Whether the reset still applies.
    armed: bool,
}

impl Drop for LingerReset<'_> {
    fn drop(&mut self) {
        if self.armed {
            let _ = self.socket.set_linger(0);
        }
    }
}

/// Ask a running server to stop (`TaskServer.stop`).
///
/// Connects a REQ socket to `tcp://localhost:<port>`, sends `terminate`, and
/// waits up to 30 s for the answer, polling every second. Any answer other
/// than `ok` is the reference's `Problem in termination` failure.
///
/// # Errors
///
/// [`ServerError::Failure`] with
/// `TaskServer: no answer to the Terminate request (the task server thread is dead)`
/// when the deadline expires, `Problem in termination` for a non-`ok`
/// answer, or a transport/parse error.
pub fn stop(port: u32) -> Result<(), ServerError> {
    let timing = Timing::default();
    stop_with(port, timing.stop_deadline, timing.stop_poll_ms)
}

/// [`stop`] with an explicit deadline and poll timeout, for tests.
fn stop_with(port: u32, deadline: Duration, poll_ms: i64) -> Result<(), ServerError> {
    // `TaskServer.stop` starts with `debug "STOP"`.
    if debug_enabled() {
        write_debug("STOP")?;
    }
    let context = Context::new();
    let req_socket = context.socket(SocketType::REQ)?;
    // The OCaml binding passes this integer straight to ZMQ_LINGER, i.e.
    // milliseconds (TaskServer.ml:151).
    req_socket.set_linger(1_000_000)?;
    let mut linger_reset = LingerReset {
        socket: &req_socket,
        armed: true,
    };
    req_socket.connect(&format!("tcp://localhost:{port}"))?;
    req_socket.send(Message::Terminate.encode(), 0)?;

    let deadline_at = Instant::now() + deadline;
    let reply = loop {
        if req_socket.poll(zmq::POLLIN, poll_ms)? > 0 {
            break req_socket.recv_bytes(0)?;
        }
        if Instant::now() > deadline_at {
            return Err(ServerError::Failure(
                "TaskServer: no answer to the Terminate request (the task server thread is dead)"
                    .to_string(),
            ));
        }
    };

    match Message::parse(&reply)? {
        Message::Ok => {}
        _ => {
            return Err(ServerError::Failure("Problem in termination".to_string()));
        }
    }

    // Reference: `set_linger_period req_socket 1_000` before close. The
    // reply has been received, so nothing is pending; the guard is disarmed
    // because the reference keeps this value.
    linger_reset.armed = false;
    req_socket.set_linger(1_000)?;
    Ok(())
}

/// The task server's mutable state (`TaskServer.t`).
#[derive(Debug)]
struct ProgramState {
    /// Queued, running and finished-but-not-deleted tasks.
    queue: QueuingSystem,
    /// Current job state name, `None` between jobs.
    state: Option<String>,
    /// PUSH endpoint advertised to TCP workers.
    address_tcp: Option<String>,
    /// PUSH endpoint advertised to inproc workers.
    address_inproc: Option<String>,
    /// Progress bar of the current job.
    progress_bar: Option<ProgressBar>,
    /// Whether the main loop keeps running (`terminate` clears it).
    running: bool,
    /// Whether `connect` is accepted.
    accepting_clients: bool,
    /// `put_data`/`get_data` hash table; values are raw bytes.
    data: HashMap<String, Vec<u8>>,
}

impl ProgramState {
    /// The initial state of `TaskServer.run` (`TaskServer.ml:851-861`).
    fn new() -> Self {
        Self {
            queue: QueuingSystem::new(),
            state: None,
            address_tcp: None,
            address_inproc: None,
            progress_bar: None,
            running: true,
            accepting_clients: false,
            data: HashMap::new(),
        }
    }
}

/// Failure of one request, classified like the reference's exception
/// handlers.
enum ReplyError {
    /// `Failure msg`; the server replies `<msg> : <raw request>`.
    Failure(String),
    /// `Assert_failure`; the server replies `<site> : <raw request>`.
    Assert(&'static str),
    /// The server cannot continue (socket or I/O error).
    Fatal(ServerError),
}

impl From<QueueError> for ReplyError {
    fn from(error: QueueError) -> Self {
        match error {
            QueueError::Failure(message) => Self::Failure(message),
            QueueError::AssertFailure(site) => Self::Assert(site),
        }
    }
}

impl From<zmq::Error> for ReplyError {
    fn from(error: zmq::Error) -> Self {
        Self::Fatal(ServerError::Zmq(error))
    }
}

impl From<io::Error> for ReplyError {
    fn from(error: io::Error) -> Self {
        Self::Fatal(ServerError::Io(error))
    }
}

/// The sending half of the inproc PAIR that feeds the PUB thread.
struct PubChannel {
    /// PAIR socket bound by [`run_with`]; only the main loop sends on it.
    socket: Socket,
    /// Reference `pub_pair_alive`: cleared when the PUB thread is gone, so
    /// a send cannot queue forever on a missing peer.
    alive: Arc<AtomicBool>,
}

impl PubChannel {
    /// Publish `state` if the PUB thread is still alive
    /// (`TaskServer.pub_send`).
    ///
    /// # Errors
    ///
    /// The 0MQ send error (fatal, as in the reference).
    fn send(&self, state: &PubState) -> Result<(), zmq::Error> {
        if self.alive.load(Ordering::Relaxed) {
            self.socket.send(state.as_str(), 0)?;
        }
        Ok(())
    }
}

/// Send one text frame (`Zmq.Socket.send`).
fn send_message(socket: &Socket, message: &Message) -> Result<(), ReplyError> {
    socket.send(message.encode(), 0)?;
    Ok(())
}

/// Send a multi-frame reply (`Zmq.Socket.send_all`).
fn send_frames(socket: &Socket, message: &Message) -> Result<(), ReplyError> {
    socket.send_multipart(message.encode_frames(), 0)?;
    Ok(())
}

/// Reply `ok`.
fn reply_ok(rep_socket: &Socket) -> Result<(), ReplyError> {
    send_message(rep_socket, &Message::Ok)
}

/// Reply `error Wrong state`.
fn reply_wrong_state(rep_socket: &Socket) -> Result<(), ReplyError> {
    reply_error(rep_socket, "Wrong state")
}

/// Reply `error <text>` (`TaskServer.error`).
fn reply_error(rep_socket: &Socket, text: &str) -> Result<(), ReplyError> {
    send_message(
        rep_socket,
        &Message::Error {
            text: text.to_string(),
        },
    )
}

/// Reply `<text> : <raw request>` for a recoverable request failure
/// (`TaskServer.ml:944-955`).
///
/// `text` is the `Failure` payload or the `Assert_failure` site, and `raw` is
/// the first frame of the offending request, echoed verbatim.
fn reply_failure(rep_socket: &Socket, text: &str, raw: &[u8]) -> Result<(), ServerError> {
    rep_socket.send(
        Message::Error {
            text: format!("{text} : {}", String::from_utf8_lossy(raw)),
        }
        .encode(),
        0,
    )?;
    Ok(())
}

/// Whether `QP_TASK_DEBUG=1` was set when the process started
/// (`TaskServer.debug_env`, forced once at module initialisation).
///
/// The reference reads the variable once, at start-up; the `OnceLock`
/// reproduces that, so a later change to the environment cannot alter the
/// already-running server.
fn debug_enabled() -> bool {
    static DEBUG_ENV: OnceLock<bool> = OnceLock::new();
    *DEBUG_ENV.get_or_init(|| std::env::var("QP_TASK_DEBUG").is_ok_and(|value| value == "1"))
}

/// Write `TASK : <text>` on stderr without an extra newline
/// (`TaskServer.debug`).
fn write_debug(text: &str) -> io::Result<()> {
    let mut stderr = io::stderr();
    stderr.write_all(b"TASK : ")?;
    stderr.write_all(text.as_bytes())?;
    stderr.flush()
}

/// The per-request debug line (`TaskServer.ml:909-914`).
fn debug_counts(queue: &QueuingSystem, message: &str) -> String {
    format!(
        "q:{}  r:{}  n:{}  c:{} : {message}\n",
        queue.number_of_queued(),
        queue.number_of_running(),
        queue.number_of_tasks(),
        queue.number_of_clients(),
    )
}

/// Start a job (`TaskServer.new_job`).
fn new_job(
    program: &mut ProgramState,
    state: &str,
    address_tcp: &str,
    address_inproc: &str,
    rep_socket: &Socket,
    pub_channel: &PubChannel,
) -> Result<(), ReplyError> {
    program.state = Some(state.to_string());
    program.progress_bar = Some(ProgressBar::new(state));
    program.address_tcp = Some(address_tcp.to_string());
    program.address_inproc = Some(address_inproc.to_string());
    program.accepting_clients = true;
    reply_ok(rep_socket)?;
    pub_channel.send(&PubState::Waiting)?;
    Ok(())
}

/// Which state [`change_pub_state`] publishes.
#[derive(Clone, Copy, Debug)]
enum PubChange {
    /// `set_running`: publish the current job's state name.
    Running,
    /// `set_waiting`: publish `Waiting`.
    Waiting,
    /// `set_stopped`: publish `Stopped`.
    Stopped,
}

/// Publish a state without touching the job (`TaskServer.change_pub_state`).
fn change_pub_state(
    program: &ProgramState,
    change: PubChange,
    rep_socket: &Socket,
    pub_channel: &PubChannel,
) -> Result<(), ReplyError> {
    let state = match change {
        PubChange::Waiting => PubState::Waiting,
        PubChange::Stopped => PubState::Stopped,
        PubChange::Running => match program.state.as_deref() {
            Some(state) => PubState::Running(state.to_string()),
            None => {
                return Err(ReplyError::Failure(
                    "Trying to change pub state while no job is ready".to_string(),
                ));
            }
        },
    };
    reply_ok(rep_socket)?;
    pub_channel.send(&state)?;
    Ok(())
}

/// Reset the job after a successful `end_job` (`TaskServer.end_job`,
/// `success`).
///
/// Replies `ok` first, then clears the progress display and every field of
/// the program state, as the reference's record construction does.
fn reset_job(program: &mut ProgramState, rep_socket: &Socket) -> Result<(), ReplyError> {
    reply_ok(rep_socket)?;
    ProgressBar::clear()?;
    program.queue = QueuingSystem::new();
    program.state = None;
    program.progress_bar = None;
    program.address_tcp = None;
    program.address_inproc = None;
    program.running = true;
    program.accepting_clients = false;
    program.data = HashMap::new();
    Ok(())
}

/// Finish the job (`TaskServer.end_job`).
fn end_job(
    program: &mut ProgramState,
    msg_state: &str,
    rep_socket: &Socket,
    pub_channel: &PubChannel,
) -> Result<(), ReplyError> {
    let Some(state) = program.state.as_deref() else {
        return reply_wrong_state(rep_socket);
    };
    if msg_state == FORCE_STATE {
        pub_channel.send(&PubState::Waiting)?;
        return reset_job(program, rep_socket);
    }
    if msg_state != state {
        return reply_wrong_state(rep_socket);
    }
    pub_channel.send(&PubState::Waiting)?;
    let clients = program.queue.number_of_clients();
    if clients == 0 {
        reset_job(program, rep_socket)
    } else {
        reply_error(rep_socket, &format!("waiting for {clients} slaves..."))
    }
}

/// Register a worker and hand it the PUSH endpoint
/// (`TaskServer.connect`).
fn connect(
    program: &mut ProgramState,
    socket_type: ConnectType,
    rep_socket: &Socket,
) -> Result<(), ReplyError> {
    if !program.accepting_clients {
        return reply_wrong_state(rep_socket);
    }
    let Some(state) = program.state.clone() else {
        return reply_wrong_state(rep_socket);
    };
    let push_address = match socket_type {
        ConnectType::Tcp => program
            .address_tcp
            .clone()
            .ok_or_else(|| ReplyError::Failure("Error: No TCP address".to_string()))?,
        ConnectType::Inproc => program
            .address_inproc
            .clone()
            .ok_or_else(|| ReplyError::Failure("Error: No inproc address".to_string()))?,
        // `TaskServer.ml:312:46` (`| Message.Connect_msg.Ipc -> assert false`).
        ConnectType::Ipc => return Err(ReplyError::Assert("TaskServer.ml:312:46")),
    };
    let client_id = program.queue.add_client();
    send_message(
        rep_socket,
        &Message::ConnectReply {
            state,
            client_id,
            push_address,
        },
    )
}

/// Remove a worker (`TaskServer.disconnect`).
fn disconnect(
    program: &mut ProgramState,
    msg_state: &str,
    client_id: Id,
    rep_socket: &Socket,
) -> Result<(), ReplyError> {
    if program.state.as_deref() != Some(msg_state) {
        return reply_wrong_state(rep_socket);
    }
    program.queue.del_client(client_id)?;
    send_message(
        rep_socket,
        &Message::DisconnectReply {
            state: msg_state.to_string(),
        },
    )
}

/// Append tasks and extend the progress bar (`TaskServer.add_task`).
///
/// The reference ignores the request state entirely: `add_task` works
/// whatever state the client believes it is in, as long as a job exists
/// (the dispatch's `Some _` guard).
fn add_task(
    program: &mut ProgramState,
    tasks: &[String],
    rep_socket: &Socket,
) -> Result<(), ReplyError> {
    for task in tasks {
        program.queue.add_task(task.clone());
        if let Some(bar) = program.progress_bar.as_mut() {
            bar.increment_end();
        }
    }
    reply_ok(rep_socket)
}

/// Delete task payloads (`TaskServer.del_task`).
fn del_task(
    program: &mut ProgramState,
    msg_state: &str,
    task_ids: &[Id],
    rep_socket: &Socket,
) -> Result<(), ReplyError> {
    if program.state.as_deref() != Some(msg_state) {
        return reply_wrong_state(rep_socket);
    }
    program.queue.del_tasks(task_ids)?;
    program.accepting_clients =
        program.queue.number_of_queued() > program.queue.number_of_clients();
    let more = program.queue.number_of_tasks() > 0;
    send_message(
        rep_socket,
        &Message::DelTaskReply {
            task_ids: task_ids.to_vec(),
            more,
        },
    )
}

/// Take one task (`TaskServer.get_task`).
fn get_task(
    program: &mut ProgramState,
    msg_state: &str,
    client_id: Id,
    rep_socket: &Socket,
    pub_channel: &PubChannel,
) -> Result<(), ReplyError> {
    if program.state.as_deref() != Some(msg_state) {
        return reply_wrong_state(rep_socket);
    }
    let popped = program.queue.pop_task(client_id)?;
    program.accepting_clients =
        program.queue.number_of_queued() > program.queue.number_of_clients();
    if program.queue.number_of_queued() == 0 {
        pub_channel.send(&PubState::Waiting)?;
    } else {
        pub_channel.send(&PubState::Running(msg_state.to_string()))?;
    }
    send_message(
        rep_socket,
        &Message::GetTaskReply {
            task_id: popped.task_id,
            task: popped.task,
        },
    )
}

/// Take up to `n_tasks` tasks (`TaskServer.get_tasks`).
///
/// The reference pops with `Queuing_system.pop_task` and, when the popped id
/// has no payload (deleted while queued), discards that pop and answers
/// `0 terminate` with the *pre-pop* queue, leaving the id queued. This port
/// detects that case with [`QueuingSystem::task_payload`] before popping, so
/// the queue is never mutated and no clone is needed. The client-membership
/// assertion that `pop_task` performs is done once up front
/// ([`QueuingSystem::has_client`]), because the client set cannot change
/// during the loop.
fn get_tasks(
    program: &mut ProgramState,
    msg_state: &str,
    client_id: Id,
    n_tasks: i64,
    rep_socket: &Socket,
    pub_channel: &PubChannel,
) -> Result<(), ReplyError> {
    if program.state.as_deref() != Some(msg_state) {
        return reply_wrong_state(rep_socket);
    }
    if !program.queue.has_client(client_id) {
        return Err(QueueError::AssertFailure("Queuing_system.ml:68:2").into());
    }

    let mut result = Vec::new();
    let mut remaining = n_tasks;
    while remaining > 0 {
        let Some(task_id) = program.queue.next_queued_id() else {
            result.push((None, "terminate".to_string()));
            break;
        };
        if program.queue.task_payload(task_id).is_none() {
            // Reference: the id stays queued and the reply ends with
            // `0 terminate` (the pop result is discarded).
            result.push((None, "terminate".to_string()));
            break;
        }
        let popped = program.queue.pop_task(client_id)?;
        let (Some(task_id), Some(task)) = (popped.task_id, popped.task) else {
            unreachable!("get_tasks: payload checked immediately before pop_task");
        };
        result.push((Some(task_id), task));
        remaining -= 1;
    }

    program.accepting_clients =
        program.queue.number_of_queued() > program.queue.number_of_clients();
    if program.queue.number_of_queued() == 0 {
        pub_channel.send(&PubState::Waiting)?;
    } else {
        pub_channel.send(&PubState::Running(msg_state.to_string()))?;
    }
    send_frames(rep_socket, &Message::GetTasksReply { tasks: result })
}

/// Mark tasks finished (`TaskServer.task_done`).
fn task_done(
    program: &mut ProgramState,
    msg_state: &str,
    client_id: Id,
    task_ids: &[Id],
    rep_socket: &Socket,
) -> Result<(), ReplyError> {
    if program.state.as_deref() != Some(msg_state) {
        return reply_wrong_state(rep_socket);
    }
    program.queue.end_tasks(task_ids, client_id)?;
    if let Some(bar) = program.progress_bar.as_mut() {
        for _ in task_ids {
            bar.increment_cur();
        }
    }
    program.accepting_clients =
        program.queue.number_of_queued() > program.queue.number_of_clients();
    reply_ok(rep_socket)
}

/// Store a raw value (`TaskServer.put_data`).
///
/// The request must carry exactly one extra frame; the raw request text ends
/// at the first frame, so the reference's `Failure` reply echoes only that
/// first frame.
fn put_data(
    program: &mut ProgramState,
    msg_state: &str,
    key: &str,
    rest_of_msg: &[Vec<u8>],
    rep_socket: &Socket,
) -> Result<(), ReplyError> {
    let [value] = rest_of_msg else {
        return Err(ReplyError::Failure(
            "Badly formed put_data message".to_string(),
        ));
    };
    match program.state.as_deref() {
        // `TaskServer.ml:646:14`.
        None => Err(ReplyError::Assert("TaskServer.ml:646:14")),
        Some(state) if state != msg_state => reply_wrong_state(rep_socket),
        Some(_) => {
            program.data.insert(key.to_string(), value.clone());
            send_message(rep_socket, &Message::PutDataReply)
        }
    }
}

/// Read a raw value (`TaskServer.get_data`).
///
/// A missing key answers a single NUL byte, as the reference does.
fn get_data(
    program: &ProgramState,
    msg_state: &str,
    key: &str,
    rep_socket: &Socket,
) -> Result<(), ReplyError> {
    match program.state.as_deref() {
        // `TaskServer.ml:678:14`.
        None => Err(ReplyError::Assert("TaskServer.ml:678:14")),
        Some(state) if state != msg_state => reply_wrong_state(rep_socket),
        Some(_) => {
            let value = program.data.get(key).cloned().unwrap_or_else(|| vec![0]);
            send_frames(rep_socket, &Message::GetDataReply { value })
        }
    }
}

/// Discard every queued task (`TaskServer.abort`).
///
/// Unlike [`QueuingSystem::abort`] (which discards one task), the reference
/// adds a throwaway client, pops the queue until it is empty, marks every
/// popped task finished, deletes it, and removes the client. The whole
/// operation is done on a clone and committed only on success, matching the
/// reference's functional update.
fn abort(program: &mut ProgramState, rep_socket: &Socket) -> Result<(), ReplyError> {
    let mut queue = program.queue.clone();
    let client_id = queue.add_client();

    let mut tasks = Vec::new();
    // The reference starts its loop with `1`, so it pops once even from an
    // empty queue, then continues while tasks remain queued.
    let popped = queue.pop_task(client_id)?;
    if let Some(task_id) = popped.task_id {
        tasks.push(task_id);
    }
    while queue.number_of_queued() > 0 {
        let popped = queue.pop_task(client_id)?;
        if let Some(task_id) = popped.task_id {
            tasks.push(task_id);
        }
    }
    // The reference accumulates with `task_id :: accu`, so both folds see
    // the popped ids in reverse order; only the first failure can observe
    // the order.
    tasks.reverse();

    for task_id in &tasks {
        queue.end_task(*task_id, client_id)?;
    }
    for task_id in &tasks {
        queue.del_task(*task_id)?;
    }
    queue.del_client(client_id)?;

    reply_ok(rep_socket)?;
    program.queue = queue;
    program.accepting_clients = false;
    Ok(())
}

/// Stop the main loop (`TaskServer.terminate`).
fn terminate(program: &mut ProgramState, rep_socket: &Socket) -> Result<(), ReplyError> {
    reply_ok(rep_socket)?;
    program.address_tcp = None;
    program.address_inproc = None;
    program.running = false;
    Ok(())
}

/// The dispatch table (`TaskServer.ml:920-941`), in the reference's arm
/// order.
///
/// `rest` holds the frames after the first one of the same multipart
/// message.
fn dispatch(
    program: &mut ProgramState,
    message: &Message,
    rest: &[Vec<u8>],
    rep_socket: &Socket,
    pub_channel: &PubChannel,
    debug: bool,
) -> Result<(), ReplyError> {
    match (program.state.is_some(), message) {
        (_, Message::Terminate) => terminate(program, rep_socket),
        (_, Message::Abort) => abort(program, rep_socket),
        (_, Message::PutData { state, key, .. }) => {
            // `TaskServer.put_data` starts with `debug (PutData_msg.to_string msg)`.
            if debug {
                write_debug(&message.encode_text())?;
            }
            put_data(program, state, key, rest, rep_socket)
        }
        (_, Message::GetData { state, key, .. }) => {
            // `TaskServer.get_data` starts with `debug (GetData_msg.to_string msg)`.
            if debug {
                write_debug(&message.encode_text())?;
            }
            get_data(program, state, key, rep_socket)
        }
        (
            false,
            Message::Newjob {
                state,
                address_tcp,
                address_inproc,
            },
        ) => new_job(
            program,
            state,
            address_tcp,
            address_inproc,
            rep_socket,
            pub_channel,
        ),
        (_, Message::Newjob { .. }) => reply_error(rep_socket, "A job is already running"),
        (true, Message::Endjob { state }) => end_job(program, state, rep_socket, pub_channel),
        (true, Message::SetRunning) => {
            change_pub_state(program, PubChange::Running, rep_socket, pub_channel)
        }
        (_, Message::SetWaiting) => {
            change_pub_state(program, PubChange::Waiting, rep_socket, pub_channel)
        }
        (_, Message::SetStopped) => {
            change_pub_state(program, PubChange::Stopped, rep_socket, pub_channel)
        }
        (false, _) => reply_error(rep_socket, "No job is running"),
        (true, Message::Connect { socket }) => connect(program, *socket, rep_socket),
        (true, Message::Disconnect { state, client_id }) => {
            disconnect(program, state, *client_id, rep_socket)
        }
        (true, Message::AddTask { tasks, .. }) => add_task(program, tasks, rep_socket),
        (true, Message::DelTask { state, task_ids }) => {
            del_task(program, state, task_ids, rep_socket)
        }
        (true, Message::GetTask { state, client_id }) => {
            get_task(program, state, *client_id, rep_socket, pub_channel)
        }
        (
            true,
            Message::GetTasks {
                state,
                client_id,
                n_tasks,
            },
        ) => get_tasks(
            program,
            state,
            *client_id,
            *n_tasks,
            rep_socket,
            pub_channel,
        ),
        (
            true,
            Message::TaskDone {
                state,
                client_id,
                task_ids,
            },
        ) => task_done(program, state, *client_id, task_ids, rep_socket),
        (_, _) => reply_error(
            rep_socket,
            &format!("Invalid message : {}", message.encode_text()),
        ),
    }
}

/// Body of the PUB thread (`TaskServer.start_pub_thread`).
///
/// Connects its PAIR socket, binds the PUB socket on `port` when the caller
/// did not pre-bind it, then republishes the last state received:
///
/// 1. publish the state (initial `Waiting`);
/// 2. poll for the next state for `pub_poll_ms`;
/// 3. publish again (this is why `Stopped` is sent twice, at most one poll
///    interval apart);
/// 4. stop once the *previous* state was `Stopped`.
///
/// Any infrastructure error is returned as a `String`, recorded by the
/// caller in the liveness flag and in the outcome.
fn pub_thread_main(
    context: &Context,
    pub_socket: Option<Socket>,
    port: u32,
    timing: Timing,
    alive: &AtomicBool,
    pub_error: &OnceLock<String>,
) -> Result<(), String> {
    let pair_socket = context
        .socket(SocketType::PAIR)
        .map_err(|error| error.to_string())?;
    pair_socket
        .connect(PAIR_ENDPOINT)
        .map_err(|error| error.to_string())?;

    let pre_bound = pub_socket.is_some();
    let pub_socket = match pub_socket {
        Some(socket) => socket,
        None => context
            .socket(SocketType::PUB)
            .map_err(|error| error.to_string())?,
    };

    if !pre_bound
        && let Err(message) = bind_socket(
            &pub_socket,
            "PUB",
            port,
            timing.bind_retries,
            timing.bind_delay,
        )
    {
        // The reference marks the channel dead, records the message and
        // reports it on stderr; the calculation keeps running.
        alive.store(false, Ordering::Relaxed);
        let _ = pub_error.set(message.clone());
        let _ = writeln!(
            io::stderr(),
            "TaskServer: {message} -- the state channel is unavailable; the calculation continues, qp_run will exit with an error."
        );
        let _ = pair_socket.set_linger(0);
        let _ = pub_socket.set_linger(0);
        return Ok(());
    }

    let mut state = PubState::Waiting;
    loop {
        let readable = pair_socket
            .poll(zmq::POLLIN, timing.pub_poll_ms)
            .map_err(|error| error.to_string())?
            > 0;
        let new_state = if readable {
            let bytes = pair_socket
                .recv_bytes(zmq::DONTWAIT)
                .map_err(|error| error.to_string())?;
            PubState::parse(&String::from_utf8_lossy(&bytes))
        } else {
            state.clone()
        };
        pub_socket
            .send(new_state.as_str(), 0)
            .map_err(|error| error.to_string())?;
        let was_stopped = state == PubState::Stopped;
        state = new_state;
        if was_stopped {
            break;
        }
    }

    let _ = pair_socket.set_linger(1_000);
    let _ = pub_socket.set_linger(1_000);
    Ok(())
}

/// Start the PUB thread and report its failure through the shared flags.
fn spawn_pub_thread(
    context: &Context,
    pub_socket: Option<Socket>,
    port: u32,
    timing: Timing,
    alive: Arc<AtomicBool>,
    pub_error: Arc<OnceLock<String>>,
) -> io::Result<thread::JoinHandle<()>> {
    let context = context.clone();
    thread::Builder::new()
        .name("task-server-pub".to_string())
        .spawn(move || {
            let result = pub_thread_main(&context, pub_socket, port, timing, &alive, &pub_error);
            if let Err(message) = result {
                alive.store(false, Ordering::Relaxed);
                let _ = pub_error.set(message);
            }
            // Once this thread is gone, a PAIR send would have no reader.
            alive.store(false, Ordering::Relaxed);
        })
}

/// Main request loop of [`run_with`].
fn main_loop(
    program: &mut ProgramState,
    rep_socket: &Socket,
    pub_channel: &PubChannel,
    debug: bool,
    timing: Timing,
) -> Result<(), ServerError> {
    while program.running {
        if rep_socket.poll(zmq::POLLIN, timing.loop_poll_ms)? == 0 {
            continue;
        }

        if let Some(bar) = program.progress_bar.as_mut()
            && bar.is_dirty()
        {
            bar.display()?;
        }

        let frames = rep_socket.recv_multipart(0)?;
        let (raw, rest) = frames
            .split_first()
            .ok_or_else(|| ServerError::Failure("Badly formed message".to_string()))?;
        let message = Message::parse(raw)?;

        if debug {
            write_debug(&debug_counts(&program.queue, &message.encode_text()))?;
        }

        match dispatch(program, &message, rest, rep_socket, pub_channel, debug) {
            Ok(()) => {}
            Err(ReplyError::Fatal(error)) => return Err(error),
            Err(ReplyError::Failure(text)) => reply_failure(rep_socket, &text, raw)?,
            Err(ReplyError::Assert(site)) => reply_failure(rep_socket, site, raw)?,
        }
    }
    Ok(())
}

/// Run the task server (`TaskServer.run`).
///
/// `rep_socket` and `pub_socket` are the sockets `qp_run` pre-bound on `P`
/// and `P+1`; passing `None` selects the reference's legacy bind path, which
/// binds its own socket with the 60 × 1 s retry. `address` is the resolved
/// `tcp://<ip>:<port>` string printed at startup; the caller must resolve it
/// once before the first server thread starts (issue #393), which is why
/// [`run`] never calls [`resolve_address`] itself.
///
/// Returns after a `terminate` request has shut the loop down and the PUB
/// thread has published its final `Stopped`.
///
/// # Errors
///
/// [`ServerError`] when a bind, socket, parse or write operation fails; a
/// parse failure is fatal, exactly as `Message.of_string` is outside the
/// reference's error handler.
pub fn run(
    rep_socket: Option<Socket>,
    pub_socket: Option<Socket>,
    port: u32,
    address: &str,
) -> Result<RunOutcome, ServerError> {
    run_with(rep_socket, pub_socket, port, address, Timing::default())
}

/// [`run`] with explicit timing, for unit tests.
fn run_with(
    rep_socket: Option<Socket>,
    pub_socket: Option<Socket>,
    port: u32,
    address: &str,
    timing: Timing,
) -> Result<RunOutcome, ServerError> {
    let context = Context::new();

    // The PAIR endpoint must exist before the PUB thread connects to it.
    let pair_socket = context.socket(SocketType::PAIR)?;
    pair_socket.bind(PAIR_ENDPOINT)?;

    // The reference starts the PUB thread before binding the REP socket.
    // Binding REP first is unobservable and lets a REP failure return
    // without leaving a live PUB thread behind.
    let rep_socket = match rep_socket {
        Some(socket) => socket,
        None => {
            let socket = context.socket(SocketType::REP)?;
            socket.set_linger(1_000_000)?;
            bind_socket(&socket, "REP", port, timing.bind_retries, timing.bind_delay)
                .map_err(ServerError::Failure)?;
            socket
        }
    };

    let alive = Arc::new(AtomicBool::new(true));
    let pub_error = Arc::new(OnceLock::new());
    let pub_thread = spawn_pub_thread(
        &context,
        pub_socket,
        port + 1,
        timing,
        Arc::clone(&alive),
        Arc::clone(&pub_error),
    )?;

    let pub_channel = PubChannel {
        socket: pair_socket,
        alive,
    };

    let mut program = ProgramState::new();
    let debug = debug_enabled();

    println!("Task server running : {address}");
    io::stdout().flush()?;

    let result = main_loop(&mut program, &rep_socket, &pub_channel, debug, timing);

    // Stop the PUB thread even after a fatal request error. The reference
    // lets the server thread die and leaks the PUB thread; a Rust thread
    // must not return while a detached thread waits forever on the PAIR.
    let _ = pub_channel.send(&PubState::Stopped);
    let _ = pub_thread.join();
    drop(pub_channel);
    drop(rep_socket);
    drop(context);

    result?;
    Ok(RunOutcome {
        pub_error: pub_error.get().cloned(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Endpoint and numeric port of a socket bound to `tcp://127.0.0.1:*`.
    fn endpoint_port(socket: &Socket) -> (String, u32) {
        let endpoint = socket
            .get_last_endpoint()
            .expect("get_last_endpoint")
            .expect("endpoint is UTF-8");
        let port = endpoint
            .rsplit(':')
            .next()
            .expect("endpoint has a port")
            .parse()
            .expect("numeric port");
        (endpoint, port)
    }

    /// Fast timing for tests; every value is the reference default scaled
    /// down, never a different algorithm.
    fn test_timing() -> Timing {
        Timing {
            bind_retries: 1,
            bind_delay: Duration::from_millis(1),
            loop_poll_ms: 10,
            pub_poll_ms: 10,
            stop_deadline: Duration::from_secs(5),
            stop_poll_ms: 10,
        }
    }

    /// A REQ client with a 5 s reply timeout.
    struct TestClient {
        socket: Socket,
    }

    impl TestClient {
        fn connect(context: &Context, endpoint: &str) -> Self {
            let socket = context.socket(SocketType::REQ).unwrap();
            socket.set_linger(0).unwrap();
            socket.connect(endpoint).unwrap();
            Self { socket }
        }

        fn send(&self, text: &str) {
            self.socket.send(text, 0).unwrap();
        }

        fn send_frames(&self, frames: &[&[u8]]) {
            self.socket
                .send_multipart(frames.iter().map(|frame| frame.to_vec()), 0)
                .unwrap();
        }

        fn recv(&self) -> String {
            let mut items = [self.socket.as_poll_item(zmq::POLLIN)];
            assert!(zmq::poll(&mut items, 5_000).unwrap() > 0, "no reply");
            String::from_utf8(self.socket.recv_bytes(0).unwrap()).unwrap()
        }

        fn recv_frames(&self) -> Vec<String> {
            let mut items = [self.socket.as_poll_item(zmq::POLLIN)];
            assert!(zmq::poll(&mut items, 5_000).unwrap() > 0, "no reply");
            self.socket
                .recv_multipart(0)
                .unwrap()
                .into_iter()
                .map(|frame| String::from_utf8(frame).unwrap())
                .collect()
        }
    }

    /// Bind a REP and a PUB socket on ephemeral ports and start the server
    /// in a thread. Returns the REP endpoint and the join handle.
    fn spawn_server() -> (String, thread::JoinHandle<Result<RunOutcome, ServerError>>) {
        let context = Context::new();
        let rep = context.socket(SocketType::REP).unwrap();
        rep.bind("tcp://127.0.0.1:*").unwrap();
        let (endpoint, port) = endpoint_port(&rep);
        let pub_socket = context.socket(SocketType::PUB).unwrap();
        pub_socket.bind("tcp://127.0.0.1:*").unwrap();
        let address = endpoint.clone();
        let handle = thread::spawn(move || {
            run_with(Some(rep), Some(pub_socket), port, &address, test_timing())
        });
        (endpoint, handle)
    }

    /// Parse `"connect_reply <state> <id> <address>"`.
    fn parse_connect_reply(reply: &str) -> (String, Id) {
        let fields: Vec<&str> = reply.splitn(4, ' ').collect();
        assert_eq!(fields.len(), 4, "reply: {reply}");
        assert_eq!(fields[0], "connect_reply");
        let state = fields[1].to_string();
        let client_id = fields[2].parse::<i64>().unwrap();
        (state, Id::new(client_id))
    }

    /// Parse `"<id> <task>"`.
    fn parse_task(frame: &str) -> (Id, String) {
        let (id, task) = frame.split_once(' ').expect("id and task");
        (Id::new(id.parse().unwrap()), task.to_string())
    }

    /// Parse `"get_task_reply <id> <task>"`.
    fn parse_get_task_reply(reply: &str) -> (Id, String) {
        parse_task(
            reply
                .strip_prefix("get_task_reply ")
                .expect("get_task_reply prefix"),
        )
    }

    #[test]
    fn pub_state_round_trip() {
        assert_eq!(PubState::parse("Waiting"), PubState::Waiting);
        assert_eq!(PubState::parse("Stopped"), PubState::Stopped);
        assert_eq!(PubState::parse("scf"), PubState::Running("scf".to_string()));
        assert_eq!(PubState::parse("").as_str(), "");
        assert_eq!(PubState::Waiting.as_str(), "Waiting");
        assert_eq!(PubState::Stopped.as_str(), "Stopped");
    }

    /// Oracle: the reference formula and its two branches
    /// (`TaskServer.ml:87-132`).
    struct FakeLookup {
        hostname: Result<String, String>,
        host_error: bool,
        interface: String,
    }

    impl HostLookup for FakeLookup {
        fn hostname(&self) -> Result<String, String> {
            self.hostname.clone()
        }

        fn ipv4_for_host(&self, host: &str) -> Result<String, String> {
            if self.host_error {
                Err("not found".to_string())
            } else if host == "node" || host == "localhost" {
                Ok("10.0.0.7".to_string())
            } else {
                Err(format!("unexpected host {host}"))
            }
        }

        fn ipv4_for_interface(&self, interface: &str) -> String {
            if interface == self.interface {
                "192.168.0.5".to_string()
            } else {
                "error: no such interface".to_string()
            }
        }
    }

    fn fake(hostname: Result<String, String>, host_error: bool) -> FakeLookup {
        FakeLookup {
            hostname,
            host_error,
            interface: "eth0".to_string(),
        }
    }

    #[test]
    fn resolve_address_from_nic() {
        let lookup = fake(Ok("node".to_string()), false);
        assert_eq!(
            resolve_address(&lookup, Some("eth0")),
            Ok("192.168.0.5".to_string())
        );
        assert_eq!(
            resolve_address(&lookup, Some("wlan0")),
            Err(ResolveError::Interface("wlan0".to_string()))
        );
    }

    #[test]
    fn resolve_address_from_hostname() {
        let lookup = fake(Ok("node".to_string()), false);
        assert_eq!(resolve_address(&lookup, None), Ok("10.0.0.7".to_string()));
    }

    #[test]
    fn resolve_address_falls_back_to_localhost() {
        // The fake only resolves "node" and "localhost", so success proves
        // the fallback was used.
        let lookup = fake(Err("gethostname failed".to_string()), false);
        assert_eq!(resolve_address(&lookup, None), Ok("10.0.0.7".to_string()));
    }

    #[test]
    fn resolve_address_host_failure() {
        let lookup = fake(Ok("node".to_string()), true);
        assert_eq!(
            resolve_address(&lookup, None),
            Err(ResolveError::HostLookup)
        );
        assert_eq!(
            ResolveError::HostLookup.to_string(),
            "Unable to find IP address from host name."
        );
        assert_eq!(
            ResolveError::Interface("eth0".to_string()).to_string(),
            "Unable to use network interface eth0"
        );
    }

    /// The two accessors added to `QueuingSystem` for `get_tasks`.
    #[test]
    fn queue_accessors() {
        let mut queue = QueuingSystem::new();
        let task_id = queue.add_task("payload");
        assert_eq!(queue.task_payload(task_id), Some("payload"));
        queue.del_task(task_id).unwrap();
        assert_eq!(queue.task_payload(task_id), None);

        let client_id = queue.add_client();
        assert!(queue.has_client(client_id));
        queue.del_client(client_id).unwrap();
        assert!(!queue.has_client(client_id));
    }

    /// The full dispatch path, in the reference's order: no job, new job,
    /// duplicate, connect, add/get/task_done, get_tasks, del_task,
    /// put/get_data, the two `end_job` branches, disconnect, terminate.
    #[test]
    fn server_round_trip() {
        let (endpoint, handle) = spawn_server();
        let context = Context::new();
        let client = TestClient::connect(&context, &endpoint);

        client.send("get_task test 1");
        assert_eq!(client.recv(), "error No job is running");

        client.send("new_job test tcp://127.0.0.1:1 inproc://push");
        assert_eq!(client.recv(), "ok");

        client.send("new_job test tcp://127.0.0.1:1 inproc://push");
        assert_eq!(client.recv(), "error A job is already running");

        client.send("connect tcp");
        let (state, client_id) = parse_connect_reply(&client.recv());
        assert_eq!(state, "test");

        client.send("add_task test aa|bb");
        assert_eq!(client.recv(), "ok");

        client.send(&format!("get_task test {client_id}"));
        let (task_a, payload) = parse_get_task_reply(&client.recv());
        assert_eq!(payload, "aa");
        client.send(&format!("task_done test {client_id} {task_a}"));
        assert_eq!(client.recv(), "ok");

        client.send(&format!("get_tasks test {client_id} 5"));
        let frames = client.recv_frames();
        assert_eq!(frames[0], "get_tasks_reply ok");
        assert_eq!(frames.len(), 3, "frames: {frames:?}");
        let (task_b, payload) = parse_task(&frames[1]);
        assert_eq!(payload, "bb");
        assert_eq!(frames[2], "0 terminate");

        client.send(&format!("del_task test {task_b}"));
        // `more` counts every payload still in the map, including the one
        // still marked running, exactly as `TaskServer.del_task` does.
        assert_eq!(client.recv(), format!("del_task_reply more {task_b}"));
        client.send(&format!("task_done test {client_id} {task_b}"));
        assert_eq!(client.recv(), "ok");

        let put = format!("put_data test {client_id} key");
        client.send_frames(&[put.as_bytes(), b"value"]);
        assert_eq!(client.recv(), "put_data_reply ok");

        client.send(&format!("get_data test {client_id} key"));
        let frames = client.recv_frames();
        assert_eq!(frames[0], "get_data_reply 5");
        assert_eq!(frames[1], "value");

        client.send(&format!("get_data test {client_id} missing"));
        let frames = client.recv_frames();
        assert_eq!(frames[0], "get_data_reply 1");
        assert_eq!(frames[1].as_bytes(), b"\0");

        client.send("end_job test");
        assert_eq!(client.recv(), "error waiting for 1 slaves...");

        client.send(&format!("disconnect test {client_id}"));
        assert_eq!(client.recv(), "disconnect_reply test");

        client.send("end_job test");
        assert_eq!(client.recv(), "ok");

        client.send("terminate");
        assert_eq!(client.recv(), "ok");
        assert_eq!(handle.join().unwrap().unwrap().pub_error, None);
    }

    /// `abort` drains every queued task, not just one. Three tasks are
    /// queued, one is popped, and `abort` must remove the other two: a
    /// discard-one implementation would leave `bb` for the next `get_task`.
    #[test]
    fn abort_drains_all_queued_tasks() {
        let (endpoint, handle) = spawn_server();
        let context = Context::new();
        let client = TestClient::connect(&context, &endpoint);

        client.send("new_job test tcp://127.0.0.1:1 inproc://push");
        assert_eq!(client.recv(), "ok");
        client.send("connect tcp");
        let (_, client_id) = parse_connect_reply(&client.recv());
        client.send("add_task test aa|bb|cc");
        assert_eq!(client.recv(), "ok");

        // One task is running while `abort` drains the queue; `abort` must
        // leave it alone (the reference only pops queued tasks).
        client.send(&format!("get_task test {client_id}"));
        let (task_aa, payload) = parse_get_task_reply(&client.recv());
        assert_eq!(payload, "aa");

        client.send("abort");
        assert_eq!(client.recv(), "ok");

        client.send(&format!("get_task test {client_id}"));
        assert_eq!(client.recv(), "get_task_reply 0");

        client.send(&format!("task_done test {client_id} {task_aa}"));
        assert_eq!(client.recv(), "ok");

        // `abort` stops accepting clients.
        client.send("connect tcp");
        assert_eq!(client.recv(), "error Wrong state");

        client.send("terminate");
        assert_eq!(client.recv(), "ok");
        handle.join().unwrap().unwrap();
    }

    /// The state commands follow the dispatch order: without a job,
    /// `set_running` falls into the `None, _` arm (`TaskServer.ml:933`)
    /// while `set_waiting` and `set_stopped` are matched with a wildcard
    /// state (`TaskServer.ml:925-926`). `end_job force` resets the job
    /// whatever its state (`TaskServer.ml:250-255`).
    #[test]
    fn pub_commands_and_force_end() {
        let (endpoint, handle) = spawn_server();
        let context = Context::new();
        let client = TestClient::connect(&context, &endpoint);

        client.send("set_running");
        assert_eq!(client.recv(), "error No job is running");
        client.send("set_waiting");
        assert_eq!(client.recv(), "ok");

        client.send("new_job test tcp://127.0.0.1:1 inproc://push");
        assert_eq!(client.recv(), "ok");
        client.send("set_running");
        assert_eq!(client.recv(), "ok");

        client.send("end_job force");
        assert_eq!(client.recv(), "ok");

        // The job is gone: the `None` arm answers again.
        client.send("get_task test 1");
        assert_eq!(client.recv(), "error No job is running");

        // `set_stopped` ends the state channel; keep it last.
        client.send("set_stopped");
        assert_eq!(client.recv(), "ok");
        client.send("terminate");
        assert_eq!(client.recv(), "ok");
        handle.join().unwrap().unwrap();
    }

    /// `get_tasks` on an unknown client trips the queue's membership assert
    /// even when the queue is empty (`Queuing_system.ml:68:2`).
    #[test]
    fn get_tasks_unknown_client_asserts() {
        let (endpoint, handle) = spawn_server();
        let context = Context::new();
        let client = TestClient::connect(&context, &endpoint);

        client.send("new_job test tcp://127.0.0.1:1 inproc://push");
        assert_eq!(client.recv(), "ok");
        client.send("get_tasks test 99 1");
        assert_eq!(
            client.recv(),
            "error Queuing_system.ml:68:2 : get_tasks test 99 1"
        );

        client.send("terminate");
        assert_eq!(client.recv(), "ok");
        handle.join().unwrap().unwrap();
    }

    /// `connect ipc` hits the reference's `assert false` and the reply names
    /// the source site.
    #[test]
    fn connect_ipc_asserts() {
        let (endpoint, handle) = spawn_server();
        let context = Context::new();
        let client = TestClient::connect(&context, &endpoint);

        client.send("new_job test tcp://127.0.0.1:1 inproc://push");
        assert_eq!(client.recv(), "ok");
        client.send("connect ipc");
        assert_eq!(client.recv(), "error TaskServer.ml:312:46 : connect ipc");

        client.send("terminate");
        assert_eq!(client.recv(), "ok");
        handle.join().unwrap().unwrap();
    }

    /// `put_data` requires exactly one extra frame.
    #[test]
    fn put_data_frame_count() {
        let (endpoint, handle) = spawn_server();
        let context = Context::new();
        let client = TestClient::connect(&context, &endpoint);

        client.send("new_job test tcp://127.0.0.1:1 inproc://push");
        assert_eq!(client.recv(), "ok");
        client.send("connect tcp");
        let (_, client_id) = parse_connect_reply(&client.recv());

        let put = format!("put_data test {client_id} key");
        client.send(&put);
        assert_eq!(
            client.recv(),
            format!("error Badly formed put_data message : {put}")
        );
        client.send_frames(&[put.as_bytes(), b"a", b"b"]);
        assert_eq!(
            client.recv(),
            format!("error Badly formed put_data message : {put}")
        );

        let put = format!("put_data test {client_id} key");
        client.send_frames(&[put.as_bytes(), b"value"]);
        assert_eq!(client.recv(), "put_data_reply ok");

        client.send("terminate");
        assert_eq!(client.recv(), "ok");
        handle.join().unwrap().unwrap();
    }

    /// The state channel publishes `Waiting` and finally `Stopped`.
    #[test]
    fn pub_channel_publishes_states() {
        let context = Context::new();
        let rep = context.socket(SocketType::REP).unwrap();
        rep.bind("tcp://127.0.0.1:*").unwrap();
        let (endpoint, port) = endpoint_port(&rep);
        let pub_socket = context.socket(SocketType::PUB).unwrap();
        pub_socket.bind("tcp://127.0.0.1:*").unwrap();
        let (pub_endpoint, _) = endpoint_port(&pub_socket);

        let sub = context.socket(SocketType::SUB).unwrap();
        sub.set_subscribe(b"").unwrap();
        sub.set_linger(0).unwrap();
        sub.connect(&pub_endpoint).unwrap();
        // PUB/SUB joins asynchronously; give the subscription time to reach
        // the publisher before the first state is sent.
        thread::sleep(Duration::from_millis(100));

        let address = endpoint.clone();
        let handle = thread::spawn(move || {
            run_with(Some(rep), Some(pub_socket), port, &address, test_timing())
        });
        let client = TestClient::connect(&context, &endpoint);

        let mut states = Vec::new();
        let mut items = [sub.as_poll_item(zmq::POLLIN)];
        assert!(zmq::poll(&mut items, 5_000).unwrap() > 0, "no state");
        states.push(String::from_utf8(sub.recv_bytes(0).unwrap()).unwrap());
        assert_eq!(states[0], "Waiting");

        client.send("terminate");
        assert_eq!(client.recv(), "ok");
        handle.join().unwrap().unwrap();

        loop {
            let mut items = [sub.as_poll_item(zmq::POLLIN)];
            assert!(zmq::poll(&mut items, 5_000).unwrap() > 0, "no Stopped");
            let state = String::from_utf8(sub.recv_bytes(0).unwrap()).unwrap();
            let stopped = state == "Stopped";
            states.push(state);
            if stopped {
                break;
            }
        }
        assert!(states.iter().any(|state| state == "Stopped"));
    }

    /// A bound but silent server trips the 30 s deadline; the test uses
    /// 100 ms.
    #[test]
    fn stop_reports_dead_server() {
        let context = Context::new();
        let rep = context.socket(SocketType::REP).unwrap();
        rep.bind("tcp://127.0.0.1:*").unwrap();
        let (_, port) = endpoint_port(&rep);

        let error = stop_with(port, Duration::from_millis(100), 10).unwrap_err();
        assert_eq!(
            error.to_string(),
            "Failure(\"TaskServer: no answer to the Terminate request (the task server thread is dead)\")"
        );
    }

    /// The per-request debug line, byte for byte
    /// (`TaskServer.ml:909-914`).
    #[test]
    fn debug_counts_matches_reference_format() {
        let mut queue = QueuingSystem::new();
        queue.add_task("aa");
        queue.add_task("bb");
        let client_id = queue.add_client();
        queue.pop_task(client_id).unwrap();
        assert_eq!(
            debug_counts(&queue, "get_task test 1"),
            "q:1  r:1  n:2  c:1 : get_task test 1\n"
        );
    }

    /// `bind_socket` gives up after its retry budget with the reference's
    /// message (trailing space included). Two wildcard sockets cannot share
    /// a port, so the failure is immediate.
    #[test]
    fn bind_failure_message_matches_reference() {
        let context = Context::new();
        let blocker = context.socket(SocketType::REP).unwrap();
        blocker.bind("tcp://*:0").unwrap();
        let (_, port) = endpoint_port(&blocker);

        let socket = context.socket(SocketType::REP).unwrap();
        socket.set_linger(0).unwrap();
        let error = bind_socket(&socket, "REP", port, 1, Duration::from_millis(1)).unwrap_err();
        assert_eq!(
            error,
            format!("Unable to bind the REP socket to port : {port} ")
        );
    }

    /// A PUB bind failure is not fatal: the calculation continues, the
    /// message is recorded in the outcome, and the server still answers
    /// `terminate` (`TaskServer.ml:765-793`).
    #[test]
    fn pub_bind_failure_is_reported_but_not_fatal() {
        let context = Context::new();
        let rep = context.socket(SocketType::REP).unwrap();
        rep.bind("tcp://127.0.0.1:*").unwrap();
        let (endpoint, _) = endpoint_port(&rep);
        // Occupy the wildcard port the PUB thread is about to bind.
        let blocker = context.socket(SocketType::PUB).unwrap();
        blocker.bind("tcp://*:0").unwrap();
        let (_, blocked_port) = endpoint_port(&blocker);
        let port = blocked_port - 1;

        let address = endpoint.clone();
        let handle =
            thread::spawn(move || run_with(Some(rep), None, port, &address, test_timing()));
        let client = TestClient::connect(&context, &endpoint);
        client.send("terminate");
        assert_eq!(client.recv(), "ok");
        let outcome = handle.join().unwrap().unwrap();
        let expected = format!("Unable to bind the PUB socket to port : {blocked_port} ");
        assert_eq!(outcome.pub_error.as_deref(), Some(expected.as_str()));
    }

    /// `stop` drives a live server to a clean shutdown.
    #[test]
    fn stop_terminates_running_server() {
        let (endpoint, handle) = spawn_server();
        let port: u32 = endpoint.rsplit(':').next().unwrap().parse().unwrap();
        stop_with(port, Duration::from_secs(5), 10).unwrap();
        assert_eq!(handle.join().unwrap().unwrap().pub_error, None);
    }

    /// An unparseable request is fatal, like `Message.of_string` outside the
    /// reference's `try`.
    #[test]
    fn malformed_request_is_fatal() {
        let (endpoint, handle) = spawn_server();
        let context = Context::new();
        let client = TestClient::connect(&context, &endpoint);
        client.send("garbage");
        let error = handle.join().unwrap().unwrap_err();
        assert!(matches!(error, ServerError::Parse(_)), "{error}");
    }
}
