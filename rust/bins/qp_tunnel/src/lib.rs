//! `qp_tunnel` — Rust port of `ocaml/qp_tunnel.ml` (P4, wave 3).
//!
//! # What it does
//!
//! `qp_tunnel` proxies the four 0MQ channels of a QP2 master process to a
//! distant cluster through one ssh tunnel. It binds four local sockets and
//! connects the corresponding remote ones:
//!
//! | thread | local (bind) | remote (connect) |
//! |---|---|---|
//! | ocaml, id 1 | ROUTER `tcp://*:<lp>` | DEALER `tcp:<host>:<p>` |
//! | fortran, id 2 | ROUTER `tcp://*:<lp+2>` | DEALER `tcp:<host>:<p+2>` |
//! | pub, id 3 | PUB `tcp://*:<lp+1>` | SUB `tcp:<host>:<p+1>` (subscribed `""`) |
//! | input, id 4 | REP `tcp://*:<lp+9>` | REQ `tcp:<host>:<p+9>` (address mode only) |
//!
//! `lp` is the first free ten-port window starting at [`DEFAULT_LOCAL_PORT`]
//! (`42379`, then `+100` per retry) on the machine's advertised IPv4 address.
//! `<p>` comes from the `(EZFIO_DIR|ADDRESS)` argument: `EZFIO` mode reads the
//! first line of `<dir>/work/qp_run_address`, `ADDRESS` mode splits the
//! argument itself on `:` into `proto:host:port`.
//!
//! The input channel is stateful. Before entering its loop it answers `test`
//! on the remote peer, and with `-g` (address mode only) it downloads the
//! remote EZFIO directory, unpacks it and writes `work/qp_run_address`. In its
//! loop it answers `test`, `get_ezfio_filename` and `get_input` (a gzipped tar
//! of the local EZFIO directory, or a forward to the remote peer).
//!
//! # Units and conventions
//!
//! Ports are `u32`; 0MQ timeouts are milliseconds. Messages are raw bytes
//! (OCaml strings are byte strings). Every user-visible string is the
//! reference's, byte for byte, except where a doc comment or the list below
//! says otherwise.
//!
//! # Deviations from the reference (all measured; see
//! `NOTES.md` in this directory)
//!
//! * The four connect endpoints keep the reference's single-slash
//!   `tcp:<host>:<port>` form. `libzmq` 4.3.5 rejects that form, so only an
//!   address argument written `tcp://<host>:<port>` works end to end (its
//!   host field keeps the `//`, and the rebuilt endpoint is well formed).
//!   This is the reference's own behaviour, not a port choice.
//! * The reference leaks a socket when `create_socket` fails and then hangs
//!   in `zmq_ctx_term`; this port closes every socket and exits 2. A thread
//!   failure prints the reference's `Thread <id> killed on uncaught
//!   exception <e>` lines, but the final `Fatal error: exception` line is the
//!   lowest-id failing thread's exception instead of the reference's
//!   `Unix.Unix_error(Unix.EINTR, "zmq_ctx_term", "")`.
//! * Linger is 0 on every socket, so shutdown takes at most one poll timeout
//!   and queued messages are dropped. The reference keeps the default
//!   (infinite) linger.
//! * Blocking receives on the remote input peer are polled in one-second
//!   slices so a `SIGUSR1`/`SIGINT` can end a peer that never answers; the
//!   reference blocks forever there.
//! * The `qp_tunnel*.tar.gz` temp file is removed when the input thread
//!   exits; the reference leaves it behind.
//! * An existing plain directory (no `.version`) gets its `work/` created and
//!   then fails reading `work/qp_run_address`; the reference fails earlier
//!   with `Unix.Unix_error(Unix.EEXIST, "mkdir", "<dir>")`.
//! * Signal handlers are installed after the threads start (the reference
//!   installs them after port discovery, before the threads), so a signal in
//!   that microsecond window takes the default action.
//!
//! # Running it
//!
//! ```text
//! cargo run -p qp_tunnel -- tcp://<host>:<port>
//! cargo run -p qp_tunnel -- -g tcp://<host>:<port>
//! cargo test -j4 -p qp_tunnel
//! ```

#![deny(missing_docs)]

mod error;
mod ocaml_num;

use std::fs;
use std::io::{self, Write as _};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;
use std::thread::{self, JoinHandle};

use qp_cli::host::SystemHost;
use qp_ezfio::EZFIO_VERSION;
use qp_zmq::task_server::resolve_address_from_env;
use tempfile::NamedTempFile;
use zmq::{Context, Socket, SocketType};

pub use error::{EXIT_FATAL, OcamlException};

/// First candidate port of the local ten-port window (`qp_tunnel.ml:7`).
pub const DEFAULT_LOCAL_PORT: u32 = 42379;

/// Poll timeout of every loop and of the interruptible receives, in
/// milliseconds (`~timeout:1000` in the reference).
const POLL_TIMEOUT_MS: i64 = 1000;

/// Offset of the input channel inside the window (`localport+9`).
const INPUT_OFFSET: u32 = 9;

/// Offset of the pub channel inside the window (`localport+1`).
const PUB_OFFSET: u32 = 1;

/// Offset of the fortran channel inside the window (`localport+2`).
const FORTRAN_OFFSET: u32 = 2;

/// Process-wide stop request, set by the signal handler and by a failing
/// thread. `Tunnel::start` clears it before spawning anything.
static STOP: AtomicBool = AtomicBool::new(false);

/// The `(EZFIO_DIR|ADDRESS)` argument, classified exactly as the reference
/// does (`Sys.file_exists x && Sys.is_directory x`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Target {
    /// An existing directory: read `<dir>/work/qp_run_address`.
    Ezfio(PathBuf),
    /// Anything else: a `proto:host:port` address.
    Address(String),
}

impl Target {
    /// Classify one command-line argument.
    ///
    /// A path that exists, is a directory (following symlinks, like
    /// `Sys.is_directory`) and can be stat'ed becomes [`Target::Ezfio`];
    /// everything else — a file, a missing path, an address — becomes
    /// [`Target::Address`].
    #[must_use]
    pub fn detect(argument: &str) -> Self {
        let path = Path::new(argument);
        if path.is_dir() {
            Self::Ezfio(path.to_path_buf())
        } else {
            Self::Address(argument.to_string())
        }
    }
}

/// The validated command line plus the classified target.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Options {
    /// `-g`: download the EZFIO directory (address mode only).
    pub get_input: bool,
    /// The single anonymous argument.
    pub target: Target,
}

/// The four proxy threads of one tunnel.
///
/// Dropping the value without calling [`Tunnel::join`] requests a stop and
/// joins every thread that is still running, so a forgotten tunnel cannot
/// block the 0MQ context's destructor.
pub struct Tunnel {
    /// Kept alive until every thread has been joined.
    ///
    /// The field is a lifetime guard, never read: dropping it first (fields
    /// drop in declaration order) would terminate a 0MQ context whose
    /// sockets the threads still use.
    #[allow(dead_code)]
    context: Context,
    /// One slot per thread, in id order: 1 ocaml, 2 fortran, 3 pub, 4 input.
    threads: [Option<JoinHandle<()>>; 4],
    /// Failing threads report `(id, exception)` here.
    errors: mpsc::Receiver<(u32, OcamlException)>,
}

impl Tunnel {
    /// Resolve the advertised IPv4 address (the ported `TaskServer.ip_address`,
    /// including `QP_NIC`) and start the tunnel.
    ///
    /// `announce` receives the text blocks the reference prints, in order:
    /// `Connect to:` on the calling thread, then `Communication [ … ]`, the
    /// `-g` block, the ssh instructions and `Ready` on the input thread.
    ///
    /// # Errors
    ///
    /// [`OcamlException::Failure`] with the reference's message when the
    /// address cannot be resolved, the address file is malformed, or the
    /// address argument is malformed; [`OcamlException::SysError`] when the
    /// address file cannot be read; `Err` from `announce` is propagated.
    pub fn start<F>(options: Options, announce: F) -> Result<Self, OcamlException>
    where
        F: FnMut(&str) -> Result<(), OcamlException> + Send + 'static,
    {
        let address = resolve_address_from_env(&SystemHost)
            .map_err(|error| OcamlException::Failure(error.to_string()))?;
        Self::start_with_ip(options, &address, announce)
    }

    /// [`Tunnel::start`] with an explicit advertised address.
    ///
    /// This is the in-process test seam: `QP_NIC` cannot be set safely from a
    /// multi-threaded test process, so tests pass `"127.0.0.1"` here.
    ///
    /// # Errors
    ///
    /// The same as [`Tunnel::start`], except that the address lookup is
    /// skipped.
    pub fn start_with_ip<F>(options: Options, ip: &str, announce: F) -> Result<Self, OcamlException>
    where
        F: FnMut(&str) -> Result<(), OcamlException> + Send + 'static,
    {
        // A new tunnel must not inherit the stop request of a previous one.
        STOP.store(false, Ordering::SeqCst);

        // Order is the reference's and observable: the long address is read
        // (and, in EZFIO mode, `work/` is created) before the port window is
        // searched.
        let remote = match &options.target {
            Target::Ezfio(directory) => {
                ensure_work(directory)?;
                let long_address = read_address(directory)?;
                parse_address(&long_address)?
            }
            Target::Address(address) => parse_address(address)?,
        };

        let context = Context::new();
        let local_port = discover_local_port(&context, ip, DEFAULT_LOCAL_PORT)?;
        let (sender, errors) = mpsc::channel();
        let mut threads: [Option<JoinHandle<()>>; 4] = [None, None, None, None];

        // Thread 1 (the OCaml channel) starts before the banner, exactly as in
        // the reference, so its `bind` failure is reported after `Connect to:`.
        threads[0] = Some(spawn_relay(
            1,
            &context,
            &sender,
            RelaySpec {
                in_endpoint: format!("tcp://*:{local_port}"),
                out_endpoint: remote.endpoint(0),
                in_type: SocketType::ROUTER,
                out_type: SocketType::DEALER,
                bind_in: true,
                forward_out: true,
            },
        ));

        let mut announce = announce;
        if let Err(error) = announce(&format!("Connect to:\ntcp://{ip}:{local_port}\n")) {
            request_stop();
            join_all(&mut threads);
            return Err(error);
        }

        threads[1] = Some(spawn_relay(
            2,
            &context,
            &sender,
            RelaySpec {
                in_endpoint: format!("tcp://*:{}", local_port + FORTRAN_OFFSET),
                out_endpoint: remote.endpoint(FORTRAN_OFFSET),
                in_type: SocketType::ROUTER,
                out_type: SocketType::DEALER,
                bind_in: true,
                forward_out: true,
            },
        ));
        threads[2] = Some(spawn_relay(
            3,
            &context,
            &sender,
            RelaySpec {
                in_endpoint: remote.endpoint(PUB_OFFSET),
                out_endpoint: format!("tcp://*:{}", local_port + PUB_OFFSET),
                in_type: SocketType::SUB,
                out_type: SocketType::PUB,
                bind_in: false,
                forward_out: false,
            },
        ));

        let setup = InputSetup {
            context: context.clone(),
            target: options.target.clone(),
            get_input: options.get_input,
            ip: ip.to_string(),
            local_port,
            remote: remote.clone(),
            remote_input: matches!(&options.target, Target::Address(_)),
        };
        threads[3] = Some(spawn_input(&sender, setup, announce));

        Ok(Self {
            context,
            threads,
            errors,
        })
    }

    /// Ask every thread to leave its loop at the next check.
    ///
    /// This is what the SIGUSR1/SIGINT handler does. The reference also resets
    /// the signal to its default action inside the handler; this port keeps
    /// its handler installed, so a second signal is harmless.
    pub fn request_stop() {
        request_stop();
    }

    /// Join the four threads in the reference's order and return the
    /// lowest-id failure, if any.
    ///
    /// Reference order: input (4), fortran (2), pub (3), ocaml (1). A failing
    /// thread has already printed its `Thread <id> killed …` line.
    ///
    /// # Errors
    ///
    /// The exception of the lowest-id thread that failed.
    ///
    /// # Panics
    ///
    /// If a thread panicked: that is a port bug (the reference has no Rust
    /// panics), and hiding it would turn it into a clean exit.
    pub fn join(mut self) -> Result<(), OcamlException> {
        for index in [3, 1, 2, 0] {
            if let Some(handle) = self.threads[index].take()
                && handle.join().is_err()
            {
                panic!("the qp_tunnel thread with id {} panicked", index + 1);
            }
        }
        let mut failures: Vec<(u32, OcamlException)> = self.errors.try_iter().collect();
        failures.sort_by_key(|(id, _)| *id);
        match failures.into_iter().next() {
            Some((_, error)) => Err(error),
            None => Ok(()),
        }
    }
}

impl Drop for Tunnel {
    fn drop(&mut self) {
        request_stop();
        join_all(&mut self.threads);
    }
}

/// Join every thread slot, ignoring a panic (this runs during unwinding too).
fn join_all(threads: &mut [Option<JoinHandle<()>>; 4]) {
    for slot in threads {
        if let Some(handle) = slot.take() {
            let _ = handle.join();
        }
    }
}

/// Write `text` to stdout and flush it (the reference's `%!`).
///
/// # Errors
///
/// [`OcamlException::SysError`] when stdout cannot be written. The text is
/// Rust's message, not the reference's `Sys_error("output_string: …")`; this
/// path is unmeasured because a closed stdout never reaches the comparison.
pub fn emit(text: &str) -> Result<(), OcamlException> {
    let mut stdout = io::stdout().lock();
    stdout
        .write_all(text.as_bytes())
        .and_then(|()| stdout.flush())
        .map_err(|error| OcamlException::SysError(error.to_string()))
}

/// The ssh instructions block printed by the input thread before `Ready`.
///
/// `hostname` is `Unix.gethostname ()` (first `ssh -L` line only; the second
/// line ends with the literal `<host>`, as measured), `ip` and `local_port`
/// describe the local side, and `remote_port` is the remote window's first
/// port (an OCaml `int`, rendered in decimal, sign included). The result
/// starts with a newline, indents both `ssh` lines by one space, and ends
/// each line with ` &`.
#[must_use]
pub fn ssh_instructions(hostname: &str, ip: &str, remote_port: i64, local_port: u32) -> String {
    let r0 = remote_port;
    let r1 = remote_port + 1;
    let r2 = remote_port + 2;
    let r9 = remote_port + i64::from(INPUT_OFFSET);
    let l0 = local_port;
    let l1 = local_port + 1;
    let l2 = local_port + 2;
    let l9 = local_port + INPUT_OFFSET;
    format!(
        "\nOn remote hosts, create ssh tunnel using:\n ssh -L {r0}:{ip}:{l0} -L {r1}:{ip}:{l1} -L {r2}:{ip}:{l2} -L {r9}:{ip}:{l9} {hostname} &\nOr from this host connect to clients using:\n ssh -R {r0}:localhost:{l0} -R {r1}:localhost:{l1} -R {r2}:localhost:{l2} -R {r9}:localhost:{l9} <host> &\n"
    )
}

/// Install the `SIGUSR1` and `SIGINT` handlers that call [`Tunnel::request_stop`].
///
/// # Panics
///
/// If `libc::signal` refuses a handler, which POSIX does not allow for these
/// two signals with a valid function pointer; the message names the signal.
pub fn install_signal_handlers() {
    // SAFETY: `SIGUSR1` and `SIGINT` are valid signal numbers and
    // `handle_stop` is an `extern "C" fn(c_int)`, the signature `signal(3)`
    // expects. The handler only stores into a lock-free `AtomicBool`, which
    // is async-signal-safe.
    let previous = unsafe { libc::signal(libc::SIGUSR1, handle_stop as libc::sighandler_t) };
    assert_ne!(
        previous,
        libc::SIG_ERR,
        "libc::signal(SIGUSR1, handler) must succeed"
    );
    // SAFETY: same contract as the SIGUSR1 registration above.
    let previous = unsafe { libc::signal(libc::SIGINT, handle_stop as libc::sighandler_t) };
    assert_ne!(
        previous,
        libc::SIG_ERR,
        "libc::signal(SIGINT, handler) must succeed"
    );
}

/// The one signal handler: record the stop request.
extern "C" fn handle_stop(_signal: libc::c_int) {
    STOP.store(true, Ordering::SeqCst);
}

/// Whether a stop was requested (signal or failing thread).
fn stop_requested() -> bool {
    STOP.load(Ordering::SeqCst)
}

/// Set the stop request (the handler cannot call this because of the `&`
/// borrow rules, so the store is spelled out there).
fn request_stop() {
    STOP.store(true, Ordering::SeqCst);
}

/// Search for a ten-port window that can be bound on `ip`, starting at `base`
/// and adding 100 per failed candidate.
///
/// One fresh REP socket probes each port: the `zmq` crate exposes no
/// `unbind`, so the socket is dropped (closed) instead of unbound. The
/// reference's unbounded scan is kept.
///
/// # Errors
///
/// [`OcamlException::SysError`] when 0MQ cannot create a probe socket, which
/// only happens while the context is terminating.
pub fn discover_local_port(context: &Context, ip: &str, base: u32) -> Result<u32, OcamlException> {
    let mut candidate = base;
    loop {
        if window_is_free(context, ip, candidate)? {
            return Ok(candidate);
        }
        // The reference retries forever (`try_new_port (port_number+100)`).
        candidate = candidate.saturating_add(100);
    }
}

/// Whether all ten ports `base..=base+9` can be bound on `ip`.
fn window_is_free(context: &Context, ip: &str, base: u32) -> Result<bool, OcamlException> {
    for offset in 0..10 {
        let port = base.saturating_add(offset);
        let endpoint = format!("tcp://{ip}:{port}");
        let socket = context.socket(SocketType::REP).map_err(|error| {
            OcamlException::SysError(format!("Unable to create a REP socket: {error}"))
        })?;
        socket.set_linger(0).map_err(zmq_exception)?;
        let bound = socket.bind(&endpoint);
        drop(socket);
        if bound.is_err() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The remote side of a tunnel: the host field as split from the address
/// (it keeps a leading `//` for a `tcp://` argument) and the first port.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Remote {
    /// Host field of the split address, verbatim.
    host: String,
    /// First port of the remote window, an OCaml `int` (`i64`): the
    /// reference's `int_of_string` accepts any value in `[-2^62, 2^62-1]`,
    /// and the endpoints render it in decimal (a negative port is printed
    /// with its sign, measured).
    port: i64,
}

impl Remote {
    /// The reference's connect endpoint `tcp:<host>:<port+offset>`
    /// (single slash, exactly as printed).
    fn endpoint(&self, offset: u32) -> String {
        format!("tcp:{}:{}", self.host, self.port + i64::from(offset))
    }
}

/// Split `proto:host:port` and parse the port as the reference's
/// `int_of_string` (full OCaml semantics: `0x`/`0o`/`0b`, `_` separators,
/// a leading sign, and the `[-2^62, 2^62 - 1]` bound).
fn parse_address(long_address: &str) -> Result<Remote, OcamlException> {
    let parts: Vec<&str> = long_address.split(':').collect();
    let [_, host, port] = parts.as_slice() else {
        return Err(OcamlException::Failure(format!(
            "{long_address} : Malformed address"
        )));
    };
    let port = ocaml_num::int_of_string(port)
        .map_err(|_| OcamlException::Failure(String::from("int_of_string")))?;
    Ok(Remote {
        host: (*host).to_string(),
        port,
    })
}

/// Read `<directory>/work/qp_run_address`: first line, OCaml-trimmed.
///
/// # Errors
///
/// [`OcamlException::SysError`] when the file cannot be read (`<path>:
/// <strerror>`, with no `open_in:` prefix: that is the measured text), and
/// [`OcamlException::EndOfFile`] when it is empty.
fn read_address(directory: &Path) -> Result<String, OcamlException> {
    let path = directory.join("work").join("qp_run_address");
    let contents = fs::read_to_string(&path).map_err(|error| error::sys_error(&path, &error))?;
    if contents.is_empty() {
        return Err(OcamlException::EndOfFile);
    }
    let first_line = contents.split('\n').next().unwrap_or("");
    Ok(ocaml_trim(first_line).to_string())
}

/// `String.trim`: the five ASCII characters OCaml's `is_space` accepts.
///
/// `\x0b` (vertical tab) is deliberately not one of them, unlike Rust's
/// Unicode `str::trim`.
fn ocaml_trim(text: &str) -> &str {
    text.trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\r' | '\u{c}'))
}

/// `Qpackage.ezfio_work`: make sure `<directory>/work` exists.
///
/// The reference calls `Ezfio.set_file` + `set_work_empty false` only when
/// `work/` is missing, which creates `work/`, writes the group marker
/// `work/.version` (`EZFIO_VERSION`) and `work/empty` for a real EZFIO tree
/// (`.version` present). A plain directory fails there with `EEXIST`; this
/// port creates `work/` for it too, which is what makes `-g` work on a
/// downloaded tree that has neither `.version` nor `work/`.
///
/// # Errors
///
/// [`OcamlException::SysError`] when `work/`, `work/.version` or
/// `work/empty` cannot be created.
fn ensure_work(directory: &Path) -> Result<(), OcamlException> {
    let work = directory.join("work");
    if work.exists() {
        return Ok(());
    }
    fs::create_dir_all(&work).map_err(|error| error::sys_error(&work, &error))?;
    if directory.join(".version").exists() {
        // `Ezfio.mkdir` writes the marker before `write_bool` writes the
        // boolean, and aborts in that order on failure.
        let version = work.join(".version");
        fs::write(&version, format!("{EZFIO_VERSION}\n"))
            .map_err(|error| error::sys_error(&version, &error))?;
        let marker = work.join("empty");
        fs::write(&marker, "F\n\n").map_err(|error| error::sys_error(&marker, &error))?;
    }
    Ok(())
}

/// A relay thread specification: two sockets, which side binds, and which
/// directions flow.
///
/// `in_endpoint` is the socket the local side receives from: the ROUTER
/// threads bind it locally, the pub thread connects it to the remote PUB.
/// `out_endpoint` is the opposite peer.
struct RelaySpec {
    /// Endpoint of `in_type`; bound when `bind_in`, connected otherwise.
    in_endpoint: String,
    /// Endpoint of `out_type`; always the opposite of `bind_in`.
    out_endpoint: String,
    /// Socket the local side receives from.
    in_type: SocketType,
    /// Socket facing the remote peer.
    out_type: SocketType,
    /// Whether `in_endpoint` is bound locally (ROUTER) or connected (SUB).
    bind_in: bool,
    /// Whether `out_type` may forward locally produced messages.
    forward_out: bool,
}

/// Start one relay thread from its [`RelaySpec`].
fn spawn_relay(
    id: u32,
    context: &Context,
    errors: &mpsc::Sender<(u32, OcamlException)>,
    spec: RelaySpec,
) -> JoinHandle<()> {
    let context = context.clone();
    let errors = errors.clone();
    thread::spawn(move || run_thread(id, errors, || relay_loop(&context, &spec)))
}

/// Forward one direction of a relay until a stop is requested.
///
/// The reference's arms are kept: when both sockets are readable, the
/// remote→local message is forwarded first; the pub thread never forwards
/// remote→local (a PUB socket is never readable).
fn relay_loop(context: &Context, spec: &RelaySpec) -> Result<(), OcamlException> {
    let socket_in = create_socket(context, spec.in_type, &spec.in_endpoint, spec.bind_in)?;
    let socket_out = create_socket(context, spec.out_type, &spec.out_endpoint, !spec.bind_in)?;
    if spec.in_type == SocketType::SUB {
        socket_in.set_subscribe(b"").map_err(zmq_exception)?;
    }

    let mut items = [
        socket_in.as_poll_item(zmq::POLLIN),
        socket_out.as_poll_item(zmq::POLLIN),
    ];
    loop {
        if stop_requested() {
            break;
        }
        match zmq::poll(&mut items, POLL_TIMEOUT_MS) {
            Ok(_) => {}
            Err(zmq::Error::EINTR) => continue,
            Err(error) => return Err(zmq_exception(error)),
        }
        let in_ready = items[0].is_readable();
        let out_ready = items[1].is_readable();
        if out_ready && spec.forward_out {
            let frames = recv_multipart(&socket_out)?;
            send_multipart(&socket_in, &frames)?;
        }
        if in_ready {
            let frames = recv_multipart(&socket_in)?;
            send_multipart(&socket_out, &frames)?;
        }
    }
    Ok(())
}

/// Everything the input thread needs; owned because it runs detached.
struct InputSetup {
    /// Shared 0MQ context.
    context: Context,
    /// Classified `(EZFIO_DIR|ADDRESS)` argument.
    target: Target,
    /// `-g`.
    get_input: bool,
    /// Locally advertised IPv4 address.
    ip: String,
    /// First port of the local window.
    local_port: u32,
    /// Remote window (always known, even in EZFIO mode).
    remote: Remote,
    /// Whether the input thread opens the remote REQ socket (address mode).
    remote_input: bool,
}

/// The sink for lines the input thread announces (stdout in the binary).
type Announce = Box<dyn FnMut(&str) -> Result<(), OcamlException>>;

/// Spawn thread 4 (the stateful input channel).
fn spawn_input<F>(
    errors: &mpsc::Sender<(u32, OcamlException)>,
    setup: InputSetup,
    announce: F,
) -> JoinHandle<()>
where
    F: FnMut(&str) -> Result<(), OcamlException> + Send + 'static,
{
    let errors = errors.clone();
    thread::spawn(move || {
        let mut announce: Announce = Box::new(announce);
        run_thread(4, errors, || input_loop(&setup, &mut announce))
    })
}

/// Run one thread body and report an uncaught exception like the reference.
fn run_thread<F>(id: u32, errors: mpsc::Sender<(u32, OcamlException)>, body: F)
where
    F: FnOnce() -> Result<(), OcamlException>,
{
    if let Err(error) = body() {
        let _ = writeln!(
            io::stderr(),
            "Thread {id} killed on uncaught exception {}",
            error.render()
        );
        let _ = errors.send((id, error));
        request_stop();
    }
}

/// The input thread: setup, optional download, then the request loop.
fn input_loop(setup: &InputSetup, announce: &mut Announce) -> Result<(), OcamlException> {
    match input_body(setup, announce.as_mut()) {
        Ok(()) | Err(InputError::Stopped) => Ok(()),
        Err(InputError::Fatal(exception)) => Err(exception),
    }
}

/// Why the input thread stopped before its loop was asked to.
enum InputError {
    /// A stop was requested while waiting for the remote peer: a clean exit.
    Stopped,
    /// A reference exception.
    Fatal(OcamlException),
}

impl From<OcamlException> for InputError {
    fn from(exception: OcamlException) -> Self {
        Self::Fatal(exception)
    }
}

/// The body of [`input_loop`].
fn input_body(
    setup: &InputSetup,
    announce: &mut (dyn FnMut(&str) -> Result<(), OcamlException> + '_),
) -> Result<(), InputError> {
    let bind_endpoint = format!("tcp://*:{}", setup.local_port + INPUT_OFFSET);
    let socket_in = create_socket(&setup.context, SocketType::REP, &bind_endpoint, true)?;
    let socket_out = if setup.remote_input {
        Some(create_socket(
            &setup.context,
            SocketType::REQ,
            &setup.remote.endpoint(INPUT_OFFSET),
            false,
        )?)
    } else {
        None
    };
    let temp = tempfile::Builder::new()
        .prefix("qp_tunnel")
        .suffix(".tar.gz")
        .tempfile()
        .map_err(|error| {
            InputError::Fatal(OcamlException::SysError(format!(
                "temp_file: {}",
                error::os_error_message(&error)
            )))
        })?;
    let session = InputSession {
        setup,
        socket_out,
        temp,
    };

    // Address mode only: the reference prints nothing in EZFIO mode.
    if let Some(socket) = &session.socket_out {
        let reply = request(socket, b"test")?;
        announce(&format!(
            "Communication [ {} ]\n",
            String::from_utf8_lossy(&reply)
        ))
        .map_err(InputError::Fatal)?;
    }

    // `-g` is a no-op in EZFIO mode, exactly like the reference.
    if setup.get_input && matches!(&setup.target, Target::Address(_)) {
        session.download(announce)?;
    }

    let hostname = qp_cli::ffi::hostname()
        .map_err(|error| InputError::Fatal(OcamlException::SysError(error.to_string())))?;
    announce(&ssh_instructions(
        &hostname,
        &setup.ip,
        setup.remote.port,
        setup.local_port,
    ))
    .map_err(InputError::Fatal)?;
    announce("Ready\n").map_err(InputError::Fatal)?;

    loop {
        if stop_requested() {
            break;
        }
        match socket_in.poll(zmq::POLLIN, POLL_TIMEOUT_MS) {
            Ok(0) => {}
            Ok(_) => {
                let request_bytes = recv_interruptible(&socket_in)?;
                let reply = match request_bytes.as_slice() {
                    b"get_input" => session.get_input()?,
                    b"get_ezfio_filename" => session.get_ezfio_filename()?,
                    b"test" => b"OK".to_vec(),
                    other => {
                        let mut reply = b"Message '".to_vec();
                        reply.extend_from_slice(other);
                        reply.extend_from_slice(b"' not understood");
                        reply
                    }
                };
                send_bytes(&socket_in, &reply)?;
            }
            Err(zmq::Error::EINTR) => {}
            Err(error) => return Err(InputError::Fatal(zmq_exception(error))),
        }
    }
    Ok(())
}

/// The sockets and temp file of one input thread.
struct InputSession<'a> {
    /// Thread configuration.
    setup: &'a InputSetup,
    /// Remote REQ socket; `None` in EZFIO mode.
    socket_out: Option<Socket>,
    /// `qp_tunnel*.tar.gz`, kept alive (and deleted on drop) for the thread.
    temp: NamedTempFile,
}

impl InputSession<'_> {
    /// `get_ezfio_filename`: the EZFIO directory, or a forward to the peer.
    ///
    /// # Errors
    ///
    /// A remote 0MQ failure, or a stop while the peer is silent.
    fn get_ezfio_filename(&self) -> Result<Vec<u8>, InputError> {
        match (&self.setup.target, &self.socket_out) {
            (Target::Ezfio(directory), _) => Ok(directory.as_os_str().as_encoded_bytes().to_vec()),
            (Target::Address(_), Some(socket)) => request(socket, b"get_ezfio_filename"),
            // Unreachable by construction (address mode always creates it).
            (Target::Address(_), None) => Err(InputError::Fatal(OcamlException::Failure(
                String::from("address mode has no remote socket"),
            ))),
        }
    }

    /// `get_input`: tar the EZFIO directory, or forward to the peer.
    ///
    /// The EZFIO branch runs the reference's exact shell command
    /// (`tar --exclude="*.gz.*" -zcf <temp> <dir>`), ignores its status, then
    /// reads the whole archive: one message. The archive is read while the
    /// named temp file is still open, which is harmless on Unix.
    ///
    /// # Errors
    ///
    /// [`OcamlException::SysError`] when the temp file cannot be read or the
    /// remote peer fails; [`InputError::Stopped`] when a stop arrives first.
    fn get_input(&self) -> Result<Vec<u8>, InputError> {
        match &self.setup.target {
            Target::Ezfio(directory) => {
                let command = format!(
                    "tar --exclude=\"*.gz.*\" -zcf {} {}",
                    self.temp.path().display(),
                    directory.display()
                );
                let _ = std::process::Command::new("/bin/sh")
                    .arg("-c")
                    .arg(command)
                    .status();
                fs::read(self.temp.path())
                    .map_err(|error| InputError::Fatal(error::sys_error(self.temp.path(), &error)))
            }
            Target::Address(_) => match &self.socket_out {
                Some(socket) => request(socket, b"get_input"),
                // Unreachable by construction (address mode always creates it).
                None => Err(InputError::Fatal(OcamlException::Failure(String::from(
                    "address mode has no remote socket",
                )))),
            },
        }
    }

    /// The `-g` block: download, unpack, write `work/qp_run_address`.
    ///
    /// The printed order is the reference's, including the missing newline
    /// between `Getting input... ` and the directory name and the space before
    /// `...done`.
    ///
    /// # Errors
    ///
    /// A remote failure, a `tar` that leaves no readable file, or a
    /// [`ensure_work`] failure.
    fn download(
        &self,
        announce: &mut (dyn FnMut(&str) -> Result<(), OcamlException> + '_),
    ) -> Result<(), InputError> {
        announce("Getting input... ").map_err(InputError::Fatal)?;
        let filename = self.get_ezfio_filename()?;
        let name = String::from_utf8_lossy(&filename);
        announce(&name).map_err(InputError::Fatal)?;
        let tarball = self.get_input()?;
        let temp_path = self.temp.path().to_path_buf();
        fs::write(&temp_path, &tarball)
            .map_err(|error| InputError::Fatal(error::sys_error(&temp_path, &error)))?;
        let _ = std::process::Command::new("/bin/sh")
            .arg("-c")
            .arg(format!("tar -zxf {}", temp_path.display()))
            .status();
        let directory = PathBuf::from(name.as_ref());
        ensure_work(&directory)?;
        let address_path = directory.join("work").join("qp_run_address");
        fs::write(
            &address_path,
            format!("tcp://{}:{}\n", self.setup.ip, self.setup.local_port),
        )
        .map_err(|error| InputError::Fatal(error::sys_error(&address_path, &error)))?;
        announce(" ...done\n").map_err(InputError::Fatal)?;
        Ok(())
    }
}

/// `create_socket`: create and bind/connect, mapping every failure to the
/// reference's message.
fn create_socket(
    context: &Context,
    kind: SocketType,
    endpoint: &str,
    bind: bool,
) -> Result<Socket, OcamlException> {
    let socket = context
        .socket(kind)
        .map_err(|_| connection_failure(endpoint))?;
    let result = if bind {
        socket.bind(endpoint)
    } else {
        socket.connect(endpoint)
    };
    result.map_err(|_| connection_failure(endpoint))?;
    socket.set_linger(0).map_err(zmq_exception)?;
    Ok(socket)
}

/// `Failure("Unable to establish connection to <endpoint>.")`.
fn connection_failure(endpoint: &str) -> OcamlException {
    OcamlException::Failure(format!("Unable to establish connection to {endpoint}."))
}

/// A 0MQ error outside `create_socket`: surfaced as a `Sys_error`.
///
/// The reference lets its binding's own exception escape here; no test
/// reaches this path, so the text is not a measured reference string.
fn zmq_exception(error: zmq::Error) -> OcamlException {
    OcamlException::SysError(format!("0MQ: {error}"))
}

/// Send one request frame and receive one reply frame on a REQ socket.
fn request(socket: &Socket, bytes: &[u8]) -> Result<Vec<u8>, InputError> {
    send_bytes(socket, bytes)?;
    recv_interruptible(socket)
}

/// Send one message, retrying `EINTR`.
fn send_bytes(socket: &Socket, bytes: &[u8]) -> Result<(), OcamlException> {
    loop {
        match socket.send(bytes, 0) {
            Ok(()) => return Ok(()),
            Err(zmq::Error::EINTR) => continue,
            Err(error) => return Err(zmq_exception(error)),
        }
    }
}

/// Send a received multipart message frame by frame, retrying `EINTR`.
///
/// Frame-by-frame keeps `SNDMORE` consistent: `send_multipart` would abort
/// mid-message if a signal interrupted it, and the peer would then see a
/// truncated multipart message.
fn send_multipart(socket: &Socket, frames: &[Vec<u8>]) -> Result<(), OcamlException> {
    let last = frames.len().saturating_sub(1);
    for (index, frame) in frames.iter().enumerate() {
        let flags = if index == last { 0 } else { zmq::SNDMORE };
        loop {
            match socket.send(frame.as_slice(), flags) {
                Ok(()) => break,
                Err(zmq::Error::EINTR) => continue,
                Err(error) => return Err(zmq_exception(error)),
            }
        }
    }
    Ok(())
}

/// Receive one multipart message, retrying `EINTR`.
fn recv_multipart(socket: &Socket) -> Result<Vec<Vec<u8>>, OcamlException> {
    loop {
        match socket.recv_multipart(0) {
            Ok(frames) => return Ok(frames),
            Err(zmq::Error::EINTR) => continue,
            Err(error) => return Err(zmq_exception(error)),
        }
    }
}

/// Receive one frame, in one-second poll slices so a stop can interrupt it.
fn recv_interruptible(socket: &Socket) -> Result<Vec<u8>, InputError> {
    loop {
        match socket.poll(zmq::POLLIN, POLL_TIMEOUT_MS) {
            Ok(0) => {
                if stop_requested() {
                    return Err(InputError::Stopped);
                }
            }
            Ok(_) => loop {
                match socket.recv_bytes(0) {
                    Ok(bytes) => return Ok(bytes),
                    Err(zmq::Error::EINTR) => {
                        if stop_requested() {
                            return Err(InputError::Stopped);
                        }
                    }
                    Err(error) => return Err(InputError::Fatal(zmq_exception(error))),
                }
            },
            Err(zmq::Error::EINTR) => {}
            Err(error) => return Err(InputError::Fatal(zmq_exception(error))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The ssh block must be the captured reference text.
    #[test]
    fn ssh_block_matches_the_reference_capture() {
        let capture = include_str!("../tests/data/ssh_block.txt");
        assert_eq!(
            ssh_instructions("lcpqsv225.univ-tlse3.fr", "130.120.229.225", 43000, 42379),
            capture
        );
        assert_eq!(capture.len(), 356);
        // The second line really ends with the literal `<host>`.
        assert!(capture.ends_with("-R 43009:localhost:42388 <host> &\n"));
    }

    #[test]
    fn address_parsing_matches_the_reference() {
        assert_eq!(
            parse_address("tcp://127.0.0.1:43000"),
            Ok(Remote {
                host: "//127.0.0.1".into(),
                port: 43000,
            })
        );
        assert_eq!(
            parse_address("tcp:127.0.0.1:1"),
            Ok(Remote {
                host: "127.0.0.1".into(),
                port: 1,
            })
        );
        let malformed = |text: &str| OcamlException::Failure(format!("{text} : Malformed address"));
        assert_eq!(parse_address("tcp:1.2.3.4"), Err(malformed("tcp:1.2.3.4")));
        assert_eq!(
            parse_address("notanaddress"),
            Err(malformed("notanaddress"))
        );
        assert_eq!(parse_address("a:b:c:d"), Err(malformed("a:b:c:d")));
        assert_eq!(parse_address(""), Err(malformed("")));
        let int_of_string = || OcamlException::Failure("int_of_string".into());
        assert_eq!(parse_address("tcp:host:notaport"), Err(int_of_string()));
        assert_eq!(parse_address("tcp:host:"), Err(int_of_string()));
        // Above the OCaml `int` bound (`2^62 - 1`), measured:
        assert_eq!(
            parse_address("tcp:host:4611686018427387904"),
            Err(int_of_string())
        );
        // The full OCaml `int_of_string` grammar is accepted (measured:
        // `0x10` is port 16, `-5` is port -5), and the port is an OCaml
        // `int`, so no port-range rejection exists.
        assert_eq!(parse_address("tcp:host:0x10").map(|r| r.port), Ok(16));
        assert_eq!(parse_address("tcp:host:-5").map(|r| r.port), Ok(-5));
        assert_eq!(parse_address("tcp:host:+1").map(|r| r.port), Ok(1));
        assert_eq!(parse_address("tcp:host:1_000").map(|r| r.port), Ok(1000));
        assert_eq!(
            parse_address("tcp:host:4611686018427387903").map(|r| r.port),
            Ok(4611686018427387903)
        );
    }

    #[test]
    fn endpoints_keep_the_reference_single_slash() {
        let remote = Remote {
            host: "//127.0.0.1".into(),
            port: 43000,
        };
        assert_eq!(remote.endpoint(0), "tcp://127.0.0.1:43000");
        assert_eq!(remote.endpoint(INPUT_OFFSET), "tcp://127.0.0.1:43009");
    }

    #[test]
    fn trim_is_ocaml_string_trim() {
        assert_eq!(ocaml_trim(" \t\r\n\u{c}abc \u{c}"), "abc");
        assert_eq!(ocaml_trim(""), "");
        // OCaml does not treat the vertical tab as a space.
        assert_eq!(ocaml_trim("\u{b}abc"), "\u{b}abc");
    }

    #[test]
    fn address_file_first_line_is_trimmed() {
        let dir = tempfile::tempdir().expect("tempdir");
        let work = dir.path().join("work");
        fs::create_dir(&work).expect("work");
        let path = work.join("qp_run_address");
        fs::write(&path, "  tcp://127.0.0.1:43000 \r\nsecond line\n").expect("write");
        assert_eq!(
            read_address(dir.path()).expect("address"),
            "tcp://127.0.0.1:43000"
        );
        // No trailing newline: the whole file is the first line.
        fs::write(&path, "tcp:h:1").expect("write");
        assert_eq!(read_address(dir.path()).expect("address"), "tcp:h:1");
    }

    #[test]
    fn empty_address_file_is_end_of_file() {
        let dir = tempfile::tempdir().expect("tempdir");
        let work = dir.path().join("work");
        fs::create_dir(&work).expect("work");
        fs::write(work.join("qp_run_address"), "").expect("write");
        assert_eq!(read_address(dir.path()), Err(OcamlException::EndOfFile));
    }

    #[test]
    fn missing_address_file_names_the_path() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("work").join("qp_run_address");
        assert_eq!(
            read_address(dir.path()),
            Err(OcamlException::SysError(format!(
                "{}: No such file or directory",
                path.display()
            )))
        );
    }

    #[test]
    fn ensure_work_matches_qpackage() {
        // A real EZFIO tree: `work/`, the group marker and `empty` are created.
        let dir = tempfile::tempdir().expect("tempdir");
        fs::write(dir.path().join(".version"), "1\n").expect(".version");
        ensure_work(dir.path()).expect("work");
        assert_eq!(
            fs::read_to_string(dir.path().join("work").join(".version")).expect("version"),
            format!("{EZFIO_VERSION}\n")
        );
        assert_eq!(
            fs::read_to_string(dir.path().join("work").join("empty")).expect("empty"),
            "F\n\n"
        );
        // A plain directory: `work/` is created, no EZFIO marker is invented.
        let plain = tempfile::tempdir().expect("tempdir");
        ensure_work(plain.path()).expect("work");
        assert!(plain.path().join("work").is_dir());
        assert!(!plain.path().join("work").join("empty").exists());
        assert!(!plain.path().join("work").join(".version").exists());
        // An existing `work/` is left untouched.
        fs::write(plain.path().join("work").join("keep"), "x").expect("keep");
        ensure_work(plain.path()).expect("work again");
        assert_eq!(
            fs::read_to_string(plain.path().join("work").join("keep")).expect("keep"),
            "x"
        );
        assert!(!plain.path().join("work").join("empty").exists());
    }

    #[test]
    fn target_detection_is_is_directory() {
        let dir = tempfile::tempdir().expect("tempdir");
        assert_eq!(
            Target::detect(&dir.path().to_string_lossy()),
            Target::Ezfio(dir.path().to_path_buf())
        );
        let file = dir.path().join("f");
        fs::write(&file, "x").expect("file");
        assert_eq!(
            Target::detect(&file.to_string_lossy()),
            Target::Address(file.to_string_lossy().into_owned())
        );
        assert_eq!(
            Target::detect("tcp:127.0.0.1:1"),
            Target::Address("tcp:127.0.0.1:1".into())
        );
    }

    #[test]
    fn port_discovery_finds_a_free_window() {
        let context = Context::new();
        let base = 45000;
        let first = discover_local_port(&context, "127.0.0.1", base).expect("window");
        assert!(first >= base);
        assert_eq!((first - base) % 100, 0);
        // A window that was free is bindable end to end.
        let mut sockets = Vec::new();
        for offset in 0..10 {
            let socket = context.socket(SocketType::REP).expect("socket");
            socket
                .bind(&format!("tcp://127.0.0.1:{}", first + offset))
                .expect("every probed port must be free");
            sockets.push(socket);
        }
        drop(sockets);
    }

    #[test]
    fn port_discovery_jumps_by_100_when_the_window_is_taken() {
        let context = Context::new();
        let base = 47000;
        let first = discover_local_port(&context, "127.0.0.1", base).expect("window");
        let taken = context.socket(SocketType::REP).expect("socket");
        taken
            .bind(&format!("tcp://127.0.0.1:{first}"))
            .expect("bind the first port of the free window");
        let second = discover_local_port(&context, "127.0.0.1", base).expect("window");
        assert!(second > first, "{second} must follow {first}");
        assert_eq!((second - base) % 100, 0);
    }
}
