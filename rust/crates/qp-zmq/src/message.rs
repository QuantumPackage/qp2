//! Byte-level codec for the QP2 task-server protocol.
//!
//! The wire format is the plain-text protocol of `ocaml/Message.ml` and
//! `ocaml/Message_lexer.mll`: one **text frame** per request (plus, for
//! `put_data`, one raw payload frame), and replies that are either the same
//! shape or a text frame followed by raw payload frames (`get_data_reply`,
//! `get_tasks_reply`).
//!
//! # Encoding
//!
//! [`Message::encode`] reproduces the single-line rendering of the reference
//! `Message.to_string` byte for byte; [`Message::encode_frames`] splits it
//! into the frames the reference sends (`Zmq.Socket.send` versus
//! `Zmq.Socket.send_all`). The eight reply-only kinds have no parser in the
//! reference lexer and are encode-only here as well.
//!
//! State names, task strings and keys are free-form text; the protocol has no
//! escaping, so a field that contains a space cannot be encoded back into the
//! same token. That is a property of the protocol, not of this port.
//!
//! # Parsing
//!
//! [`Message::parse`] is a hand-written equivalent of the OCaml lexer. The
//! grammar has three token classes with non-obvious rules, all preserved:
//!
//! * `word`: a maximal run of bytes other than space, tab and line feed.
//! * `text`: a maximal run of bytes with no `|`, of length at least two,
//!   whose first byte is not a space. Because the reference class is
//!   `[^ '|']+` for the tail, a trailing space **is** part of the token
//!   (`"error x "` parses as the error text `"x "`), while a lone byte before
//!   a `|` or the end of input is skipped.
//! * `integer`: digits only, parsed as an OCaml 63-bit `int`; a larger run is
//!   `Failure "int_of_string"`.
//!
//! Unknown bytes at keyword position are skipped one at a time, so
//! `"okay"` is `Ok` and `"terminate_x"` is `Terminate`. Input that contains
//! no keyword raises `Failure "lexing: empty token"`.
//!
//! The reference protocol is stringly typed: task and client identifiers
//! share one representation (`Id.Task.t` and `Id.Client.t` are both `int`),
//! which [`Id`] mirrors. Identifiers are strictly positive; the `x > 0`
//! assertions that the reference performs while *constructing* a message are
//! reported as [`ParseError`] values here instead of aborting the process.
//!
//! # Deliberate differences
//!
//! * Tokens must be valid UTF-8. The reference keeps raw bytes; QP2 clients
//!   only send ASCII, and rejecting early is what lets the server use `str`.
//! * The `Assert_failure` location strings in [`ParseError::AssertFailure`]
//!   identify the reference source site approximately; the server only prints
//!   them in its `TASK SERVER CRASH` diagnostic.

use std::fmt;

/// Largest value of an OCaml native `int` on a 64-bit platform (`2^62 - 1`).
const OCAML_INT_MAX: i64 = (1 << 62) - 1;

/// A task or client identifier.
///
/// The reference has two module aliases (`Id.Task` and `Id.Client`) over the
/// same `int`; the wire format carries no distinction, so the port uses one
/// type. The value is always strictly positive, as enforced by
/// [`Id::new`](Id::new) and by the parser.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Id(i64);

impl Id {
    /// Build an identifier from its wire value.
    ///
    /// # Panics
    ///
    /// Panics when `value <= 0`: identifiers are strictly positive, and
    /// [`Message::parse`] rejects such values before calling this
    /// constructor.
    #[must_use]
    pub fn new(value: i64) -> Self {
        assert!(
            value > 0,
            "Id::new: identifiers are strictly positive, got {value}"
        );
        Self(value)
    }

    /// The wire value (`Id.Task.to_int` / `Id.Client.to_int`).
    #[must_use]
    pub fn get(self) -> i64 {
        self.0
    }

    /// The next identifier (`Id.increment`).
    #[must_use]
    pub(crate) fn next(self) -> Self {
        Self::new(self.0 + 1)
    }
}

impl fmt::Display for Id {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// Connection type of a `connect` request (`Message.Connect_msg.t`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConnectType {
    /// `connect tcp`.
    Tcp,
    /// `connect inproc`.
    Inproc,
    /// `connect ipc`; encodes and parses, but the reference server asserts
    /// when it dispatches it (there is no IPC push address).
    Ipc,
}

impl ConnectType {
    /// The keyword used on the wire.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Tcp => "tcp",
            Self::Inproc => "inproc",
            Self::Ipc => "ipc",
        }
    }
}

/// A protocol message: the 18 parseable request kinds plus the 8 reply kinds
/// that the reference only ever encodes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Message {
    /// `new_job <state> <tcp> <inproc>`: a job starts.
    Newjob {
        /// Program state name (the job's `zmq_state`).
        state: String,
        /// Advertised TCP endpoint of the PUSH socket.
        address_tcp: String,
        /// Advertised inproc endpoint of the PUSH socket.
        address_inproc: String,
    },
    /// `end_job <state>`: the job is finished (or `force`).
    Endjob {
        /// State the client believes it is in, or `force`.
        state: String,
    },
    /// `connect <tcp|inproc|ipc>`: a worker asks for a PUSH address.
    Connect {
        /// Transport requested by the worker.
        socket: ConnectType,
    },
    /// `connect_reply <state> <client_id> <push_address>` (reply only).
    ConnectReply {
        /// Current state name.
        state: String,
        /// Identifier assigned to the new client.
        client_id: Id,
        /// PUSH endpoint the client must connect to.
        push_address: String,
    },
    /// `disconnect <state> <client_id>`: a worker leaves.
    Disconnect {
        /// State the client believes it is in.
        state: String,
        /// Identifier given by [`Message::ConnectReply`].
        client_id: Id,
    },
    /// `disconnect_reply <state>` (reply only).
    DisconnectReply {
        /// Current state name (echoed from the request).
        state: String,
    },
    /// `add_task <state> <task>|<task>...`: append tasks to the queue.
    AddTask {
        /// State the client believes it is in.
        state: String,
        /// Task payloads; `|` separates them on the wire.
        tasks: Vec<String>,
    },
    /// `add_task_reply <task_id>` (reply only).
    ///
    /// Part of the reference's union (`Message.ml:489`) and encoded here,
    /// but never constructed: the reference builds it nowhere either, and
    /// the server never sends it. Kept for wire-grammar parity.
    #[allow(dead_code)]
    AddTaskReply {
        /// Identifier the server assigned to the last added task.
        task_id: Id,
    },
    /// `del_task <state> <id>|<id>...`: delete tasks by identifier.
    DelTask {
        /// State the client believes it is in.
        state: String,
        /// Task identifiers to delete.
        task_ids: Vec<Id>,
    },
    /// `del_task_reply <more|done> <id>|...` (reply only).
    DelTaskReply {
        /// Identifiers that were deleted (echoed from the request).
        task_ids: Vec<Id>,
        /// Whether queued tasks remain (`more`) or not (`done`).
        more: bool,
    },
    /// `get_task <state> <client_id>`: take one task.
    GetTask {
        /// State the client believes it is in.
        state: String,
        /// Identifier of the requesting worker.
        client_id: Id,
    },
    /// `get_task_reply <task_id> <task>` (reply only); `0` when no task was
    /// available or its payload was deleted.
    GetTaskReply {
        /// Identifier of the task, if one was popped.
        task_id: Option<Id>,
        /// Payload of the task, if it still exists.
        task: Option<String>,
    },
    /// `get_tasks <state> <client_id> <n>`: take up to `n` tasks.
    GetTasks {
        /// State the client believes it is in.
        state: String,
        /// Identifier of the requesting worker.
        client_id: Id,
        /// Strictly positive number of tasks requested.
        n_tasks: i64,
    },
    /// `get_tasks_reply ok` followed by one frame per task (reply only).
    GetTasksReply {
        /// `(id, payload)` pairs; `None` is the `terminate` sentinel.
        tasks: Vec<(Option<Id>, String)>,
    },
    /// `put_data <state> <client_id> <key>` plus one raw payload frame.
    PutData {
        /// State the client believes it is in.
        state: String,
        /// Identifier of the requesting worker.
        client_id: Id,
        /// Hash-table key.
        key: String,
    },
    /// `put_data_reply ok` (reply only).
    PutDataReply,
    /// `get_data <state> <client_id> <key>`: read the hash table.
    GetData {
        /// State the client believes it is in.
        state: String,
        /// Identifier of the requesting worker.
        client_id: Id,
        /// Hash-table key.
        key: String,
    },
    /// `get_data_reply <byte_len>` plus one raw value frame (reply only).
    GetDataReply {
        /// Raw value; `[0]` (a single NUL byte) when the key is absent.
        value: Vec<u8>,
    },
    /// `task_done <state> <client_id> <id>|<id>...`: tasks finished.
    TaskDone {
        /// State the client believes it is in.
        state: String,
        /// Identifier of the worker that ran the tasks.
        client_id: Id,
        /// Identifiers of the finished tasks.
        task_ids: Vec<Id>,
    },
    /// `terminate`: stop the server loop.
    Terminate,
    /// `abort`: discard queued tasks and stop accepting clients.
    Abort,
    /// `ok`: generic positive reply.
    Ok,
    /// `error <text>`: generic negative reply (also used as a request kind).
    Error {
        /// Error text, appended after `error ` verbatim.
        text: String,
    },
    /// `set_stopped`: publish `Stopped` on the state channel.
    SetStopped,
    /// `set_running`: publish the state name on the state channel.
    SetRunning,
    /// `set_waiting`: publish `Waiting` on the state channel.
    SetWaiting,
}

impl Message {
    /// Render the message on one line, exactly as `Message.to_string`.
    ///
    /// This is the single text frame of every message except
    /// [`Message::GetDataReply`] and [`Message::GetTasksReply`], whose frame
    /// layout is given by [`Message::encode_frames`].
    #[must_use]
    pub fn encode(&self) -> Vec<u8> {
        match self {
            Self::Newjob {
                state,
                address_tcp,
                address_inproc,
            } => format!("new_job {state} {address_tcp} {address_inproc}").into_bytes(),
            Self::Endjob { state } => format!("end_job {state}").into_bytes(),
            Self::Connect { socket } => format!("connect {}", socket.as_str()).into_bytes(),
            Self::ConnectReply {
                state,
                client_id,
                push_address,
            } => format!("connect_reply {state} {client_id} {push_address}").into_bytes(),
            Self::Disconnect { state, client_id } => {
                format!("disconnect {state} {client_id}").into_bytes()
            }
            Self::DisconnectReply { state } => format!("disconnect_reply {state}").into_bytes(),
            Self::AddTask { state, tasks } => {
                format!("add_task {state} {}", tasks.join("|")).into_bytes()
            }
            Self::AddTaskReply { task_id } => format!("add_task_reply {task_id}").into_bytes(),
            Self::DelTask { state, task_ids } => {
                format!("del_task {state} {}", join_ids(task_ids)).into_bytes()
            }
            Self::DelTaskReply { task_ids, more } => {
                let more = if *more { "more" } else { "done" };
                format!("del_task_reply {more} {}", join_ids(task_ids)).into_bytes()
            }
            Self::GetTask { state, client_id } => {
                format!("get_task {state} {client_id}").into_bytes()
            }
            Self::GetTaskReply { task_id, task } => match (task_id, task) {
                (Some(task_id), Some(task)) => {
                    format!("get_task_reply {task_id} {task}").into_bytes()
                }
                _ => b"get_task_reply 0".to_vec(),
            },
            Self::GetTasks {
                state,
                client_id,
                n_tasks,
            } => format!("get_tasks {state} {client_id} {n_tasks}").into_bytes(),
            Self::GetTasksReply { .. } => b"get_tasks_reply ok".to_vec(),
            Self::PutData {
                state,
                client_id,
                key,
            } => format!("put_data {state} {client_id} {key}").into_bytes(),
            Self::PutDataReply => b"put_data_reply ok".to_vec(),
            Self::GetData {
                state,
                client_id,
                key,
            } => format!("get_data {state} {client_id} {key}").into_bytes(),
            Self::GetDataReply { value } => {
                let mut frame = format!("get_data_reply {} ", value.len()).into_bytes();
                frame.extend_from_slice(value);
                frame
            }
            Self::TaskDone {
                state,
                client_id,
                task_ids,
            } => format!("task_done {state} {client_id} {}", join_ids(task_ids)).into_bytes(),
            Self::Terminate => b"terminate".to_vec(),
            Self::Abort => b"abort".to_vec(),
            Self::Ok => b"ok".to_vec(),
            Self::Error { text } => format!("error {text}").into_bytes(),
            Self::SetStopped => b"set_stopped".to_vec(),
            Self::SetRunning => b"set_running".to_vec(),
            Self::SetWaiting => b"set_waiting".to_vec(),
        }
    }

    /// Render the message as one line of text.
    ///
    /// Uses a lossy conversion for [`Message::GetDataReply`] values, which
    /// are arbitrary bytes; the server only calls this for parsed messages,
    /// whose tokens are valid UTF-8.
    #[must_use]
    pub fn encode_text(&self) -> String {
        String::from_utf8_lossy(&self.encode()).into_owned()
    }

    /// Render the wire frames of the message.
    ///
    /// Most messages are a single text frame. [`Message::GetDataReply`] is
    /// `get_data_reply <len>` followed by the raw value, and
    /// [`Message::GetTasksReply`] is `get_tasks_reply ok` followed by one
    /// `"<id> <task>"` frame per task.
    #[must_use]
    pub fn encode_frames(&self) -> Vec<Vec<u8>> {
        match self {
            Self::GetDataReply { value } => vec![
                format!("get_data_reply {}", value.len()).into_bytes(),
                value.clone(),
            ],
            Self::GetTasksReply { tasks } => {
                let mut frames = Vec::with_capacity(tasks.len() + 1);
                frames.push(b"get_tasks_reply ok".to_vec());
                for (task_id, task) in tasks {
                    match task_id {
                        Some(task_id) => frames.push(format!("{task_id} {task}").into_bytes()),
                        None => frames.push(b"0 terminate".to_vec()),
                    }
                }
                frames
            }
            _ => vec![self.encode()],
        }
    }

    /// Parse one text frame (`Message.of_string`).
    ///
    /// # Errors
    ///
    /// Returns the exception the reference lexer or one of its constructors
    /// would raise; see [`ParseError`]. Malformed input is fatal to the
    /// reference task-server thread, so callers are expected to propagate it.
    pub fn parse(raw: &[u8]) -> Result<Self, ParseError> {
        let mut lexer = Lexer::new(raw);
        let keyword = lexer.next_keyword()?;
        match keyword {
            Keyword::AddTask => {
                let state = lexer.read_word()?.to_owned();
                let tasks = lexer.read_text()?.into_iter().map(str::to_owned).collect();
                Ok(Self::AddTask { state, tasks })
            }
            Keyword::DelTask => {
                let state = lexer.read_word()?.to_owned();
                let task_ids = lexer
                    .read_ints()?
                    .into_iter()
                    .map(parse_id)
                    .collect::<Result<_, _>>()?;
                Ok(Self::DelTask { state, task_ids })
            }
            Keyword::GetTask => {
                let state = lexer.read_word()?.to_owned();
                let client_id = parse_id(lexer.read_int()?)?;
                Ok(Self::GetTask { state, client_id })
            }
            Keyword::GetTasks => {
                let state = lexer.read_word()?.to_owned();
                let client_id = parse_id(lexer.read_int()?)?;
                let n_tasks = lexer.read_int()?;
                if n_tasks <= 0 {
                    return Err(ParseError::StrictlyPositiveInt(n_tasks));
                }
                Ok(Self::GetTasks {
                    state,
                    client_id,
                    n_tasks,
                })
            }
            Keyword::TaskDone => {
                let state = lexer.read_word()?.to_owned();
                let client_id = parse_id(lexer.read_int()?)?;
                let task_ids = lexer
                    .read_ints()?
                    .into_iter()
                    .map(parse_id)
                    .collect::<Result<_, _>>()?;
                Ok(Self::TaskDone {
                    state,
                    client_id,
                    task_ids,
                })
            }
            Keyword::Disconnect => {
                let state = lexer.read_word()?.to_owned();
                let client_id = parse_id(lexer.read_int()?)?;
                Ok(Self::Disconnect { state, client_id })
            }
            Keyword::Connect => {
                let socket = lexer.read_word()?;
                let socket = match socket {
                    "tcp" => ConnectType::Tcp,
                    "inproc" => ConnectType::Inproc,
                    "ipc" => ConnectType::Ipc,
                    _ => return Err(ParseError::AssertFailure("Message.ml:75:12")),
                };
                Ok(Self::Connect { socket })
            }
            Keyword::NewJob => {
                let state = lexer.read_word()?.to_owned();
                let address_tcp = lexer.read_word()?.to_owned();
                let address_inproc = lexer.read_word()?.to_owned();
                if !address_tcp.starts_with("tcp://") {
                    return Err(ParseError::AddressInvalid);
                }
                if !address_inproc.starts_with("inproc://") {
                    return Err(ParseError::AssertFailure("Address.ml:47:4"));
                }
                Ok(Self::Newjob {
                    state,
                    address_tcp,
                    address_inproc,
                })
            }
            Keyword::EndJob => Ok(Self::Endjob {
                state: lexer.read_word()?.to_owned(),
            }),
            Keyword::GetData => {
                let state = lexer.read_word()?.to_owned();
                let client_id = parse_id(lexer.read_int()?)?;
                let key = lexer.read_word()?.to_owned();
                Ok(Self::GetData {
                    state,
                    client_id,
                    key,
                })
            }
            Keyword::PutData => {
                let state = lexer.read_word()?.to_owned();
                let client_id = parse_id(lexer.read_int()?)?;
                let key = lexer.read_word()?.to_owned();
                Ok(Self::PutData {
                    state,
                    client_id,
                    key,
                })
            }
            Keyword::Terminate => Ok(Self::Terminate),
            Keyword::Abort => Ok(Self::Abort),
            Keyword::Ok => Ok(Self::Ok),
            Keyword::Error => match lexer.read_text()?.into_iter().next() {
                Some(text) => Ok(Self::Error {
                    text: text.to_owned(),
                }),
                None => Err(ParseError::Hd),
            },
            Keyword::SetStopped => Ok(Self::SetStopped),
            Keyword::SetRunning => Ok(Self::SetRunning),
            Keyword::SetWaiting => Ok(Self::SetWaiting),
        }
    }
}

fn join_ids(ids: &[Id]) -> String {
    ids.iter()
        .map(|id| id.get())
        .map(|id| id.to_string())
        .collect::<Vec<_>>()
        .join("|")
}

fn parse_id(value: i64) -> Result<Id, ParseError> {
    if value <= 0 {
        // `Id.of_int` asserts `x > 0` (`ocaml/Id.ml:5`).
        return Err(ParseError::AssertFailure("Id.ml:5:4"));
    }
    Ok(Id::new(value))
}

/// A parse or message-construction failure of the reference protocol.
///
/// The [`Display`](fmt::Display) form is the OCaml exception rendering, so a
/// task-server crash prints the same shape as `Printexc.to_string` in the
/// reference.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ParseError {
    /// No keyword was found (`Failure "lexing: empty token"`).
    EmptyToken,
    /// EOF inside a reader (`Failure "Error in MessageLexer (2|3|4)"`).
    MessageLexer(u8),
    /// `error` with no text token (`Failure "hd"`).
    Hd,
    /// An integer did not fit an OCaml 63-bit `int` (`Failure "int_of_string"`).
    IntOfString,
    /// A `new_job` TCP address does not start with `tcp://`.
    AddressInvalid,
    /// One of the reference `assert false` sites fired; the string names the
    /// reference source location.
    AssertFailure(&'static str),
    /// `get_tasks` asked for a non-positive number of tasks.
    StrictlyPositiveInt(i64),
    /// A token is not valid UTF-8 (deviation from the byte-string reference).
    InvalidUtf8,
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyToken => f.write_str("Failure(\"lexing: empty token\")"),
            Self::MessageLexer(n) => write!(f, "Failure(\"Error in MessageLexer ({n})\")"),
            Self::Hd => f.write_str("Failure(\"hd\")"),
            Self::IntOfString => f.write_str("Failure(\"int_of_string\")"),
            Self::AddressInvalid => f.write_str("Invalid_argument(\"Address Invalid\")"),
            Self::AssertFailure(site) => write!(f, "Assert_failure(\"{site}\")"),
            Self::StrictlyPositiveInt(value) => write!(
                f,
                "Invalid_argument(\"Strictly_positive_int : (x > 0) : x={value}\")"
            ),
            Self::InvalidUtf8 => f.write_str("Failure(\"invalid UTF-8 in message token\")"),
        }
    }
}

impl std::error::Error for ParseError {}

/// Keywords of the lexer, in no particular order: matching takes the longest
/// keyword that is a prefix of the input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Keyword {
    AddTask,
    DelTask,
    GetTask,
    GetTasks,
    TaskDone,
    Disconnect,
    Connect,
    NewJob,
    EndJob,
    GetData,
    PutData,
    Terminate,
    Abort,
    Ok,
    Error,
    SetStopped,
    SetRunning,
    SetWaiting,
}

const KEYWORDS: &[(&[u8], Keyword)] = &[
    (b"add_task", Keyword::AddTask),
    (b"del_task", Keyword::DelTask),
    (b"get_task", Keyword::GetTask),
    (b"get_tasks", Keyword::GetTasks),
    (b"task_done", Keyword::TaskDone),
    (b"disconnect", Keyword::Disconnect),
    (b"connect", Keyword::Connect),
    (b"new_job", Keyword::NewJob),
    (b"end_job", Keyword::EndJob),
    (b"get_data", Keyword::GetData),
    (b"put_data", Keyword::PutData),
    (b"terminate", Keyword::Terminate),
    (b"abort", Keyword::Abort),
    (b"ok", Keyword::Ok),
    (b"error", Keyword::Error),
    (b"set_stopped", Keyword::SetStopped),
    (b"set_running", Keyword::SetRunning),
    (b"set_waiting", Keyword::SetWaiting),
];

/// Cursor over the message bytes, reproducing the four ocamllex rules
/// `kw`, `get_word`, `get_int` and `get_text`.
struct Lexer<'a> {
    input: &'a [u8],
    pos: usize,
}

impl<'a> Lexer<'a> {
    fn new(input: &'a [u8]) -> Self {
        Self { input, pos: 0 }
    }

    fn eof(&self) -> bool {
        self.pos >= self.input.len()
    }

    /// Match `kw`: longest keyword prefix, else consume one byte and retry.
    /// EOF raises the ocamllex "empty token" failure because `kw` has no eof
    /// rule.
    fn next_keyword(&mut self) -> Result<Keyword, ParseError> {
        loop {
            let rest = &self.input[self.pos..];
            if rest.is_empty() {
                return Err(ParseError::EmptyToken);
            }
            let mut best: Option<(usize, Keyword)> = None;
            for (keyword, kind) in KEYWORDS {
                if rest.starts_with(keyword) && best.is_none_or(|(len, _)| keyword.len() > len) {
                    best = Some((keyword.len(), *kind));
                }
            }
            if let Some((len, kind)) = best {
                self.pos += len;
                return Ok(kind);
            }
            self.pos += 1;
        }
    }

    /// Match `get_word`, looping over `NONE` like `read_word`.
    fn read_word(&mut self) -> Result<&'a str, ParseError> {
        loop {
            if self.eof() {
                return Err(ParseError::MessageLexer(3));
            }
            let byte = self.input[self.pos];
            if matches!(byte, b' ' | b'\t' | b'\n') {
                self.pos += 1;
                continue;
            }
            let start = self.pos;
            while self.pos < self.input.len()
                && !matches!(self.input[self.pos], b' ' | b'\t' | b'\n')
            {
                self.pos += 1;
            }
            return std::str::from_utf8(&self.input[start..self.pos])
                .map_err(|_| ParseError::InvalidUtf8);
        }
    }

    /// Match `get_int`, looping over `NONE` like `read_int`.
    fn read_int(&mut self) -> Result<i64, ParseError> {
        loop {
            if self.eof() {
                return Err(ParseError::MessageLexer(4));
            }
            if self.input[self.pos].is_ascii_digit() {
                return self.read_digits();
            }
            self.pos += 1;
        }
    }

    /// Match `get_int`, looping over `NONE` like `read_ints`; EOF returns the
    /// identifiers accumulated so far.
    fn read_ints(&mut self) -> Result<Vec<i64>, ParseError> {
        let mut ids = Vec::new();
        loop {
            if self.eof() {
                return Ok(ids);
            }
            if self.input[self.pos].is_ascii_digit() {
                ids.push(self.read_digits()?);
            } else {
                self.pos += 1;
            }
        }
    }

    fn read_digits(&mut self) -> Result<i64, ParseError> {
        let start = self.pos;
        while self.pos < self.input.len() && self.input[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        let text = std::str::from_utf8(&self.input[start..self.pos])
            .map_err(|_| ParseError::InvalidUtf8)?;
        let value: i64 = text.parse().map_err(|_| ParseError::IntOfString)?;
        if value > OCAML_INT_MAX {
            return Err(ParseError::IntOfString);
        }
        Ok(value)
    }

    /// Match `get_text`, looping like `read_text` and accumulating every
    /// `TEXT` token until EOF. `NONE` consumes one byte.
    fn read_text(&mut self) -> Result<Vec<&'a str>, ParseError> {
        let mut tokens = Vec::new();
        while !self.eof() {
            let start = self.pos;
            let byte = self.input[start];
            if byte == b' ' || byte == b'|' {
                self.pos += 1;
                continue;
            }
            let end = self.input[start..]
                .iter()
                .position(|&b| b == b'|')
                .map_or(self.input.len(), |offset| start + offset);
            if end - start >= 2 {
                tokens.push(
                    std::str::from_utf8(&self.input[start..end])
                        .map_err(|_| ParseError::InvalidUtf8)?,
                );
                self.pos = end;
            } else {
                self.pos += 1;
            }
        }
        Ok(tokens)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn round(input: &str) {
        let message = Message::parse(input.as_bytes()).expect(input);
        assert_eq!(
            String::from_utf8_lossy(&message.encode()),
            input,
            "round trip of {input:?}"
        );
    }

    #[test]
    fn canonical_requests_round_trip_byte_exact() {
        // Golden bytes captured from the reference `Message.of_string` /
        // `Message.to_string` pair (`ocaml/Message.ml`).
        for input in [
            "new_job scf tcp://192.168.1.2:43210 inproc://scf",
            "end_job scf",
            "connect tcp",
            "connect inproc",
            "connect ipc",
            "disconnect scf 3",
            "add_task scf first|second task",
            "del_task scf 1|2",
            "get_task scf 4",
            "get_tasks scf 4 3",
            "put_data scf 4 N_states_diag",
            "get_data scf 4 N_states_diag",
            "task_done scf 4 5|6",
            "terminate",
            "abort",
            "ok",
            "error boom",
            "set_stopped",
            "set_running",
            "set_waiting",
        ] {
            round(input);
        }
    }

    #[test]
    fn parse_produces_the_reference_values() {
        assert_eq!(
            Message::parse(b"new_job scf tcp://h:1 inproc://x").unwrap(),
            Message::Newjob {
                state: "scf".into(),
                address_tcp: "tcp://h:1".into(),
                address_inproc: "inproc://x".into(),
            }
        );
        assert_eq!(
            Message::parse(b"add_task st a|b c").unwrap(),
            // A one-character token cannot be a `text` (the reference lexer
            // requires at least two characters), so the task-list rule skips
            // it and continues: the reference round-trips this input to
            // `add_task st b c`. Verified with the OCaml `golden` harness.
            Message::AddTask {
                state: "st".into(),
                tasks: vec!["b c".into()],
            }
        );
        assert_eq!(
            Message::parse(b"del_task st 7|8").unwrap(),
            Message::DelTask {
                state: "st".into(),
                task_ids: vec![Id::new(7), Id::new(8)],
            }
        );
        assert_eq!(
            Message::parse(b"connect ipc").unwrap(),
            Message::Connect {
                socket: ConnectType::Ipc
            }
        );
        assert_eq!(Message::parse(b"okay").unwrap(), Message::Ok);
        assert_eq!(Message::parse(b"terminate_x").unwrap(), Message::Terminate);
        assert_eq!(
            Message::parse(b"error boom").unwrap(),
            Message::Error {
                text: "boom".into()
            }
        );
    }

    #[test]
    fn text_tokens_include_a_trailing_space_and_skip_lone_bytes() {
        // Reference probes: `Message.to_string (Message.of_string ...)`.
        assert_eq!(
            Message::parse(b"error x ").unwrap(),
            Message::Error { text: "x ".into() }
        );
        assert_eq!(
            Message::parse(b"add_task s a ").unwrap(),
            Message::AddTask {
                state: "s".into(),
                tasks: vec!["a ".into()],
            }
        );
        // A one-byte segment is skipped, and `|` is not part of text.
        assert_eq!(
            Message::parse(b"add_task s a").unwrap(),
            Message::AddTask {
                state: "s".into(),
                tasks: vec![],
            }
        );
        assert_eq!(
            Message::parse(b"add_task s a|b").unwrap(),
            Message::AddTask {
                state: "s".into(),
                tasks: vec![],
            }
        );
        assert_eq!(
            Message::parse(b"add_task s  ab  cd|ef ").unwrap(),
            Message::AddTask {
                state: "s".into(),
                tasks: vec!["ab  cd".into(), "ef ".into()],
            }
        );
        // The encoded form of an empty task list keeps the separator space.
        assert_eq!(
            Message::parse(b"add_task s a").unwrap().encode(),
            b"add_task s "
        );
    }

    #[test]
    fn integers_use_ocaml_int_semantics() {
        assert_eq!(
            Message::parse(b"get_task s 12x").unwrap(),
            Message::GetTask {
                state: "s".into(),
                client_id: Id::new(12),
            }
        );
        assert_eq!(
            Message::parse(b"del_task s 1|2 x").unwrap(),
            Message::DelTask {
                state: "s".into(),
                task_ids: vec![Id::new(1), Id::new(2)],
            }
        );
        // 2^62 is out of range for an OCaml native int.
        assert_eq!(
            Message::parse(b"get_task s 4611686018427387904"),
            Err(ParseError::IntOfString)
        );
        assert_eq!(
            Message::parse(b"get_task s 0"),
            Err(ParseError::AssertFailure("Id.ml:5:4"))
        );
        assert_eq!(
            Message::parse(b"get_tasks s 1 0"),
            Err(ParseError::StrictlyPositiveInt(0))
        );
    }

    #[test]
    fn malformed_input_matches_reference_failures() {
        assert_eq!(Message::parse(b""), Err(ParseError::EmptyToken));
        assert_eq!(Message::parse(b"   "), Err(ParseError::EmptyToken));
        assert_eq!(Message::parse(b"ab|cd"), Err(ParseError::EmptyToken));
        assert_eq!(Message::parse(b"error x"), Err(ParseError::Hd));
        assert_eq!(Message::parse(b"error "), Err(ParseError::Hd));
        assert_eq!(
            Message::parse(b"get_task"),
            Err(ParseError::MessageLexer(3))
        );
        assert_eq!(Message::parse(b"terminate 12"), Ok(Message::Terminate));
        assert_eq!(
            Message::parse(b"new_job s udp://h:1 inproc://x"),
            Err(ParseError::AddressInvalid)
        );
        assert_eq!(
            Message::parse(b"new_job s tcp://h:1 ipc://x"),
            Err(ParseError::AssertFailure("Address.ml:47:4"))
        );
        assert_eq!(
            Message::parse(b"connect foo"),
            Err(ParseError::AssertFailure("Message.ml:75:12"))
        );
        assert_eq!(
            Message::parse(b"add_task \xff\xff x"),
            Err(ParseError::InvalidUtf8)
        );
    }

    #[test]
    fn reply_encodings_match_the_reference() {
        // Golden bytes from the reference `*_msg.to_string` functions.
        let connect_reply = Message::ConnectReply {
            state: "scf".into(),
            client_id: Id::new(3),
            push_address: "tcp://192.168.1.2:43210".into(),
        };
        assert_eq!(
            connect_reply.encode(),
            b"connect_reply scf 3 tcp://192.168.1.2:43210"
        );
        assert_eq!(
            Message::DisconnectReply {
                state: "scf".into()
            }
            .encode(),
            b"disconnect_reply scf"
        );
        assert_eq!(
            Message::AddTaskReply {
                task_id: Id::new(7)
            }
            .encode(),
            b"add_task_reply 7"
        );
        assert_eq!(
            Message::DelTaskReply {
                task_ids: vec![Id::new(1), Id::new(2)],
                more: true,
            }
            .encode(),
            b"del_task_reply more 1|2"
        );
        assert_eq!(
            Message::DelTaskReply {
                task_ids: vec![Id::new(1), Id::new(2)],
                more: false,
            }
            .encode(),
            b"del_task_reply done 1|2"
        );
        assert_eq!(
            Message::GetTaskReply {
                task_id: Some(Id::new(5)),
                task: Some("first task".into()),
            }
            .encode(),
            b"get_task_reply 5 first task"
        );
        assert_eq!(
            Message::GetTaskReply {
                task_id: None,
                task: None,
            }
            .encode(),
            b"get_task_reply 0"
        );
        assert_eq!(Message::PutDataReply.encode(), b"put_data_reply ok");

        let get_tasks = Message::GetTasksReply {
            tasks: vec![
                (Some(Id::new(5)), "first task".into()),
                (None, "terminate".into()),
            ],
        };
        assert_eq!(get_tasks.encode(), b"get_tasks_reply ok");
        assert_eq!(
            get_tasks.encode_frames(),
            vec![
                b"get_tasks_reply ok".to_vec(),
                b"5 first task".to_vec(),
                b"0 terminate".to_vec(),
            ]
        );
        assert_eq!(
            Message::GetTasksReply { tasks: vec![] }.encode_frames(),
            vec![b"get_tasks_reply ok".to_vec()]
        );

        let get_data = Message::GetDataReply {
            value: b"hello".to_vec(),
        };
        assert_eq!(get_data.encode(), b"get_data_reply 5 hello");
        assert_eq!(
            get_data.encode_frames(),
            vec![b"get_data_reply 5".to_vec(), b"hello".to_vec()]
        );
    }

    #[test]
    fn replies_are_not_parsed_as_replies() {
        // The reference lexer has no reply keywords; `get_data_reply ...`
        // matches the `get_data` prefix and parses as that request.
        assert_eq!(
            Message::parse(b"get_data_reply 5 hello").unwrap(),
            Message::GetData {
                state: "_reply".into(),
                client_id: Id::new(5),
                key: "hello".into(),
            }
        );
        assert_eq!(Message::parse(b"ok_reply").unwrap(), Message::Ok,);
    }
}
