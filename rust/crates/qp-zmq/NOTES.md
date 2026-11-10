# qp-zmq — notes for the next agent (P4, wave 3)

Status of this crate, 2026-11-05:

| module | state |
|---|---|
| `message` | **done**, byte-exact codec for the 18 request kinds + 8 reply kinds, 7 tests |
| `queuing_system` | **done**, port of `ocaml/Queuing_system.ml`, 7 tests |
| `progress_bar` | **done**, port of `ocaml/Progress_bar.ml`, 7 tests |
| `port_window` | **done** (pure part only), candidate window + `qp_stop` string, 4 tests |
| `task_server` | **done**, port of `ocaml/TaskServer.ml` minus `run`'s bind scan, 19 tests |

44 tests total. `lib.rs` declares all five. What is *not* here — by plan — is
the `qp_run` binary: the bind scan over `port_window::candidate_ports`, the
`QP_NIC` call site, the SIGUSR1 handler and the slave mode belong to
`qp-cli`/`bins/qp_run` in the successor wave. `task_server::run` takes
already-bound sockets for that reason.

Build and test what exists:

```
cd rust && cargo test -p qp-zmq -j4
cd rust && cargo clippy -p qp-zmq --all-targets -j4 -- -D warnings
cd rust && cargo fmt -p qp-zmq -- --check
```

2026-11-05 result: 14 tests pass, clippy clean, fmt clean.

## Reference facts verified against `~/qp2` (do not re-derive)

- Port window comes from the pid, **there is no `$QP_PORT`**:
  base `10 * (pid mod 2823) + 32769`, then `+100` per retry, forever
  (`ocaml/qp_run.ml:39-77`). The base window is P..P+9: the real REP binds P,
  the real PUB binds P+1, throwaway REP probes bind P+2..P+9. All are
  `tcp://*:N`. The assumption is that the Fortran collector binds PULL on P+2.
- `Zmq.Socket.set_linger_period`: REP 1_000_000 µs, PUB 1_000 µs.
- Lexer: a `text` is a maximal `|`-free run of length >= 2 starting with a
  non-space; a `word` is a run of non-space/tab/LF; `get_int` lexes
  `['0'-'9']+` (a leading `-` is dropped). Task id lists skip every token that
  is not a `text`, which is why `add_task st a|b c` parses as tasks `["b c"]`
  and `add_task s a|b` as an empty list (verified with the OCaml harness in
  `$HOME/OpenCode/tmp/qp-msg-golden/golden2`).
- `Id.of_int` asserts `x > 0` (`Id.ml:5`), failure text `Id.ml:5:4`;
  `int_of_string` rejects values above `2^62-1`.
- Queue order: `queued_front` is a newest-first stack, `queued_back` is
  next-first; `pop_task` pops the **oldest** task (refills `back` with
  `rev front` when empty). `to_string` prints `front @ rev back`, i.e. newest
  first, and sorts the maps by id like `Map.bindings`.
- `number_of_queued` / `number_of_running` / `number_of_clients` are stored
  counters in the reference; here they are derived from the collections, which
  is equivalent for every constructor in `Queuing_system.ml` (there is no
  decrement path that skips the corresponding removal).

## Deviations from the reference

1. **`set_io_threads 16` is dropped.** `TaskServer.ml` sets the context to 16
   IO threads; the `zmq` 0.9.2 crate exposes no `ctx_set` at all (the binding
   would be `zmq_sys`, an undeclared transitive dependency). Adding a crate is
   out of scope for this wave. Consequence: default libzmq IO thread count (1).
2. **Non-UTF-8 tokens are rejected.** OCaml strings are bytes, so the
   reference accepts a non-UTF-8 `text`/`word`; here `Message::parse` returns
   `ParseError::NonUtf8`. Every producer in QP2 emits UTF-8 or ASCII, so the
   wire stays byte-identical for real traffic. `GetDataReply` values are
   `Vec<u8>` and are never lossily converted on encode; `encode_text` is
   documented as lossy and must only be used for parsed messages (the server's
   `Invalid message : ...` fallback).
3. `Message::encode_text` is public (the task server needs it); it was
   `pub(crate)` while the server was missing.
4. Reference `tests/test_message.ml:19` (`disconnect 3 mystate`) is stale and
   fails against the current `Message.ml`; do not use it as a golden.

## Reference and plan documents that are wrong

- `rust/CONTRACTS.md:355` gives `connect_reply <state> <client_id> <ok|error>`.
  The wire format is `connect_reply <state> <client_id> <push_address>`
  (`Message.ml:101`), e.g. `connect_reply scf 3 tcp://host:43211`.
- `rust/CONTRACTS.md:361` gives `del_task_reply <state> <client_id>`. The wire
  format is `del_task_reply <more|done> <id>|<id>` (`Message.ml:222`), e.g.
  `del_task_reply more 1|2`.
- The P4 plan sentence that `qp_run` uses `$QP_PORT` (default 43210) is wrong;
  see the port window above. `lib.rs` was corrected.

## `work/empty` — answered empirically (was an open question)

`qp_edit -c` / `EzFio` create `work/empty`. Measured on a copy of
`~/qp2/tmp/hooktest.ezfio` (the `work/` directory removed first):

```
QP_ROOT=~/qp2 PATH=~/qp2/bin:$PATH ~/qp2/ocaml/_build/qp_edit.native -c $D/ezfio
wc -c $D/ezfio/work/empty      # -> 3
od -c $D/ezfio/work/empty      # -> F \n \n
```

So the reference writes `F\n\n`: `fortran_string_of_bool false` = `"F\n"` and
`print_bool` appends another `\n`. `qp-ezfio`'s `EzFio::set` writes `F\n`
(2 bytes). The 3-byte form is the reference; `qp-cli` must decide between
matching it and documenting a 2-byte deviation. Do not silently change the
fixture.

## Environment constraint to report, not to work around

The wave gate says "zero own `unsafe`", but `qp-cli` cannot port `qp_run.ml`
faithfully without it: `Unix.localtime_r` (Date line), `getifaddrs` (`QP_NIC`),
`gethostname` + IPv4 resolution, and `sigaction` (SIGUSR1 -> `qp_stop`) have no
`std` equivalent. Plan: 4 small `unsafe` blocks with `// SAFETY:` comments,
isolated in `qp-cli`. That needs a PI decision; it is not implemented yet.

## P4 wave-3 findings: `task_server`, `progress_bar`, `port_window`

### The plan is wrong about the port base

`ocaml/qp_run.ml:72` is `10 * (pid mod 2823) + 32769`; the brief said 2821.
The OCaml is normative and the code uses 2823 (`port_window::PORT_SLOTS`).
A test pins the formula so the typo cannot come back.

### `QueuingSystem::abort` is not `TaskServer.abort`

`queuing_system.rs` has an invented `abort` that pops **one** task; its test
`abort_discards_one_task_and_removes_its_client` locks that in. The reference
`TaskServer.abort` (`TaskServer.ml:695`) drains **all** queued tasks:

```
aux [] queue 1  -> pop one; then recurse on (number_of_queued new_queue)
```

then `end_task` all, `del_task` all, `del_client`, reply `ok`, and sets
`accepting_clients = false`. `task_server::abort` implements that loop
privately. Do not "simplify" it to `QueuingSystem::abort`: that would leave
`n-1` tasks queued and change the reply of the next `get_task`. `queuing_system.rs`
was otherwise left untouched; the only additions are the read-only accessors
`has_client` and `task_payload` (needed by `del_task`/`get_tasks`).

### Linger values are milliseconds, 1:1

`Zmq.Socket.set_linger_period` forwards the OCaml int straight to
`zmq_setsockopt(ZMQ_LINGER)` (`zmq.6.0.0/src/caml_zmq_stubs.c:233`,
`src/zmq.ml:297`), and libzmq documents ZMQ_LINGER in milliseconds. So
`1_000_000` is 1 000 000 ms, not microseconds. The port uses the same numbers:
REP 1 000 000, PUB 1 000, `stop`'s REQ 1 000 000 then 1 000
(`qp_run.ml:74-75`, `TaskServer.ml:151,182`). A `LingerReset` guard zeroes the
REQ linger on the failure path: the reference leaks the socket to the GC at
exit, the Rust `Drop` would otherwise block for ~1 000 s.

### Reference details that drive the tests

- Dispatch order (`TaskServer.ml:920-941`): Terminate, Abort, PutData,
  GetData, `None,Newjob`, `_,Newjob`, `Some,Endjob`, `Some,SetRunning`,
  `_,SetWaiting`, `_,SetStopped`, `None,_`, `Some,Connect/Disconnect/AddTask/
  DelTask/GetTask/GetTasks/TaskDone`, `Invalid`. Consequences: `set_waiting`
  and `set_stopped` answer `ok` even with no job, `set_running` does not.
- The PUB thread publishes the current state **before** testing the previous
  one, so `Stopped` is sent at least twice, up to a poll apart, and
  `set_stopped` kills the state channel. Keep such a command last in a test.
- The reference `run` spawns the PUB thread *before* binding REP. The port
  binds REP first. The reordering is unobservable (`new_job` is what makes the
  state channel matter) and is documented in `task_server.rs`.
- `put_data` requires exactly one extra frame (`Failure "Badly formed put_data
  message"`), checked before the state arm; `get_data` answers one NUL byte
  when the key is absent.
- `get_tasks` must not mutate the queue when the next id's payload was
  deleted: the port peeks `next_queued_id` + `task_payload` before popping.
  `Assert_failure` sites are spelled `TaskServer.ml:<line>:<col>` with
  OCaml's **0-based** column (`312:46`, `646:14`, `678:14`).
- `del_task` takes `state` and ids only (no client id). Its reply counts
  every payload still in the map, including running ones: deleting the second
  of two tasks still present answers `del_task_reply more <id>`.
- `stop` failure strings are verbatim: `TaskServer: no answer to the
  Terminate request (the task server thread is dead)` and
  `Problem in termination`. The port uses `Instant`, not `Unix.time`, so a
  wall-clock jump cannot extend or shorten the 30 s deadline.
- `bind_socket` retries **every** `zmq::Error`, not just `EADDRINUSE`: the
  `zmq` crate does not expose `errno`, so `Unix.Unix_error` cannot be told
  apart from the other exceptions. The reference's 60 x 1 s budget is kept.

### `QP_TASK_DEBUG`

The reference forces `debug_env` once at module initialisation; the port
memoizes it in a `OnceLock`. The per-request debug line uses
`message.encode_text()` where the reference uses `Message.to_string`; for any
message the server can parse (i.e. anything that reached dispatch) the two are
byte-identical, including `get_data_reply` (lossy only for payloads the
server never prints). `stop` also emits the `STOP` line.

### Bind behaviour measured on this host (macOS, system libzmq)

- `tcp://*:0` then `tcp://*/<port>`: **ok**; a wildcard bind does not exclude
  a later specific bind on macOS.
- Two wildcard binds on one port: second fails `Address already in use`.
  This is what the bind-failure tests use.
- `bind` on port 70000 **succeeds** (silently truncated to 16 bits), so it is
  not a failure probe — do not use it as one.

### Test-only gotcha

`Message_lexer.mll:60` requires a `text` to be at least two bytes, so
`add_task test a|b` yields the single task `"b"` and `a` is dropped. All
fixtures use two-character payloads (`aa|bb|cc`).

### Formatting gate

`cargo fmt -p qp-zmq -- --check` is clean. `cargo fmt --all -- --check`
fails on `rust/crates/qp-types/src/lib.rs:86` (`pub mod gto;` ordering), which
comes from the P3 wave. It was left alone: reformatting another wave's file is
not part of this scope.
