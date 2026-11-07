# `qp_tunnel` — dual-run gate recipe (reference vs Rust)

**Status (2026-11-05):** reworked to the plan — the CLI is a bin-local
statement-by-statement port of the reference's `Command_line.ml` +
`Getopt.parse_cmdline` (opam `getopt` 20120615) in `src/cli.rs`, replacing
the earlier clap-based parse. All 51 tests are green and the gates below
pass byte-identical. The three clap deviations are closed:

* the help usage line now **re-wraps with `argv[0]`** (the reference's
  `Format` `@[<hov 4>` box; measured rule below);
* unknown / misspelled options now print the reference's
  `Getopt.Error("…")` one-liner;
* the port of a `proto:host:port` address is now the full OCaml
  `int_of_string` (`0x`/`0o`/`0b`, `_` separators, a leading sign, the
  `[-2^62, 2^62-1]` bound) and `Remote.port` is an `i64` — the reference
  renders negative ports with their sign in every endpoint.

Reference binary: `$HOME/qp2/ocaml/qp_tunnel.native`.
Rust binary: `rust/target/debug/qp_tunnel` after
`cd rust && cargo build -j4 -p qp_tunnel`.
Gate drivers and captures live in `$HOME/OpenCode/tmp/w3-p4-tools/tunnel/`.

## 1. CLI matrix (byte comparison)

Run both binaries through `./qp_tunnel` symlinks so `argv[0]` is identical
(`./qp_tunnel` → each help block is 483 B there; the in-repo capture
`tests/data/help_qp_tunnel.txt` is the 481-B block for the 9-char bare
name). All cases exit 2:

| argv (after `./qp_tunnel`) | stderr |
| --- | --- |
| *(none)*, `a b`, `-- ezfio extra` | help once (483 B) + `Failure("EZFIO_FILE or ADDRESS is missing")` |
| `-h` | help **twice** (966 B) + the missing-argument failure |
| `-h notanaddress`, `-h ezfio` | help once + the address / missing failure |
| `notanaddress`, `tcp:1.2.3.4`, `a:b:c:d` | `Failure("<arg> : Malformed address")` |
| `tcp:1.2.3.4:notaport`, `tcp:1.2.3.4:4611686018427387904` | `Failure("int_of_string")` |
| `tcp:127.0.0.1:0x10`, `tcp:127.0.0.1:-5` | dead-remote shape with decimal endpoints (`16…`, `-5…`) |
| `-z ezfio`, `-gz ezfio`, `-g1 ezfio`, `-hx ezfio` | `Getopt.Error("Unknown option -z / -z / -1 / -x")` |
| `--nosuch ezfio`, `--nosuch=1 ezfio` | `Getopt.Error("Unknown option --nosuch")` |
| `--h ezfio` | `Getopt.Error("Unknown option --h")` |
| `--get-input=1 [ezfio]`, `-h=x ezfio`, `-g=x ezfio` | `… does not accept argument` / `Getopt.Error("Unknown option -=")` |
| `-gh ezfio`, `--get-input ezfio` | help once + `Failure("ezfio : Malformed address")` (parsing continues) |
| `-- ezfio` | `--` is consumed; `ezfio` is anonymous → malformed address |

`QP_ROOT` is checked **before** parsing (`env -u QP_ROOT` fails even for
`-h`); `QP_ROOT=` (empty) is accepted and reaches the parser.

The same matrix is pinned in `tests/cli.rs` (18 e2e cases, spawned with
`CommandExt::arg0`); the dead-remote cases compare the four `Thread …`
lines as a set (thread-scheduling race in both implementations) plus the
lowest-id fatal line.

## 2. The usage-line re-wrap rule

The usage line is the only `argv[0]`-dependent part of the help. It is the
reference's `Format` `@[<hov 4>` box, and the wrap behaviour is **not**
plain greedy packing. Derived from argv0-length sweeps against renamed
copies of the reference binary (tunnel lengths 9–80, `qp_set_mo_class`
15–89; both reproduce the same rule, and the port now matches the tunnel
reference byte-for-byte at all 72 lengths):

* line 1 starts at column 2, continuation lines are indented 6;
* a word following a line of content length `L` (separator space at column
  `s = L+1`, ending at `e = s + w`) is kept on the line only when
  **both** `s ≤ 68` and `e ≤ 77`;
* when the word wraps, the separator space stays as a **trailing space** on
  the old line only when `e ≤ 77` (the reference commits the space when
  the word would fit the 77-column margin, then still wraps a short word
  that starts past column 69 — e.g. a 4-char option word ending at column
  73 wraps with a trailing space).

`wrap_usage` in `src/cli.rs` implements exactly this and is pinned against
three bare-invocation captures of renamed reference copies
(`tests/data/help_qp_tunnel.txt` 481 B / `help_medium.txt` 499 B /
`help_long.txt` 540 B). The re-wrap sweep recipe (port vs reference,
identical `argv[0]` per side) is the one used on 2026-11-05 in
`$HOME/OpenCode/tmp/w3-p4-tools/tunnel/sweep2/`.

## 3. Port parsing

`parse_address` splits `proto:host:port` on `:` (exactly three fields) and
parses the port with the bin-local `src/ocaml_num.rs`
(OCaml `int_of_string`). Measured reference behaviour, all pinned:

* `0x10` → 16, `-5` → -5, `+1` → 1, `1_000` → 1000;
* above `2^62 - 1` → `Failure("int_of_string")`, exit 2;
* any in-range port (even 70000, 4294967290) is **accepted**: the
  reference renders the decimal endpoint (`tcp:127.0.0.1:70000`), fails to
  connect in the four proxy threads, and then **hangs** in `zmq_ctx_term`
  (documented deviation below).

## 4. Dead remote (`tcp:127.0.0.1:1`)

The reference prints the four `Thread … killed on uncaught exception
Failure("Unable to establish connection to …")` lines and then hangs in
`zmq_ctx_term` (`EINTR`). The port prints the same four lines (compared as
a set), the lowest-id failure as `Fatal error: exception …`, and exits 2.

## 5. Live banner and relay round trip

`gate_live.py` (fake input REP answering `test` → `OK`) and
`gate_relay.py` (all four channels: REQ echo, pub, input protocol incl.
`get_input`/`get_ezfio_filename`) drive one binary each. **Run them
sequentially** with different remote bases — both tunnels want the local
42379 window. Compare after normalising the printed `binary=` line and the
remote base (perl `\b`-based substitution, not BSD `sed`). Both transcripts
were byte-identical on 2026-11-05 (live: banner, `Communication [ OK ]`,
ssh block, `Ready`, clean exit 0; relay: the full message matrix in the
driver output, `qp_tunnel exited properly.`).

## 6. `-g`

Covered end to end by `tests/get_input.rs` (fake REP serving
`tests/data/no_work.tar.gz`, exact stdout, unpacked file,
`work/qp_run_address`, clean signal exit). The reference oracle is
`getB.out` with a tarball that already contains `testdir/work/`; the
reference **fails** on the `no_work` tarball with `Unix.EEXIST "mkdir"`
while the port creates `work/` (documented deviation, do not gate that
path against the reference).

## Known, deliberate deviations

* **Dead remote:** the reference hangs in `zmq_ctx_term`; the port exits 2
  (the four thread lines are identical).
* **Non-UTF-8 argv:** the reference walks raw bytes; the port decodes with
  `to_string_lossy` (inherited from T4.1).
* **`-g` on a tarball without `work/`:** the reference's `mkdir` fails
  `EEXIST`; the port creates the directory.

## Gotchas that cost time

* **The wrap rule is not greedy.** A plain "pack to column 77/78" model
  fails at exactly the trailing-space configurations (a 4-char word
  ending at column 73 wraps; a 12-char word ending at 77 fits). The two
  constants are `s ≤ 68` (commit) and `e ≤ 77` (margin + trailing space);
  they were pinned at their boundaries by the sweeps. `qp_set_mo_class`'s
  `wrap_usage` has the same fixed rule (its 75-length sweep is also
  byte-identical).
* **`-h` prints the help once and parsing continues** (`set_specs`); the
  driver prints it a second time when the anonymous argument is missing,
  so `-h` alone shows two blocks.
* **`--` is consumed**, not an anonymous argument; a single `-` and an
  empty string **are** anonymous arguments.
* **`timeout` does not exist on macOS zsh** — use the python drivers or a
  background + kill for the hang cases.
* All scratch lives under `$HOME/OpenCode/tmp/w3-p4-tools/tunnel/`,
  never `/tmp`.

## Out of scope

Multi-node ssh (two clusters) is the P5 Ampere cutover gate and needs the
user's permission; nothing here opens a real ssh tunnel.
