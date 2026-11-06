# Known deviations from the OCaml tool layer

The Rust tools are ports, not rewrites: the EZFIO on-disk format, the ZMQ
wire protocol, the help text of `qp_set_mo_class`/`qp_tunnel`, the
`qp_edit -c` contract (C1–C7), the `qp_create_ezfio` output on the six
corpus fixtures and the fresh-tree `qp_edit -c` file set (275 files vs 275)
are byte-identical to the reference. This file lists **everything that is
not**, one row per difference a user can observe, with the reason and the
pointer to the detailed note. Anything not listed here is a bug.

Reference: `~/qp2` at `be4b71c1` (the read-only oracle). "Ref" below means
that build, run side by side with the port.

## A. Command line and process behaviour

| # | Tool(s) | Ref | Port | Why | Impact |
|---|---|---|---|---|---|
| A1 | `qp_create_ezfio` | takes `argv` as bytes | rejects the argument: `<tool>: command line argument <n> is not valid UTF-8` (argument 0 included) | `String::from_utf8` instead of a byte parser; chose an error over a panic or a silent loss | only non-UTF-8 command lines, which no QP2 caller produces |
| A2 | `qp_edit`, `qp_run`, `qp_set_mo_class`, `qp_print_basis`, `qp_tunnel` | takes `argv` as bytes | `String::from_utf8_lossy` (U+FFFD for each invalid byte) | same | a non-UTF-8 option value reaches the tool with replacement characters |
| A3 | `qp_create_ezfio` | `Unknown option -<byte>` prints the first **raw byte** of an unknown short option | `Unknown option -<char>` prints the decoded **character** | UTF-8 parser | unreachable through a typo, since no `set_specs` entry is non-ASCII |
| A4 | `qp_create_ezfio`, `qp_edit`, `qp_run` | `Fatal error: exception Getopt.Error("Unknown option --nope")` | clap diagnostics (`error: unexpected argument '--nope' found`) | clap parses these three; the other two use a `getopt` port | text differs, **exit code 2 identical** |
| A5 | `qp_print_basis` | dies of `SIGPIPE` (exit 141, empty stderr) when stdout closes | ignores `SIGPIPE`, reports `Fatal error: exception Sys_error("Broken pipe (os error 32)")`, exit 2 | Rust starts with `SIGPIPE` ignored | `... \| head -c 1` only; file/terminal stdout is identical |
| A6 | `qp_run` | banner prints the git stamp baked into `Git.ml` at build time | prints `unknown` unless `QP_GIT_SHA1`/`QP_GIT_DATE`/`QP_GIT_MESSAGE` are exported | baking the commit into a binary makes the build non-reproducible and a copied binary lie | banner text only |
| A7 | `qp_run` | `Sys.command` → `system(3)`: a shell killed by a signal is 255 | `/bin/sh -c` with inherited streams: a shell killed by a signal is 128+n | `std::process::Command` | `qp_run` only tests `status = 0`; invisible to every job flow |
| A8 | `qp_run` (`qp_stop`) | signal handler runs immediately | `SIGUSR1` polled: `qp_stop` starts 0–200 ms later, two signals in one interval collapse, a signal after job exit is dropped | no async-signal-safe path to the ZMQ doorbell | job teardown latency ≤ 200 ms; `qp_stop` is idempotent |
| A9 | `qp_edit` | an uncaught exception prints `Assert error in file $QP_ROOT/ocaml/<file>, line <n>, character <c>` | prints `Assert error: <message>` | Rust cannot recover an OCaml source location | banner text only; **exit code 2 preserved** |
| A10 | `qp_edit` | sexplib `Error:  (...)` prints the full type path | two-space prefix, value snippet and exit status reproduced; type-path wording approximated | no sexplib in Rust | only malformed-input messages (C7-A) |

## B. EZFIO tree

| # | Tool(s) | Ref | Port | Why | Impact |
|---|---|---|---|---|---|
| B1 | all EZFIO users | opening a directory that is not EZFIO mostly no-ops (`Sys.file_exists` only) | `<dir>/.version` must be a regular file, else `not an EZFIO directory` | fail loudly instead of silently ignoring a wrong path | a typo'd directory is an error, not a no-op |
| B2 | all EZFIO writers | `library` stamp = the configure-time ezfio library path | `library` = `qp-ezfio` | a static pure-Rust build has no such path; nothing reads the file back | machine stamp, stripped as VOLATILE by the golden harness |
| B3 | `qp_create_ezfio` | a create failure other than a missing parent is reported as the `Unix.Unix_error` it is (e.g. `ENOTDIR`) | always the `ENOENT` message | single pre-check of the parent | wrong `Unix_error` name on exotic failures |
| B4 | `qp_print_basis` | `Input_ao_basis.read` writes `ao_basis/ao_md5` and `ao_cartesian` when absent | reads them, never writes them | the print path does not open the section for writing | only the opened tree (a missing file); stdout/stderr/exit identical |
| B5 | `qp_edit -c` | `at_exit remove_backup` deletes `work/backup.tar` even in check mode | check mode creates no backup and deletes nothing | "check writes nothing" is the honest reading | a stray pre-existing `work/backup.tar` survives |
| B6 | `qp_edit` | backup written as `./.backup.tar` then `mv`-ed into `work/`; a failed `mv` leaves the temp file in the cwd | tar built in a temp file outside the tree and copied in | avoid littering the caller's cwd | no observable difference on success |
| B7 | all EZFIO writers | `work/empty` is `F\n\n` (3 bytes) | `F\n\n` since `qp-ezfio` `d5439a70` | was 2 bytes on the branch where it was measured | none (same value) |
| B8 | `qp_edit` (generated readers) | `Qpackage.get_ezfio_default` scans `$QP_ROOT/data/ezfio_defaults/`; when no file carries the `(section, field)` pair it prints `<section>/<field> not found` and raises `Not_found`, which `get` (catching only `Sys_error`) lets escape → `Fatal error: exception Not_found`, exit 2 (measured on the `be4b71c1` build); a missing directory instead raises `Sys_error` → `Info: …` per section and a silent skip (rc 0) | `defaults.rs` returns `None` and `materialise_default` falls back to the TOML `default:` in both shapes | cutover choice (plan §9 F1, user decision 2026-11-05): a data tree missing one of its 55 files must not abort the editor or silently drop a section, and the TOML values are the ones the OCaml configs carry | a 54/55 (or directory-missing) tree: the reference fatals (rc 2) or skips with `Info:` lines, the port materialises the section — at the 2026-11-05 gate the Ampere reference (`3b10f0a8`, an ancestor of the oracle, predating the symmetry section: no `Input_symmetry.ml`, no `Symmetry` task) wrote none, while the port materialised `symmetry/auto_symmetry = T` from the built-in fallback; on a 55/55 tree the two agree; pinned by `c3_c_a_missing_defaults_file_falls_back_to_the_built_in_default` |
| B9 | `qp_edit -c` | the pre-fix driver applies every section in check mode (`set` has no `apply` flag); the determinants writer re-runs `write ~force:false` and rewrites `determinants/psi_{det,coef}_qp_edit.gz` from `psi_det`/`psi_coef` whenever `n_det <= n_det_qp_edit` | check mode applies nothing: an existing pair keeps its bytes, and is read (never written) only when `n_det_qp_edit != n_det` | C2 ("no rewrite") is the normative post-fix behaviour; the pair is a scratch copy that the next interactive write regenerates | a pair that disagrees with the working arrays stays stale where the reference refreshed it; no corpus fixture has a pair; pinned by `c4_b`/`c4_c`/`c4_d` in `bins/qp_edit/tests/check_mode.rs` |

## C. `qp_edit` document text (`-w`, `$EDITOR`)

The document is user-facing text; the parse-back ignores descriptions and
keeps full precision, so none of these rows changes an EZFIO value.

| # | Item | Ref | Port | Why | Impact |
|---|---|---|---|---|---|
| C1 | nuclei block | `%16.8f %16.8f %16.8f` Angstrom, padded | `%.17g` Angstrom, single-spaced | contract C5: `%16.8f` is a 1.9e-8 bohr grid that does not round-trip | `-w` bytes differ; a `+5e-9` bohr perturbation now survives `-c` |
| C2 | descriptions containing an RST-escaped vertical bar | the OCaml unescape pass drops the backslash before emitting `Input_*.ml`, so the document shows a bare vertical bar | prints the two characters backslash + bar, as carried by the `EZFIO.cfg` doc string | the P2 config emitter does not run the OCaml unescape pass that the `Input_*.ml` emitter does | document text only; the parse-back drops descriptions; the fix belongs to the emitter |
| C3 | section tail | every section block ends with a line of three spaces | no whitespace-only line | the port trims the block separator | `-w` is 172 bytes shorter on the fresh tree; trailing whitespace only |

## D. TaskServer and the wire

| # | Item | Ref | Port | Why | Impact |
|---|---|---|---|---|---|
| D1 | ZMQ context | `set_io_threads 16` | default (1 IO thread) | `zmq` 0.9.2 exposes no `ctx_set`; the binding would be an undeclared transitive dependency | throughput only; revisit before the Ampere multi-node gate |
| D2 | non-UTF-8 message tokens | OCaml strings are bytes, so accepted | `ParseError::NonUtf8`, message rejected | `String` tokens | every QP2 producer emits ASCII/UTF-8; wire is byte-identical for real traffic |
| D3 | dead remote (`qp_tunnel`) | hangs in `zmq_ctx_term` | exits 2 | `zmq_ctx_term` with a dead peer; the port closes and reports | an unavailable peer fails instead of hanging; the four thread lines are identical |
| D4 | `qp_tunnel -g` on a tarball without `work/` | `Unix.EEXIST "mkdir"` | creates `work/` | `create_dir_all` | an unpack failure becomes a success |

## E. Build and packaging

| # | Item | Ref | Port | Why |
|---|---|---|---|---|
| E1 | toolchain | OCaml >= 5.0 + opam (never built in CI: the in-tree opampack ships 4.12) | Rust 1.92.0, pinned by `rust/rust-toolchain.toml` | the port |
| E2 | published names | `ocaml/*.native` and `element_create_db.byte` are compiled binaries | the same names are **symlinks** into `rust/target/release/` | `etc/paths.rc`, `bin/qp_export_as_tgz` and the man pages keep working unchanged |
| E3 | `external/opampack` and the opampack tarball | downloaded and extracted by `configure -i ocaml` | `configure` no longer references either | the packages are gone from the build; the tarball path (`external/qp2-dependencies/<arch>/opampack.tar.gz`) lives inside the unchecked-out `external/qp2-dependencies` submodule, so it could not be deleted in this tree |
| E4 | man pages | generated by Sphinx from `docs/source/users_guide/qp_*.rst` | not regenerated in P5 | the `.rst` sources are current (only the `qp_tunnel` synopsis omits `-h`/`[--]`, while its option list documents `-h`); the checked-in `man/*.1` files were generated in 2019 and are stale for `qp_create_ezfio`, `qp_edit`, `qp_tunnel`, `qp_plugins`, `qp_stop` and `qp_set_frozen_core` — a pre-existing docs-build drift, not a behaviour change of the port |

## Resolved during P5 (kept for the record)

* **Runtime defaults back-fill.** `qp_edit -c` on a fresh `qp_create_ezfio`
  tree created 8 files the port did not (the whole `qmcchem` section, the
  three `ao_extra_basis` logicals and its MD5). The generator's readers take
  those values from `$QP_ROOT/data/ezfio_defaults/*.ezfio_interface_default`
  at run time (`Qpackage.get_ezfio_default`); the port now does the same
  (`qp-input/src/defaults.rs`, an exact port of the scan) and visits the
  `qmcchem` module. Fresh tree: 275 files vs 275, 273 common files
  byte-identical, 2 machine stamps skipped. Details:
  `tests/golden/contracts/check-mode.md` (C3).
* **Dead clap structs.** `qp_types::command_line::{SetMoClassArgs,
  SetMoClassConfig, TunnelArgs, TunnelConfig}` duplicated the bin-local
  `getopt` ports and no binary used them; removed in the P5 audit.
* **`qp_edit` configuration after a clean build.** Rewriting
  `ocaml/Makefile` dropped the `Input_auto_generated.ml` rule, which was
  also what wrote `data/qp_edit_config/registry.toml`: a fresh `ninja`
  produced the per-module tables but no registry, and `qp_edit` failed to
  start. The Makefile now regenerates the registry with the same
  `ei_handler.py ocaml_global` call (and ignores the two obsolete OCaml
  files it still emits). `etc/paths.rc` also exports `QP_EDIT_CONFIG`,
  because the binary's `<exe_dir>/../data` fallback holds on macOS, where
  `current_exe` keeps the `ocaml/qp_edit` symlink path, but not on Linux,
  where `/proc/self/exe` resolves it into `rust/target/release/`.

## Not deviations (gated)

* EZFIO payloads: `create` on h2, h3, c2h2, b2, h2g (47/47 files each) and
  the committed hcn fixture; `qp_edit -c` on the six corpus fixtures
  (C1–C7) and on hcn (293 common files identical, 2 machine stamps
  skipped, plus the C4-only `*_qp_edit` pair).
* Help text of `qp_set_mo_class` (75 argv0 lengths) and `qp_tunnel`
  (72 lengths), including the `Format @<[<hov 4>` re-wrap rule, and the
  `-h` double print when the anonymous argument is missing.
* The ZMQ wire format (18 message kinds) against a live server and relay.
* `data/executables` and the `.native` names; the `-h` text of the five
  Rust tools matches the reference (argv[0] and line-wrap width aside).
  The 12 `qp_*.1` man pages are discussed in E4.

## Where the details live

`CONTRACTS.md` (§4 EZFIO, §5 check mode, §6 driver), `bins/*/NOTES.md`,
`crates/*/NOTES.md` (`qp-cli`, `qp-input`, `qp-zmq`, `qp-types`), the
fixture READMEs (`crates/qp-ezfio/tests/fixtures/README.md`) and
`tests/golden/contracts/check-mode.md`.
