# P4 tools — dual-run gate recipe (`qp_set_mo_class`)

**Status (2026-11-09, wave S2):** the help/error byte contract documented
below is **retired** — user-visible text is now Rust-native (user decision
2026-11-09; the hand-rolled `getopt` port in `src/cli.rs` — including
`wrap_usage` — is deleted, the CLI is clap derive). The contract that
survives is the exit codes (0/1/2), the EZFIO payloads and the data gates
(DEVIATIONS.md A4/A11; `steps/012-simplify-tool-layer/implementation.md`
S2). The byte-comparison recipe below is kept as historical evidence of
the **reference** behaviour; claims in it about the port matching the
reference byte-for-byte are pre-S2 and no longer hold.

**Status (2026-11-05, historical):** the port is implemented and the gate
below is green —
71 compared cases byte-identical on `(exit code, stdout, stderr)`, plus the
decompressed `mo_basis/mo_class` payload and the materialised `mo_label`
where relevant. This file was the recipe the parent re-runs after merging;
the
gate script is self-contained and uses only files in the repository (the
fixtures) plus a fake `QP_ROOT` the script creates in its own scratch.

```
GATE: 71 passed, 0 failed
```

## What is compared

For every case, the reference
(`$HOME/qp2/ocaml/qp_set_mo_class.native`) and the port are each run in their
own copy of the same fixture, and the gate compares:

* the **exit code**;
* the **stdout bytes** (with the one exception below);
* the **stderr bytes**;
* when `mo_basis/mo_class.gz` exists in both resulting trees, the
  **decompressed** payload bytes — the `.gz` headers themselves are never
  compared (gzip implementations and timestamps differ);
* when the case starts without `mo_class` or with a custom `mo_label`, the
  **materialised** `mo_label` plain file.

The fixture is mutated **identically for both sides** via `MOD_*` environment
prefixes on the `gcase` call (e.g. `MOD_MO_NUM='9000' gcase q-mo9000 hcn -q
ezfio` rewrites `mo_basis/mo_num`; `MOD_BK=…` rewrites the fake root's
`bitmasks_module.f90`; `MOD_RM_CLASS` removes `mo_class.gz`).

**argv[0] shim (help cases only).** The reference help embeds `argv[0]` in an
OCaml `Format` `hov` box and re-wraps the usage line with the actual name
(measured rule: a word is kept on a line only when its separator space is at
column ≤ 68 **and** it ends at column ≤ 77; on a wrap the space stays as a
trailing space only when the word would have fit the 77 margin — see the
Gotcha below; pre-S2 the port reproduced this in `cli.rs::wrap_usage`,
since deleted — the help is clap-generated and not byte-compared). The four help cases
(`cli-noargs`, `cli-h`, `cli-two-anon`, `cli-h-ezfio`) therefore cannot be
compared with each side invoked by its own path: the gate's `gcase_help`
execs the **port** through
`python3 -c 'import os,sys; os.execv(sys.argv[1],[sys.argv[2]]+sys.argv[3:])'
"$PORT" "$REF" <args>` so both sides see the reference's exact `argv[0]`.
All other cases run each binary directly.

## Cases (71)

| group | cases | reference behaviour pinned |
| --- | --- | --- |
| CLI (18) | `cli-noargs`, `cli-h`, `cli-two-anon`, `cli-h-ezfio` | help block (argv[0]-dependent wrap), exit 2 (`-h` twice on no args); `-h ezfio` prints help *and* continues the run |
| | `cli-unknown-s/-l`, `cli-h-unknown`, `cli-missing-v1/v2/v3`, `cli-empty-val`, `cli-concat-arg/-unk`, `cli-help-eq`, `cli-dashdash`, `cli-short-eq`, `cli-immediate`, `cli-long-eq` | the measured `Getopt` quirks: `--` consumed, concat rules, `-c=1` keeps `=1` literally, `-c ""` dies in `Range.parse` before `check`, … |
| ranges (14) | `range_*`, `range-notspec1/2`, `range-overmax`, `range-dup1/2` | every error class from `Range.ml` (lines 26/45/46) and the over-max / duplicate-orbital checks |
| set (7) | `set-partition`, `set-current`, `set-del`, `set-mo-num-1a/1b`, `set-mo-num-0`, `set-mo-num-neg`, `set-mo-num-abc`, `set-mo-num-10001` | partition writes; `mo_num` edge values: `1` skips the max check, `0`/negative reach the asserts / `Array.init("Array.init")`, `abc` → `int_of_string`, `10001` → `Invalid_argument("MO_number: 10001")` |
| set_file (4) | `sf-missing`, `sf-file`, `sf-parent`, `sf-no-mo-num` | `ezfio.ml:426` shapes: `nosuch` (missing) is created and then fails at `Failure("mo_basis/mo_num not found")` — the tree is left on disk; a regular file named `ezfio` → `Unix_error(EEXIST, "mkdir", "ezfio")`; `ezfio/sub` → `Unix_error(ENOENT, "mkdir", "ezfio/sub")` because `mkdir` is single-level and the parent does not exist (nothing created); h2 tree without `mo_num` → `Failure("mo_basis/mo_num not found")` |
| `-q` (24) | `q-hcn`, `q-nint-3/200/0/neg/abc`, `q-mo9000`, `q-mo10001`, `q-alpha-miss/empty/abc/0`, `q-beta-miss/neg/0`, `q-mo0`, `q-mo-abc`, `q-label-bogus`, `q-noclass`, `q-bk16`, `q-bk42`, `q-bk-missing`, `q-bk-noline` | the full `n_int` resolution: file value → `Bitlist.n_int_of_mo_num` fallback → `max` re-read; the double `MO_number` warning at 10001; `bit_kind_size` parser accept/reject; materialised `mo_class` payload + `mo_label` compared |
| QP_ROOT unset (3) | `qproot-bare`, `qproot-q`, `qproot-dashdash` | the module-init guard fires before anything else, 0 B stdout |

The authoritative case list is in the script itself (`gcase`/`u` calls); the
table above is the family-by-family summary.

## Prerequisites

* `cargo build -j4 -p qp_set_mo_class` in the worktree's `rust/`
  (`target/debug/qp_set_mo_class` is the port under test).
* `$HOME/qp2/ocaml/qp_set_mo_class.native` (the reference; running it is fine,
  it is never written).
* `rust/tests/golden/fixtures/hcn` and `h2` (repository files; the script
  copies them into its scratch and never writes into the fixture directory).
* The script creates its own fake `QP_ROOT` in its scratch
  (`gate/qroot`, with `src/bitmask/bitmasks_module.f90` carrying
  `bit_kind_size = 64`); no real QP install is needed.
* Scratch: `$HOME/OpenCode/tmp/w3-p4-tools/set-mo-class/` (run with **bash**).

## Known, deliberate deviations (reported, not failed)

* **`.version` stamp text.** The reference stamps
  `library = qp_set_mo_class`; the port stamps `library = qp-ezfio` (the
  shared `qp_ezfio` writer contract, documented at T4.1). The gate does not
  compare the stamp files.
* **Non-UTF-8 argv.** The reference (OCaml) accepts arbitrary bytes in
  `argv`; the port (clap, since wave S2) rejects a non-UTF-8 argument
  with its own `error: invalid UTF-8 was detected in one or more
  arguments` diagnostic, exit 2 (DEVIATIONS A1/A2). The pre-S2
  `String::from_utf8_lossy` behaviour is gone; not exercised by the
  gate.
* **`mo_class` shape validation.** The Rust `qp-ezfio` `StaticRegistry`
  rejects an array whose on-disk dimensions differ from the declared spec
  (defence for the array-heavy `qp_create_ezfio`). The reference client has
  no such table: `read_array` checks only per-line syntax. The port therefore
  opens the tree with a bin-local `NoRegistry` (empty `Registry` impl in
  `main.rs`) so an inconsistent tree — e.g. `mo_num` changed to 9000 over the
  55-MO `hcn` arrays — reads and writes exactly like the reference (measured:
  both proceed to the `n_int` fallback). See `NoRegistry` in `main.rs`.

## Gotchas that cost time

* **The usage line re-wraps with `argv[0]`, and the rule is not plain
  greedy.** The help is not a static substitution: `help()`
  (`Command_line.ml`) puts the usage in a `Format @<[<hov 4>` box, so the
  first usage line starts at column 2, continuation lines at column 6, and
  a word following a line of length `L` (separator space at column
  `s = L+1`, ending at `e = s + w`) is kept on the line only when **both**
  `s ≤ 68` and `e ≤ 77`; when it wraps, the separator space stays as a
  trailing space on the old line only when `e ≤ 77` (a 4-char option word
  ending at column 73 wraps *with* a trailing space, while a 12-char word
  ending at column 77 fits — a single greedy margin cannot produce both).
  Derived from an argv0-length sweep (15–89) against renamed copies of the
  reference binary; the pre-S2 port matched the reference byte-for-byte at
  all 75 lengths. It implemented the rule (`wrap_usage` in `cli.rs`, since
  deleted in wave S2.1) and
  pinned it with three captures taken from *renamed copies* of the reference
  binary (name lengths 15/31/65): `tests/data/help_bare.txt` (639 B),
  `help_medium.txt` (657 B), `help_long.txt` (697 B) (since deleted).
* **`has` checks the file, not just the group.** `ezfio.ml:54-66`:
  `has g k` = `exists(g)` **and** `Sys.file_exists(g/k)`. A group whose
  `.version` exists but the file does not is `Absent` (the `n_int` fallback
  then applies). An earlier port draft used the group-only check and failed
  the `q-nint-*` family.
* **The 10000-MO warning prints twice** for `mo_num = 10001` under `-q`
  (once per `MO_number.of_int` conversion inside `Mo_basis.read`), followed
  by `N_int > 128` and `Invalid_argument("N_int_number: 157")`.
* **`Array.init` on a negative `mo_num`** raises
  `Invalid_argument("Array.init")` (measured) — not "negative dimension".
* **The `n_int` fallback `max` check re-reads the file** when
  `determinants/n_int` is absent, so garbage in the file is re-raised
  uncaught; `mo_num = 1` skips the maximum check entirely.
* All scratch lives under `$HOME/OpenCode/tmp/w3-p4-tools/set-mo-class/`,
  never `/tmp`.

## The script

Copy into scratch (`$HOME/OpenCode/tmp/w3-p4-tools/set-mo-class/gate.sh`)
and run with `bash`. Key variables:

```bash
GATE="$HOME/OpenCode/tmp/w3-p4-tools/set-mo-class/gate"
REF="$HOME/qp2/ocaml/qp_set_mo_class.native"
PORT="<port-checkout>/rust/target/debug/qp_set_mo_class"
FIX="<port-checkout>/rust/tests/golden/fixtures"
```

The script defines `setup_fixture` (copies `hcn`/`h2` and applies the
`MOD_*` mutations identically on both sides, including the shared fake-root
bitmask file, restored after each side), `gcase` (direct invocation, stdout
compared after normalising the argv[0] path string to `PROG` for the rare
non-help case that prints a path), `gcase_help` (the `os.execv` shim above),
and `u` (QP_ROOT-unset variants via `env -u QP_ROOT`). Final line:
`GATE: 71 passed, 0 failed`.

The full script as run on 2026-11-05 is preserved in
`$HOME/OpenCode/tmp/w3-p4-tools/set-mo-class/gate.sh` (scratch, not in the
repository; this section + the case table above are the durable record).
