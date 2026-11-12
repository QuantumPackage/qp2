# `qp_create_ezfio` — reference semantics catalogue

**Status: the Rust port is NOT written.** `src/main.rs` is still the wave-3
scaffold (prints `not implemented yet (P4)`, exits 2). This file records what
`ocaml/qp_create_ezfio.ml` actually does, as measured on the reference build, so
the port starts from evidence rather than from a re-read.

Measurements below were taken with the reference oracle build
(`ocaml/qp_create_ezfio.native`, `QP_ROOT` = the reference checkout),
in `$HOME/OpenCode/tmp/w3-create-ezfio/probe`. Exit codes and messages are
verbatim. The golden oracle `golden_harness.py verify --only create` is green
6/6 with that build (log: `$HOME/OpenCode/tmp/w3-create-ezfio/ref-verify-create.log`).

## Gamess verdict

`qp_create_ezfio.ml` has **no Gamess input path**. `grep -i gamess` over the file
returns nothing; `Gamess` appears elsewhere only as the `Gto.Gamess` printf
format in `Gto.ml` / `Basis.ml`. Molecule sources are XYZ with a Z-matrix
fallback, nothing else. `Gamess.ml` stays unported, as the plan allows.

## Error surface (all exit 2)

Two distinct classes. The port must reproduce both, including the fact that the
second class is *not* prefixed.

Prefixed with `qp_create_ezfio: ` (these are the `Failure`/`Sys_error` paths
caught in `run`):

| trigger | message |
| --- | --- |
| `-b` absent | `--basis option is missing.` |
| no `FILE` | `input file is missing` |
| `-o` dir exists | `<dir> already exists` |
| unknown basis | `Basis <b> not found` |
| unknown pseudo | `Pseudo file <p> not found.` (note the trailing dot) |
| input unreadable | `cannot parse geometry from <f>: expected xyz (atom count on line 1, then 'Element x y z' lines) or a Gaussian z-matrix (<f>: No such file or directory)` |
| Z-matrix `Failure` | `cannot parse geometry from <f>: expected xyz (atom count on line 1, then 'Element x y z' lines) or a Gaussian z-matrix (parse error: Syntax error)` |

Uncaught OCaml exceptions — printed by the runtime as
`Fatal error: exception <exn>`, **no `qp_create_ezfio: ` prefix**, exit 2:

| trigger | stderr |
| --- | --- |
| XYZ count 3 but only 2 atom lines | `Fatal error: exception Invalid_argument("Array.sub")` |
| first line `zzz` / `0` / `-1` | `Fatal error: exception Element.ElementError("Element <tok> unknown")` |
| Z-matrix with unknown element | `Fatal error: exception Element.ElementError("Element Q unknown")` |
| `-c ""` | `Fatal error: exception Invalid_argument("index out of bounds")` |

Consequences for the port, to be decided against `CONTRACTS.md` rather than
guessed:

* A non-positive or non-numeric XYZ count is **not** an error: the XYZ reader
  raises `XYZError`, `Molecule.of_file` silently falls back to Z-matrix, and the
  first line is then re-read as an element symbol. That is why `zzz` reports
  `Element Zzz unknown` (capitalised) instead of a count complaint.
* `Element.ElementError` must escape the geometry handler untouched. Do not
  widen the catch to `anyhow::Error`.
* A short XYZ (count exceeds the number of atom lines) is an OCaml
  `Invalid_argument`, not a `Failure`. It is the `Array.sub` in the XYZ reader.

## Observed quirks

1. **Trailing `qp_edit -c` is unconditional and its status is ignored**
   (`ignore (Sys.command ("qp_edit -c "^dir))`), then the directory name is
   printed on stdout and the tool exits 0. Measured: with `qp_edit` absent from
   `PATH`, stderr gets `sh: qp_edit: command not found`, stdout still gets the
   dir name, `exit=0`. The 270/284-file fixtures in `rust/tests/golden/fixtures`
   owe most of their keys to that call, so the golden gate needs a `qp_edit`
   on `PATH` for the fake QP root, not just `qp_create_ezfio`.
2. **`-h` prints the usage block with `argv[0]`** as invoked and exits 0. The
   usage line lists options in the order `[-a] -b <string> [-c <int>] [-g
   <float>] [-h] [-m <int>] [-o EZFIO_DIR] [-p <string>] [-x] [--] FILE`.
3. **Bad Z-matrix atom indices are tolerated**: `H` then `H 9 0.7` is accepted
   (all-zero geometry, no diagnostic), because `Atom.create` on a missing
   neighbour is not reached. Only syntax and unknown elements fail.
4. **Ghost centres** (`-g`): for each pair `i>j` with `d(i,j) < (r_i + r_j) *
   threshold`, one atom `Element.X`, charge 0, at
   `x_i + (x_j - x_i) * r_i/(r_i + r_j)` (same for y, z — weights by covalent
   radii, not the midpoint). Ghosts are appended after the nuclei
   (`nuclei @ ghost`); `h2g` (`-g 5.0`) yields exactly one `X`.
5. **Ghost basis falls back to the H section** when building the AO list.
6. **First-wins for the principal basis** (`not (Hashtbl.mem ...)` guard);
   pseudo entries are plain `Hashtbl.add` and therefore *shadow* earlier ones.
   Basis/pseudo specs are `|`-separated, trimmed, `Element:base` or
   `i,Element:base`; local file first, then `$QP_ROOT/data/basis` (resp.
   `data/pseudo`).
7. **Layout quirks to reproduce byte for byte**: `nucleus_shell_num` is a
   run-length fold over shells, `shell_index` is 1-based per shell,
   `ao_coef`/`ao_expo` are zero-padded to `[ao_num; ao_prim_num_max]`, and
   `ao_power` is `<ao_num ; 3>`.
8. **md5 is over EZFIO round-tripped floats** (`%24.15e` then
   `float_of_string`), not the raw library literals; `Input.Ao_basis.read`
   recomputes and `write` stamps `ao_md5` *before* the other keys. Fixture
   stamps: `h2 44c24554b5fb4cafb9632b78640a58c7`, `hcn
   c4997d5a5cac61bf6d88d46396d6c608`; `MD5("None") =
   6adf97f83acf6453d4a6a4b1070f3754`.
9. **On write failure the created directory is removed recursively** (via
   `Qputils.rmdir`, which chdirs).
10. `Zmatrix.center_of_mass` prints to stdout but is never called from
    `to_xyz`; no stdout quirk to reproduce.

## Golden gate recipe (mechanism validated; the Rust binary has not been gated)

The recipe below was run to completion with the *reference* binary in the fake
root and returns `exit=0`, 6/6 `ok`, hash `477b6d162caaa101` (h2), `be478641…`
(h3), `079df26b…` (c2h2), `2c8e302a…` (b2), `187aa754…` (h2g), `9469511e…`
(hcn). Log: `$HOME/OpenCode/tmp/w3-create-ezfio/ref-verify-create-fakeroot.log`.
It therefore gates the Rust port by substitution, once the port exists.

Three things the fake root must provide, each learned the hard way:

* **`qp_create_ezfio`, bare name, no `.native`.** `find_tools` accepts either
  spelling, but the manifest replay executes the literal argv[0]
  `qp_create_ezfio`; a `.native`-only root dies with
  `FileNotFoundError: 'qp_create_ezfio'`.
* **`qp_edit` on `PATH`.** `qp_create_ezfio` itself shells out to
  `qp_edit -c <dir>` and ignores the status; the 270/284 file counts come from
  that call. Without it the replay prints `sh: qp_edit: command not found`
  and the trees will not match.
* **`qp_run` plus a `src` symlink.** `hcn` carries a `post` step
  `qp_run scf ezfio`; `qp_run` resolves the program as
  `$QP_ROOT/src/hartree_fock/scf`, so a root with only `data/` fails with
  `exited 127`.

```
cd <port-checkout>
cargo build -j4 -p qp_create_ezfio          # from rust/
FAKE=$HOME/OpenCode/tmp/w3-create-ezfio/qproot
rm -rf $FAKE && mkdir -p $FAKE/rust/target/release
ln -sfn ~/qp2/data $FAKE/data
ln -sfn ~/qp2/src  $FAKE/src                # hcn: qp_run scf
cp ~/qp2/rust/target/release/qp_edit $FAKE/rust/target/release/
cp ~/qp2/rust/target/release/qp_run  $FAKE/rust/target/release/
cp rust/target/debug/qp_create_ezfio $FAKE/rust/target/release/
python3 rust/tools/golden-harness/golden_harness.py verify \
    --qp-root $FAKE --corpus rust/tests/golden --only create
```

Caveat to keep in mind: the `create` gate 6/6 is not a pure test of the port.
`hcn` needs the reference `qp_run scf`, and every fixture needs `qp_edit -c`, so
a green run proves `qp_create_ezfio` composition with those tools, not that the
whole chain is Rust.

`ao_md5` cross-check against `md5/basis-md5.jsonl` still has to be written; no
test does it yet. The file is 27 JSON lines, one per *basis* (not per
element/basis pair), each carrying `basis`, `element`, `pseudo` (null for 23 of
them), `n_ao`, `ao_md5` and `ao_extra_md5`; the check is one AO build per line
against those three numeric/hex fields.

`create` argvs (from `rust/tests/golden/manifest.json`): `h2` `-b cc-pvdz -o
ezfio h2_1.xyz`; `h3` `-b cc-pvdz -m 2 -o ezfio h3_2.xyz`; `c2h2` `-b
cc-pvdz_ecp_bfd -p bfd -o ezfio c2h2.xyz`; `b2` `-b 6-31g -o ezfio
b2_stretched.zmt`; `h2g` `-b cc-pvdz -g 5.0 -o ezfio h2_1.xyz`; `hcn` `-b
aug-cc-pvdz -o ezfio hcn.xyz`. Expected file counts: 270 for all but `hcn`
(284).

## Corrections after the second probe round (2026-11-04)

Item 2 above is wrong on one point: `-h` does **not** exit 0 on its own.
`getopt` prints the help block when it sees `-h` and parsing *continues*.
Measured: `qp_create_ezfio -h` alone ends in `qp_create_ezfio: input file is
missing`, exit 2; `qp_create_ezfio -h -b cc-pvdz -o dir mol.xyz` prints the
help, creates `dir`, calls `qp_edit -c dir`, and exits 0. The port must print
the help inline, not return from parsing.

Item 9 is confirmed by a second measurement (Og with `-b cc-pvdz` fails with
`Basis not found for atom 1 (Og)`, and the created directory is gone), but
note that `Qputils.rmdir` is recursive: `fs::remove_dir_all(dir).ok()`
matches, a plain `remove_dir` does not.

New exception class, not in the table above:

| trigger | stderr |
| --- | --- |
| `-o` under a missing parent | `Fatal error: exception Unix.Unix_error(Unix.ENOENT, "mkdir", "/nonexistent/deep/ezfio")` |

`qp-ezfio`'s `EzFio::create` uses `create_dir_all`; the reference calls
`Unix.mkdir` once. The bin must pre-check the parent of the output directory
and synthesise the exception. Do not change `qp-ezfio`.

Also measured, and relevant to the port:

* `-b ''` and `-o ''` are **getopt** errors, not basis/output errors:
  `args.(no+1).[0] <> '-'` indexes an empty string →
  `Fatal error: exception Invalid_argument("index out of bounds")`.
* `-c ''` hits `x.[0]` in the charge parser → the same
  `Invalid_argument("index out of bounds")`.
* `-c` accepts `m1` for `-1`; that is the documented way to pass a negative
  charge.
* The help block is 1817 bytes / 45 newlines, with `argv[0]` substituted at
  byte-lines 0 and 4 only. It is reflowed to width 79 for longer `argv[0]`;
  the port reproduces only the bare-name form. Capture:
  `$HOME/OpenCode/tmp/w3-create-ezfio/probe3/help_bare.txt`
  (md5 `603aa344f8ffc05b93b79dfae8bf9a5d`).
* `pseudo/nucl_charge_remove` is the **run-length encoding** of consecutive
  equal `n_elec` values (`qp_create_ezfio.ml:300-308`), so its dimension is
  the number of runs, not the number of nuclei.
* Cargo 1.92.0 rejects `cargo -j4 build` (`unexpected argument '-j' found`);
  write `cargo build -j4`.
* The default output name is `x ^ ".ezfio"` when `rsplit2 on '.'` yields
  `(x, "xyz")` or `(x, "zmt")`, else `xyz_file ^ ".ezfio"`.

## Corrections after re-reading the reference (2026-11-05)

Three points above are wrong or incomplete. Follow this section where they
conflict; it is the result of a direct re-read and fixture check.

1. **`pseudo/nucl_charge_remove` is NOT run-length encoded.**
   `qp_create_ezfio.ml:296-304` builds one `z_core` per pseudo/nucleus and
   writes `~rank:1 ~dim:[| nucl_num |] ~data:z_core` in nucleus order.
   Fixture evidence: `c2h2/pseudo/nucl_charge_remove.gz` has header dim 4
   and data `2,2,0,0` for its four nuclei; `h2` is dim 2 with `0,0`. The
   RLE fold is `basis/nucleus_shell_num`, not this key.
2. **`ao_extra_md5` is not written by the bin.** `qp_create_ezfio.ml` and
   `Input_ao_basis.ml` never mention `ao_extra*` (grep returns nothing).
   All 27 rows of `rust/tests/golden/md5/basis-md5.jsonl` carry
   `ao_extra_md5 = 6adf97f83acf6453d4a6a4b1070f3754`, i.e. `MD5("None")`;
   the key appears after the trailing reference `qp_edit -c`. Gate (b)
   must assert that constant with this comment (or run a real `qp_edit`),
   never claim the bin produces it.
3. **There is no `QP_EDIT` in the reference.** A grep over
   `~/qp2/ocaml/*.ml` and `~/qp2/src` returns nothing. The injectable
   editor is a port addition (`$QP_EDIT` if set, else `qp_edit`, run as
   `sh -c "<editor> -c <dir>"`, status ignored) so tests can set
   `QP_EDIT=true`.

Checked and **not** a discrepancy: `Zmatrix::to_xyz_string` uses
`format_f(x, 0)`, which is `format!("{x:.6}")` with zero padding—the same
six decimals as OCaml `Printf "%f"`, so the Z-matrix round trip through
`of_xyz_string` matches for non-integral coordinates too.

## Corrections after the front-end port (2026-11-05, session 4)

The front end (pre-S2: `src/error.rs`, `src/cli.rs`, `src/main.rs`; the
shared
numeric parser was `qp_types::num` — since deleted in wave S2: the CLI is
clap derive, errors are `qp_types::error::ToolError`, numeric parsing is
Rust `FromStr`)
is implemented and was differentially tested against
`~/qp2/ocaml/qp_create_ezfio.native`: 30 argv cases, 0 differences over
(exit status, stdout bytes, stderr bytes). Script:
`$HOME/OpenCode/tmp/p4probe/compare.py`; log:
`$HOME/OpenCode/tmp/p4probe/final-test.log`.

Findings that correct or sharpen the text above:

1. **`-axm2` is a getopt error even though `-m2` works.** getopt allows
   concatenation only of *argument-less* short options, so the `m` inside
   `-axm2` is rejected before its attached `2` is read:
   `Getopt.Error("Only non-argument short-options can be concatenated (error
   with option m in -axm2)")`. `-m2` on its own is fine.
2. **`-c m+1` is `-1`, not `+1`.** The `m` prefix negates
   `int_of_string "+1"`, i.e. `-(+1)`; `-c m-1` is `+1`. Negative charges are
   written `mN`, and the sign inside is not special-cased.
3. **`-h` skips only the mandatory `--basis` check.** It does not skip the
   charge/ghost/multiplicity conversions: `-h -c abc` prints the help, then
   `qp_create_ezfio: int_of_string`, exit 2. Nor does it skip the anonymous
   `FILE` check (`-h` alone prints help then `input file is missing`).
4. **A scan error after `-h` suppresses the help.** `-h -Z` prints no help and
   fails on `-Z`. The handler merely records the flag; the help is printed
   after a *complete* scan, not when `-h` is seen.
5. **Long options with `=` hand the driver an empty string for `=`-less
   values.** `--charge=` reaches the charge parser's `x.[0]` and is the
   unprefixed `Invalid_argument("index out of bounds")` (not a getopt error);
   `--ghost=` and `--multiplicity=` are the prefixed `float_of_string` /
   `int_of_string` failures; `--basis=` is accepted as the empty basis string
   and fails later in the run. Contrast `-b ''` / `-o ''` / `-c ''`, which
   are getopt's empty-next-token `Invalid_argument`.
6. **`QP_ROOT` is checked before argv is even collected.** With the variable
   unset, `-h` still dies with the uncaught `Failure` and prints no help. An
   empty `QP_ROOT=""` passes (presence only).
7. **Deviation — non-ASCII short option.** The reference indexes bytes and
   prints the first raw byte in `Unknown option -<byte>`; the Rust port parses
   UTF-8 and prints `Unknown option -<char>`. No entry of `set_specs` is
   non-ASCII, so this is only reachable through a typo.
8. **Deviation — non-UTF-8 argv.** The reference takes bytes; the port rejects
   the argument with the prefixed `command line argument <n> is not valid
   UTF-8` (argument 0 included, no reference counterpart). Chosen over a
   panic or a lossy conversion.
9. `pseudo/nucl_charge_remove` is confirmed per nucleus, not RLE, by the
   `c2h2` fixture (`dim 4`, data `2,2,0,0`). The RLE sentence in the
   2026-11-04 corrections is superseded by the 2026-11-05 section above.

## Corrections after the INPUTS port (2026-11-05, session 6)

`src/geometry.rs` and `src/basis_read.rs` are implemented and wired into
`create_ezfio`; the writer, `ao_md5` and the trailing `qp_edit` are still
successor #6's. Findings that correct, sharpen or extend the text above:

1. **The input half is three calls, not one.** The reference order is
   geometry → `build_basis` (every fetch, then the key bindings) →
   `build_pseudo` → the output-name and `Sys.file_exists` check → the
   per-nucleus basis reads. The port keeps that order as
   `build_basis_table` / `read_pseudo` / `read_basis`; folding the first and
   third into one `read_basis(spec, …)` would put the per-atom reads before
   the pseudo fetch, so a "pseudo not found" run would report a basis error
   instead.
2. **Measured basis error classes** (an auxiliary spec entry makes attempt 1
   report them; reference and port were byte-identical on stdout, stderr and
   exit code in all ten probe cases):
   `GTO_Read_Failure("<line>")` and `Sys_error("<strerror>")` are
   **unprefixed**; a shell line shorter than two bytes is
   `Invalid_argument("String.sub / Bytes.sub")`; `Failure("Angmom should be
   …")`, `Failure("int_of_string")`, `Failure("hd")` and
   `Failure("Failed in of_prim_coef_list")` are **prefixed**; an element
   missing from the file and EOF inside a shell are both `End_of_file` and
   therefore take the `Element … not found` or fallback path, not an error.
   The earlier "other errors become a prefixed failure" is superseded by
   this table.
3. **A directory passes `fetch_channel`.** `Sys.file_exists` and `open_in`
   both accept a directory; the failure only appears on the first read
   (`Sys_error("Is a directory")` for a basis, an empty library for a
   pseudo, because `Pseudo.find` swallows every read error). `--basis=`
   with a fake root therefore gives `Basis not found for atom 1 (H)`, not a
   `Sys_error`.
4. **Deviation — pseudo text is read eagerly and as UTF-8.** The reference
   keeps the channel open and calls `input_line` at write time; the port
   reads the whole file when the table is built, keeps the text, and feeds
   `Pseudo::read_element`. Consequences: a non-UTF-8 pseudo library becomes
   `""` instead of being parsed bytewise, and a mid-file I/O error keeps
   the bytes already read (the reference would raise `End_of_file` there).
   A directory, the case reachable from the CLI, yields `""` in both.
5. **Not reproduced: the two `Filename.temp_file`/`Sys.remove` pairs.**
   `qp_create_ezfio.ml:107-111` and `194-198` create a scratch file under
   `$TMPDIR` and remove it before building each table. It changes no
   stdout, stderr, exit status or EZFIO byte; it would only differ with an
   unwritable `$TMPDIR`, where the reference dies with `Sys_error` and this
   port proceeds. Not worth two syscalls per run.
6. **`-o` ordering is deferred to #6.** The reference checks
   `Sys.file_exists ezfio_file` *after* both tables and *before* the
   per-nucleus reads. The port reads the basis first, so "existing output
   directory plus a bad basis" currently reports the basis error; the check,
   its position and the default-name rule belong to the writer task.

## The WRITE path + the float-parity resolution (2026-11-05)

The write path landed as a parent merge artifact (the create-ezfio write
subagent died 5× on context limits with zero code; the parent wrote
`src/write.rs` + `src/specs.rs` and wired `main.rs` per the digested
`write-sequence.md`). `write.rs` transcribes `write_file ()`
(`qp_create_ezfio.ml:282+`): exists-check position (between the table
builds and the per-nucleus reads, `:273`), `set_file` semantics, the
`ezfio` group, pseudo/basis/coord/charge/electron sections, the
`ao_coef`/`ao_expo` flattening (AO index fastest, per `create_expo_coef`'s
`build` — the nested range loop is deliberate, `#[allow]` documented),
the `QP_EDIT` seam for the trailing `ignore @@ Sys.command "qp_edit -c"`,
the directory line on stdout, and best-effort `remove_dir_all` on any
post-create failure.

**`set_file` details (measured, `ezfio.ml:426`):** single-level
`Unix.mkdir` — a missing parent is the *uncaught*
`Unix.Unix_error(ENOENT, "mkdir", path)`. The port pre-checks the parent
and renders the same exception; any other create failure is also rendered
as ENOENT (documented deviation — the reference's `mkdir` on a parent
that is a file would be ENOTDIR). `set_file` also `mkdir`s the `ezfio`
group directory with its `.version` (`2.0.8\n`) even when no key is
stored there — the port does the same. The reference then stamps the
group with machine state (timestamp, `$USER`, configure-time library
path); the golden harness strips the three files as VOLATILE, but the
port writes them anyway (2026-11-05): the Fortran consumers of a
created tree probe for `ezfio/creation` and `stop 1` without it
(`src/hartree_fock/IRPF90_temp/ezfio_files/ezfio.irp.f:30`) — this is
what broke the hcn `create` gate's `qp_run scf` post-step until the
stamps were added. `creation` = `date` with `LANG` removed (the
reference's `LANG= date`), `user` = `$USER` (empty line when unset,
like `echo $USER`), `library` = the fixed string `qp-ezfio` —
documented deviation: the reference stamps its configure-time ezfio
library path (`ezfio.ml:28`), a static pure-Rust build has no such
path, and nothing in the reference tool layer reads the file back.
A stamp failure raises the same caught `Failure("Unable to create new
ezfio file:\n" ^ filename)` the reference raises.

**The 1-ULP mystery, solved.** The c2h2 `create` differential showed
exactly one differing file: `nuclei/nucl_coord.gz`, H z coordinate
`3.145260163204818e+00` (reference) vs `…817e+00` (port). The port's
`qp-types::point3d::to_bohr` **divided** by `0.52917721067121`, while the
reference `Point3d.of_string` **multiplies** by
`Units.angstrom_to_bohr = 1. /. 0.52917721067121` — the correctly-rounded
reciprocal `0x3ffe3c5176338382` (`Point3d.ml:14,22`, `Units.ml:7`).
`fl(x / k)` and `fl(x * (1 / k))` differ by one ulp on some inputs
(`1.6644`: `0x4009297e291302b0` vs `…b1`), and the `%24.15e` write
(16 significant digits, `ezfio.ml:171`) turned the in-memory difference
into a written-decimal difference. Two traps that hid this:

1. The old `point3d.rs` module doc *asserted* "the reference build was
   changed to divide … do not simplify into a multiplication". False:
   `Units.ml` is unchanged since 2020-09-25, `Point3d.ml`'s last change
   (2024-04-22, CRLF fix) predates the binary by years, and the binary's
   embedded constants (`Units.cmx`, `Point3d.cmx`, compiled
   2026-10-30/11-03) are the correctly-rounded `K` and `1/K`.
2. A naive "wrong constant" search dead-ends: no fixed multiplier or
   divisor reproduces the reference outputs — `T(4.0) ≠ 4·T(1.0)`
   (exactly-representable probe inputs `0.5/1.0/2.0/4.0`) kills every
   single-operation model, because `fl(2^a·x) = 2^a·fl(x)`. The outputs
   are consistent only with a correctly-rounded parse, a correctly-rounded
   multiply-by-reciprocal, and the 16-digit write rounding on top.

Proof: a 12-point probe through the reference binary (the four exact
values + `1.6644/0.74/1.1/2.675/0.1/3.33/7.125/123.456`): the pipeline
`float_of_string` → `× (1.0/0.52917721067121)` → `%24.15e` → re-parse
reproduces all 12 written values **bit-for-bit**; the divide pipeline
matches none. After switching `to_bohr` to the multiply form, all six
corpus fixtures produce byte-identical `create` output (47/47 files each,
volatile trio excluded): h2, h3, c2h2, b2, h2g vs the reference binary
side by side; hcn create-output vs the committed fixture.

**Documented deviation (outside every gate so far):** the reference
`Point3d.to_string` prints `%16.8f %16.8f %16.8f` (`Point3d.ml:49–53`,
width-16 padded, 8 decimals); the port's `point3d::to_string` prints
17-significant-digit components, single-spaced, unpadded. No gated tool
prints coordinates through it (`create` writes raw `%24.15e` arrays; the
`check` transcripts matched), but any future tool that prints molecule
coordinates must use the padded `%16.8f` form.
