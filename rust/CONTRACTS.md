# QP2 tool-layer contracts (P0, T0.5)

**Normative for the Rust port.** This document freezes the byte-level
behaviour the Rust executables must reproduce. Every claim is either
measured on the reference build or transcribed from source; file:line
references are against `~/qp2 @ be4b71c1` unless noted. The post-fix
`qp_edit -c` target (commit `32c9899d`, `15e725aa` on `origin/qp-edit-fix`)
is normative where it differs from the reference build; the pre-fix
behaviour of `be4b71c1` is recorded as *evidence* in
`tests/golden/transcripts/` and is the comparison target only for
`verify --contract pre-fix`.

Sources: the frozen spec (section 1 below), decisions D4 and D6,
`tests/golden/` (corpus, D6 oracle) and
`tests/golden/contracts/check-mode.md` (C1–C7).

---

## 1. EZFIO on-disk format

The format is implemented once, by `ezfio-core` (`external/ezfio`
submodule, the ezfio-rs implementation); `qp-ezfio` is the QP API layer
the Rust tools link. Citations below marked *historical C
implementation* refer to the Fortran client deleted by the ezfio-rs
switch (`external/ezfio/src/libezfio_util.irp.f`); they are kept only to
explain corpus bytes that no current writer produces.

One directory per EZFIO **group**; a key is the file inside the group.
There is no index or manifest: presence of a file is the presence of a
key. A group *exists* iff `group/.version` exists
(`external/ezfio/Ocaml/ezfio.ml:54–57`; the same rule is
`fs.rs:340` `group_exists` in `ezfio-core`). A group directory
without `.version` is invalid: the OCaml writer's `create_group` →
`mkdir` calls `Unix.mkdir` on it and fails `EEXIST`
(`ezfio.ml:69–83`). Both writers create `.version` first.

### 1.1 Scalars (plain text, one value, one line)

| type | OCaml writer (`Ocaml/ezfio.ml:166–170`) | Fortran writer (historical C implementation, superseded by `ezfio-core`) |
|---|---|---|
| int / int64 | `%20d\n` → 21 B, e.g. `h2/nuclei/nucl_num` | `'(I20)'` → 20 chars + `\n`, 21 B |
| real / double | `%24.15e\n` → 25 B, lowercase `e`, one digit before the point: `   1.000000000000000e-15\n` | `'(E24.15)'` → 25 B, uppercase `E`, zero before the point: `  -0.928871749993754E+02\n` |
| logical | `"T\n"`/`"F\n"` + an *extra* `\n` (`fortran_string_of_bool`) → 3 B: `h2/ao_basis/ao_cartesian` = `F\n\n` | `'(L1)'` → 2 B: `hcn/ao_basis/use_cgtos` = `F\n` |
| string | `%s\n`, no padding: `h2g/ao_basis/ao_basis` = `cc-pvdz\n` | `'(A)'` of the declared `character*(N)` → blank-padded to N + `\n`: `hcn/mo_basis/mo_label` = 64 + `\n` = 65 B; `hcn/ao_basis/ao_md5` = 32 + `\n` = 33 B |

Writer attribution in the corpus: `hcn/hartree_fock/energy` is the only
Fortran-written scalar real (uppercase `E`); the two 2-B logicals are
the only Fortran logicals; the 4 Fortran-written arrays are
`hcn/mo_basis/mo_{coef,occ,symmetry,class}.gz` (all other arrays are
OCaml-written). Both writers' files must be readable by the port.

**Reader rules (normative).** Read the file as text, take the first
line, `trim`, then parse: `i64::from_str`, `f64::from_str`, `T`/`t`/
`F`/`f` for booleans, the rest for strings. Do **not** assume field
widths (20/24/32 vary by writer); do not reject `F\n\n` or `F\n`.
String arrays: the OCaml reader trims each element
(`ezfio.ml:363–377`), the original C reader read the line into a fixed
`character*(32)` buffer then the declared width (`ezfio_read_array_ch`,
historical C implementation) — so blanks are padding, never data, and a
line longer than 32 chars is silently truncated by Fortran. `ezfio-core`,
the implementation the port links, trims each character record on both
ends, which is what the OCaml reader did.

### 1.2 Arrays (gzip-compressed, `group/key.gz`)

Decompressed payload (this is the D6 oracle; the gzip container bytes —
mtime, OS byte, level — are not a contract):

```
line 1:  rank, right-justified in 3 chars + \n        OCaml %3d / Fortran (I3)
line 2:  dims, one per dimension                     (see below)
lines 3+: one element per line, column-major order (first index fastest)
```

* **Rank line** is 3 chars + `\n` in both writers (`"  2"`).
* **Dims line** differs by writer and both forms occur in the corpus:
  OCaml `Printf.fprintf "%20d "` per dim → 21·rank chars
  (`h2/nuclei/nucl_coord.gz`, rank 2, = 42); Fortran `'(30(I20,X))'` →
  20·rank + (rank−1) chars, gfortran elides the trailing `X`
  (`hcn/mo_basis/mo_coef.gz`, rank 2, = 41; rank 1 = 20).
  **Readers must split on whitespace and parse the integers; they must
  accept 20 or 21 chars per dim.** Dims may name runtime symbols in the
  EZFIO registry — the port resolves them dynamically, as OCaml does
  (`ezfio.ml:347–356`).
* **Body element widths:**

| type | OCaml | Fortran |
|---|---|---|
| int / int64 | 20 chars + `\n` (`%20d`) | 25 B: 24-char `I20` field + `\n` (`write(buffer(i)(1:24),'(I20)')`, historical C implementation template) |
| real | 24 chars + `\n` (`%24.15e`) | 25 B: 24-char `E24.15` field + `\n` |
| logical | `"T\n"\n` | 25 B: 24-char `L1` field + `\n` (absent from corpus) |
| string | raw + `\n` | declared `len(dat(1))` chars + `\n` (`ezfio_write_array_ch` in the historical C implementation) |

Observed over the whole corpus: body line lengths {20: 1247, 24: 6194,
1: 17 (`nucl_label`, declared width 1), 32: 55 (`mo_class`, declared
width 32)}; dims line lengths {42, 21, 63: OCaml; 20, 41: Fortran}.

* **Numeric text forms** are not just a case difference: OCaml emits
  `d.dddd…e±XX`, Fortran `±0.dddd…E±XX`; both are 24 chars and both
  parse with `f64::from_str` (which accepts `e`/`E`; no `D` exponents
  occur — `Gto.read_one` even rewrites `D`→`e` on input, `Gto.ml:118`).
  A Rust writer must reproduce the **OCaml** form for every file it
  writes (the golden `create` comparison is byte-for-byte on the
  decompressed payload and the fixtures it owns were OCaml-written);
  it must *accept* the Fortran form when reading.
* **Float value parity (`qp_create_ezfio`):** the reals written are the
   in-memory OCaml doubles, and the coordinate conversion is a
   **multiplication** by the correctly-rounded reciprocal
   `1. /. 0.52917721067121` (`Units.ml:7`;
   `Point3d.ml:11–14,22`: `l.(i) *. Units.angstrom_to_bohr`) — **not** a
   division by `0.52917721067121`. The two forms differ by one ulp on
   some inputs (corpus: c2h2's `1.6644` Å z coordinate lands on
   `…818e+00` with the multiply form and `…817e+00` with division, after
   the 16-significant-digit `%24.15e` write), and the `create` gate
   compares bytes. Verified 2026-11-05 with a 12-point probe through
   the reference binary: parse (correctly rounded) →
   `× (1/0.52917721067121)` → `%24.15e` reproduces all 12 reference
   values bit-for-bit; the divide pipeline matches none of them.
* `determinants/psi_det` is an `integer*8` array
  (`src/determinants/EZFIO.cfg:86–90`), not sexp text.
* **Writer determinism:** the OCaml writer pipes through a shell
  `gzip` (`ezfio.ml:405–426`), the original C writer through a
  `.<name><PID>.gz` FIFO + `gzip` + `rename` (historical C
  implementation, superseded by `ezfio-core`). The port must emit
  **deterministic** gzip (mtime 0); container bytes are never compared.

### 1.3 `.version`

Text `2.0.8` + `\n`. OCaml writes `%s\n` → **6 B** (222 files in the
corpus); Fortran writes `'(A)'` of a 32-char `libezfio_version` →
**33 B** (exactly `hcn/{hartree_fock,json,mo_basis}/.version`, created
by the Fortran SCF). The port writes the 6-B form (it is the writer of
the groups it creates) and reads either. `.version` must be written
before any key of the group (section order above).

### 1.4 Volatile keys and scratch

* `ezfio/creation` = `LANG= date` output, `ezfio/user` = `$USER`,
  `ezfio/library` = absolute path of the EZFIO source tree, compiled in
  at generation time (`Ocaml/ezfio.ml:28, 430–435`). Written once when
  the directory is created; never rewritten. The corpus and the
  snapshots exclude them (machine state, no contract), **but the port
  must still write all three on create**: Fortran refuses to run on a
  directory without `ezfio/creation`
  (`src/ezfio_files/ezfio.irp.f:28–33`, `stop 1`).
* `work/` is scratch (`qpstop`, `qpkill`, `qp_edit`'s
  `backup.tar`); excluded from snapshots.
* **Units:** coordinates are bohr (`qp_create_ezfio` converts
  Angstrom→bohr; `-a` means the input is already bohr).
* **Ordering:** arrays are column-major (first index fastest).

### 1.5 What is *not* frozen

gzip container bytes, file mtimes, directory entry order,
`ezfio/{creation,user,library}`, `work/` content, and the
`ezfio_files/ezfio_convention` int (20250211 default, 20210101 with
`-q`; its value *is* part of the `create` fixture).

---

## 2. CLI surface (normative, from source — **not** the man pages)

### 2.1 Binaries, PATH, program lookup

Six executables: `qp_run`, `qp_edit`, `qp_create_ezfio`,
`qp_set_mo_class`, `qp_print_basis` and
`element_create_db`, built by
`scripts/compilation/build_rust_tools.sh` into
`rust/target/release/`. PATH after sourcing:
`$QP_ROOT/bin:$QP_ROOT/rust/target/release` (`etc/paths.rc`).
`qp_run` resolves PROGRAM through `data/executables` — 113 lines
`name $QP_ROOT/path`, `$QP_ROOT` substituted, sorted by name
(`Qpackage.ml:52–73`) — and prints a `Possible choices:` list before
failing. `QP_ROOT` unset → `Qpackage.root` fails
`"QP_ROOT environment variable is not set.\nPlease source the quantum_package.rc file."`
(exit 2 for every tool, since it is reached before argument parsing in
`qp_run` and on first use elsewhere).

### 2.2 Parsing rules (shared)

`Command_line` (`ocaml/Command_line.ml:135–200`) over the vendored
`getopt` (`external/opampack/.../lib/getopt/getopt.ml:25–95`):

* `set_specs` **prepends `-h/--help`**, parses, and *if `-h` was seen*
  prints help (stdout) and **skips the mandatory-option check**; the
  anonymous-argument check still runs in every `main`, and most tools
  call `Command_line.help()` again before failing — **hence `-h` prints
  the help twice and still exits 2** (byte counts below: no-arg help
  × 2).
* Long options `--name=arg` / `--name`; short `-xARG` immediate, else
  the next argv element **only if it does not start with `-`**;
  argument-less short options concatenate (`-cc-pvdz` is legal for
  create); `--` ends option parsing; an empty `""` argv element counts
  as an anonymous argument.
* Only Mandatory **options** are checked centrally: the error is
  `--<long> option is missing.` (`Command_line.Error`, caught only by
  `qp_create_ezfio`). Unknown options raise `Getopt.Error` which **no
  tool catches** → `Fatal error: exception Getopt.Error("Unknown option --nope")`,
  exit 2.
* Help → stdout; errors → stderr. OCaml uncaught exception → exit 2.
  `qp_edit` alone catches `Failure`/`Invalid_argument` **inside** its
  `run` (`qp_edit.ml:590–610`): banner of `=-…`, the message, backup
  restored, **exit 1** (failures raised *before* the `try` — e.g.
  missing EZFIO_DIR — are uncaught and exit 2).

### 2.3 Per-tool spec (as in source)

* **`qp_create_ezfio`** (`qp_create_ezfio.ml:700–825`): `-o EZFIO_DIR`,
  `-b <basis>` (mandatory; may be `base | Elem:base | i,Elem:base`,
  per-atom), `-a` (input in bohr), `-c <int>` (charge; negative written
  `m1`), `-g <float>` (ghost atoms), `-m <int>` (multiplicity, def 1),
  `-p <pseudo>`, `-q` (old 20210101 conventions), `-x` (cartesian, 6d),
  anonymous `FILE` (mandatory). Default output name: the input with
  `.xyz`/`.zmt` → `.ezfio`, else `<file>.ezfio`. Refuses an existing
  directory (`<dir> already exists`), takes the **first** anonymous
  arg (extras silently ignored), prints the directory name on stdout,
  ends with `qp_edit -c <dir>` (exit code ignored), removes the
  directory on write failure.
* **`qp_edit`** (`qp_edit.ml:525–615`): `-c`, `-w <file>`,
  `-r <file>` (mutually exclusive: `read and write options are
  incompatible`), `-n <int>`, `-s <range>`, anonymous `EZFIO_DIR`
  (exactly one). Non-`-c` runs write `work/backup.tar` first and
  restore it on failure.
* **`qp_run`** (thin launcher, single-node): `-p <prefix>`
  (`gdb` → `gdb --args `), anonymous `PROGRAM EZFIO_DIR` (exactly
  two). Order of checks: EZFIO dir exists
  (`the EZFIO directory does not exist: <d>`) → program in
  `data/executables` (else `Possible choices:` + `program not found in
  data/executables: <p>`) → prints `Date:`, banner, `Git Commit/Date/
  SHA1`, `EZFIO Dir` → `qp_edit -c` (non-zero → `the input is
  inconsistent (qp_edit -c failed)`) → run the program with the
  inherited environment; non-zero job exit prints
  `Program exited with code %d.` and **propagates the code**; finally
  `Wall time:`. SIGUSR1 → `qp_stop`. The launcher creates no
  `work/` entries and sets no environment variables for the job
  (multi-node removed 2026-11-12).
* **`qp_set_mo_class`** (`qp_set_mo_class.ml:189–240`):
  `-c <range>` core, `-i` inact, `-a` act, `-v` virt, `-d` del,
  `-q` query, anonymous `EZFIO_DIR`.
* **`qp_print_basis`** (`qp_print_basis.ml`): **no** `Command_line`,
  reads `Sys.argv.(1)` directly — no `-h`, no options, extra args
  ignored. Prints the basis RST, then the MO RST; either missing →
  `Error reading EZFIO file` / `Error reading basis set` / `Error
  reading the mo set`, exit 2.
* **`element_create_db`** (`element_create_db.ml`): no arguments;
  writes `$QP_ROOT/data/list_element.txt`, 110 lines
  `"%3d %3s %s %f\n"`, line 0 `  0   X Ghost 0.000000`; silent, exit 0.

### 2.4 Measured exit codes (reference build, bare argv[0])

| invocation | stdout | stderr (first line) | exit |
|---|---|---|---|
| `qp_create_ezfio` | — | `qp_create_ezfio: --basis option is missing.` | 2 |
| `qp_create_ezfio -b cc-pvdz` | — | `qp_create_ezfio: input file is missing` | 2 |
| `qp_create_ezfio -h` | 1817 B help | `qp_create_ezfio: input file is missing` | 2 |
| `qp_create_ezfio --nope` | — | `Fatal error: exception Getopt.Error("Unknown option --nope")` | 2 |
| `qp_edit` | 717 B help | `Fatal error: exception Failure("EZFIO_DIR is missing")` | 2 |
| `qp_edit -h` | 1434 B (help × 2) | same | 2 |
| `qp_run` | 2518 B help (list of programs) | `Fatal error: exception Failure("Inconsistent command line")` | 2 |
| `qp_run -h` | 5036 B | same | 2 |
| `qp_run scf /nonexistent` | — | `Fatal error: exception Failure("EZFIO directory /nonexistent not found")` | 2 |
| `qp_run nosuchprog <valid dir>` | 1659 B `Possible choices:` | `Fatal error: exception Failure("Executable nosuchprog not found")` | 2 |
| `qp_run scf <bad input>` | — | (job output) | job's exit code, after `Program exited with code %d.` |
| `qp_set_mo_class` | 639 B help | `Fatal error: exception Failure("EZFIO_FILE is missing")` | 2 |
| `qp_set_mo_class -h` | 1278 B | same | 2 |
| `qp_print_basis` | — | `Fatal error: exception Invalid_argument("index out of bounds")` | 2 |
| `qp_print_basis /nonexistent` | — | `Fatal error: exception Failure("Error reading EZFIO file")` | 2 |
| `qp_print_basis <dir without MOs>` | basis RST | `Fatal error: exception Failure("Error reading the mo set")` | 2 |
| `qp_print_basis <dir with MOs>` | basis + MO RST | — | 0 |
| `element_create_db` | — | — | 0 |

**Wave S2 note (2026-11-09).** The table above stays the reference
measurement. The port no longer reproduces its help behaviour or its
`Fatal error: exception` stderr lines: `-h`/`--help` is clap's built-in
(prints the help once to stdout and exits 0, no continuation into the
argument checks), a missing or doubled positional is a clap
required-argument / unexpected-argument diagnostic (exit 2), and every
error/notice string is plain Rust text (`<tool>: <message>`). The byte
counts and the `Fatal error: exception` first lines are reference
evidence only; the port contract is the exit column (DEVIATIONS.md
A4/A9/A10/A11).

### 2.5 Man-page drift (the man pages are wrong; do not copy them)

| man page | drift vs source |
|---|---|
| `qp_create_ezfio.1` | documents `-d/--dummy=<float>` (source: `-g/--ghost`); documents `-x/--cart` (source: `--cartesian`); omits `-q/--qmc` |
| `qp_edit.1` | omits `-w`, `-r`, `-n`, `-s` |
| `qp_print_basis.1` | **does not exist** |
| `element_create_db.1` | **does not exist** |
| `qp_run.1`, `qp_set_mo_class.1` | consistent with source |

---

## 3. D4 — the md5 contract

**Storage.** `ao_basis/ao_md5`, `ao_extra_basis/ao_extra_md5`,
`mo_basis/ao_md5`: 32 lowercase hex chars. Fortran treats the values
as opaque (no `ezfio_get_*md5` readers in tracked Fortran: existence
only). The OCaml tool layer compares them; the port must too.

**Canonical object.** `Basis.t = (Gto.t * Nucl_number.t) list`
(`Basis.ml:6`), where `Gto.t = { sym : Angmom.t; lc :
(GaussianPrimitive.t * AO_coef.t) list }` (`Gto.ml:11–15`) and
`GaussianPrimitive.t = { sym : Angmom.t; expo : AO_expo.t }`
(`GaussianPrimitive.ml:6–10`). All `[@@deriving sexp]`;
`Angmom.t` is the atom `S|P|D|F|G|H|I|J|K|L` (`Angmom.ml:4`),
`Nucl_number.t` a bare int, `AO_coef.t`/`AO_expo.t` bare float atoms
(`Qptypes.ml:168–192`). The *long* basis is first collapsed by
`Long_basis.to_basis` (`Long_basis.ml:23–39`): a Cartesian shell keeps
only the component whose `Angmom.Xyz` equals the shell's first one
(this is what makes the digest independent of the 5-component
representation of a P shell).

**Serialization.** `Sexp.to_string` (sexplib) with three traps:
adjacent list elements print with **no separator** — `((a)(b))`,
`(sym S)` — while sibling atoms are separated by one space; a float
atom prints `%.15G` and falls back to `%.17G` when the shorter form
does not round-trip (`sexplib0/float.ml:19–23`, `sexp_conv.ml:34`); a
symbol atom is **quoted only if** it is empty, contains any of
`"()\;`, contains a byte in `0x00..=0x20` or `0x7f..=0xff`, or
contains an **adjacent** `#|`/`|#` sequence (`a#b|c` prints
unquoted) — escapes are a backslash + three decimal digits per byte
(verified against a sexplib0 v0.14 oracle; implemented and pinned in
`crates/qp-types/sexplib_compat.rs`). These rules are what the Rust
reimplementation must match byte for byte; the P0 corpus pins them.

**Pipeline** (`To_md5.ml`): `sexp_of_t t |> Sexp.to_string
|> Digest.string |> Digest.to_hex |> MD5.of_string` (MD5 of the UTF-8
sexp string, lowercase hex). `MD5("None") =
6adf97f83acf6453d4a6a4b1070f3754` (the digest written when a section
cannot be read).

**Reference values** (see `tests/golden/md5/basis-md5.jsonl`, 27
bases, last element of each basis file; `generate` reproduces it):
H2/cc-pvdz `44c24554b5fb4cafb9632b78640a58c7`; HCN/aug-cc-pvdz
`c4997d5a5cac61bf6d88d46396d6c608` (the `hcn` fixture's `ao_md5`).
The two fixture values are digests of the **fixture EZFIOs** — their
basis values went through the lossy `E24.15`/`%24.15e`
writer/reader round trip — and are NOT reproducible from the original
basis-library files at full precision. Do not assert them against a
library-derived recomputation; the 27-row corpus is the golden test.
Expected by-construction collisions: `6-31g_star`/`6-31g_star_star`
(KRYPTON) and `def2-svp`/`def2-tzvpp` (RADON).

**Port notes (P3, measured).** `Element::of_string`'s error message
carries the **normalised** symbol (input `qq` →
`Element Qq unknown`), not the raw input.

**Read/write semantics (pre-fix OCaml, evidence).**
`Input_ao_basis.read` / `Input_ao_extra_basis.read` recompute the
digest on **every** read and **write** the stamp (catch-all writes
`MD5("None")`); `write` stamps before the other keys
(`Input_ao_basis.ml:166–268`). `Input_mo_basis.read_ao_md5`
(`Input_mo_basis.ml:56–77`) writes `mo_basis/ao_md5` when absent, then
fails `The current MOs don't correspond to the current AOs.` on
mismatch; `read_mo_label` writes `None` when absent; `read_mo_num`
fails `More alpha electrons than MOs`.

**Port contract.** Compute the identical digest (P3 test: all 27 rows).
**Do not inherit the read-side echo write** — the post-`32c9899d`
contract (C2) forbids any write during `qp_edit -c`; stamping moves to
explicit write paths only. Comparison in `read_ao_md5` is hard failure
(exact parity first, D4); the trust-on-first-use re-stamp (warn once,
rewrite both stamps, continue) is the documented fallback, off by
default.

---

## 4. Check-mode contract (T0.3)

Normative target: `tests/golden/contracts/check-mode.md` (rules C1–C7:
exit 0, no rewrite, missing defaults materialised — 7 files for `hcn`,
the `*_qp_edit` scratch pair not written, no coordinate quantisation,
idempotence, malformed input still reported). The pre-fix transcript of
`be4b71c1` (9 files created for `hcn`, incl. the scratch pair) is
evidence, not the target. `verify --contract {post-fix,pre-fix}` gates
against the right one.

Verified on the Rust port (2026-11-05, `w3/p4-qp-edit` merged): the
post-fix target holds — a check-mode run creates 7 files for `hcn`,
all 7 byte-identical; the other five fixtures create 0. Additional
measured facts (reference build): `-h` is not an exit — it prints the
help on stdout and continues; without `EZFIO_DIR` the help is printed
twice (1434 B) and the process exits 2 with `Fatal error: exception
Failure("EZFIO_DIR is missing")` (wave S2: the port's `-h` is clap's
built-in — exit 0, no continuation — and the missing-`EZFIO_DIR` error
is clap's diagnostic). The reference's `Input_bitmasks` is `section = false` and the
driver never visits it; the port does not model it as a handler at all (the
determinants handler carries the shared derivation itself). The
reference's `at_exit` removes
`work/backup.tar` unconditionally (a pre-existing archive is deleted
even by `-c`); the port skips both create and remove in check mode
("check writes nothing") — documented deviation. The check-mode engine
refuses to open a directory without `.version` (loud failure, not a
silent no-op). `determinants_by_hand` must stay the last handler: a
following section's title line would leak into its tail parse.

---

## 5. What the plan got wrong (corrections)

* §4 "character padded/truncated to **24**" holds for the **OCaml**
  writer only; Fortran uses the declared `character*(N)` width (1 for
  `nucl_label`, 32 for `mo_class`), and the Fortran *reader* truncates
  lines to 32 chars. Reader rule: split on whitespace, trim, accept
  any width.
* T0.5 says "CLI surface/exit codes **per man page**"; the man pages
  are wrong (§2.5) — the tables here are from source and re-measured.
* §4 "scalar one value, one line": a logical scalar written by the
  OCaml layer is **3 bytes** (`F\n\n`), not one line of `F`; Fortran
  writes 2 (`F\n`). Both must read, and the port's writer must match
  the 3-byte form where it owns the file.
* T4.2 "interactive REPL": there is no REPL, prompt, line reader or
  echo in `qp_edit`. The "interactive" mode renders the section to a
  `Filename.temp_file "qp_edit_" ".rst"`, launches `$EDITOR` (fallback
  `vi`), re-reads and applies it; `-c` skips temp file, editor and
  writes (`Rst_string` is the identity).
* `contracts/check-mode.md` C7 claimed any bad section gives
  `Info: Read error in <section>` + exit 0: wrong. Measured: a sexplib
  conversion failure prints `Error:  (...)` and exits 0 (no `Info:`);
  an escaped `Invalid_argument` prints `Info:`, the retry raises, the
  driver banners and exits 1. (C7 reworded accordingly.)
