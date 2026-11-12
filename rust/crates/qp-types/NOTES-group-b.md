# Group B port notes — `qp-types` (wave 2 / P3)

Merge notes for the parent. Everything below is in
`rust/crates/qp-types/`; nothing is committed (commit is denied in this
sandbox). The OCaml reference is the port worktree (`ocaml/`, contract commit
`be4b71c1`).

## 1. What landed

| Rust file | OCaml source | Lines | Unit tests |
|---|---|---|---|
| `src/bitlist.rs` | `Bitlist.ml` | 507 | 8 |
| `src/determinant.rs` | `Determinant.ml` | 374 | 7 |
| `src/excitation.rs` | `MO_class.ml`, `Range.ml`, `Excitation.ml` | 706 | 8 |
| `src/angmom.rs` | `Angmom.ml` | 669 | 8 |
| `src/zmatrix.rs` | `Zmatrix.ml` | 1102 | 11 |
| `src/pseudo.rs` | `Pseudo.ml`, element table | 1770 | 10 |
| `src/gamess.rs` | `Gamess.ml` (dead in QP2) | 868 | 7 |
| `src/gaussian_primitive.rs` | `GaussianPrimitive.ml` + `ocaml_fmt` | 313 | 6 |
| `src/command_line.rs` | `Command_line.ml` + 5 tool specs | 776 | 8 |

Wave w13 S1 (2026-11-09) removed `src/excitation.rs` (test-only, no
production reference) and `src/gamess.rs` (no production reference), and
shrank the other rows by the P3.5 de-duplication onto `angmom` /
`gaussian_primitive` / `element`.

Integration tests (16):

| File | Tests | What it pins |
|---|---|---|
| `tests/test_group_b_bitlist.rs` | 1 | full `ocaml/tests/test_bitlist.ml` transcript (minus the `Bit` line, see M3) |
| `tests/test_group_b_determinants.rs` | 2 | `test_determinants.ml` transcript + electron-count errors |
| `tests/test_group_b_excitation.rs` | 1 | `test_excitation.ml` transcript, 7 lines |
| `tests/test_group_b_pseudo_md5.rs` | 3 | OCaml-oracle goldens (display/sexp/md5), empty potentials, c2h2 EZFIO arrays |
| `tests/test_group_b_command_line.rs` | 9 | messages, option spellings, `-h`, exit code 2, `Range.of_string` |

Fixtures: `tests/fixtures/pseudo/` holds decompressed copies of the c2h2
`pseudo/*` EZFIO payloads plus `bfd`/`ccecp` source libraries; provenance
(sha256, regeneration commands) is in
`tests/fixtures/pseudo/README.md`. Nothing under `data/` or
`rust/tests/golden/` was modified.

## 2. Build and test

```sh
cd <port-checkout>/rust
cargo test -j4 -p qp-types
cargo clippy -j4 -p qp-types --all-targets -- -D warnings
cargo fmt -p qp-types -- --check
```

All 89 tests pass (73 unit + 16 integration); clippy and fmt are clean.
The crate compiles with the sibling group's modules **absent**: no group B
module imports `crate::element`, `crate::md5`, or `crate::mo_class`.

## 3. Merge checklist (in order)

**M1 — `ElementSymbol` → group A's `Element`.** Group B needed an element
type for `Pseudo`, `Zmatrix` and `Gamess`, and it had to compile alone.
`pseudo.rs` therefore defines a self-contained `ElementSymbol(usize)` plus
`pub(crate) ELEMENTS: [ElementData; 119]` (symbol, long name, Z, mass;
sorted by Z with `ELEMENTS[i].z == i`; `X`/`Ghost`/mass 0 first; generated
from `Element.ml` and checked). Public surface used by group B:

```rust
ElementSymbol::of_string(&str) -> Result<Self, ElementError>  // OCaml Element.of_string
.symbol() -> &'static str     .long_name() -> &'static str
.z() -> u8                    .to_charge() -> f64
.mass() -> f64                Display = symbol
```

`ElementError::Unknown { normalized }` prints `Element {normalized} unknown`.
At merge, alias `ElementSymbol` to group A's `Element` (or delete the table
and re-import), keeping `Display` = symbol and `mass()`/`to_charge()` names.
`zmatrix.rs` and `gamess.rs` import `crate::pseudo::ElementSymbol` and
`crate::pseudo::ElementError`; those two `use` lines change to
`crate::element::…`. `ElementSymbol::of_charge` was **not** ported (group A
owns it).

**M2 — `MoClass` → group A's `MO_class`.** The task list assigned
`MO_class.ml` to no one; `Excitation.ml` cannot exist without it, so
`excitation.rs` contains `MoClass` (five constructors
`create_core/_inactive/_active/_virtual/_deleted`, `of_string`, `numbers`,
`name`, `is_empty`, `to_bitlist`) and `pub type Hole = MoClass;` /
`pub type Particle = MoClass;`. If group A also ported `MO_class`, delete
one copy and re-export the other; the `Hole`/`Particle` aliases must stay
visible from `excitation`.

**M3 — `Bit` (`Bit.ml`) is not ported.** `Bit.ml` is used only by
`qp_set_mo_class.ml` (P4), not by any group B module. The ported
`test_bitlist.ml` transcript drops its one `Bit.to_string Bit.One` line;
the rest of the transcript is byte-identical.

**M4 — dependencies.** `rust/crates/qp-types/Cargo.toml` gained
`md5 = { workspace = true }` (the original task requirement, `Pseudo::to_md5`,
was removed in wave w13 S1; the dependency is still used by `md5.rs`) and
`clap = { workspace = true }` (`Command_line`; the task mandated clap
derive while saying "no other new dependencies" — this is the one
conflict, flagged here). Both workspace entries already existed; only
`Cargo.lock` and the crate manifest changed. If the parent prefers
`qp-types` to stay dependency-free, `command_line.rs` is self-contained
and can move to a future `qp-cli` crate unchanged.

**M5 — `lib.rs`.** Alphabetical `pub mod` lines for the nine modules only.
The crate-level `//!` doc was rewritten (the scaffold said "Status:
scaffold") and now has a module table and a units section; **expect a
conflict with group A there** — keep one crate doc and merge the two
module tables, and keep `#![deny(missing_docs)]`. Group B does not
reference group A's modules.

**M6 — test fixtures.** `tests/fixtures/pseudo/` is new and self-contained
(decompressed text, no flate2). `tests/fixtures/` did not exist before.

**M7 — no name collisions.** `ocaml_fmt` lives inside
`gaussian_primitive.rs` (`gaussian_primitive::ocaml_fmt`); the crate-level
`ocaml_fmt` used by `Pseudo` is the same module. If a sibling needs
`%G`/`%e`/`%f`, import it rather than duplicating.

## 4. Preserved OCaml quirks (do not "fix" without a parity decision)

**Bitlist.** `of_string ?zero` never reads `zero` (kept for call-site and
test parity). `to_mo_number_list` uses `(i lsr 6) lsl 6 + 1` as the word
offset, so words 0..63 all report MO 1..64 (invisible for `n_int = 1`).
Empty input gives `EmptyString` instead of the reference `Out_of_memory`
(`(0-1) lsr 6 + 1` words). `to_int64` fails with `"N_int > 1"`.

**Determinant.** `alpha_beta` splits on `len / 2` and silently drops the
last word when the count is odd. The constructors check the word count
before the alpha/beta electron counts, in that order, with the reference
messages. `to_string`/`Display` ignores the MO count exactly as the
reference's unused `~mo_num` argument did.

**Excitation / MO classes.** *(module removed in wave w13 S1 — test-only, no production reference; the following documents reference behaviour)* The particle-in-inactive check prints
`"Particles can not be in virtual MOs"` (reference bug). `Display` pads
the class name to nine characters only when the list is non-empty
(`"Core     : (1, 2)"` vs `"Core"`). The check order is core-hole,
core-particle, deleted-hole, deleted-particle, virtual-hole,
inactive-particle. `Range.of_string` splits a piece on its *first* `-`,
treats `""` as nothing, and requires `start <= finish`.

**Angmom.** `Xyz::of_string` feeds every non-`xyz` character to
`int_of_string`, even before any letter: `"12"` is the zero triplet,
`"a"` fails. The `Null` state means `"-2"` alone does not fail while
`"x-2"` does. (The reference's `Spd` shell representation was removed in wave w13 S1: no production reference.)

**Zmatrix.** `in_range` only adds/subtracts (NaN unchanged, infinite
input loops); angles wrap to `(-180, 180]`, dihedrals to `(-360, 360]`.
A negative bond length is parsed as a *label* (the reference catches every
exception). Two-token lines split on `' '` only, and a two-token line is always a
label (`Coord`). A three-token `Second` line ignores its middle id token.
`line_of_string` returns tokens in reference order (1,3,5,7,2).
`to_xyz`/`to_xyz_string` place atom 1 at the origin, atom 2 on `+z`, with
the reference (non-standard) angle orientation.

**Pseudo.** Element lookup keeps the *first* line whose first token
matches, skipping blank lines implicitly; a missing element is a valid
empty potential. Tokens split on single spaces (tabs are not separators);
extra header tokens are ignored. `r_power` is `column - 2`. `Display`
formats are `%20f %8d %20f` (local) and `%20f %8d %20f   |p><p|`
(non-local), headers `%20s`; empty blocks are filtered. `to_sexp` is the
sexplib layout `((element X)(n_elec N)(local(...))(non_local(...)))` with
no separators between adjacent list elements; floats use `%.15G` with the
`%.17G` round-trip fallback; `to_md5` is the lowercase md5 of that text.

**Gamess.** The module is dead in QP2 (`Sym` does not exist, no callers).
`string_of_drt` ignores its symmetry argument; `string_of_basis` is a
partial application in the reference and a plain function here.
Constructors set `Data::sym` to `"D4h\n"`/`"C4v\n"`/`"C1\n"` because
`Sym.to_data` is unknowable; the newline is required by the writer.
Integer division truncates toward zero.

**Command_line.** `-h` is prepended by `set_specs`, prints the help and
**skips the mandatory-option check**; the anonymous check still runs, so
`qp_create_ezfio -h` prints help and fails with `input file is missing`
(exit 2) and the other tools print the help twice. `qp_create_ezfio`
takes only the first anonymous argument. `qp_edit` requires exactly one
and reports `"EZFIO_DIR is missing"` otherwise; `qp_set_mo_class` says
`"EZFIO_FILE is missing"` (reference wording), `qp_run` says
`"Inconsistent command line"`. `-c m1` is -1. `qp_edit -s` uses
`Range.of_string` with no 10 000 cap.

## 5. Deliberate deviations (all documented at the definition)

* OCaml `int_of_string`/`float_of_string` accept `0x`/`0o`/`0b`, `_`
  separators and hexadecimal floats; the ports accept decimal digits with
  optional sign and `_` (floats: Rust `f64::from_str`). QP2 inputs are
  plain decimal. Applies to `Angmom`, `Zmatrix`, `Pseudo`, `Command_line`.
* `MO_class`/`Bitlist` MO indices are capped at `1..=10_000`; the
  reference only *warns* above 10 000 and rejects `< 1`. This also stops
  `[1-100000000]` from allocating 10^8 entries before failing.
* `N_int` is capped at `1..=128` (the EZFIO wrapper's range); OCaml would
  abort in `N_int_number.of_int`.
* `Zmatrix` reports `AtomOutOfBounds` instead of an unchecked array read.
* `Pseudo` validates `expo < 0` (`NegativeExpo`) and `r_power ∉ -2..=8`
  (`RPowerOutOfRange`); the reference reaches these through generated
  wrappers/`assert`.
* `Gamess` file I/O (`read_mos`, `read_until_found`,
  `read_natural_mos`, `read_canonical_mos`, `string_of_vec`) is not
  ported (needs `Str` + files); `Guess::Canonical/Natural` carry the
  `$VEC` text instead of a file name.
* CLI: `clap` replaces `getopt` + `Format`. Help text is not
  byte-identical; unknown options and missing values are `clap`
  diagnostics, still exit 2. `-b ""` (empty *next* argv element) crashes
  the reference `getopt` with `Invalid_argument("index out of bounds")`;
  `clap` accepts it as an empty basis. `-c ""`/`--charge=` likewise give
  `int_of_string` instead of an uncaught `Invalid_argument`. `--name=ARG`
  and `--name` still differ where `getopt` does (a next argv element
  starting with `-` is never consumed as a value; `-xARG` and clustered
  flag shorts are). Exit code 2 for every usage error is the measured
  value (`rust/CONTRACTS.md` §2.4), not the "1" in the task text.
* `Pseudo::read_element` takes the whole file as a `&str` and does no
  file I/O; `Pseudo`'s primitive fields are private with validating
  constructors and accessors.

## 6. Pseudo md5: what was asked vs what is possible

The task asked for a golden test of `Pseudo.to_md5` against
`tests/golden/fixtures/*/ao_extra_basis/ao_extra_md5`. **That cannot
work**: `Input_ao_extra_basis.ml` never mentions `Pseudo`; it hashes the
extra *basis* (`Basis.to_md5`), and every row of
`tests/golden/fixtures/basis-md5.jsonl` has
`ao_extra_md5 = MD5("None") = 6adf97f83acf6453d4a6a4b1070f3754`.
`Pseudo.to_md5` additionally has **no caller anywhere in the OCaml
tree** (it exists through `include To_md5`). What the test does instead:

1. pins display + sexp + md5 for H, He, C, Cu and Mg/ccecp, produced by
   the OCaml oracle (`$HOME/OpenCode/tmp/w2b/pseudo_oracle_dune/`), e.g.
   H `6db102f4312533980cca7bb54ff8dafb`, He `b9ea355a55af1007569a8ef99d3f539f`,
   C `0597b7fc1faaf3cd8fd14d166dfa2cdc`, Cu `524ec2d54625d52a022866cbf4c5436e`,
   Mg/ccecp `a1f192d2ca123bc6655d2c82fb2d27ac`. The Cu display string is
   byte-identical to what `ocaml/tests/test_pseudo.ml` prints, so that
   reference test is covered too;
2. pins the empty potential for X and H (`201bdb5b0e8d7950ef2be710b83f4503`,
   `0087fa8e946abf5bbb7ed2227b435200`);
3. parses the `bfd` C and H blocks and reconstructs the c2h2 EZFIO arrays
   (`pseudo_{dz_k,v_k,n_k,dz_kl,v_kl,n_kl,klocmax,kmax,lmax}`,
   `nucl_charge_remove`) from the golden fixture, checking every value and
   dimension. This is the strongest available end-to-end check.

## 7. Not done (and why)

* `Bit.ml` — unused by this group (M3), belongs to `qp_set_mo_class`.
* `Gamess` file readers — `Str` regex + I/O; the module is dead code.
* The CLI **executables**: help rendering, `QP_ROOT`/PATH lookup,
  filesystem and EZFIO checks, `qp_edit` read/write incompatibility,
  banner + backup/restore, exit-1 path. P4 owns them; the structs carry
  the raw strings and `EXIT_USAGE`/`EXIT_FAILURE`.
* `Element` (`of_charge`, `of_z`, the full type) and any `md5`/`mo_class`
  module — group A.
* No fuzzing or property tests beyond the analytic cases; no benchmark
  (no performance claim in the plan for these pure types).

## 8. Oracle provenance (scratch, not committed)

* `$HOME/OpenCode/tmp/w2b/types_oracle_dune/`: verbatim `Bitlist.ml`,
  `Determinant.ml`, `Range.ml`, `MO_class.ml`, `Excitation.ml`,
  `String_ext.ml` + stubs; output `oracle_transcript.txt` (27 lines).
* `$HOME/OpenCode/tmp/w2b/pseudo_oracle_dune/`: `Pseudo.ml` + stubs;
  takes `[file [element]]`, output `oracle_all.txt`, `mg_oracle.txt`.
* `$HOME/OpenCode/tmp/w2b/getopt_oracle/`: the vendored `getopt` parser
  with the `qp_create_ezfio` spec, used to fix every CLI spelling in
  §5 of `tests/test_group_b_command_line.rs`.
* `$HOME/OpenCode/tmp/w2b/gen_elements.py`: generated the 119-row
  element table from `Element.ml`; `cargo test` checks `z == index`.
