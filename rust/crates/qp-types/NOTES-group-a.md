# Group A — core types (plan P3), implementation notes

Status: module group A complete, all 27 md5 golden rows reproduced, 87 tests
green, `clippy -D warnings` clean, `cargo fmt --check` clean. Work is
**uncommitted** on branch `w2/p3-types-a`.

## What was built, and where

Everything lives in `rust/crates/qp-types/`. The crate is a library; there is
no binary.

| file | LOC | ports |
|---|---:|---|
| `src/lib.rs` | 66 | crate docs + module list (only the `mod`/`pub use` lines are mine) |
| `src/element.rs` | 1167 | `ocaml/Element.ml` (119 elements; `ALL` and the match arms of the six lookup functions generated from an oracle dump, see "Generation") |
| `src/point3d.rs` | 248 | `ocaml/Point3d.ml`, `ocaml/Units.ml` |
| `src/charge.rs` | 116 | `ocaml/Charge.ml` |
| `src/atom.rs` | 195 | `ocaml/Atom.ml` |
| `src/molecule.rs` | 606 | `ocaml/Molecule.ml`, `ocaml/Multiplicity.ml` |
| `src/gto.rs` | 851 | `ocaml/Gto.ml`, `ocaml/Angmom.ml` (`Xyz` part), `ocaml/GaussianPrimitive.ml` |
| `src/basis.rs` | 433 | `ocaml/Basis.ml`, `ocaml/Long_basis.ml` |
| `src/basis_library.rs` | 255 | no OCaml counterpart: the `qp_create_ezfio` + `Input_ao_basis.read` sequence packaged as one call |
| `src/sexplib_compat.rs` | 348 | sexplib0's `Sexp.to_string` / float printing (not in this repo; reversed from sexplib0 and checked against an oracle) |
| `src/md5.rs` | 56 | `ocaml/To_md5.ml`, `ocaml/MD5.ml` |

Tests: `src` unit tests (55) plus four integration files,
`tests/test_group_a_md5.rs` (504), `test_group_a_molecule.rs` (296),
`test_group_a_sexp.rs` (255), `test_group_a_element.rs` (203) — 32 more tests.

`Cargo.toml` gained one line, `md5 = { workspace = true }`; the workspace
already pins `md5 = "0.7"`, and `Cargo.lock` gained the single `md5 0.7.0`
package. No other dependency was added. `serde` was already in the scaffold's
dependency list and is unused by group A; I left it because the list is shared
with group B and the merge decides.

## The contract that had to be satisfied

`Basis.to_md5 = MD5(Sexp.to_string (sexp_of_t b))`, byte for byte, for the
27 rows of `rust/tests/golden/md5/basis-md5.jsonl` (`CONTRACTS.md` §3):

* the hashed object is the section read **twice**, tagged `Nucl_number` 1 and
  2, with every exponent and coefficient put through the EZFIO
  `%24.15e` → `float_of_string` round trip (`LibraryBasis`);
* sibling **lists** print with no separator, sibling **atoms** with one space;
* a float atom prints `%.15G` and falls back to `%.17G` only when the short
  form does not parse back to the same `f64`;
* atoms that need quoting are escaped the way `Sexp.to_string` does it.

The 27-row test asserts `ao_md5` and `n_ao` per row, the two by-construction
collisions, a two-nucleus `cc-pvdz` fluorine row against the OCaml test's
literal (`227bd952313a812f889a6aa83194c161`), and that the corpus itself still
has 27 rows. `basis_file_sha256` and `pseudo` are parsed so that a schema
change is noticed, but not asserted (no SHA-256 dependency; the pseudo does
not enter `ao_md5` for these rows).

## Build, test, measure

From `rust/` (toolchain 1.92.0, arm64 macOS / Apple M1, `--offline` because
the crate cache is warm):

```sh
cargo build  -j4 -p qp-types --offline
cargo test   -j4 -p qp-types --offline
cargo clippy -j4 -p qp-types --all-targets --offline -- -D warnings
cargo fmt    -p qp-types -- --check
```

The md5 test needs the installed basis library. It looks at `$QP_BASIS_DIR`,
then `$QP_ROOT/data/basis`, then `$HOME/qp2/data/basis`, and **panics with the
list of paths it tried** if none exists — a missing library must not look like
a green run.

Last run on 2026-11-05 (Apple M1, rustc 1.92.0):

| command | result |
|---|---|
| `cargo test -j4 -p qp-types --offline` | **87 passed**, 0 failed (55 + 8 + 4 + 10 + 10), 0 doc-tests |
| `cargo clippy … -D warnings` | clean |
| `cargo fmt -- --check` | clean |
| md5 test binary, `--release` | 4 tests in **0.13 s** (wall 0.32 s) |
| md5 test binary, debug | 4 tests in **0.70 s** |

The release number includes reading the 27 basis-library files twice each
(`LibraryBasis` calls `read_element` twice, and each call rewinds and scans
from offset 0). Total object hashed over the 27 rows: **178 076 bytes**; the
largest row is `ano-rcc/CURIUM`, 64 118 bytes of S-expression for 608
Cartesian functions. Sum of `n_ao` over the corpus: 2856.

## Cost

Nothing here is superlinear except the distance matrix, and nothing allocates
per iteration in a loop that runs over primitives more than a handful of
times.

| routine | time | memory |
|---|---|---|
| `Element::*`, `Charge`, `Point3d`, `Atom::*` | O(1) (match arms over the 119 variants) | O(1), except the strings they return |
| `Molecule::of_xyz_string` | O(lines × tokens), one `Vec` of nuclei | O(n_atoms) |
| `Molecule::distance_matrix` | O(n²) | O(n²) (`n²` `f64`s) |
| `Gto::read_one` | O(primitives), one line `String` + one token `Vec` per line | O(shell) |
| `Basis::read` / `find` | O(bytes of file) | O(shells) / O(longest line) |
| `LongBasis::of_basis` | O(n_ao) = Σ (l+1)(l+2)/2 | O(n_ao) |
| `Basis::to_sexp` / `to_md5` | O(output bytes) | O(output bytes): the printed form is materialised before hashing, exactly as OCaml does |
| `md5_hex` | O(len) | 16-byte state + one 32-byte hex `String` |

`Vec::with_capacity` is used where the final size is known (`of_basis`,
`AngmomXyz::of_angmom`, the `to_string_*` line buffers). In `Gto::read_one`
the capacity request is capped at 64 because `n` comes from the file: a
fabricated `S 100000000` must not reserve gigabytes before the missing lines
are detected.

## Tests, and what pins what

* **md5** (`test_group_a_md5.rs`): all 27 rows; the two collisions; the
  `cc-pvdz` fluorine reference value; deterministic re-read.
* **molecule**: port of `test_molecule.ml` (`test_atom`, the three water
  names and their errors), the post-fix `%.17g` `to_string` block, the
  distance matrix against an independently computed `sqrt`, a CRLF+tab xyz
  file, `to_xyz` → `of_xyz_file` round trip, typed `XyzError` cases.
* **sexp**: printer spacing, the `must_escape` predicate (including
  `"a#b|c"` staying unquoted and non-ASCII becoming decimal escapes), `%e`
  with a two-digit exponent, `%f`, `%G` trailing-zero removal, the
  `%.15G`/`%.17G` fallbacks, subnormal and `-0` pins.
* **element**: table lengths, `of_charge` round trip, `of_string` case forms
  (`NEON` → `Ne`), radius/mass pins against the reference dump, error texts.
* **unit tests**: always a closed-form or table-derived reference, plus an
  edge case (empty input, single shell, zero charge, `Og`).

## Deliberate deviations from the reference

These are the only places where the Rust behaviour is not the OCaml
behaviour; each is documented at the item.

1. `Charge::to_int` uses Rust's `as i32` (saturating, NaN → 0) where OCaml's
   `int_of_float` is unspecified out of range. The value is a small integer
   in every caller (display, electron count).
2. `Gto::read_one` returns `GtoError::ShortLine` where OCaml's `String.sub`
   raises `Invalid_argument` on a line shorter than two bytes.
3. `GtoError::Io` where OCaml lets `End_of_file` escape from inside a shell.
4. `Molecule::of_xyz_string*` and `of_xyz_file*` return
   `anyhow::Result`, `MultiplicityError`/`XyzError` where OCaml raises
   (`assert`, `failwith`, `Invalid_argument`).
5. `Multiplicity::of_int` returns `Option`, so an impossible value is
   rejected at construction.
6. `Element::of_charge` takes the already-truncated `i32`
   (`Charge::to_int`), because the OCaml argument is a `Charge` whose only
   use is `int_of_float`.
7. `<Type>::to_string` inherent methods are `Display` impls (`Charge`,
   `Element`, `Multiplicity`, `AngmomXyz`, `Molecule`, `LongBasis`), so
   `clippy::inherent_to_string_shadow_display` cannot fire and `format!`
   and `.to_string()` cannot disagree.
8. `float::of_string` differs from OCaml's: no `_` separators, no `0x` hex
   floats. Nothing in the tool layer emits them, and a file containing one
   is reported as an input error.

## Reference quirks reproduced on purpose (do not "fix")

* `Molecule.name` prints the elements in **reverse** order: `"OH2, Singlet"`,
  `"OH2 (+), Doublet"`, `"OH2 (2-), Singlet"`.
* `Basis.find` compares the **untrimmed** line with `Element::of_string`,
  ignores only `UnknownElement`, propagates EOF, and `find(Og)` returns
  immediately leaving the position at 0.
* `AngmomXyz::of_angmom` iterates `x` and `y` **downwards**, so the first
  monomial of a P shell is `z`, not `x`; `LongBasis::to_basis` keeps the
  entry equal to that first monomial, which is why the digest does not depend
  on the Cartesian count.
* `Gto::read_one`: symbol is bytes 0..2, shell size from `buffer[2..]`, `D`
  normalised to `e` in the coefficient only, blank line / `$END` / EOF all
  mean end of section.
* `Basis.to_string_general` emits the atom separator **before** the next
  nucleus header (the OCaml conses header-then-separator and reverses).
* Angstrom → bohr is a **division** by `bohr_to_angstrom` in `Point3d`, while
  `Element` radii **multiply** by `angstrom_to_bohr`. Both are post-fix
  behaviour and deliberately not unified.
* `Molecule.to_xyz` writes the molecule name as the title line;
  `of_xyz_file` ignores that line, as the reference does.
* Coordinates are printed with `%.17g`, charges with `%.17g`, Gamess blocks
  with `%3d`/`%16f` and a blank line between nuclei.
* `AngmomXyz::of_string("s")` returns the zero monomial (the empty
  monomial), matching `Angmom.Xyz.of_string`.

## Not done, and why

* `Molecule::of_zmt_file` / `Molecule::of_file` — they need `Zmatrix`, which
  is group B's module. `of_xyz_string*` / `of_xyz_file*` are complete.
* `Pseudo` (group B) and the `Spd` shell representation: out of scope for
  group A, and this crate compiles without them.
* EZFIO input/output: `LibraryBasis` models the `ao_expo`/`ao_coef` round
  trip but no file is read or written; that belongs to the EZFIO layer.
* SHA-256 verification of `basis_file_sha256`, and `pseudo`-aware digests:
  no dependency, and the corpus rows do not need it (the prototype that
  generated the corpus reproduces the four ECP rows without the pseudo).
* Non-UTF-8 basis-library files: `BufRead` lines go through `String`, so a
  file that is not UTF-8 is an I/O error rather than a parse error.
* The `md5` unit-test JSON parser is hand-written (~150 lines) to avoid a
  `serde_json` dependency for one file format.

## Where the plan / `CONTRACTS.md` was wrong or silent

1. **The escaping rule is missing from `CONTRACTS.md` §3**, which lists only
   two "traps". The third one decides the digest of any atom containing
   punctuation, and had to be reversed from sexplib0 v0.17.0 and confirmed
   against a compiled sexplib0 v0.14 oracle. The actual predicate is: quote
   iff the atom is empty, contains any of `"()\;`, any byte in
   `0x00..=0x20` or `0x7f..=0xff`, or an **adjacent** `#|` / `|#` pair — an
   atom like `"a#b|c"` or `"a|b#c"` prints **unquoted**. An earlier
   reading of the plan ("a `#` anywhere before a `|`") is wrong.
   Escapes are `\"`, `\\`, `\n`, `\t`, `\r`, `\b`, and otherwise `\` plus
   three **decimal** digits per byte (`é` → `"\195\169"`).
2. **`Element::of_string`'s error payload is the normalised symbol**, not the
   input: a shadowed match arm in `Element.ml:21` makes `of_string "qq"`
   report `"Element Qq unknown"`. The plan described it as echoing the input.
3. **The two fixture values in `CONTRACTS.md` §3** (`H2/cc-pvdz
   44c24554…`, `HCN/aug-cc-pvdz c4997d5a…`) are md5s of *fixture EZFIO
   files*, not of anything reproducible from the installed basis library;
   they are listed next to the 27 corpus rows, which invites asserting them.
   They are not asserted. The reproducible reference is the corpus plus the
   `cc-pvdz` fluorine literal from `test_basis.ml`.
4. The plan's module split gives `Angmom`, `GaussianPrimitive` and
   `Multiplicity` to group B while group A's types embed them. Group A
   therefore defines minimal local versions in `gto.rs` / `molecule.rs`
   (`Angmom`, `AngmomXyz`, `GaussianPrimitive`, `Multiplicity`) and a plain
   `pub type NuclNumber = i32`. The merge should replace these with group
   B's modules; nothing in group A refers to any group B file.

## Generation

`src/element.rs`'s `ALL` table, its `Display`/`long_name` texts and the
charge, covalent-radius, van der Waals-radius and mass match arms were
generated from a dump of the reference build
(`~/OpenCode/tmp/p3a/gen/element.tsv`, produced by a small OCaml driver
linking the real `Element` module), then checked by the independent pins in
`tests/test_group_a_element.rs`. The generator is scratch, not in the repo;
the table is the artifact and the tests pin it.
