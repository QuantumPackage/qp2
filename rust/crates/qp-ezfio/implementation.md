# P1 — `qp-ezfio`, native Rust EZFIO format

Status: **T1.1 and T1.2 done**, plan DoD met locally. Port worktree, branch
`w1/p1-ezfio`. **Nothing is committed**:
this session's permission policy denies `git commit` (see the commands at the
end).

## What was built, and where

`rust/crates/qp-ezfio/` only; `rust/Cargo.toml` is untouched (workspace deps
were pre-declared). `rust/Cargo.lock` changed because `tempfile` moved from a
dev-dependency to a normal one (the atomic write path) — generated, not
hand-edited.

| file | content |
|---|---|
| `src/lib.rs` | crate docs: format, precision/units, registry semantics, example; `#![deny(missing_docs)]`, `#![forbid(unsafe_code)]` |
| `src/array.rs` | `Array<T>`: extents + column-major (first index fastest) data, checked construction, index helpers |
| `src/layout.rs` | `ByteLayout::{Ocaml, Fortran}`, `ElementType`, `FieldKind::{Scalar, ArrayElement}`, `Element` for `f64`/`i64`/`bool`/`String`, real formatters/parsers |
| `src/ezfio.rs` | `EzFio` handle: open/create, `.version`, groups/keys/has/remove, scalar and gzip-array get/set, atomic writes, shape checks, `DimResolver` |
| `src/registry.rs` | `KeySpec`, `Registry` trait, `StaticRegistry` (binary search), `DimExpr` parser/evaluator |
| `src/keys.rs` | generated table: **332 attributes** from the 43 tracked `src/*/EZFIO.cfg` |
| `src/error.rs` | `Error` (`NotSet`, `Kind`, `TypeMismatch`, `Shape`, `Malformed`, `Rank`, `Name`, `DimensionTooLarge`, `Io`) with file/line context |
| `tools/gen_keys.py` | regenerates `keys.rs`; `--check` and `--cross-check` modes |
| `tests/fixtures/` | 29 reference-produced payloads, both generators, `README.md` (producers, commands, deviations) |
| `tests/common/mod.rs` | fixture loading, gzip install, temp tree, read-back of written bytes |
| `tests/golden.rs` | 13 read + byte-exact write tests per element type/producer, signed zero, reference-bug rejection |
| `tests/roundtrip.rs` | 7 round-trip/interop/runtime-dimension tests |

Registry policy: advisory. Unknown keys (plugins) are unchecked; known keys
reject a contradicting Rust type and a contradicting shape; a dimension
expression that references an attribute not set yet resolves to `None` and the
shape check is skipped (needed for `pseudo/pseudo_dz_k` written with `(3,0)`
before its scalars exist).

## Build and test

```sh
cd <port-checkout>/rust
cargo test -j4 -p qp-ezfio
cargo clippy -j4 --all-targets -p qp-ezfio -- -D warnings
cargo fmt -p qp-ezfio -- --check
cargo test -j4 --workspace
```

Results (2026-11-05, rustc/cargo 1.92.0, macmini):

- `qp-ezfio`: **38 unit + 13 golden + 7 round-trip + 1 doctest = 59 passed, 0
  failed**, ~0.1 s of test time;
- clippy `-D warnings`: clean; `cargo fmt --check`: clean;
- `cargo test --workspace`: green (qp-cli, qp-input, qp-types, qp-zmq unaffected).

Fixture provenance was re-verified: rebuilding both corpora from
`tests/fixtures/src` with gfortran 16.2.0 and OCaml 5.3.0 reproduced **all 29
committed payloads byte for byte**.

## Complexity and cost

Let `N` be the number of array elements and `B` the payload size (bytes).

- `get_array`: time O(B + N), memory O(N · size_of::<T>()) for the data plus two
  64 KiB buffers and one reused line; the gzip stream is decoded incrementally.
- `set_array`: time O(B), extra memory O(1); gzip deflate dominates. `mtime(0)`
  and no embedded file name make the bytes deterministic for a given `flate2`
  version (D6 compares decompressed payloads anyway).
- Scalars: O(size); `set` writes a temp file in the group directory then
  renames, so a reader never sees a half-written file.
- Registry lookup: O(log 332); each `DimExpr` is parsed per call (O(expression
  length)) and each named reference costs one small file read.
- Known cost: each `f64` formatting allocates two `String`s (`scientific` then
  the field). Integer/logical/string fields do not allocate per element beyond
  the reused line. No benchmark was run because the plan asks for no
  performance claim.

## Not done, and why

- **No commits.** `git commit` is denied by this session's policy; the PI must
  run the commands below. Staging was also reset.
- **Full generated registry is not here.** P1 gets a static table generated
  from the tracked configs; P2 owns the codegen/dual-emission registry.
- **No `qmcchem` keys.** `~/qp2` has an untracked `src/qmcchem` symlink, so the
  installed reference exposes 2 extra attributes; only tracked configs were
  imported (`tools/gen_keys.py --cross-check ... --allow-missing qmcchem`
  passes).
- **Non-finite reals have no fixture**, only measured unit-test values
  (`nan`/`inf`/`-inf` for OCaml, `NaN`/`Infinity`/`-Infinity` for Fortran).
- **No key-name escaping beyond validation.** Names are rejected if empty, `.`,
  `..`, leading `.`, or containing `/`, `\`, NUL; EZFIO has no escaping layer.

## What plan §4 gets wrong

Measured against `libezfio_util.irp.f` and `ocaml/ezfio.ml`:
`tests/fixtures/README.md` has the byte-level table.

1. "Each element in a 24-character field … logical `L1` (` T ` / ` F `)
   right-justified" — the Fortran array writer emits `T` + 23 blanks (left
   justified); the OCaml writer emits `T` and, for arrays, an extra `\n` after
   every element (a bug its own reader cannot parse). See the
   `reference-bugs/ocaml/logicals.txt` fixture.
2. "character padded/truncated to 24" — Fortran pads to the declared
   `character*(n)` (8 or 32 here), OCaml writes raw. `n` is a property of the
   Fortran variable, invisible in the file; the crate uses the registry width
   under the Fortran layout and writes raw otherwise.
3. "`E24.15`" is 15 significant digits and `%24.15e` is 16; the plan does not
   say `real` is **lossy**, nor that Fortran→OCaml conversion cannot recover
   the dropped digit. It is documented in the crate and pinned by tests.
4. Scalar vs array slot: Fortran scalar writers are `I20`/`E24.15`/`L1`/`A`
   with no blank padding; only array elements fill the 24-column slot. Modelled
   by `FieldKind`.
5. §4's dims line is right for Fortran (`30(I20,X)` emits no trailing blank:
   the trailing `X` is not processed after the last list item), while the
   OCaml client writes a trailing space after every extent. The reader accepts
   both; the writer distinguishes them.
6. §4 lists scalars as `double`/`logical`/`character`; `integer`/`integer*8`
   scalars exist too and both use `I20`. `integer*8` and `integer` share the
   same on-disk representation, so the crate stores both in `i64`.

## PI commands (commit denied here)

```sh
cd <port-checkout>
git add rust/Cargo.lock rust/crates/qp-ezfio/Cargo.toml \
        rust/crates/qp-ezfio/src rust/crates/qp-ezfio/tools
git commit -m "rust(w1/p1-ezfio): native EZFIO reader/writer with typed registry"
git add rust/crates/qp-ezfio/tests/fixtures
git commit -m "rust(w1/p1-ezfio): reference-produced golden payloads and generators"
git add rust/crates/qp-ezfio/tests
git commit -m "rust(w1/p1-ezfio): golden byte, round-trip and interop tests"
git add rust/crates/qp-ezfio/implementation.md
git commit -m "rust(w1/p1-ezfio): implementation notes"
```

Merge check: the branch changes only `rust/crates/qp-ezfio/**` plus the
generated `rust/Cargo.lock`; the only cross-agent risk is a `Cargo.lock`
conflict, resolved by regenerating it after the merge.
