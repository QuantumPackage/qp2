# qp-ezfio — QP2 tool-layer API over `ezfio-core`

Status: current. The data plane was delegated to `ezfio-core` on 2026-11-10
(branch `ezfio-unify`, spec `~/OpenCode/tmp/ezfio-unify/spec.md`, core commit
`72af154c`, the gitlink the Fortran/Python/OCaml layers of the same build
use). The P1 notes, when this crate had its own reader and writer, are
summarised under **History** below.

## What it computes

`qp-ezfio` is the typed Rust API the QP2 tools use to read and write an EZFIO
directory: one directory per **group**, one file per attribute, arrays
gzip-compressed. The crate owns:

- the `EzFio` handle: `open`/`create`, groups, keys, `has`, `remove`, the
  `.version` marker, and the registry's runtime dimension resolution;
- the typed registry (`KeySpec`, `Registry`, `StaticRegistry`, `DimExpr`),
  which checks the **element type on reads and writes** and the **declared
  shape on writes**, for known keys only (unknown keys are a plugin's
  business and stay unchecked);
- the `Array<T>` container: extents plus column-major (first index fastest)
  data, matching the on-disk order.

It owns **no byte code**. Every payload is formatted and parsed by
`ezfio_core::ops`; reads and writes preserve the reference tolerance
documented in `CONTRACTS.md` §1 (upper/lower-case exponents, the implied
three-digit `E24.15` overflow form, `T`/`.TRUE.` logicals). `ByteLayout`,
`with_layout`, the field formatters and the payload parsers are gone.

**One writer.** Every write passes `WriterDialect::Ocaml` (the private
`WRITE_DIALECT` in `src/ezfio.rs`, never exposed on the handle): the tools
replace the reference OCaml client, and the golden harness compares EZFIO
state reference-vs-Rust byte for byte. Core's `Canonical` dialect (Python
`%24.15E`, trimmed characters, `T\n` scalar logicals) is *not* what this
crate emits; before the swap the byte-parity gate proved the old
`ByteLayout::Ocaml` writer and `WriterDialect::Ocaml` emitted identical
payloads for every type, rank 1–3.

## Where it lives

| file | content |
|---|---|
| `src/lib.rs` | crate docs: the format as specification, delegation, precision/units, registry semantics, example; `#![deny(missing_docs)]`, `#![forbid(unsafe_code)]` |
| `src/array.rs` | `Array<T>`: extents + column-major data, checked construction, index helpers |
| `src/layout.rs` | `ElementType` (`Character`/`Integer`/`Integer64`/`Logical`/`Real`), the `Element` trait and its four `ezfio_core` conversions (`Scalar` and `ArrayData`, both ways). No formatter, no parser |
| `src/ezfio.rs` | `EzFio` handle: open/create, `.version`, groups/keys/has/remove, scalar and array get/set over `ezfio_core::ops`, registry checks, `DimResolver` |
| `src/registry.rs` | `KeySpec`, `Registry` trait, `StaticRegistry` (binary search), `DimExpr` parser/evaluator |
| `src/keys.rs` | generated table: **332 attributes** from the 43 tracked `src/*/EZFIO.cfg` |
| `src/error.rs` | `Error` (`Io`, `Core`, `NotEzFio`, `NotSet`, `Kind`, `Malformed`, `TypeMismatch`, `Shape`, `InconsistentArray`, `IndexOutOfBounds`, `Dim`, `Rank`, `DimensionTooLarge`, `InvalidName`) |
| `tools/gen_keys.py` | regenerates `keys.rs`; `--check` and `--cross-check` modes |
| `tests/fixtures/` | 29 reference-produced payloads, both generators, `README.md` (producers, commands, deviations) |
| `tests/common/mod.rs` | fixture loading, gzip install, temp tree, read-back of written bytes |
| `tests/golden.rs` | 14 read + byte-exact write tests per element type; Fortran payloads guard the reader and their write side is compared with the OCaml twin (deviation E5) |
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

Results (2026-11-11, rustc/cargo 1.92.0, laptop):

- `qp-ezfio`: **38 unit + 14 golden + 7 round-trip + 1 doctest = 60 passed, 0
  failed**, ~0.2 s of test time (the pre-swap count was 59: the eight
  `tests/parity.rs` cases are gone, the golden suite gained one);
- clippy `-D warnings`: clean; `cargo fmt --check`: clean;
- `cargo test --workspace`: green (qp-cli, qp-input, qp-types, qp-zmq
  unaffected).

Fixture provenance was re-verified before the swap: rebuilding both corpora
from `tests/fixtures/src` with gfortran 16.2.0 and OCaml 5.3.0 reproduced
**all 29 committed payloads byte for byte**.

## Complexity and cost

Let `N` be the number of array elements and `B` the payload size (bytes).
All gzip work is inside `ezfio_core`.

- `get_array`: time O(B + N). `ops::array_info` decompresses the whole stream
  to parse the two header lines, then `ops::read_array` decompresses it again
  for the elements, so the payload is inflated **twice**; peak memory is one
  array plus one decompressed text buffer (the first is dropped before the
  second call). A concurrent write between the calls makes the shapes
  disagree and the read fails with the core's `BadDims` — no tool writes an
  attribute it is reading.
- `set_array`: time O(B), extra memory one owned `ArrayData` copy plus the
  core's full text image (~25 bytes per real element, ~21 per integer)
  before deflate. Container bytes are deterministic (mtime 0).
- Scalars: O(size); the core writes a temp file in the group directory and
  renames it, so a reader never sees a half-written field.
- `get`/`get_array`/`DimResolver::resolve` open the payload once to surface an
  OS-level failure as `Error::Io` (the core folds "cannot open" and "cannot
  parse" into one `IoError` and drops the `io::Error`); the core re-opens it.
- Registry lookup: O(log 332); each `DimExpr` is parsed per call (O(expression
  length)) and each named reference costs one small file read.
- No allocation per element beyond the core's own string building; the crate
  adds none. No benchmark was run because the spec asks for no performance
  claim, and the delegation moves the hot path into `ezfio-core`.

## Behaviour deltas of the delegation

Each is deliberate, documented at the code site, and covered by a test.

- **Malformed payloads** report `Error::Core("Error reading data in <path>.")`:
  the core has one error surface and does not carry the offending line number
  the local parser used to put in `Error::Malformed`. `Malformed` survives for
  the crate's own shape guards (rank 0, negative extent).
- **Character records are trimmed on both ends** (`ezfio_core`'s
  `str::trim`, which is what the reference OCaml reader did). The previous
  local reader kept leading blanks; `"  spaced"` now reads back `"spaced"`.
  Writers still store characters verbatim (OCaml behaviour).
- **Logical word forms** `true`/`false` are no longer accepted; `T`/`t`/
  `.TRUE.`/`F`/`f`/`.FALSE.` are. No EZFIO writer emits the words.
- **Records past the declared element count are ignored** instead of rejected
  as trailing garbage; a short file is still an error.
- **Rank is capped at 10** (`MAX_RANK`, the core's `dims(10)` ABI buffer):
  a payload with a deeper rank cannot be read back. The previous local reader
  accepted up to 31. EZFIO has no array above rank 3.

## Not done, and why

- **No commits.** `git commit` is denied by this session's policy; the PI
  runs it after review. Staging is prepared, not committed.
- **Full generated registry is not here.** P1 gets a static table generated
  from the tracked configs; P2 owns the codegen/dual-emission registry.
- **No `qmcchem` keys.** `~/qp2` has an untracked `src/qmcchem` symlink, so the
  installed reference exposes 2 extra attributes; only tracked configs were
  imported (`tools/gen_keys.py --cross-check ... --allow-missing qmcchem`
  passes).
- **Non-finite reals have no fixture**, only measured unit-test values
  (`nan`/`inf`/`-inf`, which both the OCaml client and the core emit).
- **No key-name escaping beyond validation.** Names are rejected if empty, `.`,
  `..`, leading `.`, or containing `/`, `\`, NUL; EZFIO has no escaping layer.
- **The core is not patched.** The double decompression of `get_array` and the
  lost line number are `ezfio-core` API shapes (`array_info` returns shape
  only, `EzfioError` carries text); fixing them means changing the submodule,
  which this task does not do.

## What the original P1 plan (§4) got wrong

Measured against the original C client (`libezfio_util.irp.f`, deleted from
`external/ezfio` by the ezfio-rs switch) and `ocaml/ezfio.ml`;
`tests/fixtures/README.md` has the byte-level table.

1. "Each element in a 24-character field … logical `L1` (` T ` / ` F `)
   right-justified" — the C array writer emits `T` + 23 blanks (left
   justified); the OCaml writer emits `T` and, for arrays, an extra `\n` after
   every element (a bug its own reader cannot parse). See the
   `reference-bugs/ocaml/logicals.txt` fixture.
2. "character padded/truncated to 24" — the C client pads to the declared
   `character*(n)` (8 or 32 here), OCaml writes raw. `n` is a property of the
   Fortran variable, invisible in the file; this crate writes raw (OCaml).
3. "`E24.15`" is 15 significant digits and `%24.15e` is 16; the plan does not
   say `real` is **lossy**, nor that Fortran→OCaml conversion cannot recover
   the dropped digit. It is documented in the crate and pinned by tests.
4. Scalar vs array slot: the C scalar writers are `I20`/`E24.15`/`L1`/`A`
   with no blank padding; only array elements fill the 24-column slot. The
   single writer keeps the two forms apart.
5. §4's dims line is right for the C client (`30(I20,X)` emits no trailing
   blank: the trailing `X` is not processed after the last list item), while
   the OCaml client writes a trailing space after every extent. The reader
   accepts both; the writer emits the OCaml form.
6. §4 lists scalars as `double`/`logical`/`character`; `integer`/`integer*8`
   scalars exist too and both use `I20`. `integer*8` and `integer` share the
   same on-disk representation, so the crate stores both in `i64`.

## History

P1 (2026-11-05, branch `w1/p1-ezfio`, merged) implemented the format in this
crate: `ByteLayout::{Ocaml, Fortran}`, the field formatters and parsers, the
gzip payload code and `tests/parity.rs`. All of it was deleted by the
`ezfio-unify` swap except `Array<T>`, the registry and the handle. The
Fortran-layout writer had no bin user (only two test call sites), so its
removal changed no tool output; `ByteLayout::Ocaml` was already the default
every bin ran with.
