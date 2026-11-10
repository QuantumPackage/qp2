# Reference EZFIO payloads

Golden payloads for `qp-ezfio`, produced by the **reference clients** and
committed as text. They are the oracle for `tests/golden.rs`; nothing here was
written by this crate.

Reference: the on-disk contract of `src/ezfio.rs` and decision D6 (compare
**decompressed** payloads: a `.gz` byte stream depends on the compressor
version, level and mtime, so comparing it is meaningless). The files below are
therefore the decompressed contents of the array files, and the plain contents
of the scalar files. `tests/common/mod.rs` gzips an array payload with `flate2`
before installing it in a temporary tree, and never routes a reference payload
through the crate's writer.

## Producers

| directory | producer | version | format |
|---|---|---|---|
| `fortran/` | `external/ezfio/src/libezfio_util.irp.f` (IRPF90), linked into `libezfio_irp.a` | GNU Fortran (Homebrew GCC) **16.2.0** | `I20`, `E24.15`, `L1`, `A` |
| `ocaml/` | `external/ezfio/Ocaml/ezfio.ml` | OCaml **5.3.0** (`~/qp2/external/opampack/opamroot/qp2-ocaml5`) | `%20d`, `%24.15e`, raw strings |

Both were built from `~/qp2` @ `be4b71c1` (`symmetry`). The OCaml object files
carry magic `Caml1999Y035`, which matches the 5.3.0 compiler; the compiler
prints `5.3.0` from `ocamlopt -version`.

## Regenerating

The generators are `src/gen_fortran.f90` and `src/gen_ocaml.ml`. Both write one
throw-away EZFIO tree; the fixtures are copied from it with `gunzip`, never
re-compressed.

```sh
export PATH="$HOME/qp2/external/opampack/opamroot/qp2-ocaml5/bin:$PATH"

# Fortran: `ezfio_mkdir` needs the directory to exist.
cd tests/fixtures/src
gfortran -I"$HOME/qp2/lib" gen_fortran.f90 "$HOME/qp2/lib/libezfio_irp.a" \
         -lz -o gen_fortran
mkdir -p "$HOME/OpenCode/tmp/ftree" && ./gen_fortran "$HOME/OpenCode/tmp/ftree"

# OCaml: `Ezfio.set_file` creates the directory, so it must NOT exist.
ocamlfind ocamlopt -package unix -linkpkg -I . \
  "$HOME/qp2/external/ezfio/Ocaml/ezfio.ml" gen_ocaml.ml -o gen_ocaml
./gen_ocaml "$HOME/OpenCode/tmp/otree"
```

(`-I .` is needed: without it the vendored compiler does not resolve the
`ezfio.cmx` it just built for `gen_ocaml.ml`.) Then, for every committed
`<producer>/<group>/<key>.txt`, compare it with `gunzip -c
<tree>/<group>/<key>.gz`; scalar keys compare with the plain file of the same
name. The corpus in this directory was regenerated this way on 2026-11-05 with
the versions above: **all 29 payloads matched byte for byte**.

## Layout

```
fortran/                # bytes from libezfio_util.irp.f
  nuclei/nucl_label     # registry-known character*(32), dims (nucl_num)
  oracle/...            # one key per element type / shape
  oracle3/short_padded  # padded scalar, read-only case
ocaml/                  # bytes from ezfio.ml
  nuclei/nucl_label     # same values, written raw (no blank padding)
  oracle/...
reference-bugs/
  ocaml/logicals.txt    # `write_bool_array` output; see below
src/                    # the two generators
```

Array keys under `oracle/` (dims follow the rank line):

| key | type | dims | purpose |
|---|---|---|---|
| `reals` | `double` | `(6)` | signed zeros, 16-digit literal that does not round-trip |
| `redge` | `double` | `(6)` | two/three-digit exponents, smallest subnormal, `pi` |
| `r2d` | `double` | `(2,3)` | rank 2, Fortran column-major flattening |
| `ints` | `integer` | `(4)` | `I20`, negative value |
| `i3d` | `integer` | `(1,2,3)` | rank 3 |
| `i8s` | `integer*8` | `(3)` | values beyond 32 bits |
| `logicals` | `logical` | `(3)` | `L1` in the 24-column array slot |
| `chars` | `character*(8)` | `(3)` | blank padding; leading blanks are data |
| `empty` | `double` | `(3,0)` | zero extent, header only |
| `rscalar`, `iscalar`, `lscalar`, `cscalar` | scalars | | plain text, one field |
| `nuclei/nucl_label` | `character*(32)` | `(3)` | registry-known width |
| `oracle3/short_padded` | `character*(32)` scalar | | Fortran writes 32 columns, readers trim |

## Deviations between the two writers (and from plan §4)

Measured, not inferred; each is pinned by a test.

1. **`real` is lossy and the two writers lose different digits.** `%24.15e`
   keeps 16 significant digits; `E24.15` keeps 15. `1234567.8901234567`
   decodes as `1234567.890123457` from the OCaml payload and
   `1234567.89012346` from the Fortran one. A Fortran payload therefore cannot
   be rewritten into an OCaml payload that is byte-identical to one written
   from the original 16-digit literal (`tests/roundtrip.rs` asserts every line
   except the documented one).
2. **Fortran scalar files are not padded to the 24-column array slot.**
   `ezfio_write_in` / `_do` / `_lo` / `_ch` write `I20`, `E24.15`, `L1` and `A`
   only; the array writers first fill `buffer(i)(1:24)`, which leaves four
   blanks after an `I20` integer and 23 after an `L1` logical. Plan §4's
   "every element in a 24-character field" describes arrays only. The crate
   models this with `FieldKind::{Scalar, ArrayElement}`.
3. **Characters are not padded/truncated to 24.** The OCaml client writes the
   string verbatim; the Fortran client pads/truncates to the declared
   `character*(n)` length (8 or 32 here), which the *variable* knows and the
   file does not. The crate pads to the registry width under
   `ByteLayout::Fortran` and only for registry-known keys; unknown keys are
   written raw under both layouts, so `fortran/oracle/chars` is read but its
   padding is not reproduced (the test pins our raw bytes explicitly).
4. **The OCaml logical array writer is broken.** `print_bool` writes `"T\n"`
   through a `%s\n` format, so every element is followed by a blank line:
   `T\n\nF\n\nT\n\n`. The reference reader itself fails on the second element
   (`fortran_bool_of_string ""`). The crate never reproduces this; the buggy
   file is kept under `reference-bugs/` and
   `reference_ocaml_logical_array_is_rejected` asserts that we reject it
   instead of silently reading one element.
5. **Dims line.** Both clients write the rank as `%3d\n`, then each extent as
   `%20d`; the OCaml client terminates each extent with a space (line ends in
   `"6 \n"`), the Fortran `30(I20,X)` does not (`"...0"` then newline). There
   is **no extra blank line** after the dims line in either layout — the
   `\n` before the first element is the one that ends the dims line.
6. **Logicals in the reference bug file are the only non-layout difference**
   listed here; when read normally (Fortran payloads, or our own writes) a
   logical is one field per line, `T`/`F`, padded on the right to the slot.

## Known reference-only attributes

`~/qp2` carries an untracked `src/qmcchem` symlink, so the installed reference
configs expose two keys that the tracked tree does not
(`qmcchem/ci_threshold`, `qmcchem/do_frozen_core`). They are deliberately not
imported into `src/keys.rs`; `tools/gen_keys.py --cross-check` reports them as
the only accepted difference.
