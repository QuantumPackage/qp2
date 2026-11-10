# `qp-input` — porting notes for the `qp_edit` handlers

**Status: P4 complete.** The generic engine, the seven handlers, the driver
CLI and the C1–C7 contract tests are implemented on this branch and pass the
golden gate (`tools/golden-harness/golden_harness.py verify --only check
--contract post-fix`). This file records what the OCaml does, where the port
deliberately differs, and which reference quirks are preserved on purpose.

Reference: the **post-fix** tree, `qp_edit-fix` at `15e725aa`. The pre-fix
`~/qp2` tree differs in three ways that matter here (no `set ?apply`, no
backup on failure, `%16.8f` Angstrom quantisation); do not use it as the
reference for behaviour.

## Module map

| module | role |
|---|---|
| `schema.rs`, `loader.rs` | generated TOML format, validation, discovery |
| `types.rs` | declared-type table, `Value`, domain checks |
| `sexp.rs` | sexplib-subset reader and the `(key value)` record reader |
| `ocaml_compat.rs` | OCaml printers: `string_of_float`, `valid_float_lexem`, `MO_label` |
| `engine.rs` | `get`/`set`, section bounding, the retry path, `--check-config` helpers |
| `handlers/{ao,ao_extra,mo}_basis.rs` | the three basis sections |
| `handlers/{electrons,nuclei_by_hand,determinants_by_hand}.rs` | the three small sections |
| `handlers/qptypes.rs` | local port of `Qptypes.ml` validation (`MO_number`, `Nucl_number`, …) |
| `bins/qp_edit/src/main.rs` | driver: `TASKS`, editor, backup, banner, exit codes |
| `bins/qp_edit/src/inspect.rs` | `--list` / `--describe` / `--check-config` (not in the reference) |

## The config artifact

`data/qp_edit_config/` is gitignored (`data/.gitignore:4`) and is **not** on
the branch; regenerate it before building the tests or running the gate. The
wave-3 generator is `$HOME/OpenCode/tmp/w3p4/gen_config.py`. Result: 39 files
(`registry.toml` + 38 modules), `SCHEMA_VERSION = 1`, 353 fields (189
editable), 3 Fortran-formula defaults.

Discover it as `--config <dir>` > `$QP_EDIT_CONFIG` >
`<exe_dir>/../data/qp_edit_config`.

Known generator bugs (the handlers hardcode their groups/keys from the OCaml,
so the gate is unaffected, but the artifact is wrong for a future generic
consumer):

* `ao_extra_basis.toml` says `ezfio_dir = ["ao_basis"]` and duplicates the
  `ao_basis` fields. `Input_ao_extra_basis.ml` uses the **`ao_extra_basis`
  group with `ao_extra_*` keys** (`Ezfio.get_ao_extra_basis_ao_extra_num`,
  etc.). The Rust handler uses the correct group/keys, not the TOML ones.

## Driver semantics (`ocaml/qp_edit.ml` post-fix)

`get s` = `make_header s ^ Rst_string.to_string (to_rst x)` when `read ()`
returns `Some x`; `""` when it returns `None` (the section vanishes from the
document); `""` plus `Info: <Sys_error>` on stderr when the read raises
`Sys_error`.

`set ?apply str s`, faithfully ported:

1. find `header` in `str`; no match ⇒ **silent no-op**;
2. `index_begin = idx + len header`;
3. `index_end` = next `"=="` searched **from** `index_begin + len header + 1`
   (this skips the header's own `====` underline), else `len str`;
4. `Rst_string.of_string` then the section's `of_rst`;
5. on any exception: print `Info: Read error in <title>`, call `of_rst` a
   second time and **ignore the result**; a second exception escapes (C7);
   (wave S2.4: the port's retry is silent — no `Info:` line — and the
   second escape is printed as `qp_edit: <title>: <message>`, exit 1);
6. write only when `apply`; check mode passes `apply:false`.

`run` order: directory-exists check (`<dir> does not exists`, typo
preserved), make `work/`, delete `work/qpstop` and `work/qpkill`, AO reorder,
AO-extra reorder, `-n` (`update_ndet`), `-s` (`extract_states`), build the
document, run `$EDITOR` (default `vi`), re-read and `set` every task.

Section order is the OCaml `tasks` list, **not** the registry order:
`nuclei_by_hand, ao_basis, ao_extra_basis, electrons, mo_two_e_ints,
two_body_rdm, ao_two_e_ints, davidson_keywords, mrci, scf_utils, mo_one_e_ints,
pseudo, ezfio_files, dressing, utils, extra_nuclei, ao_one_e_ints,
becke_numerical_grid, symmetry, perturbation, bitmask, mu_of_r,
utils_trust_region, mo_optimization_utils, hamiltonian, density_for_dft,
nuclei, casscf_cipsi, dft_keywords, cipsi, basis, trexio, determinants,
utils_cc, mol_properties, mo_basis, determinants_by_hand`.

`determinants_by_hand` **must be last**. Its `of_rst` treats every line after
its own `\nDeterminants` header as coefficient/determinant data. A section
that followed would leave its title line in the body (the bound cuts at the
first `==` of the next underline, not before the title), and the determinant
parser would misread it. As the last task it receives the tail of the
document, exactly as in the reference. `TASKS` in
`bins/qp_edit/src/main.rs` is the reference list, and a unit test pins the
membership and the last position.

Failure handling: `Failure`/`Invalid_argument` print the 80-column `=-`
banner and the message, restore `work/backup.tar`, exit 1; `Assert_failure`
prints the banner and exits 2; a missing `EZFIO_DIR` prints the help then
`Failure("EZFIO_DIR is missing")` (help twice, exit 2); `-r` with `-w` is
`Failure "read and write options are incompatible"` (exit 1, after the
backup). `-h` prints the help and **continues** — it is not an exit.

## Generic `get`/`set`

`get` materialises absent keys from the TOML `default` before rendering
(C3) — a real EZFIO write, allowed in check mode. Values are printed with the
OCaml printers, never with Rust's `Display`: `string_of_float` is `%.12g`
plus the OCaml dot rule, integers are 32-bit, choices keep the canonical
spelling from the type table, arrays are Fortran (column-major) order.

`set` parses the section body as a sexp of `(key value)` pairs (every line
containing `=`, `=` replaced by a space, each line parenthesised), checks the
domain with `types::TypeSpec`, and writes with `ByteLayout::Ocaml`. A
conversion failure (the sexplib `Error:  (...)` of the reference) yields
`None` (section skipped, exit 0); anything else takes the retry path above.
Wave S2.4: the port prints the conversion failure as
`qp_edit: parse error in <section>: <message>` (exit 0) and the second
escape as `qp_edit: <section>: <message>` (exit 1) — plain text, no
sexplib reporter, no `Info:` line.

**The read path never writes the `*_md5` stamps** (D4/C2). The OCaml
`Input_ao_basis.read` writes the computed stamp on success and `md5("None")`
on any exception; the port computes the same `Some`/`None` and drops the
write. `Input_mo_basis.read_ao_md5` writes `mo_basis/ao_md5` when absent and
then fails on mismatch; the port compares without writing and keeps the
failure.

## The seven handlers

### `ao_basis` / `ao_extra_basis`

Same structure. `read` is wrapped in `try … with _ -> None` in the reference,
so **any** exception (including a domain assertion) renders the section as
`""`, not as a fatal error. `ao_prim_num_max` has a Fortran-formula default
the config cannot evaluate; the handler supplies `max` of `ao_prim_num`.

`to_rst` is byte-exact and pinned by a test. The `set` arm of the reference
is `() (* TODO *)`, so both handlers' `set` returns `Ok(())`: only `get`
matters for the corpus.

### `mo_basis`

`mo_coef_to_string` prints blocks of **four** MOs (`if i+4 < mo_num`, then a
final 1–4 block), blocks joined by `"\n\n"`, each a header
`"  #%15d %15d %15d %15d"` and one `"  %3d %15.10f …"` line per AO. The
reference function is named `print_five` and has a `| 5 ->` branch, but the
caller only reaches it with 1–4: the branch is dead. `mo_num = 0` hits
`assert false`.

Read side materialises `mo_label` (`"None"`), `mo_class` and `mo_occ` when
absent (C3); `read_mo_occ` is gated on `has_mo_basis_mo_label`, not on
`mo_occ`, and clamps into `[0,2]`. `set` is `Ok(())` like the reference.

### `electrons`

Multiplicity from `Multiplicity.of_alpha_beta`, one line
`"Spin multiplicity is %s."`. `elec_num` comes from the formula default; the
handler supplies `alpha + beta`.

### `nuclei_by_hand`

xyz text in **Angstrom**, `"  %-3s  %3s %s"` per atom with `%.17g` charge
and coordinates. Coordinates are stored in bohr and converted by **division**
by `Units.bohr_to_angstrom` (not multiplication by the reciprocal: the two
differ in the last bit, and goldens are byte-compared). C5 forbids the
pre-fix `%16.8f` Angstrom quantisation; the port keeps `%.17g` and a
`+5e-9` bohr perturbation survives `-c` and `-w`/`-r`.

`set` writes the molecule: `nucl_num`, `nucl_label`, `nucl_charge`,
`nucl_coord`, in that order (group `nuclei`). Point charges live in the
generic `Nuclei` section (`is_periodic`, `point_charges`) and are not touched
here.

### `determinants_by_hand`

The riskiest section and the last task. `read ~full:false` is what `get`
calls; field order is the materialisation order (`n_det_qp_edit`, `n_det`,
`n_int`, `bit_kind`, `expected_s2`, `psi_coef`, `psi_det`, `n_states`,
`state_average_weight`, `mo_label`).

* `read_only = (n_det_qp_edit <> n_det)`; when false, `read_psi_coef` /
  `read_psi_det` create the canonical HF pair — the seven files C3 expects on
  `hcn`; when true they read the `_qp_edit` pair instead.
* `n_int` is computed as `(mo_num - 1)/64 + 1`, not taken from the data.
* `write ~force:false` skips the `psi_*` writes when `n_det >
  n_det_qp_edit`, and is the only writer of the `_qp_edit` pair. Check mode
  never writes, so the pair must not appear (C4). The pre-fix corpus contains
  `psi_{det,coef}_qp_edit.gz`; the post-fix goldens do not.
* `of_rst` uses its own header filter that additionally requires
  `line[0] = ' '`, and consumes the **rest of the document** — hence the
  ordering rule above.
* Coefficients are tab-joined per state with `Float.to_string` and `"0."`
  past the end; each determinant is two spaces plus the
  `Determinant.to_string` lines, prefixed by two more.

## Deliberate deviations

1. **Check mode creates no backup** and does not run the reference's
   `at_exit remove_backup`. The reference `at_exit` always removes
   `work/backup.tar`, even in `-c`, so a stray pre-existing archive would be
   deleted by a check; the port leaves the tree alone (C2 is about writes,
   but "check mode writes nothing" is the honest reading).
2. **`EzFio::open` requires `.version`.** The reference mostly no-ops on a
   directory that is not EZFIO. A `qp_edit` run on such a directory now fails
   loudly instead of silently doing nothing.
3. **clap diagnostics** for an unknown option differ from OCaml's `Getopt`
   text; the exit code is the same (2). Documented, not a behavioural bug.
4. **`-h` with no `EZFIO_DIR`** prints the help twice (once for `-h`, once
   from the missing-argument branch) and exits 2, matching the reference's
   "print and continue" `help ()`.
5. **`version` is not a subcommand.** `--version` is a usage error (the
   reference has no version command either).
6. **Inspection modes** (`--list`, `--describe --module`, `--check-config`,
   `--config`) are a maintenance interface on the same binary, dispatched
   before `EditArgs`. They have no reference counterpart; they replace the
   deleted prototype binary.
7. **Assert banner text.** The reference prints
   `Assert error in file $QP_ROOT/ocaml/<file>, line <n>, character <c>`;
   Rust cannot recover the OCaml source location, so the banner prints the
   message only. Exit code 2 is preserved.
8. **`create_backup`** builds the tar in a temporary file outside the tree
   and copies it to `<dir>/work/backup.tar`, instead of writing
   `.backup.tar` in the current directory and `mv`-ing it. Observable
   difference: a failed `mv` no longer leaves `.backup.tar` in the user's
   cwd. Errors are swallowed in both.
9. **`main` returns `ExitCode`** rather than `exit`-ing, so the temporary
   document and the backup are cleaned up deterministically (the reference
   relies on `at_exit`).
10. **type-table storage** keeps values in their declared width: EZFIO
    integers are 32-bit (`i64` on the Rust side, range-checked), reals are
    written with `%24.15e`, logicals as `T`/`F`.

## Preserved quirks (do not "fix")

* `update_ndet` has `if n_det_new > n_det_new` — always false, so the
  "Requested n_det should be less than %d" guard never fires and `-n` larger
  than the stored count runs the copy loop out of bounds (`Invalid_argument`
  from the array index). The port drops the dead guard and keeps the
  unconditional `Reducing n_det to <n>` on stdout.
* `read_nucl_num` returns the stored `Nucl_number.get_max()`.
* `MO_class.of_string` accepts only the bare keyword.
* `MO_number` / `Nucl_number` / `N_int_number` / `Det_number.of_int` assert
  `x >= min`, warn above 10 000 (50e9 for determinants), and default `max`
  to the stored EZFIO value; `Bit_kind ∈ {1,2,4,8}`; `MD5` asserts length 32
  and all bytes `< b'g'`.
* `MO_label.of_string` prints the rejected string to stdout before failing
  (multi-line `Failure` payload).
* (pre-S2) sexplib's `Error:  (...)` messages were approximated, not
  byte-identical: the snippets and the two-space prefix were reproduced,
  the type-path wording was not. Wave S2.4 deleted the reporter: the
  messages are plain Rust text and C7-A pins the full line.
* A section whose header was deleted by the user is a silent no-op.

## Tests and gate

* `crates/qp-input`: 119 unit tests (printers against independently computed
  reference strings, sexp reader, every handler's `get`, the `set` retry
  path; the sexp `fail_msg`/`snippet` tests are gone since wave S2.4).
* `bins/qp_edit`: 7 unit tests (wave S2: the help byte count and banner
  tests are gone — the help is clap-generated), 17 C1–C7 contract tests
  (including the S2.4 exact-line pins), 19 config tests, 3 interactive
  tests (no-op editor, edit reaches EZFIO, failed read-back
  restores the backup). 655 workspace tests in total.
* Gate (debug binary copied to `$FAKE/ocaml/qp_edit`):

  ```
  python3 rust/tools/golden-harness/golden_harness.py verify \
      --only check --contract post-fix --qp-root $FAKE
  ```

  Expected: the five small fixtures `0 created, 0 compared byte for byte`,
  `hcn 7 created, 7 compared byte for byte`, no failures, empty transcripts,
  exit 0.

## What is not done

* `ao_basis`, `ao_extra_basis` and `mo_basis` `set` are no-ops because the
  reference `set` arm is `() (* TODO *)`. Editing those sections in the
  interactive editor therefore does nothing; the sections are read-only.
* The config-artifact bugs above are not fixed: they need a change in
  `scripts/ezfio_interface/`, not in this crate.
