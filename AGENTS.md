# Quantum Package 2 — agent guide

Repo root = this directory. `VERSION` = 2.3.1. Determinant-driven
quantum-chemistry suite: **Fortran + IRPF90** code generation for the compute
part, **Rust** for the CLI / EZFIO bindings, **EZFIO** for the data model,
**GNU make** for the build (default), **Ninja** as the opt-in alternative
(`./configure -c <cfg> -b ninja`).

Authoritative in-tree references (prefer them over guessing):

- `docs/source/programmers_guide/programming.rst` — architecture, modules, algorithms
- `docs/source/programmers_guide/conventions.rst` — enforced style rules
- `docs/source/programmers_guide/ezfio.rst` — `EZFIO.cfg` syntax
- `docs/source/programmers_guide/index_providers.rst` — generated list of every provider
- `docs/source/users_guide/` — one page per `qp_*` command
- `INSTALL.rst`, `src/README.rst` (how DFT functionals are wired)
- `RELEASE_NOTES.org` — what changed in each released version

Online: <https://quantum-package.readthedocs.io> (manual) and
<https://quantumpackage.github.io/qp2> (project site).

## 1. Rules that prevent damage

1. **Never run `./configure` unless the user asks for it** — not even `-h`.
   Before parsing its options it executes `git submodule init/update` and
   `git pull` in `external/qp2-dependencies`. That needs the network, can
   hang, and if interrupted mid-checkout it leaves `external/{irpf90,ezfio}`
   empty. The help text is the comment block at the top of `configure`.
2. `external/{ezfio,irpf90,qp2-dependencies}` are **git submodules**. Never edit them.
3. In a module directory, the only files to create/edit are:
   `NEED`, `README.rst`, `EZFIO.cfg`, `*.irp.f`, `*.template.f`, `*.include.f`,
   `*.ezfio_config`, and the per-module `FC`/`FCFLAGS`/`CFLAGS`/`LIB`.
   Everything else there — symlinks to needed modules, `.gitignore`,
   `ezfio_interface.irp.f`, `IRPF90_temp/`, the per-module `Makefile`
   (make mode) or `build.ninja` (ninja mode) — is generated.
   Edit the source, then build; never patch generated output.
4. Never edit generated files: the build manifest (`config/Makefile` in
   make mode, `config/build.ninja` in ninja mode), the per-module
   `Makefile`/`build.ninja`, `irpf90.make`, `config/qp_create_ninja.pickle`,
   `data/executables`,
   `data/qp_edit_config/`, `data/list_element.txt`, `rust/target/`, `lib/`,
   `share/`.
5. **Ask the user before any git commit.**
6. Match `conventions.rst`: no tabs, ≤80 columns, single-quoted Fortran
   strings, `BEGIN_DOC`/`END_DOC` on every provider/subroutine/function,
   LaTeX between `$…$` in doc blocks, program name == file basename.

## 2. Layout

```
src/<module>/          77 core modules (see §4 for the anatomy)
plugins/local/<name>/  29 in-tree plugins; installed = symlinked into src/
rust/                  the Rust tool layer (six binaries, built into
                       `rust/target/release/` on PATH);
                       `scripts/compilation/build_rust_tools.sh` drives it
scripts/compilation/   qp_create_ninja (build graph), build_rust_tools.sh, qp_path.py, read_compilation_cfg.py
scripts/module/        module_handler.py (NEED genealogy, clean/tidy), create_executables_list.sh
scripts/ezfio_interface/ ei_handler.py, qp_edit_config.py
scripts/utility/       qp_path.py and misc
bin/                   hand-written qp_* commands — NOTE: `bin/` is in .gitignore,
                       a new file there needs `git add -f`
data/                  basis sets, pseudos; generated `data/executables` (name → binary path)
etc/*.rc               fragments sourced by the tracked `quantum_package.rc`
config/*.cfg           compiler/build configs: gfortran (+_mkl, _openblas,
                       _armpl, _avx, _avx2, _macos, _debug), ifort 2019/2021
                       (avx/avx2/sse4/rome/xHost/debug), nvfortran,
                       cray_{gfortran,intel}, bull, flang, travis
tests/                 bats suite (`tests/bats/`, inputs in `tests/input/`)
docs/, man/            Sphinx sources, man pages
lib/, include/, share/ tracked but empty placeholders that the build fills
external/              git submodules: ezfio, irpf90, qp2-dependencies
                       (dev-stable also carries armadillo-code, OpenOrbitalOptimizer)
```

## 3. Build and run

```bash
source quantum_package.rc      # or: ./bin/qpsh   (interactive QP shell)
make                           # root: build every module, then the Rust tool layer
make                           # inside a module dir: build that module + Rust tool layer
make -f $QP_ROOT/config/Makefile module_<name>   # build one module from anywhere
```

The build system is GNU make by default. The opt-in ninja alternative is
`./configure -c config/<cfg> -b ninja`; from then on the same three
invocations work with `ninja` instead of `make` (`ninja -f
$QP_ROOT/config/build.ninja module_<name>`).

- `./configure -c config/<cfg> [-g nvidia|intel|none] [-b make|ninja]`
  writes the build manifest — `config/Makefile` plus one `Makefile` per
  module in make mode (default), or `config/build.ninja` with `-b ninja` —
  `etc/00.qp_root.rc` (which exports `QP_BUILD_SYSTEM`), and the
  `src/gpu_arch` symlink. Only run it when the user asks.
- The root `Makefile` is a tracked dispatcher: it re-runs the build-graph
  generator (`update`) and delegates to `config/Makefile` in make mode, or
  to `ninja -f config/build.ninja` in ninja mode; in make mode it runs the
  build inside a shell that has sourced `quantum_package.rc` (the per-module
  compiler commands rely on that environment, e.g. `CPATH`/`C_INCLUDE_PATH`
  pointing at `include/` for the C headers).
- Before all that, `configure` checks: `irpf90`, the selected build tool
  (`make` or `ninja`), `cmake`, TREXIO (`-ltrexio`), zlib, the python
  packages `docopt` and `resultsFile`, SimpleGPU (`-lsimple_gpu-cpu`), and
  BLAS/LAPACK through the compilation config. Missing dependencies are
  installed with `./configure -i <pkg>` (or `-i all`), where `<pkg>` is one of:
  armadillo, bats, bse, cmake, docopt, libcint, ninja, openorb,
  qmckl, simple_gpu, trexio, zlib (plus the variants `trexio-nohdf5`
  and `qmckl-intel`). `irpf90` and `ezfio` are *submodules*, not `-i` targets —
  the tab-completion in `etc/configure.rc` misleadingly lists them.
  (Note: INSTALL.rst no longer offers the `conda env create -f qp2.yml` path,
  since `qp2.yml` is not tracked in the repo.)
- **Intel `ifx` is not recommended** (verified by the maintainer,
  2026-10-19); the supported Intel configs are the `ifort_2019_*` /
  `ifort_2021_*` ones.
- `-g amd` in the script points at `plugins/local/gpu_amd`, which is *not*
  in-tree (it ships with a plugin repository); the in-tree GPU backends are
  `gpu_x86`, `gpu_nvidia`, `gpu_intel`.
- Every build first re-runs `qp_create_ninja update` (the `update` target
  of the root `Makefile`; the ninja graph does the same through its
  `dummy_target` edge), so new files, new modules and edited `NEED` files
  are picked up automatically. Nothing is ever edited by hand in the build
  files.
- A program `program foo` living in `src/m/foo.irp.f` is built as
  `src/m/IRPF90_temp/foo`; `data/executables` maps `foo` → that path, and
  `qp_run foo <dir.ezfio>` resolves it there.
- Run one program: `qp_run <program> <ezfio_dir>`; in `qpsh` the argument
  is omitted (current file). `USR1` triggers a clean stop (`qp_stop`).

## 4. Anatomy of a module

A module is a directory in `src/` containing:

- `NEED` — one module name per line: at build time each needed module is
  **symlinked** into this directory, and all its `*.irp.f` are compiled in
  (so grep inside a module can show code from other modules: use `git grep`).
- `README.rst` — included in the Sphinx docs.
- `EZFIO.cfg` — control parameters (see §6); generates
  `ezfio_interface.irp.f`, the EZFIO config/default tables, and the
  `data/qp_edit_config/` editor entries.
- `*.irp.f` — IRPF90 sources, including `program` blocks. File basename ==
  executable name, so `program fci` must live in `fci.irp.f`.
- `*.template.f`, `*.include.f` — IRPF90 codelet templates (see §5).
- `*.ezfio_config` — hand-written EZFIO group description (arrays with
  computed shape, e.g. `src/determinants/spindeterminants.ezfio_config`).
- `FC`, `FCFLAGS`, `LIB` — optional extra compiler/linker flags for this
  module (e.g. `src/trexio/LIB` = `-ltrexio`).

Root modules — nothing NEEDs them, they are the build entry points:
`dummy` (aggregates all the core modules, used for a development build),
`casscf_cipsi`, `mrci`, `mp2`, `mo_optimization`, `trexio`,
`generators_full_tc`, `ao_extra_basis`, `utils_periodic`.
Program modules such as `fci` and `ccsd` are not roots: they are NEEDed
only by `dummy` and `tools` (on `master`, `tools` needs `fci` only, so
`ccsd` is still a root there).
`module_handler.py print_descendant <module>` prints the genealogy.

## 5. IRPF90 in QP2

```fortran
BEGIN_PROVIDER [ double precision, my_var, (0:mo_num) ]
  implicit none
  BEGIN_DOC
  ! What this is, in LaTeX if useful: $E = \sum_i c_i$.
  END_DOC
  PROVIDE mo_num
  ...
END_PROVIDER
```

- An entity is computed lazily, exactly when first used. Force with
  `PROVIDE name` (mandatory at the top of a routine before using an array
  provider); an entity never referenced is never computed.
- `TOUCH x` after assigning `x` by hand (outside its provider): marks
  everything that depends on `x` as stale and re-provides the entities in
  the current scope. `SOFT_TOUCH` does the same but defers the recomputation
  (use it when nothing else in the routine uses them). Typical pattern:
  `read_wf = .True.` / `SOFT_TOUCH read_wf`.
- `FREE x` releases the memory of a provider (used for big buffers).
- `ASSERT(cond)` is only compiled in if `irpf90 -a` is passed (i.e. if the
  compilation config has `--assert`); otherwise it is silently dropped.
- `IRP_IF MACOS` is the conditional block in use (older branches also use
  `IRP_IF SET_NESTED`; the MPI conditional blocks existed until the
  2026-11 single-node cleanup).
- `BEGIN_TEMPLATE … SUBST [ A, B ]` + `a1; b1 ;;` lines + `END_TEMPLATE`
  unrolls a block with `$A`/`$B` substituted (see `src/cipsi/selection.irp.f`).
  Separately, a `*.template.f` file is a codelet template whose `$subroutine`
  placeholders are instantiated by the `.irp.f` that includes it
  (see `src/determinants/h_apply.template.f`).
- IRPF90 is invoked serialized (one call per module) as
  `irpf90 -I <needed modules> $IRPF90_FLAGS` — by the per-module `Makefile`
  in make mode (where `$IRPF90_FLAGS` is the config's flags minus the
  `--ninja` token), or by the ninja graph in ninja mode; **`IRPF90_FLAGS`
  comes from the compilation config** (`config/*.cfg`): `--ninja
  --align=32` there, while `*debug*.cfg` / `travis.cfg` add `--assert`.
  Don't hardcode the flags: read the config in use.
  Useful tools inside a module:
  `irpf90_entities` (table `file : type :: name (dims)`), `irpman <entity>`,
  `irpf90 -t <entity>` (what breaks when it is touched), `tags`.

## 6. EZFIO

An EZFIO "file" is a **directory**: one subdirectory per EZFIO *group*
(= module name), one plain-text file per parameter, plus `.version` and
`work/` (large scratch; delete before archiving). Everything is atomic units
(Bohr, Hartree).

`EZFIO.cfg` entry:

```
[threshold_generators]
type: Threshold
doc: Threshold on the generators
interface: ezfio,provider,ocaml
default: 0.999
size: (determinants.n_det)        # arrays only; note the dot: module.value
```

- `interface` may be `ezfio` (Fortran I/O API `ezfio_get_<module>_<param>` /
  `ezfio_set_…`), `provider` (an IRPF90 provider with the parameter name),
  `ocaml` (editable in `qp_edit`; `default:` required for scalars).
- Types are the generator's fancy types (historically OCaml's); list them
  with `ei_handler.py list_supported_types`.
- Adding a parameter and building is enough: `ezfio_interface.irp.f` and the
  `data/qp_edit_config/` tables are regenerated. Never commit them.

Command line: `qp set <module> <param> <value>`, `qp get …`, `qp has …`,
`qp set_file <dir>`, `qp edit <dir>` (interactive editor), `qp reset`.
The current file is `$EZFIO_FILE`.

Parameters that drive most calculations:

| module | parameter | default | meaning |
|---|---|---|---|
| determinants | `n_det_max` | 1e6 | stop CIPSI above this |
| determinants | `n_states` | 1 | number of states |
| determinants | `read_wf` | False | start from the stored wave function |
| determinants | `s2_eig` | True | enforce a spin eigenstate |
| determinants | `threshold_generators` | 0.999 | generator prefilter |
| determinants | `selection_factor` | 1. | selection size multiplier |
| perturbation | `do_pt2` | True | compute PT2 / use stochastic CIPSI |
| perturbation | `pt2_max` | 1e-4 | stop CIPSI below this \|PT2\| |
| perturbation | `pt2_relative_error` | 2e-3 | stochastic PT2 error target |
| perturbation | `h0_type` | EN | denominator (EN/CFG/HF) |
| davidson_keywords | `threshold_davidson` | 1e-10 | eigensolver convergence |
| davidson_keywords | `n_states_diag` | 4 | vectors in Davidson (doubled on retry) |
| davidson_keywords | `n_det_max_full` | 1000 | below this, full diagonalization/LAPACK |
| davidson_keywords | `disk_based_davidson` | True | store vectors on disk |
| scf_utils | `thresh_scf` | 1e-10 | SCF convergence |
| scf_utils | `scf_algorithm` | DIIS | DIIS or damping |
| scf_utils | `mo_guess_type` | Huckel | initial MOs |
| dressing | `thresh_dressed_ci` | 1e-5 | dressed-CI convergence |

MO classes (`Core`, `Inactive`, `Active`, `Virtual`, `Deleted`) are set with
`qp_set_mo_class`; `qp_set_frozen_core` is the frozen-core shortcut
(`-c "[1-n]" -a "[n+1-mo_num]"`).

## 7. Data model worth knowing

- Determinants: `integer(bit_kind) psi_det(N_int,2,psi_det_size)` with
  `bit_kind = 8`, `bit_kind_size = 64`, `N_int = (mo_num-1)/64+1`,
  dimension 2 = alpha/beta bit strings; `psi_coef(psi_det_size, N_states)`;
  `psi_det_size ≥ N_det` (arrays are allocated once and over-sized).
  Bit tricks (`popcnt`, `iand`, `ieor`, `shiftl`, `not`) come from the
  `bitmasks` module — always `use bitmasks`.
- `psi_det_sorted`, `psi_det_alpha_unique`/`psi_det_beta_unique`, and the
  bilinear-form arrays `psi_bilinear_matrix_{rows,columns,values,order}`
  exist to speed up H·u and the selection; selection code builds/uses them.
- Integrals live in open-addressing hash maps (`src/utils/map_module.f90`)
  keyed by a single integer: `two_e_integrals_index` folds four indices
  through triangular pair indices `p = min(i,k) + max(i,k)(max(i,k)-1)/2`,
  exploiting the 8-fold permutation symmetry. Sets: `ao_integrals_map`,
  `mo_integrals_map` (+ `_erf` variants for range separation), built by the
  `*_two_e_ints` modules; `all_mo_integrals` is the usual trigger provider.
  The MO map has a small window cache and an optional Cholesky path
  (`cholesky` module, `do_mo_cholesky`) that stores 3-index tensors.
- H·u and the connected/selected space come from `h_apply` codelets
  (`src/determinants/h_apply.template.f`) driven by the selector modules
  (`selectors_full`, `selectors_cassd`, `dressing`).

## 8. Programs and algorithms

Compute programs are run through `qp_run <name> <ezfio>`. Their stdout is a
human summary; the numbers to consume are in the run's JSON document,
`<EZFIO>/json/NNNNN.json` (§9).

- `scf` (RHF/ROHF/UHF, DIIS/damping in `scf_utils`); `ks_scf` (Kohn-Sham),
  `rs_ks_scf` (range-separated). DFT dispatch is by string providers
  `exchange_functional` / `correlation_functional` → `energy_x`, `energy_c`,
  `potential_{x,c}_{alpha,beta}_ao`; a new functional is one file
  `src/functionals/<name>.irp.f` providing `energy_x_<name>` and friends
  (recipe in `src/README.rst` + `programmers_guide/new_ks.rst`).
- `fci` = selected CI / CIPSI driver: `run_cipsi` (deterministic PT2) or
  `run_stochastic_cipsi` when `do_pt2` with stochastic PT2; selection is
  parallelized over OpenMP tasks, stopping on `n_det_max`, `pt2_max`,
  `variance_max`, `correlation_energy_ratio_max`. `pt2` recomputes PT2 on a
  stored wave function; `dressing`/`davidson_dressed` add the dressed
  (PT2-renormalized) Hamiltonian; variants for CSF and TC (bi-orthogonal)
  spaces exist.
- Excited states: `cis` / `cisd` / `casscf` guesses, then `fci` with
  `n_states > 1` (`state_following`, `s2_eig`).
- Post-HF and orbital work: `mrci`, `ccsd`, `lccsd`, `mp2`, `casscf`,
  `orb_opt` / `optimization` (trust region, `utils_trust_region`),
  `four_idx_transform`, `save_natorb*`, `save_one_e_dm`, `rotate_mos`,
  `swap_mos`, `sort_wf`, `truncate_wf`, `print_*`, `properties`
  (`mol_properties`: dipole, spin density, …).
- Interfacing: `fcidump`, `fcidump_pyscf`, `molden`, `export_trexio_*`,
  `import_trexio_*`, `qp_convert_output_to_ezfio`.
- Parallelism: single machine only, OpenMP threads inside each
  executable. The distributed layer (task server, in-process and remote
  slaves, `qp_run -s`, `qp_tunnel`, the `is_*_slave` provider) and the
  multi-node build configurations and cluster launchers were removed in
  2026-11; `mpi_rank` / `mpi_size` / `mpi_master` survive as constants
  in `src/mpi/`. Thread counts come from `etc/local.rc`
  (`QP_NTHREADS_DAVIDSON`, `QP_NTHREADS_PT2`, `OMP_NUM_THREADS`).

## 9. Output: human summary on stdout, JSON record

Two streams with different readers. Only the second is a contract.

**stdout is the human summary.** Short banner, a table that grows by one
row per iteration, the extrapolation table at the end, the answer in a box.
It is ASCII-only (0x20-0x7E) and aimed at a person watching a run: no
script should parse it, and a number meant for a machine never goes there.
`human_quiet` (`src/utils/human.irp.f`) silences the legacy provenance
dump through the `write_double` / `write_int` / `write_bool` helpers of
`src/ezfio_files/output.irp.f`; `QP_VERBOSE=1` restores it for debugging.
No machine-readable line is printed: scripts read the JSON document, never
stdout (`scripts/get_fci_conv.sh` and `scripts/qp_geom_opt.py` go through
`scripts/utility/qp_json.py`).

**The JSON document is the record, and it holds more than stdout.**
`<EZFIO>/json/NNNNN.json` (first unused number, so a restart never touches
an earlier run), `"schema": "qp3/1"`, opened once by
`json_init('<program>')` (`scf`, `ks_scf`, `rs_ks_scf`, `cis`, `cisd`,
`fci`, `pt2`, `mrci`, `casscf`, `cipsi_orb_opt`) and filled by
`call json_add('/results/<method>/<key>', value)`.
`src/json/README.rst` is the reference, `src/json/DESIGN.md` the design.

- `provenance`: `ezfio_file`, `host`, `pid`, `started`; `units`: hartree,
  bohr, second, declared once; `run`: `kind`, `timestamp`, `elapsed`.
- `input`: the scalar EZFIO inputs the calculation actually read
  (integer, double precision, logical, fixed-length string), nested as
  `input/<ezfio_group>/<ezfio_name>`. Array-valued inputs and `integer*8`
  are not hooked. Emitted by the EZFIO provider generator
  (`scripts/ezfio_interface/`), not by the call site.
- `results`: grouped by method, not by program. `scf` logs `iteration`,
  `energy`, `delta_energy`, `max_error_diis`, `level_shift`, `dim_diis`,
  `converged`, plus the per-iteration `history/*`; `cipsi` logs
  `algorithm`, `n_det`, `n_cfg` (when `s2_eig`), per-state
  `states/{energy,s2,pt2,rpt2,variance,extrapolated_energy}` plus
  `pt2_err`, `rpt2_err`, `variance_err`, and the per-iteration series under
  `history/`; `davidson` logs one `calls/NNNN/{algorithm,n_state,epoch,
  itertot,converged,wall_time,shape,iteration,energy,residual,s2}` block
  per solve; `symmetry` logs `group`, `n_irrep` and per-state
  `states/<s>/{irrep,purity}`. The history convention (one index-aligned
  array per quantity, `shape[n,nh]` with the state index fastest, and the
  `(nh)` scalars such as `n_det`, `epoch`, `memory`) is stated once in
  `src/json/DESIGN.md` §10.7; reshape the 2-D arrays with `order='F'`.
  `ks`, `casscf` and `mo_opt` have their own objects. stdout carries the
  same last-row numbers and one extrapolated energy, but the full
  per-iteration history survives only here.

Each write replaces the whole document (`.tmp.<pid>` + `fsync` +
`rename`), so the file is valid JSON after `kill -9` and always complete;
doubles round-trip bit-exactly; a non-finite value becomes `null`, meaning
"not available". `QP_JSON_MIN_INTERVAL` (default 1 s, `0` = synchronous
write after each `json_add`) bounds the staleness. After a restart the
highest-numbered file is the complete state; older files are history and
are never merged. `scripts/utility/qp_json.py` (`files`, `load`,
`load_all`, `load_last`) implements that rule and is what
`qp_extract_cipsi_data.py` and `qp_exc_energy.py` use.

**The `input` hooks live in generated per-module files.**
`src/<module>/ezfio_interface.irp.f` is written by `ei_handler.py` from
`EZFIO.cfg`; it is not versioned (`.gitignore`d) and carries the
`call json_add('/input/...')` lines. The ninja edge
`build_ezfio_interface` takes the generator scripts as implicit inputs
(`scripts/compilation/qp_create_ninja`), so `ninja` regenerates them
whenever `EZFIO.cfg` or the generator changes. If a tree shows an `input`
object with only `ezfio_files.ezfio_convention`, its interfaces predate the
generator: `rm src/*/ezfio_interface.irp.f && ninja` and run again.

## 10. Recipes

- **New keyword**: add a block to the module `EZFIO.cfg`, then use it as a
  provider (or `ezfio_get_<module>_<param>`). No registration step.
- **New program**: new `src/<module>/<name>.irp.f` with `program <name>`;
  build and it appears in `data/executables`.
- **New subroutine/provider in an existing module**: add it to an existing
  `*.irp.f`; nothing to declare. Keep the file's theme.
- **New module**: `mkdir src/<name>`, add `NEED`, `README.rst`, sources;
  the next build (`make`, or `ninja` in a ninja-configured tree) finds it.
  `module_handler.py print_descendant <name>` shows what it pulls in.
- **New plugin**: `qp_plugins create -n <name> [-r <repo>] <needed modules…>`;
  `qp_plugins install <name>` symlinks `plugins/<repo>/<name>` into `src/`.
  Never edit files through the symlink in `src/`.
- **Docs**: module docs come from `README.rst`; the provider index is
  regenerated with `cd docs/source && python3 auto_generate.py`. Building the
  docs is finicky: `auto_generate.py` reads every module's `IRPF90_man/`
  (built artifacts), so the tree must be fully `ninja`-built first — with only
  some built it aborts with "you did not run ninja at the root". Recent Sphinx
  (9.x, pinned in `docs/requirements.txt`) needs **Python ≥ 3.12** (the conda
  base here is 3.10 → use `uv venv -p 3.12`). Also note: IRPF90's generated
  provider pages emit C-domain directives whose bare identifier the modern
  C-domain parser rejects (verified against Sphinx 9.1.0 on 2026-10-20).
- **User-visible change**: add a line under the current `* Version` heading of
  `RELEASE_NOTES.org` (org-mode, one section per release).
- **Debugging a number**: `qp set` the keyword, rerun the program; every
  program is restartable from the EZFIO state. For reproducibility, scripts
  should build a fresh EZFIO (`qp create_ezfio -b <basis> <file.xyz>`) rather
  than reuse a modified one.

## 11. Tests and CI

- `cd tests && qp_test [-a] [TEST]` (bats). `-a` runs all installed modules,
  otherwise only the current directory; `TEST` is a bats number/name.
  Helpers: `tests/bats/common.bats.sh` (`eq` for float compare, `test_exe`,
  `run_only_test`), inputs in `tests/input/`.
- Modules may carry their own `NN.<name>.bats` (e.g. `src/hartree_fock/10.hf.bats`).
- CI: `.github/workflows/compilation.yml` (a ninja job and a make job, the
  make one exercising the default build; both run `qp_test -a`),
  `.github/workflows/configuration.yml` (exercises `./configure -i <pkg>` for
  each dependency on a bare ubuntu-24.04) and `.drone.yml`.

## 12. Pitfalls checklist

- Executable name = **file** basename, not the `program` name.
- `bin/` is gitignored: `git add -f bin/qp_newcommand`.
- A Fortran run is a function of the EZFIO state, not of the input file you
  passed; diff the EZFIO directory when results surprise you.
- Grepping inside `src/<module>/` also scans symlinked modules — use `git grep`.
- `src/gpu_arch` is a configure-created symlink and `gpu` NEEDs it; if it is
  missing, `module_handler` aborts on an unresolved NEED (create it with
  `./configure -g none`, or by hand: `ln -s plugins/local/gpu_x86 src/gpu_arch`).
- IRPF90 output lines are long and contain `! file: line` comments:
  `-ffree-line-length-none` is mandatory, and never put `!` inside a string
  literal.
- `PROVIDE` is needed before using array providers inside a routine; forgetting
  it gives wrong or stale values, not a compile error.
- `TOUCH` is for values assigned by side effect; using it in an iteration
  without a convergence test silently produces stale dynamics.
- `src/<module>/ezfio_interface.irp.f` is generated output: it holds the
  `json_add('/input/...')` hooks and is rebuilt whenever `EZFIO.cfg` or
  `scripts/ezfio_interface/` changes (implicit ninja inputs). An `input`
  block holding only `ezfio_files.ezfio_convention` means the interfaces
  are stale: `rm src/*/ezfio_interface.irp.f && ninja` (§9).

## 13. Code-model RAG (`make rag` / `ninja rag` / `qp rag`)

`make rag` (repo root; `ninja rag` in a ninja-configured tree) builds a
machine-readable model of the whole source tree — one `irpf90 -J` dump
per build edge of the configured build manifest (`config/Makefile` in
make mode, `config/build.ninja` in ninja mode; 25 modules, run with the
exact `$IRPF90_FLAGS` of the build) — merged into one corpus with a
SQLite FTS5/BM25 index. It lands in
`$QP_ROOT/share/rag/` (gitignored): `per-module/<M>.json` × 25,
`qp2.json` (merged model), `qp2_rag.db` (index), `merge-report.json`
(counts, collisions, applied renames). It is not part of the default
build; the build tool rebuilds it automatically when any module source
or the tool changes, and a rebuild takes a few minutes.

Requirement: an irpf90 with `-J` (≥ 2.3.3) — the in-tree
`external/irpf90` pin carries it (it is on `$PATH` through
`etc/irpf90.rc`, so `make rag` just works). `RAG_IRPF90=/path/to/bin/irpf90`
is a dev override to run the RAG against a different irpf90; with an
irpf90 that lacks `-J` the target fails with a clear error.

Interrogation (read-only, works in a `qpsh` shell):

```bash
qp rag query "PT2 selection" -k 5
qp rag query "bitmask" --kind routine
qp rag graph <entity> --depth 3
qp rag stats
qp rag build                    # == make rag / ninja rag (explicit rebuild)
```

Underneath, `qp rag …` calls `bin/qp_rag` →
`python3 scripts/rag/qp_rag.py …` (bin/ is gitignored; `git add -f` when
committing). Names that collide across modules are suffixed `__<module>`
in the corpus (e.g. `run__fci`); the `file` field of each record is the
ground truth for where it lives.
