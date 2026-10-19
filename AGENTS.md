# Quantum Package 2 — agent guide

Repo root = this directory. `VERSION` = 2.3.1. Determinant-driven
quantum-chemistry suite: **Fortran + IRPF90** code generation for the compute
part, **OCaml** for the CLI / EZFIO bindings / ZMQ task server, **EZFIO** for
the data model, **Ninja** for the build.

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
   `ezfio_interface.irp.f`, `IRPF90_temp/`, `build.ninja` — is generated.
   Edit the source, then build; never patch generated output.
4. Never edit generated files: root `build.ninja`, `config/build.ninja`,
   `config/qp_create_ninja.pickle`, `data/executables`, `ocaml/*.native`,
   `ocaml/qp_edit.ml`, `ocaml/Input_*.ml`, `ocaml/Qptypes.ml`, `ocaml/ezfio.ml`,
   `lib/`, `share/`.
5. **Ask the user before any git commit.**
6. Match `conventions.rst`: no tabs, ≤80 columns, single-quoted Fortran
   strings, `BEGIN_DOC`/`END_DOC` on every provider/subroutine/function,
   LaTeX between `$…$` in doc blocks, program name == file basename.

## 2. Layout

```
src/<module>/          75 core modules (see §4 for the anatomy)
plugins/local/<name>/  28 in-tree plugins; installed = symlinked into src/
ocaml/                 OCaml: qp_run, qp_create_ezfio, qp_edit, TaskServer, EZFIO bindings
scripts/compilation/   qp_create_ninja (build graph), qp_path.py, read_compilation_cfg.py
scripts/module/        module_handler.py (NEED genealogy, clean/tidy), create_executables_list.sh
scripts/ezfio_interface/ ei_handler.py, qp_edit_template
scripts/utility/       qp_path.py and misc
bin/                   hand-written qp_* commands — NOTE: `bin/` is in .gitignore,
                       a new file there needs `git add -f`
data/                  basis sets, pseudos; generated `data/executables` (name → binary path)
etc/*.rc               fragments sourced by the tracked `quantum_package.rc`
config/*.cfg           compiler/build configs: gfortran (+_mpi, _mkl, _openblas,
                       _armpl, _avx, _avx2, _macos, _debug), ifort 2019/2021
                       (avx/avx2/sse4/rome/xHost/debug/_mpi), nvfortran,
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
ninja                          # root: build every module, then the OCaml tools
ninja                          # inside a module dir: build that module + OCaml tools
ninja -f $QP_ROOT/config/build.ninja module_<name>   # build one module from anywhere
make                           # == source quantum_package.rc; ninja
```

- `./configure -c config/<cfg> [-g nvidia|intel|none]` writes the root
  `build.ninja`, `config/build.ninja`, `etc/00.qp_root.rc`, and the
  `src/gpu_arch` symlink. Only run it when the user asks.
- Before all that, `configure` checks: `irpf90`, `ninja`, `cmake`, ZeroMQ
  (`-lzmq`) + f77zmq, TREXIO (`-ltrexio`), GMP, zlib, `ocamlc`, the python
  packages `docopt` and `resultsFile`, SimpleGPU (`-lsimple_gpu-cpu`), and
  BLAS/LAPACK through the compilation config. Missing dependencies are
  installed with `./configure -i <pkg>` (or `-i all`), where `<pkg>` is one of:
  armadillo, bats, bse, cmake, docopt, f77zmq, gmp, ninja, ocaml, openorb,
  qmckl, simple_gpu, trexio, zeromq, zlib (plus the variants `trexio-nohdf5`
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
- Every build first re-runs `qp_create_ninja update` (ninja target
  `dummy_target`), so new files, new modules and edited `NEED` files are
  picked up automatically. Nothing is ever edited by hand in the build files.
- A program `program foo` living in `src/m/foo.irp.f` is built as
  `src/m/IRPF90_temp/foo`; `data/executables` maps `foo` → that path, and
  `qp_run foo <dir.ezfio>` resolves it there.
- Run one program: `qp_run <program> <ezfio_dir>`; in `qpsh` the argument
  is omitted (current file). `qp_run -s <program> <ezfio>` starts a remote
  slave; `USR1` triggers a clean stop (`qp_stop`).

## 4. Anatomy of a module

A module is a directory in `src/` containing:

- `NEED` — one module name per line: at build time each needed module is
  **symlinked** into this directory, and all its `*.irp.f` are compiled in
  (so grep inside a module can show code from other modules: use `git grep`).
- `README.rst` — included in the Sphinx docs.
- `EZFIO.cfg` — control parameters (see §6); generates
  `ezfio_interface.irp.f` and the OCaml bindings.
- `*.irp.f` — IRPF90 sources, including `program` blocks. File basename ==
  executable name, so `program fci` must live in `fci.irp.f`.
- `*.template.f`, `*.include.f` — IRPF90 codelet templates (see §5).
- `*.ezfio_config` — hand-written EZFIO group description (arrays with
  computed shape, e.g. `src/determinants/spindeterminants.ezfio_config`).
- `FC`, `FCFLAGS`, `LIB` — optional extra compiler/linker flags for this
  module (e.g. `src/zmq/LIB` = `-lf77zmq -lzmq`, `src/trexio/LIB` = `-ltrexio`).

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
- `IRP_IF MPI` / `IRP_ELSE` / `IRP_ENDIF`, `IRP_IF MPI_DEBUG`, `IRP_IF ZMQ_PUSH`,
  `IRP_IF MACOS` are the conditional blocks in use (older branches also use
  `IRP_IF SET_NESTED`).
- `BEGIN_TEMPLATE … SUBST [ A, B ]` + `a1; b1 ;;` lines + `END_TEMPLATE`
  unrolls a block with `$A`/`$B` substituted (see `src/cipsi/selection.irp.f`).
  Separately, a `*.template.f` file is a codelet template whose `$subroutine`
  placeholders are instantiated by the `.irp.f` that includes it
  (see `src/determinants/h_apply.template.f`).
- IRPF90 is invoked by ninja (pool `depth = 1`, i.e. serialized) as
  `irpf90 -I <needed modules> $IRPF90_FLAGS`; **`IRPF90_FLAGS` comes from the
  compilation config** (`config/*.cfg`): `--ninja --align=32` there, plus
  `-DMPI` for `gfortran_mpi.cfg`, while `*debug*.cfg` / `travis.cfg` add
  `--assert`. Don't hardcode the flags: read the config in use.
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
  `ocaml` (bindings + `default:` required for scalars).
- Types come from the OCaml side; list them with `ei_handler.py list_supported_types`.
- Adding a parameter and building is enough: `ezfio_interface.irp.f` and the
  `qp_edit`/`Input_*` OCaml files are regenerated. Never commit them.

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

Compute programs are run through `qp_run <name> <ezfio>`:

- `scf` (RHF/ROHF/UHF, DIIS/damping in `scf_utils`); `ks_scf` (Kohn-Sham),
  `rs_ks_scf` (range-separated). DFT dispatch is by string providers
  `exchange_functional` / `correlation_functional` → `energy_x`, `energy_c`,
  `potential_{x,c}_{alpha,beta}_ao`; a new functional is one file
  `src/functionals/<name>.irp.f` providing `energy_x_<name>` and friends
  (recipe in `src/README.rst` + `programmers_guide/new_ks.rst`).
- `fci` = selected CI / CIPSI driver: `run_cipsi` (deterministic PT2) or
  `run_stochastic_cipsi` when `do_pt2` with stochastic PT2; selection is
  distributed over ZMQ tasks, stopping on `n_det_max`, `pt2_max`,
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
- Parallelism: `qp_run` starts an OCaml `TaskServer`; masters spawn
  in-process slaves and accept remote ones (`qp_run -s`, `qp_tunnel`).
  The ZMQ address is `$EZFIO/work/qp_run_address`; slave code paths are
  selected by `is_zmq_slave`. MPI builds use `config/gfortran_mpi.cfg`
  (`-DMPI`) plus `qp_mpirun` / `qp_srun`. Thread counts come from `etc/local.rc`
  (`QP_NTHREADS_DAVIDSON`, `QP_NTHREADS_PT2`, `OMP_NUM_THREADS`, `QP_NIC`).

## 9. Recipes

- **New keyword**: add a block to the module `EZFIO.cfg`, then use it as a
  provider (or `ezfio_get_<module>_<param>`). No registration step.
- **New program**: new `src/<module>/<name>.irp.f` with `program <name>`;
  build and it appears in `data/executables`.
- **New subroutine/provider in an existing module**: add it to an existing
  `*.irp.f`; nothing to declare. Keep the file's theme.
- **New module**: `mkdir src/<name>`, add `NEED`, `README.rst`, sources;
  the next `ninja` finds it. `module_handler.py print_descendant <name>`
  shows what it pulls in.
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

## 10. Tests and CI

- `cd tests && qp_test [-a] [TEST]` (bats). `-a` runs all installed modules,
  otherwise only the current directory; `TEST` is a bats number/name.
  Helpers: `tests/bats/common.bats.sh` (`eq` for float compare, `test_exe`,
  `run_only_test`), inputs in `tests/input/`.
- Modules may carry their own `NN.<name>.bats` (e.g. `src/hartree_fock/10.hf.bats`).
- CI: `.github/workflows/compilation.yml` (build + `qp_test -a`),
  `.github/workflows/configuration.yml` (exercises `./configure -i <pkg>` for
  each dependency on a bare ubuntu-24.04) and `.drone.yml`.

## 11. Pitfalls checklist

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
