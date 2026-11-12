# A selectable libcint engine for the AO integrals

Design note. Target branch: `libcint-ao-ints`, based on `origin/dev-stable`
(`d2b3bfd5`). This document is the spec for PR1, updated to describe the
implementation and the numbers it produced.

Author: A. Scemama. Written 2026-11.

---

## 1. Goal

Quantum Package 2 computes its AO integrals with hand-written Fortran. Add a
second engine, **[libcint](https://github.com/sunqm/libcint)**, and let the
user choose between the two with one EZFIO keyword. libcint is the default.

The choice must be a *runtime* choice in a single binary: the same `qp_run`
executable serves both, so that a calculation can be benchmarked both ways
without rebuilding.

PR1 covers the integrals that dominate everything else:

* one-electron: overlap, kinetic energy, nucleus-electron attraction;
* two-electron: the plain $1/r_{12}$ Coulomb operator `ao_two_e_integral`
  and its range-separated counterpart `ao_two_e_integral_erf`.

Everything else (pseudopotentials, point charges, Schwartz screening, the
$|\chi_i||\chi_j|$ overlap used for screening, per-atom NAI, cGTO/periodic
paths) keeps using the QP2 code.

## 2. Non-goals

* No change to the EZFIO data model besides one new logical.
* No change to the integral consumers (Cholesky, map, HF, DFT, CIPSI).
* No GPU work, no 3- and 4-centre one-electron operators, no
  gradients, no relativistic/spinor integrals.
* No removal of the QP2 code. The QP2 engine stays the reference
  implementation and the fallback for every unsupported case.

## 3. User-visible behaviour

New EZFIO parameter in the `ao_basis` group, next to `use_cgtos`:

```ini
[use_libcint]
type: logical
doc: If true, compute the AO integrals with libcint
interface: ezfio
default: True
```

* `use_libcint = True` (default): libcint computes overlap, kinetic, NAI and
  the two-electron integrals.
* `use_libcint = False`: nothing changes with respect to today.
* When the basis or the Hamiltonian is not representable in libcint, the
  corresponding integral falls back to the built-in code with a message on
  the master rank; the EZFIO flag is left alone. §7.8 lists the cases.
* `qp edit` shows the keyword; old EZFIO files that lack it get `True`
  written on first read, exactly like `use_cgtos` gets `False`.

Because the default flips, `RELEASE_NOTES.org` gets a line under the current
version. The default was flipped after the validation of §13 passed on
macOS (the platform this branch was developed and measured on); the Linux
side is exercised by the `./configure -i libcint` job added to
`.github/workflows/configuration.yml` (§5), and by the bats test in the
normal CI build.

## 4. Why libcint, and how it is pinned

libcint is a small, dependency-free C library of Gaussian integrals
(GPL-3.0, same licence family as QP2's GPL-3.0). It is the integral engine
of PySCF and is well tested for Cartesian and spherical Gaussians up to high
angular momentum.

Pin **v6.2.1** as a git submodule (`external/libcint`), like `external/ezfio`
and `external/irpf90`. No system libcint is assumed or accepted, so the ABI
(`FINT` width, presence of the `cintXe_*_cart` wrappers) is the one we
validated.

Build it **static**, with `I8=off`:

```
cmake -DCMAKE_INSTALL_PREFIX=$QP_ROOT \
      -DBUILD_SHARED_LIBS=OFF \
      -DENABLE_STATIC=ON \
      -DWITH_CINT2_INTERFACE=ON \
      -DI8=off \
      -DCMAKE_DISABLE_FIND_PACKAGE_QUADMATH=ON \
      -DCMAKE_BUILD_TYPE=Release
```

* `ENABLE_STATIC=ON` forces `BUILD_SHARED_LIBS=0` in libcint's own
  `CMakeLists.txt` (lines 167-173), so `libcint.a` is what lands in
  `$QP_ROOT/lib`.
* `I8=off` is a hard ABI requirement: `FINT` must be `int` (32-bit). With
  `I8=on` every `FINT` argument silently changes width and the Fortran
  interfaces are wrong.
* `WITH_CINT2_INTERFACE=ON` exposes the `cintXe_*_cart` entry points we call
  (`cint1e_ovlp_cart`, `cint1e_kin_cart`, `cint1e_nuc_cart`, `cint2e_cart`).
  It is already the default, but state it explicitly.
* `QUADMATH` is only linked when found; disabling the search keeps the static
  archive free of a `-lquadmath` dependency on compilers that ship it.
* Static linking avoids `LD_LIBRARY_PATH`/`DYLD_LIBRARY_PATH` problems on
  macOS and in CI.

libcint needs no BLAS or LAPACK: the `dgemm` that appears in `nm -u
lib/libcint.a` (`_CINTdgemm_NN1`, `_CINTdgemm_NT`, `_CINTdgemm_TN`) is
defined inside the archive itself (`src/fblas.c`), and the only undefined
symbols left are libc and libm ones. So `-lcint -lm` is the whole
dependency, and `-lcint` does not have to sit after the BLAS/LAPACK flags
that the compilation config already emits.

`libcint.a` is not compiled as position-independent code; static objects are
linked into the executables, never `dlopen`ed, so this is fine (verified by
the macOS build of this branch).

## 5. Provisioning

`configure` already runs `git submodule init/update` before anything else
(lines 18-20), so adding `external/libcint` to `.gitmodules` is enough to get
the sources.

Install path, mirroring the `trexio`/`simple_gpu` cases:

```bash
elif [[ ${PACKAGE} = libcint ]] ; then
        execute << EOF
          cd "\${QP_ROOT}"/external/libcint
          rm -rf build && mkdir build && cd build
          cmake -DCMAKE_INSTALL_PREFIX="\${QP_ROOT}" \
                -DBUILD_SHARED_LIBS=OFF -DENABLE_STATIC=ON \
                -DWITH_CINT2_INTERFACE=ON -DI8=off \
                -DCMAKE_DISABLE_FIND_PACKAGE_QUADMATH=ON \
                -DCMAKE_BUILD_TYPE=Release ..
          (cmake --build . -j 8 || cmake --build .) && cmake --build . --target install
EOF
```

Plus:

* `libcint` added to the `PACKAGES="..."` list used by `-i all` (line 226).
* A **mandatory dependency check** next to the TREXIO/QMCKL checks:

  ```bash
  LIBCINT=$(find_lib -L"${QP_ROOT}"/lib -lcint -lm)
  if [[ ${LIBCINT} = $(not_found) ]] ; then
          error "libcint (libcint) is not installed."
          fail
  fi
  ```

  The explicit `-L${QP_ROOT}/lib` is needed because `configure` does not
  source `etc/paths.rc`, so `LIBRARY_PATH` is not set yet.
* `etc/configure.rc` gets `libcint` in the `-i` completion list.
* `.github/workflows/configuration.yml` gets a `./configure -i libcint` step
  before the final `./configure -c config/gfortran_debug.cfg`, so the
  submodule build is exercised on a bare ubuntu-24.04; the compilation
  workflow installs it the same way as the other dependencies, since the
  `find_lib` check above is mandatory.
* `INSTALL.rst` and `AGENTS.md` §3 list of `-i` packages get `libcint`.

The new module gets a one-line `LIB` file:

```
-L$QP_ROOT/lib -lcint -lm
```

which `qp_create_ninja` appends to the end of the link command
(`scripts/compilation/qp_create_ninja`, lines 108-133). `-lcint` is a static
archive and `-lm` follows it, so the linker resolves `pow`, `erf`, `sqrt`
references from inside the archive.

## 6. Module layout

New module `src/libcint_ao_integrals/`:

| file | contents |
|---|---|
| `NEED` | `ao_basis`, `hamiltonian` |
| `LIB` | `-L$QP_ROOT/lib -lcint -lm` |
| `README.rst` | module doc, included in Sphinx |
| `libcint_interface.f90` | plain module: `bind(C)` interfaces, slot constants, ABI check, and the `libcint_block_cache` module (§9.2) |
| `libcint_data.irp.f` | providers that build `atm`, `bas`, `env`, `libcint_env_lr`, `libcint_scale`; validity checks |
| `libcint_one_e.irp.f` | `libcint_ao_one_e(op, mat)` and the `libcint_ao_overlap/kinetic/nuc` wrappers |
| `libcint_two_e.irp.f` | `libcint_ao_two_e_integral(i,j,k,l,lr)`, uses `libcint_block_cache` |
| `src/ao_two_e_ints/test_libcint.irp.f` | `program test_libcint`, prints the integrals (§13) |
| `src/ao_two_e_ints/10.libcint.bats` | engine-vs-engine regression test (§13) |

`src/ao_basis/libcint.irp.f` holds the `use_libcint` provider (copy of the
`use_cgtos` provider; the `MPI_BCAST` block that followed the master read at
the base was removed with the multi-node layer, 2026-11). It lives in
`ao_basis` because that is where the EZFIO variable is declared, and because
both integral modules already NEED `ao_basis`.

### NEED graph

```
ao_one_e_ints      -> ... libcint_ao_integrals
ao_two_e_ints      -> ... libcint_ao_integrals
libcint_ao_integrals -> ao_basis hamiltonian
ao_basis           -> nuclei basis
hamiltonian        -> ezfio_files nuclei
```

`test_libcint` lives in `src/ao_two_e_ints/` rather than in the new module:
it calls `ao_overlap`, `ao_kinetic_integrals`, `ao_integrals_n_e` and
`ao_two_e_integral`, so putting it in `libcint_ao_integrals` would create a
cycle (`libcint_ao_integrals` -> `ao_two_e_ints` -> `libcint_ao_integrals`).

No cycle: nothing in the `ao_basis`/`nuclei`/`basis`/`ezfio_files`/
`hamiltonian` closure refers to `ao_overlap` or `ao_two_e_integral`
(checked with `git grep`). `hamiltonian` is NEEDed only to read `mu_erf`
for the range-separated environment.

## 7. Mapping QP2 data onto libcint

### 7.1 Conventions

* All libcint pointers (`atm[PTR_COORD]`, `bas[PTR_EXP]`, `bas[PTR_COEFF]`,
  `PTR_RANGE_OMEGA`) are **0-based offsets into `env`**; the Fortran index is
  `offset + 1`.
* `atm` and `bas` are flat 1-based Fortran arrays of `c_int`, `natm*6` and
  `nbas*8` long. Shell `s` (1-based) uses `bas(8*(s-1)+1 : 8*s)`. A slot
  inside shell `s` is read as `bas(SLOT + BAS_SLOTS*(s-1) + 1)`, and a slot
  of atom `i` as `atm(SLOT + ATM_SLOTS*(i-1) + 1)`: the trailing `+1` turns
  the 0-based C slot offset into the 1-based Fortran index.
* QP2 is Cartesian internally even when `ao_cartesian = False` (the spherical
  transform happens in `src/ao_basis/spherical_to_cartesian.irp.f`). So we
  always call the `_cart` entry points.
* QP2's `(ij|kl)` has `i,j` on electron 1 and `k,l` on electron 2, matching
  `cint2e_cart(out, shls={i,j,k,l}, ...)`.
* libcint output layout: `out[ip + jp*ni + kp*ni*nj + lp*ni*nj*nk]` with
  `ip` the fastest index.

### 7.2 `atm`: one entry per nucleus

`natm = nucl_num`.

| slot | value |
|---|---|
| `CHARGE_OF` | `nint(nucl_charge(i))` (see the fractional-charge note below) |
| `PTR_COORD` | `20 + 3*(i-1)` |
| `NUC_MOD_OF` | `POINT_NUC` (= 1) |
| `PTR_ZETA` | `0` |

and `env(ptr+1:ptr+3) = nucl_coord(i,:)`.

`CINTg1e_nuc` multiplies by `-abs(atm[CHARGE_OF])` (`g1e.c:233`), so the
attractive sign is automatic, and `CINTnuc_mod` only reads `PTR_ZETA` for
`GAUSSIAN_NUC`: `POINT_NUC` is the explicit point-charge model and gives
the same kernel as the default `NUC_MOD_OF = 0`. The `abs` is also why a
negative `nucl_charge` cannot go to libcint: QP multiplies by
`-nucl_charge(i)`, so the charge would turn from repulsive to attractive.
libcint v6.2.1 also supports `FRAC_CHARGE_NUC`; that would be the clean way
to support non-integer `nucl_charge`, but PR1 takes the simpler route: if
any charge is negative or `abs(nucl_charge - nint(nucl_charge)) > 1e-12`,
the NAI operator falls back to QP2 with a warning (2e integrals do not
depend on the charges, so they stay on libcint).

### 7.3 `bas`: one entry per QP2 shell

`nbas = shell_num`. For shell `s`, with `l = shell_ang_mom(s)`,
`i0 = ao_first_of_shell(s)` and `np = shell_prim_num(s)`:

| slot | value |
|---|---|
| `ATOM_OF` | `basis_nucleus_index(s) - 1` |
| `ANG_OF` | `l` |
| `NPRIM_OF` | `np` |
| `NCTR_OF` | `1` |
| `KAPPA_OF` | `0` |
| `PTR_EXP` | `oexp` |
| `PTR_COEFF` | `ocoef` |

with

```
env(oexp+1 : oexp+np)   = ao_expo_ordered(i0, 1:np)
env(ocoef+1: ocoef+np)  = ao_coef_normalized_ordered(i0, 1:np) / fac(l)
```

`fac(l)` is libcint's `CINTcommon_fac_sp(l)` (§7.6). Only the coefficients
of the first |AO| of the shell are stored; the per-|AO| normalization is
applied after the block, see §7.7.

One libcint shell per QP2 shell works because a QP2 shell is exactly a
contraction with a fixed angular momentum and a fixed set of primitives;
`ao_expo`/`ao_coef` are stored per AO only for historical reasons, and the
components of a shell share them (this is checked at initialisation, §7.8).

The exponents are passed in ascending order (`*_ordered`), the order QP2's
own integrals use. libcint does not require sorted exponents
(`optimizer.c:289`, "Exponents are not required to be sorted"), but using
the sorted arrays keeps the two engines bit-comparable on the same data.

`NCTR_OF = 1` means one contraction per shell. libcint's `PTR_COEFF` block is
`NCTR_OF*NPRIM_OF` long; with `NCTR_OF = 1` it is just the primitive
coefficients.

### 7.4 `env`: one array, two variants

Layout, in doubles:

```
1..20                  header (PTR_EXPCUTOFF, PTR_RANGE_OMEGA, ...)
21 .. 20+3*natm        nuclear coordinates
next                   exponents   (nprim_total)
next                   coefficients (nprim_total)
```

`nenv = 20 + 3*natm + 2*nprim_total + 8`.

Header slots that matter:

* `env(1)` = `PTR_EXPCUTOFF` = 0. Zero means "use libcint's default"
  (`60`, i.e. essentially no exponential cutoff). We keep QP2's
  threshold-free behaviour: callers screen with `ao_two_e_integral_zero`, and
  `ao_two_e_integral` itself must return the exact integral.
* `env(9)` = `PTR_RANGE_OMEGA` = 0 for the plain Coulomb operator,
  `mu_erf` for the long-range operator.

Two read-only `env` arrays are built once, `libcint_env` (omega = 0) and
`libcint_env_lr` (omega = `mu_erf`), and shared by all threads. libcint only
*reads* `env[PTR_EXPCUTOFF]` (`g1e.c:37`, `g2e.c:57`) and
`env[PTR_RANGE_OMEGA]`; it never writes to it. Building both up front avoids
writing `PTR_RANGE_OMEGA` from inside an OpenMP region.

### 7.5 Shell ↔ AO correspondence

libcint orders the Cartesian components of a shell as
`xx, xy, xz, yy, yz, zz, ...`, which is exactly `ao_power_index` order:
`ao_power_index(nx,ny,nz) = (l-nx)(l-nx+1)/2 + nz + 1` gives `x=1, y=2, z=3`
for `l=1` and `xx=1, xy=2, xz=3, yy=4, yz=5, zz=6` for `l=2`.
QP2 fills the AOs of a shell in that order (`aos.irp.f:58-71`, comment in
`spherical_to_cartesian.irp.f:832`). Therefore

```
component ip (1..ncart(l)) of shell s  <->  AO  ao_first_of_shell(s) + ip - 1
```

No power table is needed at run time. This was verified numerically: the
d self-overlap diagonal came out `{3,1,1,1,1,1}/16 * (pi/2a)^{3/2}`, i.e.
`xx = 0.3691314831`, `xy = xz = yy = yz = zz = 0.1230438277` for `a = 1`.

### 7.6 The `CINTcommon_fac_sp(l)` factor

This is the one non-obvious piece of the mapping, so it is written out in
full.

libcint multiplies every integral by `prod_shells CINTcommon_fac_sp(l)`
(`g2e.c:55-56` for two-electron, `cint1e.c:121` for one-electron;
`g1e.c` sets `envs->common_factor = 1`, so nothing else is folded in):

```c
double CINTcommon_fac_sp(FINT l)
{
        switch (l) {
                case 0: return 0.282094791773878143;
                case 1: return 0.488602511902919921;
                default: return 1;
        }
}
```

The factor depends on the shell (`l`), not on the Cartesian component.
Consequence: if we handed libcint the plain `ao_coef_normalized`, every
integral would be off by `f(l_i) f(l_j) f(l_k) f(l_l)`, i.e. by
`0.07957747` per s shell and `0.23873241` per p shell.

Two ways to cancel it:

1. divide the coefficients we put in `env` by `f(l)` (§7.3); libcint then
   returns QP2's raw integrals, and the per-AO scale of §7.7 is the only
   correction left;
2. divide each output block by `f(l_i) f(l_j)` (1e) or
   `f(l_i) f(l_j) f(l_k) f(l_l)` (2e).

Option 1 is chosen: one division per coefficient at initialisation, and one
`env` serves both the 1e and 2e engines, so the factors cannot drift apart.

**Evidence** (probe `test_1e8.c` in `~/OpenCode/tmp/libcint-check`, one s
shell with `a = 1`, `c = 1`, one unit charge at the origin):

| quantity | libcint | raw analytic value | ratio |
|---|---|---|---|
| `<s|s>` | 0.1566642672 | 1.9687012432 | 0.0795774715 = f(0)^2 |
| `<s|T|s>` | 0.2349964007 | 2.9530518649 | 0.0795774715 |
| `<s|-1/r|s>` | -0.2500000000 | -3.1415926536 | 0.0795774715 |
| `(ss|ss)` | 0.0276945914 | 4.3733... | 0.0063325739 = f(0)^4 |

and for p and d shells the same probe gives `f(1)^2 = 0.2387324146` and
`f(2)^2 = 1` exactly, uniformly over all components of the shell.

### 7.7 Per-AO scaling `r_i`

After the `f(l)` division, libcint computes the integrals with the
contraction coefficients of the first |AO| of each shell. QP2's integrals
are in terms of
`ao_coef_normalized(i,j) = ao_coef(i,j) * ao_coef_normalization_factor(i) *
(primitive norm)`, where the normalization factor is per AO, not per shell:

```
ao_coef_normalization_factor(i) = 1 / sqrt(<chi_i|chi_i>)
```

which differs from component to component inside a shell (`d_xx` and `d_xy`
have different self-overlaps). libcint cannot represent per-component
coefficients, so we apply them *after* the block is computed:

```
r_i = ao_coef_normalization_factor(i) / ao_coef_normalization_factor(i0)
```

with `i0 = ao_first_of_shell(shell(i))`, and

```
1e:  integral_qp2(i,j)         = r_i * r_j * integral_libcint(i,j)
2e:  integral_qp2(i,j,k,l)     = r_i * r_j * r_k * r_l
                               * integral_libcint(i,j,k,l)
```

Multiplying the block by `r_i` replaces the first-|AO| coefficient by the
one of |AO| `i`: `ao_coef_normalized(i0,j) * r_i = ao_coef_normalized(i,j)`,
one factor per index (two for 1e, four for 2e).

This is exact when the raw contraction coefficients are the same for all
components of a shell, i.e. when

```
ao_coef_normalized(i,j) / ao_coef_normalization_factor(i)
```

is independent of the component `i` inside a shell. That holds when
`prim_normalization_factor` depends on `l` only (`basis.irp.f:32-44`,
`powA = (l,0,0)`; the only supported convention, 20250211) and for
`primitives_normalized = False`. It does **not** hold for a basis whose
primitive normalization is computed per component with the powers of that
component (a TREXIO file written by such a code; the old GAMESS EZFIO
convention did this and is now rejected in `aos.irp.f`). §7.8 turns this
into a runtime check rather than a version test.

### 7.8 Validity checks

`libcint_data` builds the tables and sets `libcint_available`. It is
`.False.`, with a warning printed once on the master rank, when

* `use_cgtos` is true (complex AOs);
* `is_periodic` is true (periodic path is out of scope);
* any shell is not "well formed": for every shell `s`, every AO `i` in `s`
  must be centered on `basis_nucleus_index(s)` (`ao_nucl(i)`), carry the
  Cartesian powers libcint expects at that position (`ao_power(i,:)`, whose
  `ao_power_index` must be `i - ao_first_of_shell(s) + 1`, and whose powers
  must add up to `l`), satisfy `ao_prim_num(i) == shell_prim_num(s)`, have
  the same `ao_expo_ordered(i,1:np)`, and the same
  `ao_coef_normalized_ordered(i,j)/ao_coef_normalization_factor(i)`
  (relative tolerance 1e-12). A failure here catches any file whose
  per-component normalization factors, centers or component order would
  permute the libcint blocks;
* `ao_l_max > 12` (libcint's Rys-quadrature limit);
* the basis is empty (`ao_num`, `shell_num` or `nucl_num` below 1), the
  `shell_prim_num` do not add up to `prim_num`, or the shells do not add
  up to `ao_num` Cartesian components;
* a shell has a negative angular momentum, more primitives than
  `ao_prim_num_max`, a non-positive exponent, or an out-of-range
  `basis_nucleus_index`.

`libcint_check_abi()` is called first and stops if `c_int` or `c_double`
is not the default `integer` or `double precision` kind; it cannot see an
`I8=on` library (§15).

`ao_integrals_n_e` additionally falls back when a nuclear charge is negative
or not an integer (§7.2): that test is the separate provider
`libcint_nuc_available`, which is `libcint_available` plus "every
`nucl_charge` is integral and non-negative", so the two-electron integrals
stay on libcint for fractional charges. Everything else falls back wholesale
by testing a single logical:

```
use_libcint_effective = use_libcint .and. libcint_available
```

which the dispatch sites (§10) consult. When `libcint_available` is true it
also `PROVIDE`s the `atm`/`bas`/`env` tables and `libcint_scale`, so the
first integral call does not trigger a provider from inside a parallel
region.

## 8. One-electron engine

Three operators, same pattern: a subroutine filling a passed matrix.

```fortran
subroutine libcint_ao_one_e(op, mat)
  integer, intent(in)           :: op        ! LIB_OP_OVLP/KIN/NUC
  double precision, intent(out) :: mat(ao_num,ao_num)
```

It is wrapped by three named subroutines, `libcint_ao_overlap(mat)`,
`libcint_ao_kinetic(mat)` and `libcint_ao_nuc(mat)`, which is what the
providers call: the `LIB_OP_*` constants then never leave
`libcint_ao_integrals`.

The three operators share the shell-pair loop but not the call: the
`select case (op)` sits before the loops, and each branch contains its own
copy of the loop, calling `cint1e_ovlp_cart`, `cint1e_kin_cart` or
`cint1e_nuc_cart` directly. Nothing is tested per shell pair, and the
duplication is confined to the twenty-line loop body.

```fortran
select case (op)
case (LIB_OP_OVLP)
  !$OMP PARALLEL ...
  do sj = 1, shell_num
    do si = 1, sj
      info = cint1e_ovlp_cart(blk, shls, atm, natm, bas, nbas, &
          env, c_null_ptr)
      ...
    enddo
  enddo
  !$OMP END PARALLEL
case (LIB_OP_KIN)
  ! same loop, cint1e_kin_cart
case (LIB_OP_NUC)
  ! same loop, cint1e_nuc_cart
case default
  stop 'libcint_ao_one_e: unknown operator'
end select
```

(An earlier draft bound the entry point to a procedure pointer, but IRPF90's
declaration parser does not recognise `procedure(...)`: the pointer is left
implicitly typed and the call does not compile, so the pointer would have to
live in a separate `.f90` module. Duplicating the loop is simpler.)

Loop over `sj = 1..shell_num`, `si = 1..sj` (the three operators are
symmetric), OpenMP over `sj`. The block buffer is allocated once per thread
(`ncart(lmax)**2` doubles) inside the parallel region and declared
`PRIVATE`, so no `THREADPRIVATE` is needed. For each shell pair compute the
block

```
info = cint1e_ovlp_cart(blk, shls, atm, natm, bas, nbas, env, c_null_ptr)
```

(the `FINT` return value is ignored), then

```
mat(i,j) = r_i*r_j*blk(ip + jp*ni)      ! i = ao_first_of_shell(si)+ip-1
mat(j,i) = mat(i,j)
```

The block buffer is private to each OpenMP thread, so no `THREADPRIVATE`
is needed for the 1e engine.

Why a subroutine and not a provider: `ao_overlap`, `ao_kinetic_integrals`
and `ao_integrals_n_e` are three `ao_num^2` arrays; the QP2 path needs them
anyway, and duplicating them in a "libcint provider" would double the
memory for nothing. The providers keep their `BEGIN_PROVIDER` declaration
and call the subroutine.

Operator mapping:

| QP2 provider | libcint | notes |
|---|---|---|
| `ao_overlap` | `cint1e_ovlp_cart` | |
| `ao_kinetic_integrals` | `cint1e_kin_cart` | libcint's `1/2 p.p`, same sign as QP2's `-1/2 (d2x+d2y+d2z)` |
| `ao_integrals_n_e` | `cint1e_nuc_cart` | sums over all `natm` internally; pseudo and point-charge terms are added afterwards by the existing code |

Verified numerically for an s function (`a=1`): `<s|T|s> = 0.2349964007`
and `<s|-1/r|s> = -0.25`, both equal to the analytic values times `f(0)^2`.

### 8.1 Splitting `ao_overlap`

`ao_overlap.irp.f` currently computes four arrays in one chained provider:
`ao_overlap`, `ao_overlap_x`, `ao_overlap_y`, `ao_overlap_z`, plus a
disk-read branch and a write branch.

New shape:

```
BEGIN_PROVIDER [double precision, ao_overlap_qp2, (ao_num,ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_x,  (ao_num,ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_y,  (ao_num,ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_z,  (ao_num,ao_num)]
  ! existing body, unchanged: use_cgtos copy or the QP2 loop,
  ! writing ao_overlap_qp2 instead of ao_overlap
END_PROVIDER

BEGIN_PROVIDER [double precision, ao_overlap, (ao_num,ao_num)]
  if (read_ao_integrals_overlap) then
    call ezfio_get_ao_one_e_ints_ao_integrals_overlap(ao_overlap)
  else if (use_libcint_effective) then
    call libcint_ao_overlap(ao_overlap)
  else
    PROVIDE ao_overlap_qp2
    ao_overlap = ao_overlap_qp2
  endif
  if (write_ao_integrals_overlap) call ezfio_set_...
END_PROVIDER
```

`ao_overlap_qp2` is lazy: it is only built when the QP2 path is taken, or
when the validation program asks for it (§13). The libcint path therefore
does not pay for it. `ao_overlap_x/_y/_z` have no consumer anywhere in the
tree (`git grep`), so letting them be unbuilt on the libcint path is
harmless; on the QP2 path they are still filled exactly as today.

`ao_overlap_abs` (the `|chi_i| |chi_j|` screening matrix), `ao_overlap_imag`
and `ao_overlap_complex` keep the QP2 code.

### 8.2 Kinetic and NAI providers

Both keep their disk-read and write branches and gain one `else if`:

```fortran
  else if (use_libcint_effective) then
    call libcint_ao_kinetic(ao_kinetic_integrals)
  else
    ... existing QP2 code ...
```

and `libcint_ao_nuc(ao_integrals_n_e)`, guarded by `libcint_nuc_available`
instead of `use_libcint_effective`, for the NAI.

`ao_deriv2_x/_y/_z` become unused on the libcint path (no consumer
elsewhere), so they are simply not computed; the QP2 path is untouched.

For `ao_integrals_n_e` the pseudo/point-charge additions at lines 103-108
stay where they are: libcint returns the bare nuclear attraction, exactly
like the QP2 branch does before those lines.

## 9. Two-electron engine

### 9.1 Element-wise API

The consumers of `ao_two_e_integral` are element-wise functions; there is no
block API today. We keep the signature:

```fortran
double precision function libcint_ao_two_e_integral(i,j,k,l,lr)
```

The long-range flag is an **argument**, not a read of `use_only_lr`: that
variable is declared in `ao_two_e_ints`, and referencing it from
`libcint_ao_integrals` would create a NEED cycle. The dispatcher passes
`use_only_lr`.

Per call:

1. map `i,j,k,l` to shells `si,sj,sk,sl` (`ao_shell`);
2. if the per-thread cache holds the block for `(si,sj,sk,sl)` and the same
   range-separation flag, return the cached element;
3. otherwise call

   ```
   call cint2e_cart(blk, shls, atm, natm, bas, nbas, env, c_null_ptr)
   ```

   (a `CINTOpt` is not built in PR1 — `CINT2e_drv` handles `opt == NULL`,
   `cache == NULL` and `dims == NULL`, verified in v6.2.1 `cint2e.c:797-819`
   and with the `test_1e8.c` probe);
4. cache the block (replace the previous one) and return
   `r_i*r_j*r_k*r_l*blk(ip + jp*ni + kp*ni*nj + lp*ni*nj*nk)`.

`do_schwartz_accel` is bypassed: with libcint the dispatcher returns before
that branch. `ao_two_e_integral_zero` is unchanged and still screens the
callers.

### 9.2 The per-thread block cache

One entry per thread, keyed by the four shell indices plus the LR flag,
living in `libcint_block_cache` — a plain Fortran module appended to
:file:`libcint_interface.f90`, not an IRPF90 source:

```fortran
integer, save :: cache_shls(4) = 0
integer, save :: cache_dim(4)  = 0
logical, save :: cache_lr     = .False.
real(c_double), save, allocatable, target :: cache_blk(:)
!$OMP THREADPRIVATE(cache_shls, cache_dim, cache_lr, cache_blk)
```

The reason for a `.f90` module is ordering, not taste: IRPF90 inserts the
`PROVIDE` calls it generates at the top of the routine body, before any
directive written in that body, while OpenMP requires `THREADPRIVATE` in
the specification part of the program unit. A module has no generated
statements, so the directive stays where the compiler expects it. The
module takes `atm`, `bas` and the two `env` arrays as arguments (it cannot
`use` the IRPF90 providers) and returns the block as a `pointer` plus the
four Cartesian dimensions; `libcint_ao_two_e_integral` provides the tables
in its own scope.

The block size is `ncart(l_si)*ncart(l_sj)*ncart(l_sk)*ncart(l_sl)`; the
buffer grows to the largest block ever requested by that thread, so the
worst-case memory is `(2*lmax+1)^4` doubles per thread (1296 doubles = 10 kB
for `lmax = 2`, 10000 doubles = 80 kB for `lmax = 4`).

Lifetime: the buffer is allocated by the first quartet a thread computes and
is never freed during the run, because the next quartet reuses it. OpenMP
frees the storage of a `THREADPRIVATE` copy like that of a static variable
of the base language, at an unspecified point (OpenMP 5.0, §2.19.2), so an
explicit `deallocate` would only release it early and force the next quartet
to allocate again. Nothing has to be reset between parallel regions either:
`cache_shls` starts zeroed on a fresh thread copy, which matches no shell
index, so the first call misses and recomputes.

`SAVE` + `!$OMP THREADPRIVATE` on an allocatable module variable was
validated with gfortran `-fopenmp` (the `tp.f90` probe in
`~/OpenCode/tmp/libcint-check` also checks that concurrent threads do not
share the buffer). Before this PR QP2 had no `THREADPRIVATE` anywhere
(there is no `COMMON`/module mutable state in the integral code), so this
is the first use; it is confined to this module.

A single entry is enough for the access patterns that matter:

* `compute_ao_two_e_integrals(j,k,l,...)` loops over all `i` with `j,k,l`
  fixed: every shell quartet is computed once and reused `ncart(l_i)` times.
* the Cholesky direct path (`cholesky.irp.f:145` and `:418`) walks a list of
  shell pairs in lexicographic order, so consecutive calls share three of
  the four shells and the block is reused `ncart` times before a new one is
  needed.

A hash of blocks was rejected: the number of shell quartets grows as
`shell_num^4`, so any cache that could actually hold the working set of the
Cholesky path would cost hundreds of megabytes per thread.

### 9.3 Range separation

`use_only_lr` selects the long-range operator `erf(mu_erf*r12)/r12` in the
QP2 dispatcher (line 54). With libcint, `ao_two_e_integral` routes directly
to `libcint_ao_two_e_integral(i,j,k,l,use_only_lr)`, which picks
`libcint_env_lr` (where `env(PTR_RANGE_OMEGA+1) = mu_erf`) instead of
`libcint_env`; `ao_two_e_integral_erf` routes to
`libcint_ao_two_e_integral(i,j,k,l,.True.)` before its own primitive-count
dispatch. The cache key includes the flag, so the two operators cannot be
confused.

`mu_erf = 0` is special: libcint reads `PTR_RANGE_OMEGA = 0` as "no range
separation", i.e. the plain Coulomb operator (`g2e.c:4443-4455`), whereas
`erf(0*r12)/r12` is zero. `libcint_ao_two_e_integral` therefore returns
zero before touching the cache when `lr` is set and `mu_erf = 0`, which
covers both `ao_two_e_integral_erf` and `ao_two_e_integral` with
`use_only_lr`.

Nothing else changes: `ao_two_e_integral_erf_schwartz` is built from
`ao_two_e_integral_erf`, so it follows the same engine, and
`ao_two_e_integral_schwartz_accel_erf` keeps the QP2 code (the latter is
only reached when libcint is off, because the libcint branch returns
first). The short-range `erfc` operator is not exposed through these
dispatchers; if it is ever wired, it needs a third `env` with a negative
`PTR_RANGE_OMEGA`.

## 10. Dispatch and fallbacks

`ao_two_e_integral` (dispatcher at `two_e_integrals.irp.f:20`):

```fortran
  if (use_cgtos) then
    ao_two_e_integral = ao_two_e_integral_cgtos(i,j,k,l)
  else if (use_libcint_effective) then
    ao_two_e_integral = libcint_ao_two_e_integral(i,j,k,l,use_only_lr)
  else if (use_only_lr) then
    ao_two_e_integral = ao_two_e_integral_erf(i,j,k,l)
  else if (do_schwartz_accel(i,j,k,l)) then
    ...
```

`use_cgtos` first: complex AOs can never go through libcint, and a file that
sets both flags must keep working (with a warning from `libcint_available`).

`ao_two_e_integral_erf` (`two_e_integrals_erf.irp.f:20`) has the same
branch at its top, before its usual primitive-count dispatch:

```fortran
  if (use_libcint_effective) then
    ao_two_e_integral_erf = libcint_ao_two_e_integral(i,j,k,l,.True.)
    return
  endif
```

Summary of the paths:

| case | overlap | kinetic | NAI | 2e |
|---|---|---|---|---|
| `read_ao_integrals_*` | disk | disk | disk | — |
| `use_cgtos` | QP2 (cGTO) | QP2 (cGTO) | QP2 (cGTO) | QP2 (cGTO) |
| `use_libcint` + available | libcint | libcint | libcint | libcint |
| non-integer or negative charges | libcint | libcint | **QP2** | libcint |
| `use_libcint` + unavailable | QP2 | QP2 | QP2 | QP2 |
| `use_libcint = False` | QP2 | QP2 | QP2 | QP2 |

The 2e column covers both `ao_two_e_integral` and `ao_two_e_integral_erf`;
on libcint, the `erf` one simply passes `lr = .True.`.

Staying on QP2 in all cases: `ao_overlap_abs`, `ao_overlap_x/_y/_z` (via
`ao_overlap_qp2`), `ao_sphe_overlap`, `ao_deriv2_x/_y/_z`,
`ao_integrals_n_e_per_atom`, `ao_pseudo_integrals`, `ao_integrals_pt_chrg`,
`ao_two_e_integral_general`, `ao_two_e_integral_schwartz_accel`,
`ao_two_e_integral_schwartz`, `ao_two_e_integral_schwartz_accel_erf`, the
map builder and Cholesky drivers.

## 11. Thread safety

* `atm`, `bas`, `env`, `env_lr` and the shell tables are built once, before
  any integral call, and are read-only afterwards. libcint writes nothing
  into them (it writes only into the `out` buffer).
* Each `cintXe_*` call gets a buffer owned by the calling thread: a
  per-thread allocation in `libcint_ao_one_e`, the `THREADPRIVATE` array in
  `libcint_block_cache` for the two-electron blocks.
* The 2e cache key and block are both `THREADPRIVATE`, so two threads never
  share a block, and a cached block cannot be returned to the wrong
  quartet.
* The libcint return status (`info`) is unused, but it is in the `PRIVATE`
  list of every parallel region in `libcint_ao_one_e`: under
  `DEFAULT(SHARED)` the assignment would be a write/write race.
* No `CINTOpt` is shared or built in PR1, so there is no hidden mutable
  state to protect.
* The `use_libcint` provider is called outside any parallel region
  (providers are evaluated lazily at first use, from serial code in practice;
  the master-read pattern of `use_cgtos` is copied verbatim, the `MPI_BCAST`
  that followed it removed with the multi-node layer, 2026-11).

## 12. Exact edit list

New files:

```
external/libcint                          (submodule, v6.2.1)
src/libcint_ao_integrals/NEED
src/libcint_ao_integrals/LIB
src/libcint_ao_integrals/README.rst
src/libcint_ao_integrals/libcint_interface.f90
src/libcint_ao_integrals/libcint_data.irp.f
src/libcint_ao_integrals/libcint_one_e.irp.f
src/libcint_ao_integrals/libcint_two_e.irp.f
src/ao_two_e_ints/test_libcint.irp.f   (validation program)
src/ao_two_e_ints/10.libcint.bats      (engine-vs-engine bats test)
src/ao_basis/libcint.irp.f              (use_libcint provider)
docs/libcint_design.md                  (this note; may move to docs/source)
```

Edited files:

| file | change |
|---|---|
| `.gitmodules` | `[submodule "external/libcint"]` |
| `configure` | `libcint` install case, `-i all` list, `find_lib -lcint` check |
| `etc/configure.rc` | `libcint` in the `-i` completion list |
| `.github/workflows/configuration.yml` | `./configure -i libcint` step |
| `.github/workflows/compilation.yml` | `./configure -i libcint` step (the check is mandatory) |
| `ci/macos-ci.sh` | `libcint` in the explicit dependency list |
| `INSTALL.rst`, `AGENTS.md` | mention `libcint` among the `-i` packages |
| `RELEASE_NOTES.org` | one line under the current version |
| `src/ao_basis/EZFIO.cfg` | `[use_libcint]` block after `[use_cgtos]` (line 58) |
| `src/ao_one_e_ints/NEED` | `+= libcint_ao_integrals` |
| `src/ao_two_e_ints/NEED` | `+= libcint_ao_integrals` |
| `src/ao_one_e_ints/ao_overlap.irp.f` | split provider (§8.1) |
| `src/ao_one_e_ints/kin_ao_ints.irp.f` | one `else if` branch in `ao_kinetic_integrals` |
| `src/ao_one_e_ints/pot_ao_ints.irp.f` | one `else if` branch in `ao_integrals_n_e` |
| `src/ao_two_e_ints/two_e_integrals.irp.f` | two lines in the dispatcher |
| `src/ao_two_e_ints/two_e_integrals_erf.irp.f` | two lines at the top of `ao_two_e_integral_erf` (§9.3) |

Nothing else. In particular no change to `cholesky.irp.f`,
`map_integrals.irp.f`, the screening functions, or the GPU module.

## 13. Validation

### 13.1 Unit level: `program test_libcint`

`src/ao_two_e_ints/test_libcint.irp.f` prints, with `ES24.16` formatting,
one tagged line per integral:

```
INFO ao_num           35
INFO shell_num        18
INFO nucl_num         5
INFO ao_l_max         2
INFO use_libcint      T
INFO libcint_avail    T
INFO libcint_nuc_avail T
INFO mu_erf           5.0000000000000000E-01
OVLP i j value
KIN  i j value
NAI  i j value
ERI  i j k l value
ERF  i j k l value
INFO max|ovlp-libcint - ovlp-qp2|   3.6765535060823140E-11
INFO sum(iijj)   5.0341697488988416E+02
```

`ERI` covers all `i <= j`, `k <= l` quartets (no 8-fold symmetry assumed);
`ERF` is the same list for the long-range operator at the EZFIO `mu_erf`.
The `INFO` lines are not compared between engines: they record the engine
selection, `libcint_available`, `libcint_nuc_available`, and two cheap
scalars (`max|ao_overlap - ao_overlap_qp2|`, `sum(i,i,j,j)`) that a wrong
`r_i` or a wrong `f(l)` shifts immediately.

`10.libcint.bats` runs the program twice on the same EZFIO, flipping
`use_libcint` between the runs, and compares the two `.dat` files keyed by
tag and indices: `|a-b| <= 1e-10` for `OVLP`, `1e-11` for the other tags,
with a `1e-12` relative escape. It creates its EZFIO from
`tests/input/{ch4,f}.xyz` and deletes it, plus its output files, on exit.

The `OVLP` tolerance is deliberately loose: QP2's `overlap_gaussian_xyz`
returns a hard-coded `1e-10` for each Cartesian factor when the Gaussian
product prefactor underflows (`fact_p < 1d-20`,
`src/utils/one_e_integration.irp.f:57-66`). That needs two *tight*
primitives on *different* centres; in CH4/cc-pVDZ the carbon 1s
($\alpha = 6665$) against the hydrogen 1s ($\alpha = 13.01$) at 2.06 bohr
has a Gaussian damping factor $\exp(-55) \approx 10^{-24}$ and a
$(\pi/p)^{3/2} \approx 10^{-5}$ prefactor, so `fact_p` lands around
$10^{-29}$ and the floor fires. libcint has no such cutoff. Measured max
$|$libcint $-$ QP2$|$:

| system | `ao_l_max` | `OVLP` | `KIN` | `NAI` | `ERI` | `ERF` |
|---|---|---|---|---|---|---|
| CH4 / STO-3G, $\mu=0.5$ | 1 | 8.9e-16 | 5.3e-15 | 2.4e-14 | 6.9e-15 | 3.4e-15 |
| CH4 / cc-pVDZ, $\mu=0.5$ | 2 | 3.68e-11 | 1.07e-14 | 9.24e-14 | 1.59e-14 | 2.72e-14 |
| CH4 / cc-pVDZ, $\mu=1.4$ | 2 | 3.68e-11 | 1.07e-14 | 9.24e-14 | 1.59e-14 | 1.67e-14 |
| F / cc-pVTZ, $\mu=0.5$ | 3 | 1.11e-15 | 1.07e-14 | 7.1e-15 | 1.40e-14 | 4.52e-14 |

The overlap floor is visible as the only entry that does not reach machine
precision: on CH4/cc-pVDZ, 28 of the 1225 overlaps differ by more than
`1e-13`, the largest being element (1,16), where both engines give
8.011119e-02 and the difference is 3.68e-11. F/cc-pVTZ has a single
nucleus, so no inter-atomic primitive pair can underflow and the overlap
agrees to 1.1e-15; STO-3G's tightest cross pair (C 1s at 71.6 against H 1s
at 3.4) stays above the cutoff, hence 8.9e-16.

The F/cc-pVTZ case is the one that matters for the mapping: it has `d` and
`f` shells (10 shells, 35 Cartesian functions), so a wrong Cartesian
ordering, a wrong `f(l)` or a wrong `r_i` cannot pass unnoticed.

Range separation is exercised by the `ERF` tag itself: the operator is
recomputed at `mu_erf = 1.4` and compared with the same run at 0.5. The
two sets of integrals differ by up to 8.611e-01, while each set matches the
QP2 value to 1.7e-14 and 2.7e-14 respectively: libcint really reads
`env(PTR_RANGE_OMEGA)`.

### 13.2 End to end

`src/hartree_fock/10.hf.bats` pins RHF/ROHF/UHF energies and the
point-charge test (`-92.79920682236470`), but it expects pre-existing EZFIO
directories (`h2_1.ezfio`, `hcn.ezfio`, ...) that no target in the
repository creates, so it cannot be run as-is. The end-to-end validation
was done with a small driver that creates a fresh EZFIO per molecule and
runs `scf` twice, flipping `use_libcint`:

| system | $E_\text{libcint} - E_\text{QP2}$ |
|---|---|
| H2 / STO-3G | 3.0e-15 |
| CH4 / cc-pVDZ, Cholesky (`do_ao_cholesky True`) | 3.3e-11 |
| CH4 / cc-pVDZ, map (`do_ao_cholesky False`) | 3.3e-11 |
| N2 / cc-pVTZ | 7.2e-10 |
| CO2 / cc-pVDZ | 1.0e-9 |
| H2S / cc-pVDZ | 1.4e-10 |
| F / cc-pVTZ (`-m 2`) | 6.5e-10 |
| HCN / def2-SVP + 2 point charges | 5.4e-10 |

The `1e-10 .. 1e-9` spread is the SCF convergence threshold
(`thresh_scf = 1e-10`) amplified by DIIS, not an integral error: CO2, H2S
and F reproduce the energies pinned in `10.hf.bats` to better than `1e-9`
with libcint, and the integral-level comparisons of §13.1 are three to
four orders of magnitude tighter. The point-charge run reproduces
`-92.79920682236470` to 5.4e-10 with libcint and to 8.4e-11 with QP2;
the NAI is on libcint there because `ao_integrals_pt_chrg` stays on QP2
and is simply added to it.

Range separation through `use_only_lr = True`, SCF: H2/STO-3G agrees to
1.0e-14 at $\mu=0.5$ and 4.0e-15 at $\mu=1.4$; CH4/cc-pVDZ to 1.1e-10 and
7.0e-11, the same DIIS explanation. A first version of that script forgot
to export `EZFIO_FILE`, so both legs ran libcint and the two energies were
identical to the last digit; `ezfio set` with an inline
`EZFIO_FILE=... assignment` writes nothing.

Threading: `test_libcint` at `OMP_NUM_THREADS=1` and `4` produces
bit-identical integral files, so the per-thread block cache does not leak
between threads and the summation order per integral is unchanged.

Fractional nuclear charges: a temporary program that sets
`nucl_charge(1) = 5.5` in memory (the EZFIO reader rounds a hand-edited
`5.0e-01` to zero -- a pre-existing quirk, it reproduces with
`use_libcint False`) makes libcint print its warning, report
`libcint_nuc_available F` and `use_libcint_effective T`, and return the
QP2 NAI matrix: `sum(ao_integrals_n_e) = -6.1513526647196954E+02`,
identical digits to QP2, while the 2e integrals stay on libcint.

## 14. Performance expectations

Correctness first, but the reason for the exercise is speed. PR1 is expected
to be:

* **faster** than QP2 for 1e integrals and for contracted bases with many
  primitives (libcint does the Rys quadrature once per primitive quartet and
  expands all Cartesian components at once);
* **comparable** for 2e integrals in the map/J-build paths, where the
  element-wise calls have strong shell-block locality and the single-entry
  cache turns `ncart` consecutive calls into one block;
* **no better, possibly slower** in the Cholesky direct path, where the
  element order is driven by the `Lset`/`Dset` pair lists and the block reuse
  is only `ncart`.

If the Cholesky benchmark regresses, the fix is a block API in a follow-up
(§16), not a change of PR1.

Benchmark protocol: `test_libcint` with `use_libcint` on/off for the 1e part,
and an `scf` run on CH4/cc-pVDZ with `do_ao_cholesky` both true and false
for the 2e part, on a quiet machine, interleaved A/B/B/A, reporting CPU time
(see the `benchmarking-noisy-machines` skill).

**Status: not benchmarked.** This PR validates correctness only; the
ordering above is a design argument, not a measurement. That matters
because it flips a default: if the Cholesky path turns out slower, the
choices are to run the protocol and ship `use_libcint = False` as the
default, or to accept the regression knowingly. A silent slowdown on the
path most production runs take is the one outcome to avoid.

## 15. Risks and open questions

1. **`r_i` scaling is the only place where a basis could be silently
   mis-integrated.** Mitigated by the §7.8 check and by the unit test; the
   failure mode would otherwise be a wrong energy, not a crash.
2. **`FINT` width.** If someone builds QP2 against an `I8=on` libcint, every
   call is silently wrong (no type error in `bind(C)` for assumed-size
   arrays). `libcint_check_abi()` only compares the C and Fortran kinds of
   `int` and `double`, so it cannot see the width of `FINT`, and the
   `find_lib` check cannot either. Mitigation: the module `LIB` and the
   `configure` case are the only supported path, and the `README.rst` says
   so; a real `sizeof(FINT)` check needs a tiny C shim (§16).
3. **`ao_expo`/`ao_coef` per AO vs per shell.** Checked at initialisation
   (§7.8) but the check costs one pass over the AOs. If it ever fails for a
   basis produced by `qp_create_ezfio`, we fall back to QP2 rather than
   producing wrong numbers.
4. **`CINTOpt` absent.** `cint2e_cart(..., NULL)` works (verified) but
   skips libcint's own screening and its `non0coeff` index. That is the
   right behaviour for QP2 (which does its own screening) but it is a
   performance item, not a correctness one.
5. **`natm = 0`.** An empty molecule cannot reach libcint: `libcint_available`
   is false when `nucl_num < 1`, so the dispatch never fires and
   `libcint_natm = nucl_num` is never passed to C. If that invariant ever
   changes, a one-line `max(nucl_num, 1)` guard with a zeroed first atom
   would make the tables well formed.
6. **Threading.** `THREADPRIVATE` is new in this code base. If a compiler
   rejects `SAVE` + `THREADPRIVATE` on an allocatable (older NVHPC?), the
   fallback is a fixed-size `cache_blk(0:1296)` dimensioned
   `(2*LIB_MAX_L+1)**4`, which is fine for `lmax <= 4`.
7. **`use_libcint` default and reproducibility.** Two QP2 versions
   disagreeing on the default is a reproducibility hazard. The EZFIO file
   carries the flag, and `qp edit` shows it, so a stored calculation
   reproduces itself; the release note must say the default changed.
8. **Energies move at the `1e-10` level, and libcint is the more accurate
   side.** The CH4/cc-pVDZ overlaps differ by up to `3.7e-11` because QP2
   substitutes `1e-10` for underflowed Cartesian overlaps (§13.1). Any
   stored energy or property computed with QP2 differs from the libcint one
   by roughly that amount times the integral's magnitude and the density.
   The unit test enforces `1e-10`, not machine precision, for that reason;
   tightening it would mean fixing `overlap_gaussian_xyz`, which is out of
   scope here.

## 16. Out of scope, and the obvious follow-ups

* Block-level 2e API (`cint2e_cart` for a whole shell quartet, consumed by
  a batched Cholesky/J-build) — the main performance follow-up.
* `CINTOpt` construction and libcint's own screening.
* `FRAC_CHARGE_NUC` for non-integer nuclear charges.
* Spherical (`_sph`) entry points, if QP2 ever wants spherical AOs natively.
* Periodic/cGTO integrals, ECP integrals, 3- and 4-centre one-electron
  operators, gradients.
* A runtime `I8` sanity check via a C shim.

---

## Appendix A — probe reproduction

`~/OpenCode/tmp/libcint-check/` contains a scratch build of libcint v6.2.1
and two probes:

* `test_shell.c`: one shell per atom, self-overlap for `l = 0,1,2`, which
  establishes `f(0)^2`, `f(1)^2`, `f(2)^2` and the Cartesian ordering
  `xx,xy,xz,yy,yz,zz`.
* `test_1e8.c`: the 8-argument `cint1e_*_cart` / `cint2e_cart` signatures
  with a trailing `CINTOpt *opt` (hand-declared, because the installed
  `cint.h` only exports the 10-argument drivers), and the s-only analytic
  values of `<s|s>`, `<s|T|s>`, `<s|-1/r|s>`, `(ss|ss)`.

The probes are not part of the PR; they are here to document where the
numbers in §7.6 and §7.7 come from.
