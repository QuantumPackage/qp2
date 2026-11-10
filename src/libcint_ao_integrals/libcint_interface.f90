!> C bindings to the Cartesian integral kernels of libcint.
!>
!> Only the entry points used by |QP| are declared here, and only in their
!> `_cart` flavour: the |AO| basis is Cartesian at this level, the spherical
!> transformation being applied later by
!> :file:`src/ao_basis/spherical_to_cartesian.irp.f`.
!>
!> The library must be built with `I8=off`, so that `FINT` is a 32-bit
!> integer and matches the default `integer` kind of Fortran, and with
!> `WITH_CINT2_INTERFACE=ON`, so that the `cintXe_*_cart` wrappers are
!> exported. With `I8=on` the arguments silently change width and every
!> integral returned through this interface is garbage.
module libcint_interface

  use iso_c_binding

  implicit none

  private

  !> Slot offsets in the flat `atm` array (see `cint.h`).
  integer(c_int), parameter, public :: &
      CHARGE_OF  = 0, &
      PTR_COORD  = 1, &
      NUC_MOD_OF = 2, &
      PTR_ZETA   = 3

  !> Slot offsets in the flat `bas` array (see `cint.h`).
  integer(c_int), parameter, public :: &
      ATOM_OF   = 0, &
      ANG_OF    = 1, &
      NPRIM_OF  = 2, &
      NCTR_OF   = 3, &
      KAPPA_OF  = 4, &
      PTR_EXP   = 5, &
      PTR_COEFF = 6

  !> Number of integers per atom in `atm` and per shell in `bas`.
  integer(c_int), parameter, public :: ATM_SLOTS = 6, BAS_SLOTS = 8

  !> Nuclear model code for a point charge.
  integer(c_int), parameter, public :: POINT_NUC = 1

  !> Offsets in the flat `env` array (see `cint.h`).
  integer(c_int), parameter, public :: &
      PTR_EXPCUTOFF   = 0, &
      PTR_RANGE_OMEGA = 8, &
      PTR_ENV_START   = 20

  !> Largest angular momentum libcint can integrate with Rys quadrature.
  integer, parameter, public :: libcint_l_max_supported = 12

  !> Identity of the one-electron operators handled by
  !> :f:func:`libcint_ao_one_e`. Callers outside this directory should use
  !> the named wrappers `libcint_ao_overlap`, `libcint_ao_kinetic` and
  !> `libcint_ao_nuc`.
  integer, parameter, public :: LIB_OP_OVLP = 1, LIB_OP_KIN = 2, LIB_OP_NUC = 3

  public :: libcint_ncart, libcint_check_abi

  !> libcint multiplies every integral by `CINTcommon_fac_sp(l)` for each
  !> shell. Declared here to document the ABI; the value is folded into the
  !> contraction coefficients at initialisation, see `libcint_data.irp.f`.
  public :: CINTcommon_fac_sp

  interface

    !> Overlap integrals, all Cartesian components of the two shells.
    function cint1e_ovlp_cart(out, shls, atm, natm, bas, nbas, env, opt) &
        bind(C, name="cint1e_ovlp_cart") result(res)
      use iso_c_binding
      real(c_double), intent(out) :: out(*)
      integer(c_int), intent(in)  :: shls(*)
      integer(c_int), intent(in)  :: atm(*)
      integer(c_int), value       :: natm
      integer(c_int), intent(in)  :: bas(*)
      integer(c_int), value       :: nbas
      real(c_double), intent(in)  :: env(*)
      type(c_ptr), value          :: opt
      integer(c_int)              :: res
    end function cint1e_ovlp_cart

    !> Kinetic energy integrals.
    function cint1e_kin_cart(out, shls, atm, natm, bas, nbas, env, opt) &
        bind(C, name="cint1e_kin_cart") result(res)
      use iso_c_binding
      real(c_double), intent(out) :: out(*)
      integer(c_int), intent(in)  :: shls(*)
      integer(c_int), intent(in)  :: atm(*)
      integer(c_int), value       :: natm
      integer(c_int), intent(in)  :: bas(*)
      integer(c_int), value       :: nbas
      real(c_double), intent(in)  :: env(*)
      type(c_ptr), value          :: opt
      integer(c_int)              :: res
    end function cint1e_kin_cart

    !> Nucleus-electron attraction, summed over every nucleus in `atm`.
    function cint1e_nuc_cart(out, shls, atm, natm, bas, nbas, env, opt) &
        bind(C, name="cint1e_nuc_cart") result(res)
      use iso_c_binding
      real(c_double), intent(out) :: out(*)
      integer(c_int), intent(in)  :: shls(*)
      integer(c_int), intent(in)  :: atm(*)
      integer(c_int), value       :: natm
      integer(c_int), intent(in)  :: bas(*)
      integer(c_int), value       :: nbas
      real(c_double), intent(in)  :: env(*)
      type(c_ptr), value          :: opt
      integer(c_int)              :: res
    end function cint1e_nuc_cart

    !> Two-electron repulsion integrals of a whole shell quartet.
    function cint2e_cart(out, shls, atm, natm, bas, nbas, env, opt) &
        bind(C, name="cint2e_cart") result(res)
      use iso_c_binding
      real(c_double), intent(out) :: out(*)
      integer(c_int), intent(in)  :: shls(*)
      integer(c_int), intent(in)  :: atm(*)
      integer(c_int), value       :: natm
      integer(c_int), intent(in)  :: bas(*)
      integer(c_int), value       :: nbas
      real(c_double), intent(in)  :: env(*)
      type(c_ptr), value          :: opt
      integer(c_int)              :: res
    end function cint2e_cart

    !> Factor applied to every integral of a shell of angular momentum `l`.
    function CINTcommon_fac_sp(l) bind(C, name="CINTcommon_fac_sp") result(res)
      use iso_c_binding
      integer(c_int), value :: l
      real(c_double)        :: res
    end function CINTcommon_fac_sp

  end interface

  public :: cint1e_ovlp_cart, cint1e_kin_cart, cint1e_nuc_cart, cint2e_cart

contains

  !> Number of Cartesian functions in a shell of angular momentum `l`.
  !>
  !> |QP| stores Cartesian |AOs|, so this is also the number of |AOs| of the
  !> shell, and consecutive shells start `libcint_ncart(l)` |AOs| apart.
  integer function libcint_ncart(l)
    implicit none
    integer, intent(in) :: l
    libcint_ncart = ((l+1)*(l+2))/2
  end function libcint_ncart

  !> Check that the Fortran kinds match the C ABI assumed by libcint.
  !>
  !> Called once at initialisation. It cannot detect a library built with
  !> `I8=on` (that changes `FINT`, not `double`), which is why
  !> :file:`INSTALL.rst` states that `libcint.a` must come from
  !> `./configure -i libcint`.
  subroutine libcint_check_abi()
    implicit none
    if (c_int /= kind(0)) then
      stop 'libcint: C int and default integer kinds differ'
    endif
    if (c_double /= kind(0.d0)) then
      stop 'libcint: C double and double precision kinds differ'
    endif
  end subroutine libcint_check_abi

end module libcint_interface

!> Thread-local cache for the two-electron shell-quartet blocks.
!>
!> `cint2e_cart` returns the integrals of a whole shell quartet, while |QP|
!> consumes them one at a time. The callers walk the quartets in an order
!> that reuses the same block many times in a row, so a one-entry cache per
!> thread removes almost all of the recomputation.
!>
!> The cache is a plain Fortran module, not an IRPF90 source, because OpenMP
!> requires `THREADPRIVATE` in the specification part of the routine, while
!> IRPF90 inserts the calls that provide the tables above any directive
!> found in the body.
!>
!> Lifetime: the block buffer is allocated by the first quartet a thread
!> computes and only grows, up to `(2*lmax+1)^4` doubles per thread (80 kB
!> for `lmax = 4`). It is deliberately not freed during the run, because the
!> next quartet reuses it: OpenMP frees the storage of a `THREADPRIVATE`
!> copy like that of a static variable of the base language, at an
!> unspecified point (OpenMP 5.0, 2.19.2), so an explicit `deallocate` would
!> only release it early and make the next quartet allocate again. If a
!> runtime ever reinitialises a thread's copy between parallel regions, the
!> zero-initialised `cache_shls` misses and the block is recomputed.
module libcint_block_cache

  use iso_c_binding
  use libcint_interface, only: cint2e_cart, libcint_ncart, ANG_OF, &
      BAS_SLOTS

  implicit none

  private

  public :: libcint_two_e_block

  !> Shell quartet currently cached by this thread, its dimensions, and
  !> whether it was computed with the long-range operator.
  integer, save :: cache_shls(4) = 0
  integer, save :: cache_dim(4)  = 0
  logical, save :: cache_lr     = .False.
  real(c_double), save, allocatable, target :: cache_blk(:)
  !$OMP THREADPRIVATE(cache_shls, cache_dim, cache_lr, cache_blk)

contains

  !> Block of the shell quartet `(si, sj, sk, sl)`, computed if this thread
  !> has not already computed it.
  !>
  !> `dims` receives the number of Cartesian components of the four shells,
  !> in the order `si`, `sj`, `sk`, `sl`. `atm`, `bas`, `env` are the
  !> libcint tables of :file:`libcint_data.irp.f`, and `env_lr` its copy
  !> carrying the range-separation parameter; they are arguments of this
  !> routine so that the module does not depend on the IRPF90 providers.
  subroutine libcint_two_e_block(si, sj, sk, sl, lr, atm, natm, bas, nbas, &
      env, env_lr, blk, dims)

    implicit none

    integer, intent(in)  :: si, sj, sk, sl
    logical, intent(in)  :: lr
    integer, intent(in)  :: natm, nbas
    integer, intent(in)  :: atm(*), bas(*)
    real(c_double), intent(in) :: env(*), env_lr(*)
    real(c_double), pointer, intent(out) :: blk(:)
    integer, intent(out) :: dims(4)

    integer(c_int) :: shls(4)
    integer(c_int) :: info
    integer        :: ni, nj, nk, nl, nsz

    if (cache_shls(1) == si .and. cache_shls(2) == sj .and. &
        cache_shls(3) == sk .and. cache_shls(4) == sl .and. &
        (cache_lr .eqv. lr)) then
      blk => cache_blk
      dims = cache_dim
      return
    endif

    ! `bas` is the Fortran (1-based) view of the libcint table: the +1 turns
    ! the C slot offset into the Fortran element index
    ni = libcint_ncart(bas(ANG_OF + BAS_SLOTS*(si-1) + 1))
    nj = libcint_ncart(bas(ANG_OF + BAS_SLOTS*(sj-1) + 1))
    nk = libcint_ncart(bas(ANG_OF + BAS_SLOTS*(sk-1) + 1))
    nl = libcint_ncart(bas(ANG_OF + BAS_SLOTS*(sl-1) + 1))
    nsz = ni*nj*nk*nl

    if (allocated(cache_blk)) then
      ! keep the capacity: mixed s/p/d bases alternate between block sizes
      ! on every quartet, and a larger buffer already holds a smaller block
      if (size(cache_blk) < nsz) then
        deallocate(cache_blk)
        allocate(cache_blk(nsz))
      endif
    else
      allocate(cache_blk(nsz))
    endif

    shls = (/ si-1, sj-1, sk-1, sl-1 /)
    ! `cint2e_cart` returns the "block not empty" flag, which is unused
    if (lr) then
      info = cint2e_cart(cache_blk, shls, atm, natm, bas, nbas, env_lr, &
          c_null_ptr)
    else
      info = cint2e_cart(cache_blk, shls, atm, natm, bas, nbas, env, &
          c_null_ptr)
    endif

    cache_shls = (/ si, sj, sk, sl /)
    cache_dim  = (/ ni, nj, nk, nl /)
    cache_lr   = lr

    blk => cache_blk
    dims = cache_dim

  end subroutine libcint_two_e_block

end module libcint_block_cache
