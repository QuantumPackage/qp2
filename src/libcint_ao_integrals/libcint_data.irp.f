use libcint_interface

!> Number of doubles needed by the `env` array.
BEGIN_PROVIDER [ integer, libcint_nenv ]
  implicit none
  BEGIN_DOC
  ! Size of the :c:data:`env` array passed to libcint: 20 header slots,
  ! 3 coordinates per nucleus, and the exponents and coefficients of every
  ! shell.
  END_DOC
  PROVIDE nucl_num prim_num
  libcint_nenv = PTR_ENV_START + 3*nucl_num + 2*prim_num + 8
END_PROVIDER

!> Number of atoms in the libcint `atm` table.
BEGIN_PROVIDER [ integer, libcint_natm ]
  implicit none
  BEGIN_DOC
  ! Number of entries in the :c:data:`libcint_atm` table, i.e. the number
  ! of nuclei.
  END_DOC
  PROVIDE nucl_num
  libcint_natm = nucl_num
END_PROVIDER

!> Number of shells in the libcint `bas` table.
BEGIN_PROVIDER [ integer, libcint_nbas ]
  implicit none
  BEGIN_DOC
  ! Number of entries in the :c:data:`libcint_bas` table, i.e. the number
  ! of |QP| shells.
  END_DOC
  PROVIDE shell_num
  libcint_nbas = shell_num
END_PROVIDER

!> libcint atom table, basis table and environment.
!>
!> All three arrays are built together, from the same data, and are
!> read-only afterwards: libcint never writes into `atm`, `bas` or `env`.
BEGIN_PROVIDER [ integer, libcint_atm, (ATM_SLOTS*nucl_num) ]
&BEGIN_PROVIDER [ integer, libcint_bas, (BAS_SLOTS*shell_num) ]
&BEGIN_PROVIDER [ double precision, libcint_env, (libcint_nenv) ]

  implicit none

  BEGIN_DOC
  ! libcint description of the |AO| basis.
  !
  ! - `atm` has one point-charge entry per nucleus;
  ! - `bas` has one shell per |QP| shell, with `NCTR_OF = 1`;
  ! - `env` holds the coordinates, the exponents and the contraction
  !   coefficients, the latter divided by `CINTcommon_fac_sp(l)` so that
  !   libcint returns the raw integrals in the |QP| convention. Only the
  !   coefficients of the first |AO| of each shell are stored; the
  !   per-|AO| normalization is carried by :f:func:`libcint_scale`.
  END_DOC

  integer          :: i, k, s, l, np, i0, ptr, ip
  double precision :: fac

  PROVIDE ao_expo_ordered ao_coef_normalized_ordered
  PROVIDE ao_coef_normalization_factor ao_first_of_shell shell_prim_num
  PROVIDE shell_ang_mom basis_nucleus_index nucl_coord

  libcint_atm = 0
  libcint_bas = 0
  libcint_env = 0.d0

  ! Nuclei
  do i = 1, nucl_num
    libcint_atm(ATM_SLOTS*(i-1) + CHARGE_OF  + 1) = nint(nucl_charge(i))
    libcint_atm(ATM_SLOTS*(i-1) + PTR_COORD  + 1) = PTR_ENV_START + 3*(i-1)
    libcint_atm(ATM_SLOTS*(i-1) + NUC_MOD_OF + 1) = POINT_NUC
    libcint_atm(ATM_SLOTS*(i-1) + PTR_ZETA   + 1) = 0
    do k = 1, 3
      libcint_env(PTR_ENV_START + 3*(i-1) + k) = nucl_coord(i,k)
    enddo
  enddo

  ! Exponents: one contiguous block per shell
  ptr = PTR_ENV_START + 3*nucl_num
  do s = 1, shell_num
    l  = shell_ang_mom(s)
    np = shell_prim_num(s)
    i0 = ao_first_of_shell(s)
    libcint_bas(BAS_SLOTS*(s-1) + ATOM_OF   + 1) = basis_nucleus_index(s) - 1
    libcint_bas(BAS_SLOTS*(s-1) + ANG_OF    + 1) = l
    libcint_bas(BAS_SLOTS*(s-1) + NPRIM_OF  + 1) = np
    libcint_bas(BAS_SLOTS*(s-1) + NCTR_OF   + 1) = 1
    libcint_bas(BAS_SLOTS*(s-1) + KAPPA_OF  + 1) = 0
    libcint_bas(BAS_SLOTS*(s-1) + PTR_EXP   + 1) = ptr
    do ip = 1, np
      libcint_env(ptr+ip) = ao_expo_ordered(i0,ip)
    enddo
    ptr = ptr + np
  enddo

  ! Coefficients: same layout, divided by the factor libcint will apply
  do s = 1, shell_num
    l  = shell_ang_mom(s)
    np = shell_prim_num(s)
    i0 = ao_first_of_shell(s)
    libcint_bas(BAS_SLOTS*(s-1) + PTR_COEFF + 1) = ptr
    fac = CINTcommon_fac_sp(l)
    do ip = 1, np
      libcint_env(ptr+ip) = ao_coef_normalized_ordered(i0,ip) / fac
    enddo
    ptr = ptr + np
  enddo

  ! Plain Coulomb operator: a zero exponent cutoff means the libcint default
  libcint_env(PTR_EXPCUTOFF + 1) = 0.d0

END_PROVIDER

!> `libcint_env` with the range-separation parameter set.
!>
!> A copy rather than a shared array with a per-call write, so that the
!> long-range and the plain Coulomb operators can be used concurrently from
!> different threads.
BEGIN_PROVIDER [ double precision, libcint_env_lr, (libcint_nenv) ]
  implicit none
  BEGIN_DOC
  ! Environment for the :math:`\mathrm{erf}(\mu r_{12})/r_{12}` operator.
  END_DOC
  PROVIDE libcint_env mu_erf
  libcint_env_lr = libcint_env
  libcint_env_lr(PTR_RANGE_OMEGA + 1) = mu_erf
END_PROVIDER

!> Per-|AO| correction applied after each libcint block.
!>
!> libcint has one coefficient set per shell, whereas |QP| normalizes each
!> |AO| separately (`d_xx` and `d_xy` do not have the same norm). The
!> coefficients given to libcint are those of the first |AO| of the shell,
!> and the remaining per-|AO| factor is applied to the integral.
BEGIN_PROVIDER [ double precision, libcint_scale, (ao_num) ]
  implicit none
  BEGIN_DOC
  ! Ratio of the |AO| normalization factor to the one of the first |AO| of
  ! the shell. Equal to 1 for `s` and `p` shells.
  END_DOC
  integer :: i, i0
  PROVIDE ao_coef_normalization_factor ao_shell ao_first_of_shell
  do i = 1, ao_num
    i0 = ao_first_of_shell(ao_shell(i))
    libcint_scale(i) = ao_coef_normalization_factor(i) &
                     / ao_coef_normalization_factor(i0)
  enddo
END_PROVIDER

!> True when libcint can represent the current basis and Hamiltonian.
!>
!> When true, the data providers are built, so that any caller testing
!> :c:data:`use_libcint_effective` has the tables ready before entering a
!> parallel region.
BEGIN_PROVIDER [ logical, libcint_available ]
  implicit none

  BEGIN_DOC
  ! True if the |AO| basis and the operators it is used with can be
  ! represented in libcint. If false, every integral falls back to the
  ! built-in |QP| code.
  END_DOC

  integer          :: i, j, ip, s, l, np, i0, ncart_total
  integer, external :: ao_power_index
  double precision :: ref, tol
  logical          :: ok

  PROVIDE use_cgtos is_periodic ao_l_max ao_num shell_num nucl_num
  PROVIDE prim_num shell_prim_num shell_ang_mom basis_nucleus_index
  PROVIDE ao_first_of_shell ao_prim_num ao_prim_num_max
  PROVIDE ao_expo_ordered ao_coef_normalized_ordered
  PROVIDE ao_coef_normalization_factor ao_power ao_nucl

  call libcint_check_abi()

  libcint_available = .False.
  ok = .True.

  if (ao_num < 1 .or. shell_num < 1 .or. nucl_num < 1) ok = .False.
  if (use_cgtos) ok = .False.
  if (is_periodic) ok = .False.
  if (ao_l_max > libcint_l_max_supported) ok = .False.
  if (sum(shell_prim_num) /= prim_num) ok = .False.

  ! The |AOs| must be exactly the Cartesian components of the shells
  ncart_total = 0
  if (ok) then
    do s = 1, shell_num
      if (shell_ang_mom(s) < 0) then
        ok = .False.
        exit
      endif
      ncart_total = ncart_total + libcint_ncart(shell_ang_mom(s))
    enddo
    if (ncart_total /= ao_num) ok = .False.
  endif

  ! Every component of a shell must share the same contraction: |QP|
  ! stores the coefficients per |AO|, libcint has one set per shell. This
  ! rejects a basis whose primitive normalization depends on the powers of
  ! the component (e.g. a TREXIO file written by a code that normalizes
  ! each Cartesian component separately).
  if (ok) then
    do s = 1, shell_num
      l  = shell_ang_mom(s)
      np = shell_prim_num(s)
      i0 = ao_first_of_shell(s)
      if (np < 1 .or. np > ao_prim_num_max) then
        ok = .False.
        exit
      endif
      if (i0 < 1 .or. i0 > ao_num) then
        ok = .False.
        exit
      endif
      if (basis_nucleus_index(s) < 1 &
          .or. basis_nucleus_index(s) > nucl_num) then
        ok = .False.
        exit
      endif
      do ip = 0, libcint_ncart(l)-1
        i = i0 + ip
        if (i > ao_num) then
          ok = .False.
          exit
        endif
        ! libcint finds the center through bas(ATOM_OF) and enumerates the
        ! components in its own order (x, y, z, xx, xy, ...). QP stores both
        ! per AO, so a file whose centers or powers disagree with the shell
        ! would give permuted blocks without any error.
        if (ao_nucl(i) /= basis_nucleus_index(s)) then
          ok = .False.
          exit
        endif
        if (ao_power(i,1) + ao_power(i,2) + ao_power(i,3) /= l) then
          ok = .False.
          exit
        endif
        if (ao_power_index(ao_power(i,1), ao_power(i,2), &
                ao_power(i,3)) /= ip + 1) then
          ok = .False.
          exit
        endif
        if (ao_prim_num(i) /= np) then
          ok = .False.
          exit
        endif
        do j = 1, np
          if (ao_expo_ordered(i0,j) <= 0.d0) then
            ok = .False.
            exit
          endif
          tol = 1.d-12 * max(1.d0, abs(ao_expo_ordered(i0,j)))
          if (abs(ao_expo_ordered(i,j) - ao_expo_ordered(i0,j)) > tol) then
            ok = .False.
            exit
          endif
          ref = ao_coef_normalized_ordered(i0,j) &
              / ao_coef_normalization_factor(i0)
          tol = 1.d-12 * max(1.d0, abs(ref))
          if (abs(ao_coef_normalized_ordered(i,j) &
                  / ao_coef_normalization_factor(i) - ref) > tol) then
            ok = .False.
            exit
          endif
        enddo
        if (.not. ok) exit
      enddo
      if (.not. ok) exit
    enddo
  endif

  if (.not.ok) then
    if (mpi_master) then
      write(6,'(A)') ' libcint: basis not representable, '// &
                     'using the built-in AO integrals'
    endif
    return
  endif

  libcint_available = .True.
  PROVIDE libcint_scale libcint_atm libcint_env_lr

END_PROVIDER

!> True when the nucleus-electron operator can be computed with libcint.
!>
!> libcint has a point-charge model and a fractional-charge model; the
!> latter is not wired yet, so a non-integer or negative nuclear charge
!> sends this operator back to the built-in code. The two-electron
!> integrals do not depend on the charges and stay on libcint.
BEGIN_PROVIDER [ logical, libcint_nuc_available ]
  implicit none
  integer :: i
  PROVIDE use_libcint nucl_num nucl_charge
  libcint_nuc_available = .False.
  if (use_libcint) then
    PROVIDE libcint_available
    libcint_nuc_available = libcint_available
  endif
  if (libcint_nuc_available) then
    do i = 1, nucl_num
      ! libcint replaces the charge by -abs(CHARGE_OF) (g1e.c:233), so a
      ! negative charge would attract where the built-in kernel repels.
      if (nucl_charge(i) < 0.d0) then
        libcint_nuc_available = .False.
        if (mpi_master) then
          write(6,'(A)') ' libcint: negative nuclear charge, '// &
                         'using the built-in nucleus-electron integrals'
        endif
        exit
      endif
      if (abs(nucl_charge(i) - dble(nint(nucl_charge(i)))) > 1.d-12) then
        libcint_nuc_available = .False.
        if (mpi_master) then
          write(6,'(A)') ' libcint: non-integer nuclear charge, '// &
                         'using the built-in nucleus-electron integrals'
        endif
        exit
      endif
    enddo
  endif
END_PROVIDER

!> Single switch consulted by every dispatch site.
BEGIN_PROVIDER [ logical, use_libcint_effective ]
  implicit none
  BEGIN_DOC
  ! :c:data:`use_libcint` restricted to the cases libcint can handle.
  END_DOC
  PROVIDE use_libcint
  if (use_libcint) then
    use_libcint_effective = libcint_available
  else
    use_libcint_effective = .False.
  endif
END_PROVIDER
