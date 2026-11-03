!-----------------------------------------------------------------------
! Exact AO and MO representation matrices.
!
! The only data that exists is the group elements' rotation matrices in the
! *input* frame (:c:data:`pg_rot_mat`), the atom permutation (:c:data:`pg_perm`) and the AO
! basis as EZFIO stores it. Everything else is rebuilt on demand:
!
!   M_g[q,p] = T_g[q,p] * B_p/B_q     (AO basis, block diagonal by shell)
!   X_g      = C^T S M_g C            (MO basis)
!
! with T_g = Sym^l(R_g) on the Cartesian monomials of the shell and B_p the
! AO factor, ``ao_coef_normalized(p,j) = A_j B_p``. B_p
! is read from :c:data:`ao_coef_normalized` itself (:c:func:`pg_ao_ratio`), never from
! :c:data:`ao_coef_normalization_factor`: the two agree in the modern EZFIO
! convention but not in the pre-2025 one, where the primitive factor also
! depends on the powers and :c:data:`ao_coef_normalization_factor` is only the last
! of the two multiplicative pieces. No matrix is stored across calls:
! :math:`|G|` can be 128 and ao_num 500, which would be 256 MB of dense matrices.
!
! The exact representation is M_g[q,p] = T_g[q,p] B_p/B_q, derived from
!   chi_p = sum_j A_j B_p v_p g_j ,  A_j = ao_coef(p,j)*prim_normalization(j)
! which holds exactly in the modern and in the legacy EZFIO conventions but
! *not* in the pre-2025 GAMESS one (there A_j depends on the AO powers and
! a d shell is not rotation covariant at all). :c:data:`pg_basis_ierr` records that
! case, and the j-independence test of :c:func:`pg_check_ao_factor` is what detects
! it: it uses ao_coef_normalized ratios only, never ao_coef/ao_coef_normalized.
!
! The authoritative test of all of this is G-iso, M_g^T S M_g = S, which
! needs no MOs.
!-----------------------------------------------------------------------

BEGIN_PROVIDER [ integer, pg_shell_target, (shell_num, pg_n_elem) ]
&BEGIN_PROVIDER [ integer, pg_basis_ierr ]
&BEGIN_PROVIDER [ character*(512), pg_basis_msg ]
  implicit none
  BEGIN_DOC
  ! ``pg_shell_target(s,g)`` is the unique shell on the atom
  ! ``pg_perm(basis_nucleus_index(s),g)`` whose contraction equals the one of
  ! shell ``s``. :c:data:`pg_basis_ierr` is 0 when the AO basis is invariant under the
  ! point group and the AO coefficients factorize as A_j*B_p inside every
  ! shell; otherwise it is 1 (no unique image shell), 2 (ambiguous image
  ! shell) or 3 (the primitive factor depends on the AO powers) and
  ! :c:data:`pg_basis_msg` says which atom, shell, angular momentum and element.
  !
  ! Reporting is left to the caller (:c:data:`pg_ao_ok` and the gate report): a
  ! provider that stops cannot be reported on.
  END_DOC
  integer :: s, sp, g, a, ap, n_match, s_match
  logical, external :: pg_shell_same

  pg_basis_ierr = 0
  pg_basis_msg = ''
  pg_shell_target(1:shell_num,1:pg_n_elem) = 0

  ! The primitive normalization must not depend on the AO powers: this is
  ! the factorization ao_coef_normalized(p,j) = A_j*B_p, without which no
  ! shell is rotation covariant (pre-2025 GAMESS convention, d shells).
  call pg_check_ao_factor(pg_basis_ierr, pg_basis_msg)
  if (pg_basis_ierr /= 0) return

  ! Every shell must have one and only one image shell on the image atom.
  do g = 1, pg_n_elem
    do s = 1, shell_num
      a = basis_nucleus_index(s)
      ap = pg_perm(a,g)
      n_match = 0
      s_match = 0
      do sp = 1, shell_num
        if (basis_nucleus_index(sp) /= ap) cycle
        if (.not. pg_shell_same(s,sp)) cycle
        n_match = n_match + 1
        s_match = sp
      enddo
      if (n_match == 0) then
        write(pg_basis_msg,'(A,I0,A,I0,A,I0,A,I0,A)') &
             'shell ', s, ' (atom ', a, ', l = ', &
             shell_ang_mom(s), ') has no counterpart on atom ', ap, &
             ' with the same contraction'
        pg_basis_ierr = 1
        return
      endif
      if (n_match > 1) then
        write(pg_basis_msg,'(A,I0,A,I0,A,I0,A)') &
             'shell ', s, ' (atom ', a, ', l = ', shell_ang_mom(s), &
             ') has several counterparts on atom with the same contraction'
        pg_basis_ierr = 2
        return
      endif
      pg_shell_target(s,g) = s_match
    enddo
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! Does the shell's primitive normalization depend on the AO powers?
!-----------------------------------------------------------------------
subroutine pg_check_ao_factor(ierr, msg)
  implicit none
  BEGIN_DOC
  ! Verifies, for every shell, that the ratios
  ! ``ao_coef_normalized(p,j)/ao_coef_normalized(p0,j)`` are independent of the
  ! primitive ``j`` (``p0`` = first AO of the shell). Cross-multiplied so that a
  ! zero coefficient needs no special case. This is the j-independence test
  ! of the representation: no double factorial, no hardcoded convention.
  END_DOC
  integer, intent(out) :: ierr
  character*(512), intent(out) :: msg
  integer, parameter :: nmax = 55
  integer :: s, p, p0, j, jref, np, nprim
  double precision :: c_ref(nmax), ratio, c0_ref, dev

  ierr = 0
  msg = ''
  do s = 1, shell_num
    p0 = ao_first_of_shell(s)
    nprim = ao_prim_num(p0)
    np = (shell_ang_mom(s)+1)*(shell_ang_mom(s)+2)/2

    ! Reference primitive: the one carrying the largest coefficient, so that
    ! the ratios below divide by the largest number available.
    jref = 1
    do j = 2, nprim
      if (abs(ao_coef_normalized(p0,j)) > abs(ao_coef_normalized(p0,jref))) jref = j
    enddo
    c0_ref = ao_coef_normalized(p0,jref)
    if (abs(c0_ref) == 0.d0) cycle
    do p = 1, np
      c_ref(p) = ao_coef_normalized(p0+p-1,jref)
    enddo

    do j = 1, nprim
      if (abs(ao_coef_normalized(p0,j)) == 0.d0) cycle
      do p = 1, np
        ! ratio(j) = C[p,j]/C[p0,j] must equal ratio(jref)
        dev = abs(ao_coef_normalized(p0+p-1,j)*c0_ref &
                - c_ref(p)*ao_coef_normalized(p0,j))
        ratio = abs(ao_coef_normalized(p0+p-1,j)*c0_ref) &
              + abs(c_ref(p)*ao_coef_normalized(p0,j))
        if (dev > 1.d-10*ratio) then
          write(msg,'(A,I0,A,I0,A,I0,A)') &
               'the AO coefficients do not factorize as A_j*B_p in shell ', &
               s, ' (atom ', basis_nucleus_index(s), ', l = ', &
               shell_ang_mom(s), '): the primitive normalization depends on ' &
               //'the AO powers (pre-2025 EZFIO convention); regenerate the ' &
               //'EZFIO with the current QP2'
          ierr = 3
          return
        endif
      enddo
    enddo
  enddo
end subroutine pg_check_ao_factor

!-----------------------------------------------------------------------
! Same contraction?
!-----------------------------------------------------------------------
logical function pg_shell_same(s1, s2)
  implicit none
  BEGIN_DOC
  ! True when the two shells have the same angular momentum and the same
  ! contraction: same number of primitives, same *sorted* exponents and the
  ! primitive coefficients that go with them, to 1e-10.
  ! ``ao_coef`` (raw), never :c:data:`ao_coef_normalized`: the latter legitimately
  ! differs between AOs of the same shell when their powers differ.
  END_DOC
  integer, intent(in) :: s1, s2
  double precision :: e1(ao_prim_num_max), c1(ao_prim_num_max)
  double precision :: e2(ao_prim_num_max), c2(ao_prim_num_max)
  integer :: n1, n2, j

  pg_shell_same = .false.
  if (shell_ang_mom(s1) /= shell_ang_mom(s2)) return
  n1 = ao_prim_num(ao_first_of_shell(s1))
  n2 = ao_prim_num(ao_first_of_shell(s2))
  if (n1 /= n2) return
  call pg_sorted_contraction(s1, e1, c1, n1)
  call pg_sorted_contraction(s2, e2, c2, n2)
  do j = 1, n1
    if (abs(e1(j)-e2(j)) > 1.d-10*max(1.d0, abs(e1(j)), abs(e2(j)))) return
    if (abs(c1(j)-c2(j)) > 1.d-10*max(1.d0, abs(c1(j)), abs(c2(j)))) return
  enddo
  pg_shell_same = .true.
end function pg_shell_same

!-----------------------------------------------------------------------
! Sorted (exponent, coefficient) list of a shell.
!-----------------------------------------------------------------------
subroutine pg_sorted_contraction(s, ex, cf, n)
  implicit none
  BEGIN_DOC
  ! Exponents of shell ``s`` in increasing order, each with its own raw
  ! coefficient. The pair must be sorted together: a basis that lists its
  ! primitives in another order is still the same contraction.
  END_DOC
  integer, intent(in) :: s
  integer, intent(out) :: n
  double precision, intent(out) :: ex(ao_prim_num_max)
  double precision, intent(out) :: cf(ao_prim_num_max)
  double precision :: te, tc
  integer :: p0, j, k, jmin

  p0 = ao_first_of_shell(s)
  n = ao_prim_num(p0)
  do j = 1, n
    ex(j) = ao_expo(p0,j)
    cf(j) = ao_coef(p0,j)
  enddo
  do j = 1, n-1
    jmin = j
    do k = j+1, n
      if (ex(k) < ex(jmin)) jmin = k
    enddo
    if (jmin /= j) then
      te = ex(j)
      tc = cf(j)
      ex(j) = ex(jmin)
      cf(j) = cf(jmin)
      ex(jmin) = te
      cf(jmin) = tc
    endif
  enddo
end subroutine pg_sorted_contraction

!-----------------------------------------------------------------------
! AO representation matrix of one group element.
!-----------------------------------------------------------------------
subroutine pg_ao_rep(ig, m)
  implicit none
  BEGIN_DOC
  ! ``m`` = M_g in the AO basis, a dense (ao_num x ao_num) matrix, block
  ! diagonal by shell: shell ``s`` maps onto ``pg_shell_target(s,g)`` and its
  ! block is Sym^l(R_g) rescaled by B_p/B_q.
  END_DOC
  integer, intent(in) :: ig
  double precision, intent(out) :: m(ao_num,ao_num)
  call pg_ao_rep_ex(ig, m)
end subroutine pg_ao_rep

!-----------------------------------------------------------------------
! AO factor ratio B_p/B_q.
!-----------------------------------------------------------------------
double precision function pg_ao_ratio(p, q)
  implicit none
  BEGIN_DOC
  ! ``B_p/B_q`` from ``ao_coef_normalized(p,j) = A_j B_p``: the ``j`` dependence
  ! cancels, so any primitive gives the same number. ``jref`` is the primitive
  ! carrying the largest ``ao_coef_normalized(q,j)`` so that the division is by
  ! the largest coefficient available (:c:func:`pg_check_ao_factor` has already
  ! proved that the ratio does not depend on ``j``).
  !
  ! This is exact in both EZFIO conventions. In the pre-2025 one
  ! :c:data:`ao_coef_normalization_factor` is NOT B_p: the primitive factor there
  ! carries a power-dependent part, so using it would make a d shell not
  ! rotation covariant by an O(1) factor (found by G-iso on ammonia/6-31G*).
  END_DOC
  integer, intent(in) :: p, q
  integer :: j, jref, nprim
  double precision :: den

  nprim = ao_prim_num(q)
  jref = 1
  do j = 2, nprim
    if (abs(ao_coef_normalized(q,j)) > abs(ao_coef_normalized(q,jref))) jref = j
  enddo
  den = ao_coef_normalized(q,jref)
  if (den == 0.d0) then
    ! Identically zero AO: cannot happen in a basis EZFIO stores. Return 0 so
    ! that the block vanishes and G-iso fails loudly in :c:func:`pg_report_ao`.
    pg_ao_ratio = 0.d0
  else
    pg_ao_ratio = ao_coef_normalized(p,jref)/den
  endif
end function pg_ao_ratio

!-----------------------------------------------------------------------
! Builder of the AO representation matrix (see :c:func:`pg_ao_rep` above).
!-----------------------------------------------------------------------
subroutine pg_ao_rep_ex(ig, m)
  implicit none
  BEGIN_DOC
  ! Block-diagonal builder: the block of shell ``s`` is Sym^l(R_g) rescaled by
  ! B_p/B_q and written at the rows of shell ``pg_shell_target(s,g)``.
  END_DOC
  integer, intent(in) :: ig
  double precision, intent(out) :: m(ao_num,ao_num)
  integer, parameter :: lmax = 9
  integer :: s, sp, p0, q0, np, nq, l, p, q
  double precision, external :: pg_ao_ratio
  integer :: powp(3, (lmax+1)*(lmax+2)/2), powq(3, (lmax+1)*(lmax+2)/2)
  double precision :: t((lmax+1)*(lmax+2)/2, (lmax+1)*(lmax+2)/2)

  ! ``m`` is intent(out) and zeroed first, so an early return still leaves a
  ! defined matrix behind.
  m(1:ao_num,1:ao_num) = 0.d0

  if (pg_basis_ierr /= 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A,A)') 'symmetry: the AO basis is not symmetric: ', &
           trim(pg_basis_msg)
      pg_detection_failed = .true.
    endif
    return
  endif

  do s = 1, shell_num
    sp = pg_shell_target(s,ig)
    if (sp == 0) cycle
    l = shell_ang_mom(s)
    if (l > lmax) then
      if (.not. pg_detection_failed) then
        write(6,'(A,I0,A,I0)') 'symmetry: shell ', s, &
             ' has angular momentum above the supported maximum ', lmax
        pg_detection_failed = .true.
      endif
      return
    endif
    np = (l+1)*(l+2)/2
    nq = np
    p0 = ao_first_of_shell(s)
    q0 = ao_first_of_shell(sp)
    do p = 1, np
      powp(1:3,p) = ao_power(p0+p-1,1:3)
      powq(1:3,p) = ao_power(q0+p-1,1:3)
    enddo
    call pg_monomial_rep(l, pg_rot_mat(1,1,ig), np, powp, nq, powq, t)
    do p = 1, np
      do q = 1, nq
        if (t(q,p) == 0.d0) cycle
        m(q0+q-1,p0+p-1) = t(q,p) &
             *pg_ao_ratio(p0+p-1, q0+q-1)
      enddo
    enddo
  enddo
end subroutine pg_ao_rep_ex

!-----------------------------------------------------------------------
! Sym^l of a 3x3 rotation, on the Cartesian monomials of a shell.
!-----------------------------------------------------------------------
subroutine pg_monomial_rep(l, r, np, powp, nq, powq, t)
  implicit none
  BEGIN_DOC
  ! Builds ``t`` such that ``U_g v_p = sum_q t(q,p) v_q`` with
  ! ``v_p(r) = r_x^a r_y^b r_z^c`` and ``U_g f(r) = f(R^-1 r)``.
  ! ``R(j,i)`` appears with the *row* index running over x,y,z: expanding
  ! ``(sum_j R(j,i) r_j)^p_i``, the coefficient of ``r_j`` is ``R(j,i)``.
  !
  ! The product is accumulated as p_i successive multiplications by the
  ! linear form ``sum_j R(j,i) r_j``, which needs no factorials and no
  ! hardcoded degree: it is the same expansion as the multinomial formula,
  ! only cheaper to get right.
  !
  ! ``powp``, ``powq`` and ``t`` are declared with their *full* extent (nmax = 55
  ! monomials for lmax = 9), not with (3,np) and (nq,np): Fortran associates
  ! an explicit-shape dummy with the actual argument by sequence, so a dummy
  ! (3,np) would stride by np while the caller's array strides by nmax, and
  ! every column above the first would land at the wrong offset. The extents
  ! here must stay equal to those of :c:func:`pg_ao_rep`.
  END_DOC
  integer, parameter :: nmax = 55
  integer, intent(in) :: l, np, nq
  double precision, intent(in) :: r(3,3)
  integer, intent(in) :: powp(3,nmax), powq(3,nmax)
  double precision, intent(out) :: t(nmax,nmax)
  double precision :: c(0:l,0:l,0:l), d(0:l,0:l,0:l)
  integer :: p, q, i, k, a, b, cc, deg

  t(1:nq,1:np) = 0.d0
  do p = 1, np
    c(0:l,0:l,0:l) = 0.d0
    c(0,0,0) = 1.d0
    deg = 0
    do i = 1, 3
      do k = 1, powp(i,p)
        d(0:l,0:l,0:l) = 0.d0
        do a = 0, deg
          do b = 0, deg - a
            cc = deg - a - b
            if (c(a,b,cc) == 0.d0) cycle
            d(a+1,b,cc) = d(a+1,b,cc) + c(a,b,cc)*r(1,i)
            d(a,b+1,cc) = d(a,b+1,cc) + c(a,b,cc)*r(2,i)
            d(a,b,cc+1) = d(a,b,cc+1) + c(a,b,cc)*r(3,i)
          enddo
        enddo
        c(0:l,0:l,0:l) = d(0:l,0:l,0:l)
        deg = deg + 1
      enddo
    enddo
    do q = 1, nq
      t(q,p) = c(powq(1,q),powq(2,q),powq(3,q))
    enddo
  enddo
end subroutine pg_monomial_rep

!-----------------------------------------------------------------------
! MO representation matrix of one group element.
!-----------------------------------------------------------------------
subroutine pg_get_mo_rep(ig, x)
  implicit none
  BEGIN_DOC
  ! ``x`` = C^T S M_g C, the matrix of group element ``ig`` in the MO basis.
  ! Nothing is kept between calls: the caller loops on
  ! the elements and holds one ``x`` at a time.
  !
  ! ``matmul``, not ``dgemm``: the matrices are at most a few hundred wide in
  ! practice and correctness matters more than the constant. Three dgemms
  ! replace it verbatim if a production run ever makes this the hot spot.
  END_DOC
  integer, intent(in) :: ig
  double precision, intent(out) :: x(mo_num,mo_num)
  double precision, allocatable :: m(:,:), cs(:,:), w(:,:)

  if (mo_num == 0) then
    return
  endif
  allocate(m(ao_num,ao_num), cs(mo_num,ao_num), w(mo_num,ao_num))
  call pg_ao_rep(ig, m)
  cs(1:mo_num,1:ao_num) = matmul(transpose(mo_coef), ao_overlap)
  w(1:mo_num,1:ao_num) = matmul(cs, m)
  x(1:mo_num,1:mo_num) = matmul(w, mo_coef)
  deallocate(m, cs, w)
end subroutine pg_get_mo_rep
