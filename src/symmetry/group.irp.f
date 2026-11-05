! Group elements, classes and identification.
!
! The whole finite set is named: 204 groups, checked one by one against the
! closed-form construction of the reference oracle (refsym2.py + ident.py,
! 204/204). There is no Abelian-only shortcut: the enumeration, the classes
! and the identification are general, :math:`|G|` <= 128.
!
! Placement notes:
!  * pg_n_max and pg_n_axes_max belong to the group description but are
!    needed to build pg_rot, so they are computed here as providers that
!    do not depend on pg_rot. pg_rot uses a local scan with the same
!    routine.
!  * The linear-molecule branch is analytic: pg_is_linear (frame.irp.f)
!    flags it and pg_build_linear builds C31v or C31v x {E,sigma_h}
!    directly, because the (permutation, determinant) key is blind to
!    every rotation about a collinear axis.
!  * pg_has_i, pg_n_sigma_h, pg_n_sigma_vtot, pg_n_c2_perp are providers
!    so that the driver can report the measured invariants. There is no
!    sigma_v/sigma_d split: it is a labelling device decided from the
!    axis slot in columns.irp.f, and it is not needed to name a group.

! Largest group order in the data set. D32d and D32h have :math:`|G|` = 128
! (checked over the 211 files; Ih is 120), so the prototype's 120 silently
! truncates two groups. This is the module-interface bound; the providers
! size their scratch arrays with allocate() rather
! than an explicit bound, because IRPF90 copies a file-scope parameter used
! as a provider declaration bound into the generated unit, where it clashes
! with this module-level copy.
integer, parameter :: pg_n_elem_max = 160

!-----------------------------------------------------------------------
! Group construction: the linear molecules take an analytic branch.
!-----------------------------------------------------------------------
subroutine pg_build_group(nat, xyz, znum, ctr, evecs, frame, linear, tol, &
     nmax, rot, perm, n_rot, ierr)
  implicit none
  BEGIN_DOC
  ! Build the group of the molecule and store it as matrices in the INPUT
  ! frame plus the induced atom permutations.
  !
  ! ``frame`` is :c:data:`pg_rot` (rows = canonical axes in the input frame),
  ! used by the linear branch only. ``evecs`` is :c:data:`pg_inertia_evec`,
  ! used by the generic branch only. This dispatcher exists because IRPF90
  ! does not look inside a subroutine body: the two callers
  ! (:c:data:`pg_n_elem`, :c:data:`pg_rot_mat`) name :c:data:`pg_is_linear` and
  ! :c:data:`pg_rot` themselves so that the providers are ordered.
  END_DOC
  integer, intent(in) :: nat, nmax
  double precision, intent(in) :: xyz(nat,3), ctr(3), evecs(3,3), frame(3,3), tol
  integer, intent(in) :: znum(nat)
  logical, intent(in) :: linear
  double precision, intent(out) :: rot(3,3,nmax)
  integer, intent(out) :: perm(nat,nmax)
  integer, intent(out) :: n_rot, ierr

  if (linear) then
    call pg_build_linear(nat, xyz, znum, ctr, frame, tol, nmax, &
         rot, perm, n_rot, ierr)
  else
    call pg_enumerate(nat, xyz, znum, ctr, evecs, tol, nmax, &
         rot, perm, n_rot, ierr)
  endif
end subroutine pg_build_group

!-----------------------------------------------------------------------
! Canonical -> input frame: R = frame^T R_c frame.
!-----------------------------------------------------------------------
subroutine pg_linear_to_input(frame, rc, r)
  implicit none
  BEGIN_DOC
  ! Rows of ``frame`` are the canonical axes expressed in the input frame,
  ! so ``s = frame (x - center)`` are the canonical coordinates and the
  ! input-frame matrix of a canonical operator ``R_c`` is
  ! ``R = frame^T R_c frame``.
  END_DOC
  double precision, intent(in) :: frame(3,3), rc(3,3)
  double precision, intent(out) :: r(3,3)

  r(1:3,1:3) = matmul(transpose(frame(1:3,1:3)), matmul(rc(1:3,1:3), &
       frame(1:3,1:3)))
end subroutine pg_linear_to_input

!-----------------------------------------------------------------------
! Append one canonical element to the linear group.
!-----------------------------------------------------------------------
subroutine pg_linear_add(nat, xyz, znum, ctr, frame, rc, tol, nmax, &
     rot, perm, n_rot, ierr)
  implicit none
  BEGIN_DOC
  ! Transform ``rc`` to the input frame, check that it maps the molecule
  ! onto itself, and append it. ``ierr`` = 1 when ``nmax`` is exhausted, 2
  ! when the element does not realise a permutation.
  END_DOC
  integer, intent(in) :: nat, nmax
  double precision, intent(in) :: xyz(nat,3), ctr(3), frame(3,3), rc(3,3), tol
  integer, intent(in) :: znum(nat)
  double precision, intent(inout) :: rot(3,3,nmax)
  integer, intent(inout) :: perm(nat,nmax)
  integer, intent(inout) :: n_rot, ierr
  double precision :: r(3,3)
  logical :: ok

  if (n_rot >= nmax) then
    ierr = 1
    return
  endif
  call pg_linear_to_input(frame, rc, r)
  call pg_maps_to_self(nat, xyz, znum, ctr, r, tol, perm(1,n_rot+1), ok)
  if (.not. ok) then
    ierr = 2
    return
  endif
  n_rot = n_rot + 1
  rot(1:3,1:3,n_rot) = r(1:3,1:3)
end subroutine pg_linear_add

!-----------------------------------------------------------------------
! Group of a linear molecule, built analytically in the canonical frame.
!-----------------------------------------------------------------------
subroutine pg_build_linear(nat, xyz, znum, ctr, frame, tol, nmax, &
     rot, perm, n_rot, ierr)
  implicit none
  BEGIN_DOC
  ! The finite group used for a linear molecule is C31v, extended to
  ! C31v x {E, sigma_h} when sigma_h maps the molecule onto itself. The
  ! order is 62 or 124.
  !
  ! Why not enumerate: for a collinear geometry every rotation about the
  ! axis and every reflection in a plane containing it induces the identity
  ! permutation, so the (permutation, determinant) key of
  ! :c:func:`pg_add_element` collapses the 31 sigma_v planes with sigma_h and
  ! the enumeration keeps only {E, one sigma_v}. The matrices are therefore
  ! built analytically from the canonical frame of :c:data:`pg_rot`: every
  ! angle is an exact multiple of 2 pi / 31, so the inverse, product and
  ! conjugation lookups by matrix (1e-8) are stable.
  !
  ! C31v, not Cinfv, and C31v x {E,sigma_h}, not C31v x Ci: the tables
  ! shipped are C31v and D31h. With the inversion the improper elements
  ! would be i C31^k = S62^{31+2k}, k = 1..15: fifteen columns named
  ! S62^33 .. S62^61, a notation no other table uses, and the inversion
  ! itself would be S62^31, an unreduced symbol (gcd(31,62) = 31) that the
  ! column parser rejects. The horizontal plane keeps sigma_h and S31^k.
  ! The rows of the D31h table are the restrictions of the Dinfh irreps to
  ! this group, so Sigma_g+, Pi_u, ... keep their chemical meaning.
  !
  ! ``ierr`` = 1 when the group does not fit in ``nmax``, 2 when a built
  ! element fails :c:func:`pg_maps_to_self`.
  END_DOC
  integer, intent(in) :: nat, nmax
  double precision, intent(in) :: xyz(nat,3), ctr(3), frame(3,3), tol
  integer, intent(in) :: znum(nat)
  double precision, intent(out) :: rot(3,3,nmax)
  integer, intent(out) :: perm(nat,nmax)
  integer, intent(out) :: n_rot, ierr
  double precision, parameter :: pi = 3.1415926535897932d0
  integer, parameter :: n = 31
  double precision :: rc(3,3), r(3,3), p(3), c, s, ang
  integer :: k, j, i
  logical :: ok, has_h

  n_rot = 0
  ierr = 0

  ! Is there an inversion centre? sigma_h (z -> -z) is central in the point
  ! group of a linear molecule, so it either maps the molecule onto itself
  ! (Dinfh) or it does not (Cinfv). Tested once, on the exact matrix.
  rc(1:3,1:3) = 0.d0
  rc(1,1) = 1.d0
  rc(2,2) = 1.d0
  rc(3,3) = -1.d0
  call pg_linear_to_input(frame, rc, r)
  call pg_maps_to_self(nat, xyz, znum, ctr, r, tol, perm(1,1), ok)
  has_h = ok

  ! E, then C31^k.
  rc(1:3,1:3) = 0.d0
  rc(1,1) = 1.d0
  rc(2,2) = 1.d0
  rc(3,3) = 1.d0
  call pg_linear_add(nat, xyz, znum, ctr, frame, rc, tol, nmax, &
       rot, perm, n_rot, ierr)
  if (ierr /= 0) return
  do k = 1, n - 1
    ang = 2.d0*pi*dble(k)/dble(n)
    c = cos(ang)
    s = sin(ang)
    rc(1:3,1:3) = 0.d0
    rc(1,1) = c
    rc(1,2) = -s
    rc(2,1) = s
    rc(2,2) = c
    rc(3,3) = 1.d0
    call pg_linear_add(nat, xyz, znum, ctr, frame, rc, tol, nmax, &
         rot, perm, n_rot, ierr)
    if (ierr /= 0) return
  enddo

  ! sigma_v_j: the plane through z with normal p = z x m_j, where m_j is the
  ! in-plane direction at 2 pi j / 31 from x. p = (-sin, cos, 0) and the
  ! reflection is 1 - 2 p p^T.
  do j = 0, n - 1
    ang = 2.d0*pi*dble(j)/dble(n)
    c = cos(ang)
    s = sin(ang)
    p(1) = -s
    p(2) = c
    p(3) = 0.d0
    do i = 1, 3
      do k = 1, 3
        rc(i,k) = -2.d0*p(i)*p(k)
      enddo
      rc(i,i) = rc(i,i) + 1.d0
    enddo
    call pg_linear_add(nat, xyz, znum, ctr, frame, rc, tol, nmax, &
         rot, perm, n_rot, ierr)
    if (ierr /= 0) return
  enddo
  if (.not. has_h) return

  ! sigma_h, S31^k = sigma_h C31^k, C2'_j = sigma_h sigma_v_j (a rotation
  ! by pi about m_j: 2 m m^T - 1).
  rc(1:3,1:3) = 0.d0
  rc(1,1) = 1.d0
  rc(2,2) = 1.d0
  rc(3,3) = -1.d0
  call pg_linear_add(nat, xyz, znum, ctr, frame, rc, tol, nmax, &
       rot, perm, n_rot, ierr)
  if (ierr /= 0) return
  do k = 1, n - 1
    ang = 2.d0*pi*dble(k)/dble(n)
    c = cos(ang)
    s = sin(ang)
    rc(1:3,1:3) = 0.d0
    rc(1,1) = c
    rc(1,2) = -s
    rc(2,1) = s
    rc(2,2) = c
    rc(3,3) = -1.d0
    call pg_linear_add(nat, xyz, znum, ctr, frame, rc, tol, nmax, &
         rot, perm, n_rot, ierr)
    if (ierr /= 0) return
  enddo
  do j = 0, n - 1
    ang = 2.d0*pi*dble(j)/dble(n)
    p(1) = cos(ang)
    p(2) = sin(ang)
    p(3) = 0.d0
    do i = 1, 3
      do k = 1, 3
        rc(i,k) = 2.d0*p(i)*p(k)
      enddo
      rc(i,i) = rc(i,i) - 1.d0
    enddo
    call pg_linear_add(nat, xyz, znum, ctr, frame, rc, tol, nmax, &
         rot, perm, n_rot, ierr)
    if (ierr /= 0) return
  enddo
end subroutine pg_build_linear

!-----------------------------------------------------------------------
! Enumerate every group element.
!  * identity is element 1
!  * candidate directions: the axes measured in frame.irp.f
!  * proper rotations of order 32..2, reflections, inversion
!  * closure under matrix multiplication, permutations composed
!-----------------------------------------------------------------------
subroutine pg_enumerate(nat, xyz, znum, ctr, evecs, tol, nmax, &
     rot, perm, n_rot, ierr)
  implicit none
  BEGIN_DOC
  ! Build the full group from the geometry. ``ierr`` = 1 when ``nmax`` is too
  ! small, 2 when the candidate list overflows, 3 on a non-orthogonal
  ! product (should never happen).
  END_DOC
  integer, parameter :: ncandmax = 2000
  integer, intent(in) :: nat, nmax
  double precision, intent(in) :: xyz(nat,3), ctr(3), evecs(3,3), tol
  integer, intent(in) :: znum(nat)
  double precision, intent(out) :: rot(3,3,nmax)
  integer, intent(out) :: perm(nat,nmax)
  integer, intent(out) :: n_rot, ierr
  double precision, parameter :: pi = 3.1415926535897932d0
  double precision :: cand(3,ncandmax), r(3,3)
  integer :: perm_i(nat), pnew(nat)
  integer :: ncand, i, j, k, n, nold
  logical :: ok, changed
  double precision :: a

  n_rot = 1
  ierr = 0
  rot(1:3,1:3,1) = 0.d0
  rot(1,1,1) = 1.d0
  rot(2,2,1) = 1.d0
  rot(3,3,1) = 1.d0
  do i = 1, nat
    perm(i,1) = i
  enddo

  call pg_build_candidates(nat, xyz, znum, ctr, evecs, ncandmax, &
       cand, ncand, ierr)
  if (ierr /= 0) return

  ! Proper rotations. 32 is the highest order in the data set (C32 in
  ! D32h/D32d); a cap of 12 mislabels every molecule with a C_n axis, n>12.
  do k = 1, ncand
    do n = 32, 2, -1
      call pg_rotation(cand(1,k), 2.d0*pi/dble(n), r)
      call pg_maps_to_self(nat, xyz, znum, ctr, r, tol, perm_i, ok)
      if (ok) then
        call pg_add_element(r, perm_i, nat, nmax, rot, perm, n_rot, ierr)
        if (ierr /= 0) return
      endif
    enddo
  enddo

  ! Reflections.
  do k = 1, ncand
    call pg_reflection(cand(1,k), r)
    call pg_maps_to_self(nat, xyz, znum, ctr, r, tol, perm_i, ok)
    if (ok) then
      call pg_add_element(r, perm_i, nat, nmax, rot, perm, n_rot, ierr)
      if (ierr /= 0) return
    endif
  enddo

  ! Inversion.
  r(1:3,1:3) = 0.d0
  r(1,1) = -1.d0
  r(2,2) = -1.d0
  r(3,3) = -1.d0
  call pg_maps_to_self(nat, xyz, znum, ctr, r, tol, perm_i, ok)
  if (ok) then
    call pg_add_element(r, perm_i, nat, nmax, rot, perm, n_rot, ierr)
    if (ierr /= 0) return
  endif

  ! Closure. Products are composed from the stored permutations, never
  ! re-matched against the geometry.
  changed = .true.
  do while (changed)
    changed = .false.
    nold = n_rot
    do i = 1, nold
      do j = 1, nold
        r(1:3,1:3) = matmul(rot(1:3,1:3,i), rot(1:3,1:3,j))
        do k = 1, nat
          pnew(k) = perm(perm(k,j),i)
        enddo
        call pg_add_element(r, pnew, nat, nmax, rot, perm, n_rot, ierr)
        if (ierr /= 0) return
        if (n_rot > nold) then
          changed = .true.
          nold = n_rot
        endif
      enddo
    enddo
  enddo

  ! Orthogonality sanity check.
  do k = 1, n_rot
    a = maxval(abs(matmul(transpose(rot(1:3,1:3,k)), rot(1:3,1:3,k))))
    if (abs(a - 1.d0) > 1.d-8) then
      ierr = 3
      return
    endif
  enddo
end subroutine pg_enumerate

!-----------------------------------------------------------------------
! Determinant of a 3x3 matrix.
!-----------------------------------------------------------------------
double precision function pg_det3(r)
  implicit none
  BEGIN_DOC
  ! Determinant of a 3x3 matrix.
  END_DOC
  double precision, intent(in) :: r(3,3)

  pg_det3 = r(1,1)*(r(2,2)*r(3,3) - r(2,3)*r(3,2)) &
          - r(1,2)*(r(2,1)*r(3,3) - r(2,3)*r(3,1)) &
          + r(1,3)*(r(2,1)*r(3,2) - r(2,2)*r(3,1))
end function pg_det3

!-----------------------------------------------------------------------
! Append an element unless it is already stored.
!-----------------------------------------------------------------------
subroutine pg_add_element(r, p, nat, nmax, rot, perm, n_rot, ierr)
  implicit none
  BEGIN_DOC
  ! Add (r,p) to the list. Identical matrices to 1e-8 are duplicates;
  ! beyond that, two elements are the same group element when they induce
  ! the same atom permutation *and* have the same determinant. If two
  ! orthogonal maps permute the atoms identically, their quotient fixes
  ! every atom, hence is the identity unless the molecule is planar, in
  ! which case it is the reflection in the molecular plane (the only
  ! orthogonal map fixing a whole plane pointwise), which differs by its
  ! determinant. Linear molecules, the other degenerate case, are
  ! rejected before enumeration.
  !
  ! The permutation test is required, not an optimisation: a candidate
  ! direction read off a geometry that is symmetric only to :c:data:`pg_tol` is
  ! tilted by O(:c:data:`pg_tol` / lever arm), so a product of two accepted
  ! elements differs from the stored element by far more than 1e-8 and
  ! would be appended again, without bound, until the array overflows.
  END_DOC
  integer, intent(in) :: nat, nmax
  double precision, intent(in) :: r(3,3)
  integer, intent(in) :: p(nat)
  double precision, intent(inout) :: rot(3,3,nmax)
  integer, intent(inout) :: perm(nat,nmax)
  integer, intent(inout) :: n_rot, ierr
  double precision :: rr(3,3)
  integer :: k, i, det
  logical :: same_perm
  double precision, external :: pg_det3

  det = nint(sign(1.d0, pg_det3(r)))
  do k = 1, n_rot
    if (maxval(abs(r(1:3,1:3) - rot(1:3,1:3,k))) <= 1.d-8) return
    rr(1:3,1:3) = rot(1:3,1:3,k)
    if (nint(sign(1.d0, pg_det3(rr))) /= det) cycle
    same_perm = .true.
    do i = 1, nat
      if (p(i) /= perm(i,k)) then
        same_perm = .false.
        exit
      endif
    enddo
    if (same_perm) return
  enddo
  if (n_rot >= nmax) then
    ierr = 1
    return
  endif
  n_rot = n_rot + 1
  rot(1:3,1:3,n_rot) = r(1:3,1:3)
  do i = 1, nat
    perm(i,n_rot) = p(i)
  enddo
end subroutine pg_add_element

!-----------------------------------------------------------------------
! (det, angle, axis) of an orthogonal 3x3 matrix. Angle in [0,2pi).
!-----------------------------------------------------------------------
subroutine pg_elem_props(r, angle, axis, det)
  implicit none
  BEGIN_DOC
  ! Angle and axis: for an improper element
  ! the angle is theta_p - pi, where theta_p is the angle of -r. The
  ! oriented angle is obtained from atan2 of the antisymmetric part, not
  ! from acos, which is why S4 gives theta_p = 3pi/2 and theta = pi/2.
  END_DOC
  double precision, intent(in) :: r(3,3)
  double precision, intent(out) :: angle, axis(3)
  integer, intent(out) :: det
  double precision :: rp(3,3), m(3,3), row(3,3), c(3,3), a(3), t(3)
  double precision :: cp, sp, tp, d, dn
  integer :: i, j
  double precision, external :: pg_det3

  a(1:3) = 0.d0
  det = nint(sign(1.d0, pg_det3(r)))

  if (det > 0) then
    rp(1:3,1:3) = r(1:3,1:3)
  else
    rp(1:3,1:3) = -r(1:3,1:3)
  endif

  ! Axis: cross products of the rows of rp - 1 (all perpendicular to the
  ! rotation axis), largest norm first.
  m(1:3,1:3) = rp(1:3,1:3)
  do i = 1, 3
    m(i,i) = m(i,i) - 1.d0
  enddo
  do i = 1, 3
    row(i,1:3) = m(i,1:3)
  enddo
  dn = 0.d0
  do i = 1, 2
    do j = i + 1, 3
      call pg_cross(row(i,1:3), row(j,1:3), t)
      d = sqrt(t(1)**2 + t(2)**2 + t(3)**2)
      if (d > dn) then
        dn = d
        a(1:3) = t(1:3)
      endif
    enddo
  enddo

  if (dn < 1.d-10) then
    ! Small angle (or identity): use the axial vector sin(theta) n.
    a(1) = 0.5d0*(rp(3,2) - rp(2,3))
    a(2) = 0.5d0*(rp(1,3) - rp(3,1))
    a(3) = 0.5d0*(rp(2,1) - rp(1,2))
    dn = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
    if (dn < 1.d-10) then
      ! theta_p = 0 (identity): any axis; use columns of rp + 1.
      do i = 1, 3
        do j = 1, 3
          c(j,i) = rp(j,i)
        enddo
        c(i,i) = c(i,i) + 1.d0
      enddo
      dn = 0.d0
      do i = 1, 3
        d = sqrt(c(1,i)**2 + c(2,i)**2 + c(3,i)**2)
        if (d > dn) then
          dn = d
          a(1:3) = c(1:3,i)
        endif
      enddo
      if (dn < 1.d-10) then
        a(1:3) = (/ 0.d0, 0.d0, 1.d0 /)
        dn = 1.d0
      endif
    endif
  endif

  axis(1:3) = a(1:3)/sqrt(a(1)**2 + a(2)**2 + a(3)**2)
  call pg_fix_sign(axis)

  cp = 0.5d0*(rp(1,1) + rp(2,2) + rp(3,3) - 1.d0)
  cp = max(-1.d0, min(1.d0, cp))
  sp = 0.5d0*((rp(3,2)-rp(2,3))*axis(1) &
            + (rp(1,3)-rp(3,1))*axis(2) &
            + (rp(2,1)-rp(1,2))*axis(3))
  tp = atan2(sp, cp)
  if (tp < 0.d0) tp = tp + 2.d0*3.1415926535897932d0
  if (det > 0) then
    angle = tp
  else
    angle = tp - 3.1415926535897932d0
  endif
  do while (angle < 0.d0)
    angle = angle + 2.d0*3.1415926535897932d0
  enddo
  do while (angle >= 2.d0*3.1415926535897932d0)
    angle = angle - 2.d0*3.1415926535897932d0
  enddo
end subroutine pg_elem_props

!-----------------------------------------------------------------------
! Is the axis of an element on a cartesian direction of the canonical frame?
!-----------------------------------------------------------------------
logical function pg_axis_along(axis, k)
  implicit none
  BEGIN_DOC
  ! True when the unit vector ``axis`` (INPUT frame) lies on the ``k``-th
  ! canonical axis (``k`` = 1, 2, 3 for x, y, z). The canonical direction in
  ! the input frame is row ``k`` of pg_rot, so this is a dot
  ! product, exactly as the identification invariants do it with row 3.
  !
  ! This is a LABELLING test and it is deliberately tight (1e-4, 0.8
  ! degree): it must answer "is this element one of the cartesian axes",
  ! not "is it roughly aligned". The closest non-axis direction in the data
  ! set is the sigma_d normal of D32d at pi/(2*32) = 2.8 degrees from x,
  ! ``cos`` = 0.9988, while an element that IS an axis has ``cos`` = 1 up to the
  ! tilt between two candidate directions built from the same geometry: for
  ! a geometry symmetric to 1e-12 that tilt is 1e-9, and it stays below 1e-7
  ! for any :c:data:`pg_tol` up to 1e-2.
  !
  ! The identification invariants use 1/2 instead: they must survive a
  ! geometry that is symmetric only within :c:data:`pg_tol`.
  END_DOC
  double precision, intent(in) :: axis(3)
  integer, intent(in) :: k
  double precision, parameter :: tol = 1.d-4
  double precision :: d

  d = axis(1)*pg_rot(k,1) + axis(2)*pg_rot(k,2) + axis(3)*pg_rot(k,3)
  pg_axis_along = abs(d) > 1.d0 - tol
end function pg_axis_along

!-----------------------------------------------------------------------
! Reduce an angle to the smallest m <= 64 and k coprime with m.
!-----------------------------------------------------------------------
subroutine pg_reduce_angle(angle, m, k, ierr)
  implicit none
  BEGIN_DOC
  ! Rational reduction of an angle: ``angle - 2 pi k/m`` < 1e-6.
  ! With this convention a reflection is (m=1,k=0) and the inversion is
  ! (m=2,k=1) (S1 = sigma, S2 = i).
  END_DOC
  double precision, intent(in) :: angle
  integer, intent(out) :: m, k, ierr
  double precision, parameter :: pi = 3.1415926535897932d0
  double precision :: d
  integer :: mm, kk, a, b, t

  ierr = 1
  m = 0
  k = 0
  ! 64 is required, not cosmetic: D32d contains S64 (k=1, m=64). A cap of 32
  ! makes its improper elements unrationalisable and the table unmatchable.
  do mm = 1, 64
    do kk = 0, mm - 1
      a = kk
      b = mm
      do while (b /= 0)
        t = mod(a, b)
        a = b
        b = t
      enddo
      if (a /= 1) cycle
      d = abs(angle - 2.d0*pi*dble(kk)/dble(mm))
      d = min(d, abs(d - 2.d0*pi))
      if (d < 1.d-6) then
        m = mm
        k = kk
        ierr = 0
        return
      endif
    enddo
  enddo
end subroutine pg_reduce_angle

!-----------------------------------------------------------------------
! Integer key encoding (det, m, k). det = +1 -> even, det = -1 -> odd.
!-----------------------------------------------------------------------
subroutine pg_encode_kind(det, m, k, kind)
  implicit none
  BEGIN_DOC
  ! Injective integer encoding of (det,m,k): kind = 2*(64*m + k) + (1-det)/2,
  ! with m <= 64 and k < m. For m <= 64 the pair (m,k) is recovered as
  ! k = mod(kind/2, 64), m = (kind/2 - k)/64.
  END_DOC
  integer, intent(in) :: det, m, k
  integer, intent(out) :: kind
  kind = 2*(64*m + k) + (1 - det)/2
end subroutine pg_encode_kind

subroutine pg_decode_kind(kind, det, m, k)
  implicit none
  BEGIN_DOC
  ! Inverse of pg_encode_kind.
  END_DOC
  integer, intent(in) :: kind
  integer, intent(out) :: det, m, k
  k = mod(kind/2, 64)
  m = (kind/2 - k)/64
  det = 1 - 2*mod(kind, 2)
end subroutine pg_decode_kind

!-----------------------------------------------------------------------
! Closed-form group order of a name, used as a cross-check.
!-----------------------------------------------------------------------
subroutine pg_expected_order(name, order, ierr)
  implicit none
  BEGIN_DOC
  ! ``C_n`` = n, ``C_nv`` = ``C_nh`` = 2n, ``S_n`` = n, ``D_n`` = 2n,
  ! ``D_nh`` = ``D_nd`` = 4n, :math:`|T|` = 12, ``Th`` = ``Td`` = :math:`|O|` = 24,
  ! ``Oh`` = 48, :math:`|I|` = 60, ``Ih`` = 120.
  END_DOC
  character*(*), intent(in) :: name
  integer, intent(out) :: order, ierr
  integer :: i, n

  ierr = 1
  order = 0
  n = len_trim(name)
  if (n == 0) return

  if (name(1:1) == 'T') then
    if (n == 1) then
      order = 12
    else if (n == 2 .and. (name(2:2) == 'h' .or. name(2:2) == 'd')) then
      order = 24
    endif
  else if (name(1:1) == 'O') then
    if (n == 1) then
      order = 24
    else if (n == 2 .and. name(2:2) == 'h') then
      order = 48
    endif
  else if (name(1:1) == 'I') then
    if (n == 1) then
      order = 60
    else if (n == 2 .and. name(2:2) == 'h') then
      order = 120
    endif
  else if (name(1:1) == 'C') then
    if (n == 1) then
      return
    else if (name(2:2) == 'i') then
      order = 2
    else if (name(2:2) == 's') then
      order = 2
    else
      i = 2
      call pg_read_int(name, i, n, order)
      if (order > 0) then
        if (i > n) then
          ! C_n: order n
        else if (i == n .and. (name(i:i) == 'h' .or. name(i:i) == 'v')) then
          order = 2*order
        else
          order = 0
        endif
      endif
    endif
  else if (name(1:1) == 'D') then
    i = 2
    call pg_read_int(name, i, n, order)
    if (order > 0) then
      if (i > n) then
        order = 2*order
      else if (i == n .and. (name(i:i) == 'h' .or. name(i:i) == 'd')) then
        order = 4*order
      else
        order = 0
      endif
    endif
  else if (name(1:1) == 'S') then
    i = 2
    call pg_read_int(name, i, n, order)
    if (order > 0 .and. i <= n) order = 0
  endif

  if (order > 0) ierr = 0
end subroutine pg_expected_order

subroutine pg_read_int(name, i, n, value)
  implicit none
  BEGIN_DOC
  ! Read a decimal integer starting at name(i); advance i past it.
  END_DOC
  character*(*), intent(in) :: name
  integer, intent(inout) :: i
  integer, intent(in) :: n
  integer, intent(out) :: value
  integer :: c

  value = 0
  do while (i <= n)
    c = ichar(name(i:i)) - 48
    if (c < 0 .or. c > 9) exit
    value = 10*value + c
    i = i + 1
  enddo
end subroutine pg_read_int

!=======================================================================
! Providers
!=======================================================================

 BEGIN_PROVIDER [ integer, pg_n_max ]
&BEGIN_PROVIDER [ integer, pg_n_axes_max ]
  implicit none
  BEGIN_DOC
  ! Highest proper-rotation order, and number of distinct lines carrying an
  ! axis of that order.
  END_DOC
  integer :: znum(nucl_num), ierr, i
  double precision :: axis(3)
  integer :: nmax, naxes

  ! Defaults before the first return: an early return still sets the built
  ! flag, so what is written here is cached until a TOUCH.
  pg_n_max = 0
  pg_n_axes_max = 0
  if (pg_detection_failed) return

  do i = 1, nucl_num
    znum(i) = nint(nucl_charge(i))
  enddo
  call pg_proper_scan(nucl_num, nucl_coord, znum, pg_center, &
       pg_inertia_eval, pg_inertia_evec, pg_tol, 2000, &
       nmax, naxes, axis, pg_is_linear, ierr)
  if (ierr /= 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A)') 'symmetry: candidate-direction scan failed'
      pg_detection_failed = .true.
    endif
    return
  endif
  pg_n_max = nmax
  pg_n_axes_max = naxes
END_PROVIDER

 BEGIN_PROVIDER [ integer, pg_n_elem ]
  implicit none
  BEGIN_DOC
  ! Order of the molecular point group.
  END_DOC
  integer :: znum(nucl_num), ierr, i, n
  double precision, allocatable :: rot(:,:,:)
  integer, allocatable :: perm(:,:)
  double precision :: frame(3,3)
  logical :: linear

  ! :c:data:`pg_n_elem` keeps its 0 default unless the enumeration succeeded: the
  ! result goes to a local first, so a failed call cannot leave a
  ! half-enumerated order behind.
  pg_n_elem = 0
  if (pg_detection_failed) return

  allocate(rot(3,3,pg_n_elem_max), perm(nucl_num,pg_n_elem_max))
  do i = 1, nucl_num
    znum(i) = nint(nucl_charge(i))
  enddo
  ! IRPF90 does not look inside pg_build_group: naming the two providers here
  ! is what orders them before this provider.
  linear = pg_is_linear
  frame(1:3,1:3) = pg_rot(1:3,1:3)
  call pg_build_group(nucl_num, nucl_coord, znum, pg_center, &
       pg_inertia_evec, frame, linear, pg_tol, pg_n_elem_max, &
       rot, perm, n, ierr)
  if (ierr == 1) then
    if (.not. pg_detection_failed) then
      write(6,'(A,I0,A)') 'symmetry: group order exceeds ', pg_n_elem_max, &
           ' (pg_n_elem_max)'
      pg_detection_failed = .true.
    endif
    return
  else if (ierr == 2 .and. linear) then
    if (.not. pg_detection_failed) then
      write(6,'(A)') 'symmetry: an element built for a linear molecule '// &
           'does not map it onto itself'
      pg_detection_failed = .true.
    endif
    return
  else if (ierr /= 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A)') 'symmetry: group enumeration failed'
      pg_detection_failed = .true.
    endif
    return
  endif
  pg_n_elem = n
END_PROVIDER

 BEGIN_PROVIDER [ double precision, pg_rot_mat, (3,3,pg_n_elem) ]
&BEGIN_PROVIDER [ integer, pg_perm, (nucl_num,pg_n_elem) ]
&BEGIN_PROVIDER [ double precision, pg_elem_angle, (pg_n_elem) ]
&BEGIN_PROVIDER [ double precision, pg_elem_axis, (3,pg_n_elem) ]
&BEGIN_PROVIDER [ integer, pg_elem_det, (pg_n_elem) ]
&BEGIN_PROVIDER [ integer, pg_elem_kind, (pg_n_elem) ]
  implicit none
  BEGIN_DOC
  ! Rotation matrices (input frame), induced atom permutations, angle,
  ! axis (plane normal for a reflection), determinant sign, and the
  ! integer (det,m,k) key of every group element.
  END_DOC
  integer :: znum(nucl_num), ierr, i, j, m, k
  double precision, allocatable :: rot(:,:,:)
  double precision :: axis(3), v(3), d, frame(3,3)
  integer, allocatable :: perm(:,:)
  integer :: n
  logical :: linear

  allocate(rot(3,3,pg_n_elem_max), perm(nucl_num,pg_n_elem_max))
  do i = 1, nucl_num
    znum(i) = nint(nucl_charge(i))
  enddo
  linear = pg_is_linear
  frame(1:3,1:3) = pg_rot(1:3,1:3)
  call pg_build_group(nucl_num, nucl_coord, znum, pg_center, &
       pg_inertia_evec, frame, linear, pg_tol, pg_n_elem_max, &
       rot, perm, n, ierr)
  if (ierr /= 0 .or. n /= pg_n_elem) then
    if (.not. pg_detection_failed) then
      write(6,'(A)') 'symmetry: inconsistent group enumeration'
      pg_detection_failed = .true.
    endif
    return
  endif

  do i = 1, n
    pg_rot_mat(1:3,1:3,i) = rot(1:3,1:3,i)
    do j = 1, nucl_num
      pg_perm(j,i) = perm(j,i)
    enddo
    ! Consistency: the permutation must realise the matrix.
    do j = 1, nucl_num
      v(1:3) = nucl_coord(j,1:3) - pg_center(1:3)
      d = (rot(1,1,i)*v(1) + rot(1,2,i)*v(2) + rot(1,3,i)*v(3) &
           + pg_center(1) - nucl_coord(perm(j,i),1))**2 &
        + (rot(2,1,i)*v(1) + rot(2,2,i)*v(2) + rot(2,3,i)*v(3) &
           + pg_center(2) - nucl_coord(perm(j,i),2))**2 &
        + (rot(3,1,i)*v(1) + rot(3,2,i)*v(2) + rot(3,3,i)*v(3) &
           + pg_center(3) - nucl_coord(perm(j,i),3))**2
      if (sqrt(d) > 10.d0*pg_tol) then
        if (.not. pg_detection_failed) then
          write(6,'(A,I0,A,I0)') 'symmetry: element ', i, &
               ' does not realise its permutation (atom ', j, ')'
          pg_detection_failed = .true.
        endif
        return
      endif
    enddo

    call pg_elem_props(rot(1:3,1:3,i), pg_elem_angle(i), axis, pg_elem_det(i))
    ! pg_elem_axis is the eigenvector of the stored matrix, and pg_rot_mat
    ! is in the input frame: keep pg_elem_axis in that same frame. The
    ! identification then compares it with pg_rot(3,:) (the canonical z
    ! expressed in the input frame); the column matching and the MO
    ! projectors must apply pg_rot if they want canonical components.
    pg_elem_axis(1:3,i) = axis(1:3)
    call pg_reduce_angle(pg_elem_angle(i), m, k, ierr)
    if (ierr /= 0) then
      if (.not. pg_detection_failed) then
        write(6,'(A,F14.8)') 'symmetry: cannot rationalise element angle ', &
             pg_elem_angle(i)
        pg_detection_failed = .true.
      endif
      return
    endif
    call pg_encode_kind(pg_elem_det(i), m, k, pg_elem_kind(i))
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! Measured symmetry breaking of the stored geometry.
! Uses the stored matrices and permutations: no rotation is rebuilt here
! and pg_maps_to_self is not consulted, so the number reports the input,
! not the matching. It is what tells a G1 failure apart: residual ~1e-16
! with a large G1 is a code bug, residual ~1e-8 is the fixture.
! Units are bohr, like nucl_coord.
!-----------------------------------------------------------------------
 BEGIN_PROVIDER [ double precision, pg_geom_residual ]
  implicit none
  BEGIN_DOC
  ! max over the elements g and the atoms i of
  ! ``R_g (r_i - center) + center - r_{pg_perm(i,g)}``, in bohr.
  END_DOC
  integer :: ig, i
  double precision :: v(3), e(3)

  pg_geom_residual = 0.d0
  do ig = 1, pg_n_elem
    do i = 1, nucl_num
      v(1:3) = nucl_coord(i,1:3) - pg_center(1:3)
      e(1) = pg_rot_mat(1,1,ig)*v(1) + pg_rot_mat(1,2,ig)*v(2) &
           + pg_rot_mat(1,3,ig)*v(3) + pg_center(1) &
           - nucl_coord(pg_perm(i,ig),1)
      e(2) = pg_rot_mat(2,1,ig)*v(1) + pg_rot_mat(2,2,ig)*v(2) &
           + pg_rot_mat(2,3,ig)*v(3) + pg_center(2) &
           - nucl_coord(pg_perm(i,ig),2)
      e(3) = pg_rot_mat(3,1,ig)*v(1) + pg_rot_mat(3,2,ig)*v(2) &
           + pg_rot_mat(3,3,ig)*v(3) + pg_center(3) &
           - nucl_coord(pg_perm(i,ig),3)
      pg_geom_residual = max(pg_geom_residual, sqrt(e(1)**2+e(2)**2+e(3)**2))
    enddo
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [ integer, pg_elem_inv, (pg_n_elem) ]
  implicit none
  BEGIN_DOC
  ! Index of the inverse of every element: the inverse atom permutation with
  ! the same determinant. Exact, no tolerance: ``h`` is the inverse of ``g`` iff
  ! its permutation is the inverse one and ``det`` agrees, and the enumeration
  ! stores permutations, never re-matches matrices.
  !
  ! The column matching needs it to close a column under inversion: a
  ! geometric class is not closed under inversion (in C3 the class {C3}
  ! has its inverse in the class {C3^2}), while a file column always is.
  END_DOC
  integer :: inv_perm(nucl_num), i, j
  double precision :: rm(3,3)
  integer, external :: pg_find_by_perm
  integer, external :: pg_find_by_mat

  ! Linear molecules: the group was built analytically, so the matrices are
  ! exact and the inverse is the transpose. The permutation route cannot be
  ! used there: every element induces the identity permutation.
  if (pg_is_linear) then
    do i = 1, pg_n_elem
      rm(1:3,1:3) = transpose(pg_rot_mat(1:3,1:3,i))
      pg_elem_inv(i) = pg_find_by_mat(rm)
      if (pg_elem_inv(i) == 0) then
        if (.not. pg_detection_failed) then
          write(6,'(A,I0)') 'symmetry: no inverse found for element ', i
          pg_detection_failed = .true.
        endif
        return
      endif
    enddo
    return
  endif

  do i = 1, pg_n_elem
    do j = 1, nucl_num
      inv_perm(pg_perm(j,i)) = j
    enddo
    pg_elem_inv(i) = pg_find_by_perm(inv_perm, pg_elem_det(i))
    if (pg_elem_inv(i) == 0) then
      if (.not. pg_detection_failed) then
        write(6,'(A)') 'symmetry: no inverse found in the group'
        pg_detection_failed = .true.
      endif
      return
    endif
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [ integer, pg_elem_mul, (pg_n_elem,pg_n_elem) ]
  implicit none
  BEGIN_DOC
  ! Multiplication table: ``pg_elem_mul(g,h)`` is the index of R_g R_h.
  !
  ! Built from the atom permutations and the determinants, never from the
  ! matrices. A stored element is only a *representative* of a symmetry
  ! operation known to :c:data:`pg_tol`: a candidate direction read off the geometry
  ! is accepted when the atoms move by less than :c:data:`pg_tol`, which for a
  ! rotation about a slightly tilted axis means a matrix accurate to
  ! :c:data:`pg_tol` / lever arm ~ 1e-8, and a product of two such matrices differs
  ! from the stored representative by the same order. A matrix comparison
  ! with a 1e-8 threshold is then a coin flip (measured on NH3: the C3 axis
  ! is tilted by 4.9e-9, 18 of the 36 products miss). The (permutation,
  ! determinant) key is exact and is the identity criterion the enumeration
  ! itself uses.
  !
  ! The table is :math:`|G|^2 integers: 64 kB at the worst case :math:`|G|` = 128. Today the
  ! only consumer is gate G2; it is also what a product of two irreps will
  ! need once the module is in the core.
  END_DOC
  integer :: g, h, k, pnew(nucl_num)
  logical :: seen(pg_n_elem)
  double precision :: rm(3,3)
  integer, external :: pg_find_by_perm
  integer, external :: pg_find_by_mat

  do g = 1, pg_n_elem
    do h = 1, pg_n_elem
      if (pg_is_linear) then
        ! The linear group is analytic: the product of two stored matrices
        ! is within 1e-15 of a stored representative, so the lookup by
        ! matrix is exact in practice. The permutation route cannot be used:
        ! all the elements of a collinear geometry induce the identity.
        rm(1:3,1:3) = matmul(pg_rot_mat(1:3,1:3,g), pg_rot_mat(1:3,1:3,h))
        pg_elem_mul(g,h) = pg_find_by_mat(rm)
      else
        do k = 1, nucl_num
          pnew(k) = pg_perm(pg_perm(k,h), g)
        enddo
        pg_elem_mul(g,h) = pg_find_by_perm(pnew, &
             pg_elem_det(g)*pg_elem_det(h))
      endif
      if (pg_elem_mul(g,h) == 0) then
        if (.not. pg_detection_failed) then
          write(6,'(A,I0,A,I0)') 'symmetry: no product element for g = ', g, &
               ', h = ', h
          pg_detection_failed = .true.
        endif
        return
      endif
    enddo
  enddo

  ! A product table that is not a Latin square is not a group table. This
  ! catches a missing element or a duplicated (permutation,det) key, i.e. the
  ! two ways the enumeration could be wrong, without testing associativity
  ! (permutation composition is associative by construction).
  do g = 1, pg_n_elem
    seen(1:pg_n_elem) = .false.
    do h = 1, pg_n_elem
      if (seen(pg_elem_mul(g,h))) then
        if (.not. pg_detection_failed) then
          write(6,'(A,I0)') 'symmetry: product table is not a Latin '// &
               'square (g = ', g
          pg_detection_failed = .true.
        endif
        return
      endif
      seen(pg_elem_mul(g,h)) = .true.
    enddo
  enddo
 END_PROVIDER

!-----------------------------------------------------------------------
! Look up a group element by (atom permutation, determinant).
!-----------------------------------------------------------------------
integer function pg_find_by_perm(p, det)
  implicit none
  BEGIN_DOC
  ! Index of the unique element that induces the atom permutation ``p`` and has
  ! determinant ``det``, or 0. This is the enumeration's own identity
  ! criterion (:c:func:`pg_add_element`): two orthogonal maps that permute the atoms
  ! identically and share their determinant are the same element. No matrix
  ! is ever compared.
  END_DOC
  integer, intent(in) :: p(nucl_num), det
  integer :: j, k
  logical :: same

  pg_find_by_perm = 0
  do j = 1, pg_n_elem
    if (pg_elem_det(j) /= det) cycle
    same = .true.
    do k = 1, nucl_num
      if (pg_perm(k,j) /= p(k)) then
        same = .false.
        exit
      endif
    enddo
    if (same) then
      pg_find_by_perm = j
      return
    endif
  enddo
end function pg_find_by_perm

!-----------------------------------------------------------------------
! Look up a group element by its matrix.
!-----------------------------------------------------------------------
integer function pg_find_by_mat(r)
  implicit none
  BEGIN_DOC
  ! Index of the stored element whose matrix equals ``r`` to 1e-8, or 0.
  !
  ! Only for the linear branch, where the matrices are analytic (angles are
  ! exact multiples of 2 pi / 31) and a transpose or a product of two stored
  ! matrices is within 1e-15 of a stored one. The generic branch must keep
  ! :c:func:`pg_find_by_perm`: there the matrices carry the tilt of a candidate
  ! direction, and a 1e-8 matrix comparison is a coin flip.
  END_DOC
  double precision, intent(in) :: r(3,3)
  integer :: j

  pg_find_by_mat = 0
  do j = 1, pg_n_elem
    if (maxval(abs(pg_rot_mat(1:3,1:3,j) - r(1:3,1:3))) <= 1.d-8) then
      pg_find_by_mat = j
      return
    endif
  enddo
end function pg_find_by_mat

 BEGIN_PROVIDER [ integer, pg_n_class ]
&BEGIN_PROVIDER [ integer, pg_class, (pg_n_elem) ]
  implicit none
  BEGIN_DOC
  ! Classes by conjugation. Element 1 is the identity, so class 1 is the
  ! identity class. Direct orbit enumeration (:math:`|G|` <= 120).
  !
  ! The conjugation is done on the atom permutations, not on the
  ! matrices: the stored matrices carry the :c:data:`pg_tol`-sized tilt of the
  ! candidate directions, so a conjugated matrix is not within 1e-8 of any
  ! stored one. The permutation (with the determinant, to separate the
  ! identity from the reflection in the molecular plane) is the exact
  ! invariant used by :c:func:`pg_add_element` too. Direct orbit enumeration,
  ! :math:`|G|` <= 128 so :math:`|G|^2` = 16384 comparisons.
  !
  ! Linear molecules are the exception: all their permutations are the
  ! identity, so the conjugation is done on the analytic matrices there
  ! (:c:func:`pg_find_by_mat`).
  END_DOC
  integer :: p(nucl_num)
  integer :: ih, c, g, k, i, j
  logical :: same
  double precision :: rm(3,3)
  integer, external :: pg_find_by_mat

  pg_n_class = 0
  if (pg_detection_failed) return

  pg_class(1:pg_n_elem) = 0
  c = 0
  do g = 1, pg_n_elem
    if (pg_class(g) /= 0) cycle
    c = c + 1
    do k = 1, pg_n_elem
      if (pg_is_linear) then
        ! The permutation of every element of a collinear geometry is the
        ! identity, so the class is the orbit of the matrix under
        ! conjugation. The matrices are analytic: R_k R_g R_k^T is within
        ! 1e-15 of a stored element.
        rm(1:3,1:3) = matmul(matmul(pg_rot_mat(1:3,1:3,k), &
             pg_rot_mat(1:3,1:3,g)), transpose(pg_rot_mat(1:3,1:3,k)))
        ih = pg_find_by_mat(rm)
      else
        ! p_h = p_k o p_g o p_{k^-1}
        do i = 1, nucl_num
          p(i) = pg_perm(pg_perm(pg_perm(i,pg_elem_inv(k)), g), k)
        enddo
        ih = 0
        do i = 1, pg_n_elem
          if (pg_elem_det(i) /= pg_elem_det(g)) cycle
          same = .true.
          do j = 1, nucl_num
            if (pg_perm(j,i) /= p(j)) then
              same = .false.
              exit
            endif
          enddo
          if (same) then
            ih = i
            exit
          endif
        enddo
      endif
      if (ih == 0) then
        if (.not. pg_detection_failed) then
          write(6,'(A)') 'symmetry: conjugation left the group'
          pg_detection_failed = .true.
        endif
        return
      endif
      pg_class(ih) = c
    enddo
  enddo
  pg_n_class = c

  do i = 1, pg_n_elem
    if (pg_class(i) == 0) then
      if (.not. pg_detection_failed) then
        write(6,'(A)') 'symmetry: an element was left without a class'
        pg_detection_failed = .true.
      endif
      return
    endif
  enddo
END_PROVIDER

 BEGIN_PROVIDER [ integer, pg_class_size, (pg_n_class) ]
&BEGIN_PROVIDER [ integer, pg_class_rep, (pg_n_class) ]
  implicit none
  BEGIN_DOC
  ! Size of every class and its lowest element index.
  END_DOC
  integer :: i, c

  if (pg_n_class < 1) return
  pg_class_size(1:pg_n_class) = 0
  pg_class_rep(1:pg_n_class) = 0
  do i = 1, pg_n_elem
    c = pg_class(i)
    pg_class_size(c) = pg_class_size(c) + 1
    if (pg_class_rep(c) == 0) pg_class_rep(c) = i
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! Measured invariants used by the identification tree.
!-----------------------------------------------------------------------
 BEGIN_PROVIDER [ logical, pg_has_i ]
&BEGIN_PROVIDER [ integer, pg_n_sigma_h ]
&BEGIN_PROVIDER [ integer, pg_n_sigma_vtot ]
&BEGIN_PROVIDER [ integer, pg_n_c2_perp ]
  implicit none
  BEGIN_DOC
  ! Measured invariants used by the identification.
  !
  ! pg_n_sigma_vtot counts *all* vertical reflections, i.e. the planes whose
  ! normal is not along z. It equals n_max in both D_nh and D_nd, which is
  ! all the tree needs (n_sigma_h is tested first, so the two never
  ! conflict). There is no sigma_v/sigma_d split here: the split is a
  ! labelling device, it is decided from the plane normal in columns.irp.f, and
  ! it cannot be derived reliably from the geometry alone (the prototype's
  ! bisector test was inverted and still "worked" because both counts were
  ! only ever summed).
  END_DOC
  integer :: i, m, k, det
  double precision :: d
  double precision :: ident(3,3)

  ! "Parallel" and "perpendicular" to the principal axis are decided on the
  ! cosine with a threshold of 1/2 (45 degrees), not with a 1e-6
  ! tolerance. The elements are accepted when the geometry maps onto
  ! itself within :c:data:`pg_tol`, so their axes carry a tilt of the order of the
  ! geometric asymmetry: a reflection of a planar water whose H is moved by
  ! 1e-2 bohr has a normal at 0.26 degrees from the molecular plane, which
  ! a 1e-6 test rejects while the geometry still satisfies
  ! :c:data:`pg_tol` = 1e-4. In the groups of interest a line is either parallel or
  ! perpendicular to z, 90 degrees apart, so 45 degrees is the widest
  ! possible margin and no intermediate case can be misassigned.
  ! (:c:data:`pg_has_i` keeps the 1e-8 test: the inversion is tested as the exact
  ! matrix -1, never built from a candidate direction.)
  double precision, parameter :: pg_axis_cos = 0.5d0

  ident(1:3,1:3) = 0.d0
  ident(1,1) = 1.d0
  ident(2,2) = 1.d0
  ident(3,3) = 1.d0

  pg_has_i = .false.
  pg_n_sigma_h = 0
  pg_n_sigma_vtot = 0
  pg_n_c2_perp = 0

  do i = 1, pg_n_elem
    d = maxval(abs(pg_rot_mat(1:3,1:3,i) + ident(1:3,1:3)))
    if (d <= 1.d-8) pg_has_i = .true.

    call pg_decode_kind(pg_elem_kind(i), det, m, k)
    ! ``d`` is the cosine between the element axis (plane normal for a
    ! reflection) and the canonical z, i.e. row 3 of pg_rot.
    d = pg_elem_axis(1,i)*pg_rot(3,1) + pg_elem_axis(2,i)*pg_rot(3,2) &
      + pg_elem_axis(3,i)*pg_rot(3,3)
    if (det < 0 .and. m == 1) then
      if (abs(d) > pg_axis_cos) then
        pg_n_sigma_h = pg_n_sigma_h + 1
      else
        pg_n_sigma_vtot = pg_n_sigma_vtot + 1
      endif
    endif
    if (det > 0 .and. m == 2) then
      if (abs(d) <= pg_axis_cos) pg_n_c2_perp = pg_n_c2_perp + 1
    endif
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! Identification, restricted to the step-2 groups.
!-----------------------------------------------------------------------
 BEGIN_PROVIDER [ character*(16), pg_group_name ]
  implicit none
  BEGIN_DOC
  ! Name of the point group from the measured invariants.
  ! Covers all 204 finite groups of the data set; gives up instead of
  ! guessing on inconsistent invariants.
  !
  ! Verified against all 209 finite files with the closed-form construction
  ! of the reference oracle (refsym2.py + ident.py): the tree returns the
  ! file name for 204 of them. The other five (C1h, C1v, D1, D1d, D1h) are
  ! degenerate aliases of Cs, Cs, C2, C2h and C2v: they are the n = 1
  ! cases where the standard naming rules collapse two names onto one
  ! group, and they cannot be produced by the tree. Their files are never
  ! opened then; do not "fix" the tree to produce them.
  !
  ! The two facts that make the tree short are
  !  * n_axes_max counts the axes carrying a rotation of order n_max, so a
  !    D<n> group has n_axes_max = 1 for n >= 3 (its perpendicular axes
  !    have order 2) and cannot be mistaken for a cubic group;
  !  * in the n_max = 1 branch the orientation of the single mirror is not
  !    usable: z is the molecular plane normal for a planar molecule but an
  !    inertia axis otherwise, so C1/Ci/Cs are decided by :math:`|G|` and by
  !    "there is exactly one reflection".
  END_DOC
  character*(16) :: name
  integer :: i, order, ierr, n_size

  pg_group_name = ''
  if (pg_detection_failed) return

  ! Cubic and icosahedral groups first: more than one axis of order n_max.
  if (pg_n_max >= 3 .and. pg_n_axes_max >= 2) then
    name = ''
    if (pg_n_elem == 12) then
      name = 'T'
    else if (pg_n_elem == 24 .and. pg_n_max == 3) then
      if (pg_has_i) then
        name = 'Th'
      else
        name = 'Td'
      endif
    else if (pg_n_elem == 24 .and. pg_n_max == 4) then
      name = 'O'
    else if (pg_n_elem == 48 .and. pg_n_max == 4) then
      name = 'Oh'
    else if (pg_n_elem == 60 .and. pg_n_max == 5) then
      name = 'I'
    else if (pg_n_elem == 120 .and. pg_n_max == 5) then
      name = 'Ih'
    endif
    if (len_trim(name) == 0) then
      if (.not. pg_detection_failed) then
        write(6,'(A,A,I0,A,I0,A,I0)') 'symmetry: cubic/icosahedral '// &
             'candidate: ', '|G| = ', pg_n_elem, ', n_max = ', pg_n_max, &
             ', n_axes_max = ', pg_n_axes_max
        pg_detection_failed = .true.
      endif
      return
    endif
  else if (pg_n_max <= 1) then
    ! No proper rotation: C1 (nothing), Ci (S2, the inversion), Cs (S1, one
    ! plane). Nothing else can have :math:`|G|` <= 2.
    if (pg_n_elem == 1) then
      name = 'C1'
    else if (pg_has_i) then
      name = 'Ci'
    else if (pg_n_sigma_h + pg_n_sigma_vtot == 1) then
      name = 'Cs'
    else
      name = ''
    endif
  else if (pg_n_sigma_h >= 1) then
    ! A horizontal plane. C<n>h has no C2 perpendicular to z, D<n>h has n.
    if (pg_n_c2_perp > 0) then
      write(name,'(A,I0,A)') 'D', pg_n_max, 'h'
    else
      write(name,'(A,I0,A)') 'C', pg_n_max, 'h'
    endif
  else if (pg_n_sigma_vtot >= 1) then
    ! Vertical planes only. There is no sigma_v/sigma_d split here: the
    ! count is the same in D<n>d and D<n>h and the split is decided from
    ! the plane normal in columns.irp.f.
    if (pg_n_c2_perp > 0) then
      write(name,'(A,I0,A)') 'D', pg_n_max, 'd'
    else
      write(name,'(A,I0,A)') 'C', pg_n_max, 'v'
    endif
  else if (pg_n_elem == pg_n_max) then
    ! No plane, no inversion (:math:`|G|` = n_max excludes i, which would double the
    ! order): the cyclic group C<n>.
    write(name,'(A,I0)') 'C', pg_n_max
  else if (pg_n_elem == 2*pg_n_max) then
    ! Pure rotations of order 2*n_max: D<n> when there are perpendicular
    ! 2-fold axes, S<2n> (improper cyclic, proper part C<n>) when there are
    ! none.
    if (pg_n_c2_perp > 0) then
      write(name,'(A,I0)') 'D', pg_n_max
    else
      write(name,'(A,I0)') 'S', 2*pg_n_max
    endif
  else
    name = ''
  endif

  if (len_trim(name) == 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
           'symmetry: cannot name the group: |G| = ', pg_n_elem, &
           ', n_max = ', pg_n_max, ', n_axes_max = ', pg_n_axes_max, &
           ', has_i = ', merge(1, 0, pg_has_i), ', n_sigma_h = ', &
           pg_n_sigma_h, ', n_sigma_vtot = ', pg_n_sigma_vtot, &
           ', n_C2_perp = ', pg_n_c2_perp
      write(6,'(A)') 'symmetry: these invariants do not match any '// &
           'point group; this is a bug, please report it with the line above.'
      pg_detection_failed = .true.
    endif
    return
  endif

  ! Closed-form order cross-check. Never guess.
  call pg_expected_order(name, order, ierr)
  if (ierr /= 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A,A)') 'symmetry: cannot compute the order of ', trim(name)
      pg_detection_failed = .true.
    endif
    return
  endif
  if (order /= pg_n_elem) then
    if (.not. pg_detection_failed) then
      write(6,'(A,A,A,A,I0,A,I0)') 'symmetry: ', trim(name), ': closed-form ', &
           'order ', order, ' but |G| = ', pg_n_elem
      write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') '  invariants: n_max = ', &
           pg_n_max, ', n_axes_max = ', pg_n_axes_max, ', n_sigma_h = ', &
           pg_n_sigma_h, ', n_sigma_vtot = ', pg_n_sigma_vtot, &
           ', n_C2_perp = ', pg_n_c2_perp
      write(6,'(A)') 'symmetry: the accepted elements do not close to '// &
           'a group of this geometry within pg_tol. Elements are '// &
           'accepted one at a time and products are not re-matched, so a '// &
           'product of two elements accepted at '// &
           'pg_tol is only a symmetry to 2*pg_tol. This happens for a '// &
           'geometry that is symmetric to within a few pg_tol only; '// &
           'lower pg_tol (fewer elements accepted) or raise it above '// &
           'the residual.'
      pg_detection_failed = .true.
    endif
    return
  endif

  ! Classes partition the group. This is the only class property step 2 can
  ! assert without the character table: the Abelian identity n_class = :math:`|G|`
  ! is NOT one of them, and the table-dependent checks live in
  ! ``columns.irp.f`` (columns partition the group, orthogonality, sum of
  ! d^2/mult).
  n_size = 0
  do i = 1, pg_n_class
    if (pg_class_size(i) < 1) then
      if (.not. pg_detection_failed) then
        write(6,'(A,A,A,I0,A)') 'symmetry: ', trim(name), ': class ', i, &
             ' is empty'
        pg_detection_failed = .true.
      endif
      return
    endif
    n_size = n_size + pg_class_size(i)
  enddo
  if (n_size /= pg_n_elem) then
    if (.not. pg_detection_failed) then
      write(6,'(A,A,A,I0,A,I0)') 'symmetry: ', trim(name), &
           ': classes cover ', n_size, ' elements, |G| = ', pg_n_elem
      pg_detection_failed = .true.
    endif
    return
  endif

  ! The measured group of a linear molecule is C31v or D31h (see
  ! pg_build_linear), but the group it stands for is Cinfv or Dinfh. The
  ! translation comes after the closed-form order check on purpose:
  ! pg_expected_order has no "inf" arithmetic, and the symbols
  ! C31v / D31h are exactly what the finite model is. pg_has_i must not be
  ! used here: this group contains no inversion for Dinfh, and it must not.
  if (trim(name) == 'C31v') then
    name = 'Cinfv'
  else if (trim(name) == 'D31h') then
    name = 'Dinfh'
  endif

  pg_group_name = name
END_PROVIDER
