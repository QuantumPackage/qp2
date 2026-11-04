! Canonical frame.
!
! Conventions fixed here (all geometric vectors are stored in the INPUT frame,
! the frame of nucl_coord; pg_rot maps input coordinates to the canonical
! frame):
!
!   pg_center       origin of the canonical frame, the centre of mass
!   pg_rot          rows are the canonical x, y and z unit vectors
!                   standard coordinates are pg_rot . (r - pg_center)
!   pg_is_planar    all nuclei within pg_tol of a plane
!   pg_plane_normal unit normal of that plane, INPUT frame, valid when
!                   pg_is_planar; oriented deterministically
!
! pg_plane_normal is given in the INPUT frame. Every other vector of the
! module (candidate directions, rotation matrices) lives in the input
! frame and pg_rot is the only object that maps between the two. A normal
! expressed in the standard frame would be the constant (0,0,1), which
! carries no information.
!
! Trap: the inertia tensor uses element_mass(int(nucl_charge(i))). There is no
! nucl_mass provider in QP2 (checked): that is what nuclei/inertia.irp.f and
! center_of_mass use.
!
! Axis convention (Cotton/Mulliken, decided by the canonical frame and used
! by every symbol predicate of the column matching): for a planar molecule
! with z in the plane,
! x is the plane normal and the molecular plane is yz. Hence sv = plane xz
! (normal y) and sd = plane yz (normal x); getting these two backwards swaps
! B2/B3 in every Abelian group.

 BEGIN_PROVIDER [ double precision, pg_center, (3) ]
  implicit none
  BEGIN_DOC
  ! Centre of mass: origin of the canonical frame.
  END_DOC
  pg_center(1:3) = center_of_mass(1:3)
END_PROVIDER

!-----------------------------------------------------------------------
! Principal moments of inertia, ascending, and their eigenvectors
! (columns of pg_inertia_evec). Recomputed with dsyev on the fly; do not
! reuse inertia_tensor_eigenvectors.
!-----------------------------------------------------------------------
 BEGIN_PROVIDER [ double precision, pg_inertia_eval, (3) ]
&BEGIN_PROVIDER [ double precision, pg_inertia_evec, (3,3) ]
  implicit none
  BEGIN_DOC
  ! Principal moments of inertia (ascending) and eigenvectors (columns).
  ! Eigenvector signs are fixed deterministically (largest component
  ! positive) so that repeated runs are byte-identical.
  END_DOC
  integer :: i
  double precision :: tensor(3,3), v(3), m

  ! Defaults first: IRPF90 sets the built flag after the body returns, so
  ! an early return caches whatever was written. A failed dsyev must hand
  ! its callers a zero block rather than an uninitialized one, and the zero
  ! block is permanent until a TOUCH.
  pg_inertia_eval(1:3) = 0.d0
  pg_inertia_evec(1:3,1:3) = 0.d0
  if (pg_detection_failed) return

  tensor(1:3,1:3) = 0.d0
  do i = 1, nucl_num
    m = element_mass(int(nucl_charge(i)))
    v(1:3) = nucl_coord(i,1:3) - center_of_mass(1:3)
    tensor(1,1) = tensor(1,1) + m*(v(2)**2 + v(3)**2)
    tensor(2,2) = tensor(2,2) + m*(v(1)**2 + v(3)**2)
    tensor(3,3) = tensor(3,3) + m*(v(1)**2 + v(2)**2)
    tensor(1,2) = tensor(1,2) - m*v(1)*v(2)
    tensor(1,3) = tensor(1,3) - m*v(1)*v(3)
    tensor(2,3) = tensor(2,3) - m*v(2)*v(3)
  enddo
  tensor(2,1) = tensor(1,2)
  tensor(3,1) = tensor(1,3)
  tensor(3,2) = tensor(2,3)

  call pg_diag3(tensor, pg_inertia_eval)
  if (pg_detection_failed) return
  pg_inertia_evec(1:3,1:3) = tensor(1:3,1:3)
  do i = 1, 3
    call pg_fix_sign(pg_inertia_evec(1:3,i))
  enddo

  ! Relative test, before any enumeration. A maxval test can never fire:
  ! maxval >= sum/3 > 1e-6*sum whenever sum > 0. Only the "all nuclei at
  ! the centre" case is rejected here; a linear molecule is legal and is
  ! detected by pg_is_linear, which
  ! uses the same minval test on a local copy.
  if (sum(pg_inertia_eval) <= 0.d0) then
    if (.not. pg_detection_failed) then
      write(6,'(A)') 'symmetry: no nucleus, or all nuclei at the centre: '// &
           'unsupported'
      pg_detection_failed = .true.
    endif
    return
  endif
END_PROVIDER

 BEGIN_PROVIDER [ logical, pg_is_linear ]
  implicit none
  BEGIN_DOC
  ! .True. for a collinear molecule. Two tests, both needed: the moment
  ! about the molecular axis is the *smallest* inertia eigenvalue and it is
  ! exactly zero for a collinear geometry (a maxval test can never fire),
  ! and every nucleus must sit within half pg_tol of the line through
  ! pg_center along pg_inertia_evec(:,1).
  !
  ! Half pg_tol, not one pg_tol: the group of the molecule is built from
  ! matrices that are exact in the canonical frame, so the only error is the
  ! distance d of a nucleus from the axis, and a sigma_v moves it by 2d. A
  ! nucleus at d = pg_tol would make every sigma_v fail pg_maps_to_self, and
  ! the run would end on "element does not realise its permutation" instead
  ! of being treated as linear.
  END_DOC
  integer :: i
  double precision :: v(3), d, ax(3)

  pg_is_linear = .false.
  if (pg_detection_failed) return
  if (sum(pg_inertia_eval) <= 0.d0) return
  if (minval(pg_inertia_eval) >= 1.d-6*sum(pg_inertia_eval)) return

  ax(1:3) = pg_inertia_evec(1:3,1)
  do i = 1, nucl_num
    v(1:3) = nucl_coord(i,1:3) - pg_center(1:3)
    d = v(1)*ax(1) + v(2)*ax(2) + v(3)*ax(3)
    d = (v(1)-d*ax(1))**2 + (v(2)-d*ax(2))**2 + (v(3)-d*ax(3))**2
    if (sqrt(d) > 0.5d0*pg_tol) return
  enddo
  pg_is_linear = .true.
 END_PROVIDER

 BEGIN_PROVIDER [ logical, pg_is_planar ]
  implicit none
  BEGIN_DOC
  ! .True. when every nucleus lies within pg_tol of a plane. The candidate
  ! normal is the inertia eigenvector with the largest eigenvalue (for a
  ! planar molecule that eigenvalue belongs to the normal).
  END_DOC
  integer :: i
  double precision :: v(3), d

  pg_is_planar = .true.
  do i = 1, nucl_num
    v(1:3) = nucl_coord(i,1:3) - pg_center(1:3)
    d = abs(v(1)*pg_inertia_evec(1,3) + v(2)*pg_inertia_evec(2,3) &
          + v(3)*pg_inertia_evec(3,3))
    if (d > pg_tol) then
      pg_is_planar = .false.
      exit
    endif
  enddo
END_PROVIDER

 BEGIN_PROVIDER [ double precision, pg_plane_normal, (3) ]
  implicit none
  BEGIN_DOC
  ! Unit normal of the molecular plane, in the INPUT frame.
  ! Valid when pg_is_planar is true; defined for every molecule.
  END_DOC
  pg_plane_normal(1:3) = pg_inertia_evec(1:3,3)
  call pg_fix_sign(pg_plane_normal)
END_PROVIDER

!-----------------------------------------------------------------------
! Canonical frame. Rows of pg_rot are x, y, z.
!-----------------------------------------------------------------------
 BEGIN_PROVIDER [ double precision, pg_rot, (3,3) ]
  implicit none
  BEGIN_DOC
  ! Rows are the canonical x, y, z unit vectors, in the INPUT frame.
  ! Standard coordinates of a point r are pg_rot . (r - pg_center).
  !
  ! Orientation follows the steps below, including the tie-breaks. Cubic and
  ! icosahedral labels are invariant under relabelling of their equivalent
  ! axes, so no attempt is made to impose a "standard cube orientation".
  END_DOC
  integer :: znum(nucl_num), ierr, n_max, n_axes
  integer :: i
  double precision :: zaxis(3), xaxis(3), yaxis(3)
  double precision :: d

  do i = 1, nucl_num
    znum(i) = nint(nucl_charge(i))
  enddo

  ! Step 3-4: z.
  call pg_proper_scan(nucl_num, nucl_coord, znum, pg_center, &
       pg_inertia_eval, pg_inertia_evec, pg_tol, 2000, &
       n_max, n_axes, zaxis, pg_is_linear, ierr)
  if (ierr /= 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A,I0,A)') 'symmetry: more than ', 2000, &
           ' candidate directions: increase pg_cand_max'
      pg_detection_failed = .true.
    endif
    return
  endif

  if (n_max >= 2) then
    ! Highest-order proper axis; already sign-oriented by pg_proper_scan.
  else if (pg_is_planar) then
    zaxis(1:3) = pg_plane_normal(1:3)
    call pg_orient_axis(nucl_num, nucl_coord, pg_center, pg_tol, zaxis)
  else
    zaxis(1:3) = pg_inertia_evec(1:3,3)
    call pg_orient_axis(nucl_num, nucl_coord, pg_center, pg_tol, zaxis)
  endif

  ! Step 5: x. Branch (a) is the Mulliken rule (x perpendicular to the
  ! molecular plane) and it must be tested first: it is the only
  ! deterministic choice for a planar D2h molecule.
  call pg_select_x(nucl_num, nucl_coord, znum, pg_center, &
       pg_inertia_evec, pg_tol, 2000, zaxis, pg_is_planar, &
       pg_plane_normal, xaxis, ierr)
  if (ierr /= 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A,I0,A)') 'symmetry: more than ', 2000, &
           ' candidate directions: increase pg_cand_max'
      pg_detection_failed = .true.
    endif
    return
  endif

  ! Orthonormalise x against z, then y = z cross x.
  d = xaxis(1)*zaxis(1) + xaxis(2)*zaxis(2) + xaxis(3)*zaxis(3)
  xaxis(1:3) = xaxis(1:3) - d*zaxis(1:3)
  d = sqrt(xaxis(1)**2 + xaxis(2)**2 + xaxis(3)**2)
  if (d < 1.d-10) then
    if (.not. pg_detection_failed) then
      write(6,'(A)') 'symmetry: cannot build a canonical x axis'
      pg_detection_failed = .true.
    endif
    return
  endif
  xaxis(1:3) = xaxis(1:3)/d
  call pg_cross(zaxis, xaxis, yaxis)
  pg_rot(1,1:3) = xaxis(1:3)
  pg_rot(2,1:3) = yaxis(1:3)
  pg_rot(3,1:3) = zaxis(1:3)
END_PROVIDER

!=======================================================================
! Helper routines. All array arguments are explicit-shape: IRPF90 does
! not generate interfaces, so assumed-shape dummies silently break
! (verified in a standalone IRPF90 harness).
!=======================================================================

!-----------------------------------------------------------------------
! Diagonalise a real symmetric 3x3 matrix. On exit a holds the
! eigenvectors as columns and e the eigenvalues in ascending order.
!-----------------------------------------------------------------------
subroutine pg_diag3(a, e)
  implicit none
  BEGIN_DOC
  ! dsyev wrapper for a 3x3 symmetric matrix (destroys ``a``).
  END_DOC
  double precision, intent(inout) :: a(3,3)
  double precision, intent(out) :: e(3)
  integer :: info, lwork
  double precision :: work(1)
  double precision, allocatable :: w(:)

  lwork = -1
  call dsyev('V', 'U', 3, a, 3, e, work, lwork, info)
  if (info /= 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A)') 'symmetry: dsyev workspace query failed'
      pg_detection_failed = .true.
    endif
    return
  endif
  lwork = int(work(1))
  allocate(w(lwork))
  call dsyev('V', 'U', 3, a, 3, e, w, lwork, info)
  deallocate(w)
  if (info /= 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A)') 'symmetry: dsyev failed in the inertia tensor'
      pg_detection_failed = .true.
    endif
    return
  endif
end subroutine pg_diag3

!-----------------------------------------------------------------------
! Flip v so that its largest-magnitude component (lowest index on ties)
! is positive. Purely deterministic sign convention.
!-----------------------------------------------------------------------
subroutine pg_fix_sign(v)
  implicit none
  BEGIN_DOC
  ! Deterministic orientation of a vector.
  END_DOC
  double precision, intent(inout) :: v(3)
  integer :: i, imax
  double precision :: amax

  imax = 1
  amax = abs(v(1))
  do i = 2, 3
    if (abs(v(i)) > amax) then
      amax = abs(v(i))
      imax = i
    endif
  enddo
  if (v(imax) < 0.d0) v(1:3) = -v(1:3)
end subroutine pg_fix_sign

!-----------------------------------------------------------------------
! c = a cross b
!-----------------------------------------------------------------------
subroutine pg_cross(a, b, c)
  implicit none
  BEGIN_DOC
  ! Cross product.
  END_DOC
  double precision, intent(in) :: a(3), b(3)
  double precision, intent(out) :: c(3)
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
end subroutine pg_cross

!-----------------------------------------------------------------------
! Rotation matrix of angle theta about axis (Rodrigues).
!-----------------------------------------------------------------------
subroutine pg_rotation(axis, theta, rot)
  implicit none
  BEGIN_DOC
  ! Proper rotation matrix, angle in radians.
  END_DOC
  double precision, intent(in) :: axis(3), theta
  double precision, intent(out) :: rot(3,3)
  double precision :: u(3), c, s, t, d

  d = sqrt(axis(1)**2 + axis(2)**2 + axis(3)**2)
  u(1:3) = axis(1:3)/d
  c = cos(theta)
  s = sin(theta)
  t = 1.d0 - c
  rot(1,1) = t*u(1)*u(1) + c
  rot(1,2) = t*u(1)*u(2) - s*u(3)
  rot(1,3) = t*u(1)*u(3) + s*u(2)
  rot(2,1) = t*u(1)*u(2) + s*u(3)
  rot(2,2) = t*u(2)*u(2) + c
  rot(2,3) = t*u(2)*u(3) - s*u(1)
  rot(3,1) = t*u(1)*u(3) - s*u(2)
  rot(3,2) = t*u(2)*u(3) + s*u(1)
  rot(3,3) = t*u(3)*u(3) + c
end subroutine pg_rotation

!-----------------------------------------------------------------------
! Reflection through the plane of normal n.
!-----------------------------------------------------------------------
subroutine pg_reflection(n, rot)
  implicit none
  BEGIN_DOC
  ! Householder reflection R = 1 - 2 n n^T.
  END_DOC
  double precision, intent(in) :: n(3)
  double precision, intent(out) :: rot(3,3)
  double precision :: u(3), d
  integer :: i, j

  d = sqrt(n(1)**2 + n(2)**2 + n(3)**2)
  u(1:3) = n(1:3)/d
  do j = 1, 3
    do i = 1, 3
      rot(i,j) = -2.d0*u(i)*u(j)
    enddo
    rot(j,j) = rot(j,j) + 1.d0
  enddo
end subroutine pg_reflection

!-----------------------------------------------------------------------
! Orient an axis so that axis . (r_i - ctr) > 0 for the lowest-index atom
! with a non-zero projection; fall back to pg_fix_sign.
!-----------------------------------------------------------------------
subroutine pg_orient_axis(nat, xyz, ctr, tol, axis)
  implicit none
  BEGIN_DOC
  ! Deterministic orientation of an axis.
  END_DOC
  integer, intent(in) :: nat
  double precision, intent(in) :: xyz(nat,3), ctr(3), tol
  double precision, intent(inout) :: axis(3)
  integer :: i
  double precision :: v(3), d

  do i = 1, nat
    v(1:3) = xyz(i,1:3) - ctr(1:3)
    d = v(1)*axis(1) + v(2)*axis(2) + v(3)*axis(3)
    if (abs(d) > tol) then
      if (d < 0.d0) axis(1:3) = -axis(1:3)
      return
    endif
  enddo
  call pg_fix_sign(axis)
end subroutine pg_orient_axis

!-----------------------------------------------------------------------
! Candidate directions, standard frame, deduplicated up to sign
! Candidate points: r_i - centre, (r_i + r_j)/2 - centre for pairs of
! equal Z, the three inertia eigenvectors, and every pairwise cross
! product of those.
!-----------------------------------------------------------------------
subroutine pg_build_candidates(nat, xyz, znum, ctr, evecs, ncandmax, &
     cand, ncand, ierr)
  implicit none
  BEGIN_DOC
  ! Build the list of candidate symmetry axes / plane normals.
  ! ``ierr`` = 1 when the list overflows ``ncandmax``.
  END_DOC
  integer, intent(in) :: nat, ncandmax
  double precision, intent(in) :: xyz(nat,3), ctr(3), evecs(3,3)
  integer, intent(in) :: znum(nat)
  double precision, intent(out) :: cand(3,ncandmax)
  integer, intent(out) :: ncand, ierr
  double precision :: base(3,ncandmax), v(3), w(3), c(3)
  double precision :: dn
  integer :: i, j, nbase

  ncand = 0
  nbase = 0
  ierr = 0

  ! Atomic positions relative to the centre.
  do i = 1, nat
    v(1:3) = xyz(i,1:3) - ctr(1:3)
    if (sqrt(v(1)**2 + v(2)**2 + v(3)**2) <= 1.d-8) cycle
    call pg_add_direction(v, ncandmax, base, nbase, ierr)
    if (ierr /= 0) return
  enddo

  ! Midpoints of pairs of atoms of the same element.
  do i = 1, nat - 1
    do j = i + 1, nat
      if (znum(j) /= znum(i)) cycle
      w(1:3) = 0.5d0*(xyz(i,1:3) + xyz(j,1:3)) - ctr(1:3)
      if (sqrt(w(1)**2 + w(2)**2 + w(3)**2) <= 1.d-8) cycle
      call pg_add_direction(w, ncandmax, base, nbase, ierr)
      if (ierr /= 0) return
    enddo
  enddo

  ! Inertia axes.
  do i = 1, 3
    call pg_add_direction(evecs(1:3,i), ncandmax, base, nbase, ierr)
    if (ierr /= 0) return
  enddo

  ! Pairwise cross products of the base directions.
  do i = 1, nbase - 1
    do j = i + 1, nbase
      call pg_cross(base(1,i), base(1,j), c)
      dn = sqrt(c(1)**2 + c(2)**2 + c(3)**2)
      if (dn <= 1.d-8) cycle
      call pg_add_direction(c, ncandmax, cand, ncand, ierr)
      if (ierr /= 0) return
    enddo
  enddo

  ! The base directions themselves are candidates too.
  do i = 1, nbase
    call pg_add_direction(base(1,i), ncandmax, cand, ncand, ierr)
    if (ierr /= 0) return
  enddo

  ! Cheap consistency check: ``ncand`` > 0 for any molecule with >= 2 atoms.
  if (ncand == 0) ierr = 2

end subroutine pg_build_candidates

!-----------------------------------------------------------------------
! Normalise v and append it to list if no stored direction is parallel or
! antiparallel to it (tolerance 1e-8 on ``cos``). Silently skips ``v`` ~ 0.
!-----------------------------------------------------------------------
subroutine pg_add_direction(v, nmax, list, n, ierr)
  implicit none
  BEGIN_DOC
  ! Append a direction to a deduplicated (up to sign) list.
  END_DOC
  double precision, intent(in) :: v(3)
  integer, intent(in) :: nmax
  double precision, intent(inout) :: list(3,nmax)
  integer, intent(inout) :: n
  integer, intent(inout) :: ierr
  double precision :: u(3), d, dn
  integer :: k

  dn = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
  if (dn <= 1.d-8) return
  u(1:3) = v(1:3)/dn
  do k = 1, n
    d = abs(u(1)*list(1,k) + u(2)*list(2,k) + u(3)*list(3,k))
    if (d > 1.d0 - 1.d-8) return
  enddo
  if (n >= nmax) then
    ierr = 1
    return
  endif
  n = n + 1
  list(1:3,n) = u(1:3)
end subroutine pg_add_direction

!-----------------------------------------------------------------------
! Does rot map the molecule onto itself? perm(i) = image of atom i.
! ok = .false. when an atom has no partner within tol, or when the
! nearest-partner map is not injective.
!-----------------------------------------------------------------------
subroutine pg_maps_to_self(nat, xyz, znum, ctr, rot, tol, perm, ok)
  implicit none
  BEGIN_DOC
  ! Test rot r_i == r_{perm(i)} for every atom.
  END_DOC
  integer, intent(in) :: nat
  double precision, intent(in) :: xyz(nat,3), ctr(3), rot(3,3), tol
  integer, intent(in) :: znum(nat)
  integer, intent(out) :: perm(nat)
  logical, intent(out) :: ok
  integer :: i, j, jmin, used(nat)
  double precision :: v(3), w(3), d, dmin

  ok = .false.
  used(1:nat) = 0
  do i = 1, nat
    v(1:3) = xyz(i,1:3) - ctr(1:3)
    w(1) = rot(1,1)*v(1) + rot(1,2)*v(2) + rot(1,3)*v(3) + ctr(1)
    w(2) = rot(2,1)*v(1) + rot(2,2)*v(2) + rot(2,3)*v(3) + ctr(2)
    w(3) = rot(3,1)*v(1) + rot(3,2)*v(2) + rot(3,3)*v(3) + ctr(3)
    jmin = 0
    dmin = huge(1.d0)
    do j = 1, nat
      if (znum(j) /= znum(i)) cycle
      d = (w(1)-xyz(j,1))**2 + (w(2)-xyz(j,2))**2 + (w(3)-xyz(j,3))**2
      if (d < dmin) then
        dmin = d
        jmin = j
      endif
    enddo
    if (jmin == 0) return
    if (sqrt(dmin) > tol) return
    if (used(jmin) /= 0) return
    used(jmin) = 1
    perm(i) = jmin
  enddo
  ok = .true.
end subroutine pg_maps_to_self

!-----------------------------------------------------------------------
! Highest proper-rotation order, number of distinct axes carrying it,
! and the selected principal axis.
!-----------------------------------------------------------------------
subroutine pg_proper_scan(nat, xyz, znum, ctr, evals, evecs, tol, &
     ncandmax, n_max, n_axes_max, axis, linear, ierr)
  implicit none
  BEGIN_DOC
  ! Enumerate proper rotations of order 2..32 among the candidate
  ! directions, keep the highest order and select the principal axis.
  ! 32, not 12: C32 occurs in D32h/D32d, the largest groups of the data set.
  ! A linear molecule is capped at 31: its improper elements are S_n^k and
  ! the shipped file names the 15 improper columns S31^k, k = 1..15. A
  ! composite order would not keep that scheme -- some powers are not S_n^k
  ! at all but reduce to a lower order (S25^5 = S5), so the columns would
  ! mix orders and "one column per power of S_n" would no longer hold.
  ! 31 is the largest odd order <= 32, and 4*31 = 124 still fits the
  ! 160-element arrays; the tables shipped are C31v and D31h.
  END_DOC
  integer, intent(in) :: nat, ncandmax
  double precision, intent(in) :: xyz(nat,3), ctr(3), evals(3), evecs(3,3), tol
  integer, intent(in) :: znum(nat)
  logical, intent(in) :: linear
  integer, intent(out) :: n_max, n_axes_max, ierr
  double precision, intent(out) :: axis(3)
  double precision, parameter :: pi = 3.1415926535897932d0
  double precision :: cand(3,ncandmax), rot(3,3), v(3)
  integer :: order(ncandmax), perm(nat), ncap
  integer :: ncand, i, j, k, n, kbest, icount, ibest
  double precision :: d, dmetric, dbest, gap, gapbest
  logical :: ok

  axis(1:3) = 0.d0
  n_max = 1
  n_axes_max = 0
  ncap = 32
  if (linear) ncap = 31
  call pg_build_candidates(nat, xyz, znum, ctr, evecs, ncandmax, &
       cand, ncand, ierr)
  if (ierr /= 0) return

  order(1:ncand) = 1
  do k = 1, ncand
    do n = ncap, 2, -1
      call pg_rotation(cand(1,k), 2.d0*pi/dble(n), rot)
      call pg_maps_to_self(nat, xyz, znum, ctr, rot, tol, perm, ok)
      if (ok) then
        order(k) = n
        exit
      endif
    enddo
    if (order(k) > n_max) n_max = order(k)
  enddo

  if (n_max == 1) return

  do k = 1, ncand
    if (order(k) == n_max) n_axes_max = n_axes_max + 1
  enddo

  ! Index of the inertia eigenvalue that is the most distinct from the
  ! other two: used by the axis tie-break below.
  kbest = 1
  gapbest = -1.d0
  do k = 1, 3
    gap = huge(1.d0)
    do j = 1, 3
      if (j /= k) gap = min(gap, abs(evals(k)-evals(j)))
    enddo
    if (gap > gapbest) then
      gapbest = gap
      kbest = k
    endif
  enddo

  ! Selection: most nuclei on the axis, then largest ``axis.W(:,kbest)``,
  ! then lowest candidate index (strict comparisons keep the first).
  ibest = 0
  icount = -1
  dbest = -1.d0
  do k = 1, ncand
    if (order(k) /= n_max) cycle
    j = 0
    do i = 1, nat
      v(1:3) = xyz(i,1:3) - ctr(1:3)
      d = v(1)*cand(1,k) + v(2)*cand(2,k) + v(3)*cand(3,k)
      d = (v(1)-d*cand(1,k))**2 + (v(2)-d*cand(2,k))**2 &
        + (v(3)-d*cand(3,k))**2
      if (sqrt(d) <= tol) j = j + 1
    enddo
    dmetric = abs(cand(1,k)*evecs(1,kbest) + cand(2,k)*evecs(2,kbest) &
                + cand(3,k)*evecs(3,kbest))
    if (j > icount .or. (j == icount .and. dmetric > dbest)) then
      ibest = k
      icount = j
      dbest = dmetric
    endif
  enddo

  if (ibest == 0) then
    ierr = 3
    return
  endif
  axis(1:3) = cand(1:3,ibest)
  call pg_orient_axis(nat, xyz, ctr, tol, axis)
end subroutine pg_proper_scan

!-----------------------------------------------------------------------
! Canonical x axis.
!-----------------------------------------------------------------------
subroutine pg_select_x(nat, xyz, znum, ctr, evecs, tol, ncandmax, &
     zaxis, planar, pnormal, xaxis, ierr)
  implicit none
  BEGIN_DOC
  ! Choose the canonical x from the branches below.
  END_DOC
  integer, intent(in) :: nat, ncandmax
  double precision, intent(in) :: xyz(nat,3), ctr(3), evecs(3,3), tol
  integer, intent(in) :: znum(nat)
  double precision, intent(in) :: zaxis(3), pnormal(3)
  logical, intent(in) :: planar
  double precision, intent(out) :: xaxis(3)
  integer, intent(out) :: ierr
  double precision, parameter :: pi = 3.1415926535897932d0
  double precision :: cand(3,ncandmax), rot(3,3), v(3)
  integer :: perm(nat)
  integer :: ncand, i, j, k, best, icount, ibest
  double precision :: d, dmetric, dbest
  logical :: ok

  ierr = 0
  xaxis(1:3) = 0.d0

  ! (a) Planar and z lies in the plane: x is the plane normal (Mulliken).
  if (planar) then
    d = zaxis(1)*pnormal(1) + zaxis(2)*pnormal(2) + zaxis(3)*pnormal(3)
    if (abs(d) < 0.5d0) then
      xaxis(1:3) = pnormal(1:3)
      return
    endif
  endif

  ! (b) A 2-fold proper axis perpendicular to z.
  call pg_build_candidates(nat, xyz, znum, ctr, evecs, ncandmax, &
       cand, ncand, ierr)
  if (ierr /= 0) return
  ibest = 0
  icount = -1
  dbest = -1.d0
  do k = 1, ncand
    d = cand(1,k)*zaxis(1) + cand(2,k)*zaxis(2) + cand(3,k)*zaxis(3)
    if (abs(d) > 1.d-6) cycle
    call pg_rotation(cand(1,k), pi, rot)
    call pg_maps_to_self(nat, xyz, znum, ctr, rot, tol, perm, ok)
    if (.not. ok) cycle
    j = 0
    do i = 1, nat
      v(1:3) = xyz(i,1:3) - ctr(1:3)
      d = v(1)*cand(1,k) + v(2)*cand(2,k) + v(3)*cand(3,k)
      d = (v(1)-d*cand(1,k))**2 + (v(2)-d*cand(2,k))**2 &
        + (v(3)-d*cand(3,k))**2
      if (sqrt(d) <= tol) j = j + 1
    enddo
    dmetric = abs(cand(1,k)*evecs(1,3) + cand(2,k)*evecs(2,3) &
                + cand(3,k)*evecs(3,3))
    if (j > icount .or. (j == icount .and. dmetric > dbest)) then
      ibest = k
      icount = j
      dbest = dmetric
    endif
  enddo
  if (ibest /= 0) then
    xaxis(1:3) = cand(1:3,ibest)
    return
  endif

  ! (c)-(d) Projection of the lowest-index atom not on the z axis.
  do i = 1, nat
    v(1:3) = xyz(i,1:3) - ctr(1:3)
    d = v(1)*zaxis(1) + v(2)*zaxis(2) + v(3)*zaxis(3)
    v(1:3) = v(1:3) - d*zaxis(1:3)
    if (sqrt(v(1)**2 + v(2)**2 + v(3)**2) > tol) then
      xaxis(1:3) = v(1:3)
      return
    endif
  enddo

  ! Fallback: deterministic vector perpendicular to z.
  best = 1
  do i = 2, 3
    if (abs(zaxis(i)) < abs(zaxis(best))) best = i
  enddo
  xaxis(1:3) = 0.d0
  xaxis(best) = 1.d0
  d = xaxis(1)*zaxis(1) + xaxis(2)*zaxis(2) + xaxis(3)*zaxis(3)
  xaxis(1:3) = xaxis(1:3) - d*zaxis(1:3)
  d = sqrt(xaxis(1)**2 + xaxis(2)**2 + xaxis(3)**2)
  if (d > 1.d-12) xaxis(1:3) = xaxis(1:3)/d
end subroutine pg_select_x
