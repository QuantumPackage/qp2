program symmetry
  implicit none
  BEGIN_DOC
  ! Step-2 driver of the symmetry module.
  !
  ! Reads the nuclei from the EZFIO directory given on the command line
  ! (through the standard :c:data:`ezfio_filename` provider), builds the canonical
  ! frame, enumerates the point group, names it, prints the classes, and
  ! prints the character table shipped in ``QP_ROOT/data/Symmetry``.
  !
  ! The EZFIO directory is the first command-line argument, as for every
  ! other QP2 program.
  !
  ! An optional environment variable ``PG_WRITE_MO_SYMMETRY=1`` writes the
  ! labels into ``mo_basis/mo_symmetry``. Nothing is
  ! written when a gate fails, so a written field is always a checked one.
  ! An environment variable is used because :c:data:`ezfio_filename` stops unless the
  ! command line has exactly one argument.
  END_DOC
  character*(256) :: arg
  logical :: write_mo

  write_mo = .false.
  call getenv('PG_WRITE_MO_SYMMETRY', arg)
  write_mo = (trim(arg) == '1')
  call pg_report(write_mo)
end program symmetry

!-----------------------------------------------------------------------
! Report. Kept in the same file as the program: IRPF90 requires one
! program per file, and the report is the whole program.
!-----------------------------------------------------------------------
subroutine pg_report(write_mo)
  implicit none
  BEGIN_DOC
  ! Prints the canonical frame, the group, its classes and the table.
  ! ``write_mo`` forwards the "-w" request of the driver to :c:func:`pg_report_mo`.
  END_DOC
  logical, intent(in) :: write_mo
  character*(1024) :: path
  character*(4096) :: plist
  double precision :: v(3), w(3), d, delta_mo
  integer :: ierr, n_mem
  integer :: i, c, a, g, det, m, k

  delta_mo = 0.d0

  write(6,'(A)') '============================================================'
  write(6,'(A)') ' symmetry driver: canonical frame, point group, character table'
  write(6,'(A)') '============================================================'
  write(6,'(A,A)') 'EZFIO directory : ', trim(ezfio_filename)
  write(6,'(A,I0)') 'nuclei          : ', nucl_num
  write(6,'(A,1X,ES12.4)') 'pg_tol (bohr)   : ', pg_tol
  write(6,'(A,1X,ES12.4)') 'pg_thresh       : ', pg_thresh
  write(6,'(A,3(1X,F16.8))') 'centre of mass (bohr):', pg_center(1:3)

  !---------------------------------------------------------------------
  ! Canonical frame.
  !---------------------------------------------------------------------
  write(6,'(A)') ''
  write(6,'(A)') 'canonical frame (rows of pg_rot, in the input frame):'
  write(6,'(A,3(1X,F16.10))') '  x =', pg_rot(1,1:3)
  write(6,'(A,3(1X,F16.10))') '  y =', pg_rot(2,1:3)
  write(6,'(A,3(1X,F16.10))') '  z =', pg_rot(3,1:3)
  write(6,'(A,L1)') 'is planar       : ', pg_is_planar
  write(6,'(A,3(1X,F16.10))') 'plane normal    : ', pg_plane_normal(1:3)

  write(6,'(A)') ''
  write(6,'(A)') 'nuclei in the canonical frame:'
  write(6,'(A)') '    i     Z              x                y                z'
  do i = 1, nucl_num
    v(1:3) = nucl_coord(i,1:3) - pg_center(1:3)
    w(1:3) = matmul(pg_rot(1:3,1:3), v(1:3))
    write(6,'(I5,1X,F6.2,3(1X,F16.10))') i, nucl_charge(i), w(1:3)
  enddo

  !---------------------------------------------------------------------
  ! Group and classes.
  !---------------------------------------------------------------------
  write(6,'(A)') ''
  write(6,'(A,A)') 'point group     : ', trim(pg_group_name)
  write(6,'(A,I0)') 'order |G|       : ', pg_n_elem
  write(6,'(A,I0)') 'n_max           : ', pg_n_max
  write(6,'(A,I0)') 'n_axes_max      : ', pg_n_axes_max
  write(6,'(A,I0)') 'classes         : ', pg_n_class
  ! Measured symmetry breaking of the stored geometry.
  write(6,'(A,ES12.4,A,ES12.4,A)') 'geometry residual : ', pg_geom_residual, &
       '  bohr  (max |R r - r_perm|; pg_tol = ', pg_tol, ')'

  write(6,'(A)') ''
  write(6,'(A)') 'classes:'
  write(6,'(A)') ' class  size   rep  det   m   k      theta'
  do c = 1, pg_n_class
    i = pg_class_rep(c)
    call pg_decode_kind(pg_elem_kind(i), det, m, k)
    write(6,'(I6,I6,I6,I5,I4,I4,1X,F12.6)') c, pg_class_size(c), i, &
         pg_elem_det(i), m, k, pg_elem_angle(i)
  enddo

  write(6,'(A)') ''
  write(6,'(A)') 'elements:'
  write(6,'(A)') '   ig  class  det   m   k      theta   axis (input frame)' &
       //'          permutation'
  do i = 1, pg_n_elem
    call pg_decode_kind(pg_elem_kind(i), det, m, k)
    call pg_perm_string(pg_perm(1,i), nucl_num, plist)
    write(6,'(I5,I6,I5,I4,I4,1X,F12.6,3(1X,F10.6),1X,A)') i, &
         pg_class(i), pg_elem_det(i), m, k, pg_elem_angle(i), &
         pg_elem_axis(1:3,i), trim(plist)
  enddo

  write(6,'(A)') ''
  write(6,'(A,L1,A,I0,A,I0,A,I0)') 'invariants      : has_i=', &
       pg_has_i, '  n_sigma_h=', pg_n_sigma_h, '  n_sigma_vtot=', &
       pg_n_sigma_vtot, '  n_C2_perp=', pg_n_c2_perp

  !---------------------------------------------------------------------
  ! AO/MO representation matrices and gates.
  !---------------------------------------------------------------------
  write(6,'(A)') ''
  call pg_report_ao(delta_mo)

  !---------------------------------------------------------------------
  ! Character table, from the installed data.
  !---------------------------------------------------------------------
  call pg_character_table_path(pg_group_name, path, ierr)
  if (ierr /= 0) then
    write(6,'(A,A)') 'character table : unavailable for ', &
         trim(pg_group_name)
    write(6,'(A,A)') 'searched        : ', trim(path)
    if (ierr == 1) then
      write(6,'(A)') 'QP_ROOT is not set, so the data directory is unknown'
    endif
    write(6,'(A)') 'install it with: make -C $QP_ROOT/src/symmetry install'
    return
  endif

  ! The table is reported through the providers of columns.irp.f, not by
  ! re-parsing the file here: this is what makes the driver a test of the
  ! column matcher and not only of the parser.
  write(6,'(A,A)') 'character table : ', trim(path)
  write(6,'(A,A,A,I4,A,I4,A,I4)') 'group           : ', &
       trim(adjustl(pg_group_name)), '   |G| = ', pg_n_elem, &
       '   irreps = ', pg_n_irrep, '   columns = ', pg_n_col
  write(6,'(A)') ''
  write(6,'(A8,2X,A4,A5,A2)',advance='no') 'irrep', 'dim', 'mult', ' |'
  do c = 1, pg_n_col
    write(6,'(1X,A8)',advance='no') trim(adjustl(pg_data_sym(c)))
  enddo
  write(6,'(A)') ''
  do a = 1, pg_n_irrep
    write(6,'(A8,2X,I4,I5,A2)',advance='no') &
         trim(adjustl(pg_irrep_name(a))), pg_irrep_dim(a), pg_irrep_mult(a), ' |'
    do c = 1, pg_n_col
      if (abs(pg_char(a,c)) < 5.d-3) then
        write(6,'(F9.3)',advance='no') 0.d0
      else
        write(6,'(F9.4)',advance='no') pg_char(a,c)
      endif
    enddo
    write(6,'(A)') ''
  enddo

  write(6,'(A)') ''
  write(6,'(A)') 'columns (a column is NOT a geometric class: it may be the union of'
  write(6,'(A)') 'several classes and of their inverse classes; the symmetry module page'
  write(6,'(A)') 'explains how a table column is matched to group elements)'
  write(6,'(A6,A10,A10,A)') 'col', 'symbol', 'elements', '  element indices'
  do c = 1, pg_n_col
    write(6,'(I6,3X,A8,I8,4X)',advance='no') c, &
         trim(adjustl(pg_data_sym(c))), pg_col_len(c)
    do g = pg_col_start(c), pg_col_start(c) + pg_col_len(c) - 1
      write(6,'(I4)',advance='no') pg_col_member(g)
    enddo
    write(6,'(A)') ''
  enddo

  ! Invariants of the matched table, printed so that a wrong match is visible
  ! without a debugger: the columns partition the group, and
  ! sum_a d_a^2/mult_a = :math:`|G|` (never sum_a d_a^2).
  n_mem = 0
  do c = 1, pg_n_col
    n_mem = n_mem + pg_col_len(c)
  enddo
  d = 0.d0
  do a = 1, pg_n_irrep
    d = d + dble(pg_irrep_dim(a))**2 / dble(pg_irrep_mult(a))
  enddo
  write(6,'(A)') ''
  write(6,'(A,I6,A,I4)') 'elements in the columns: ', n_mem, '   |G| = ', &
       pg_n_elem
  write(6,'(A,F10.4,A,I4)') 'sum_a d_a^2/mult_a     : ', d, '   |G| = ', &
       pg_n_elem

  !---------------------------------------------------------------------
  ! MO irrep labels and gates G3..G6.
  !---------------------------------------------------------------------
  call pg_report_mo(delta_mo, write_mo)
end subroutine pg_report


!-----------------------------------------------------------------------
! One-line rendering of an atom permutation, e.g. "1 3 2".
!-----------------------------------------------------------------------
subroutine pg_perm_string(perm, nat, s)
  implicit none
  BEGIN_DOC
  ! Comma-free list of ``perm(1:nat)`` in ``s``; truncated with "..." when ``s``
  ! is too short.
  END_DOC
  integer, intent(in) :: nat, perm(nat)
  character*(*), intent(out) :: s
  character*(12) :: item
  integer :: i, pos

  s = ''
  pos = 0
  do i = 1, nat
    write(item,'(I0)') perm(i)
    if (pos + len_trim(item) + 2 > len(s)) then
      s = trim(s)//' ...'
      return
    endif
    s(pos+1:pos+len_trim(item)) = trim(item)
    pos = pos + len_trim(item) + 1
  enddo
end subroutine pg_perm_string


!-----------------------------------------------------------------------
! AO representation matrices, and gates G-iso, G1, G2.
!-----------------------------------------------------------------------
subroutine pg_report_ao(delta_mo)
  implicit none
  BEGIN_DOC
  ! Prints the AO basis summary, then the gates:
  !
  ! ``G-iso``
  !     M_g^T S M_g = S. Authoritative test of the AO transform, with S the
  !     RAW ao_overlap of the EZFIO, tolerance 1e-7; needs no MOs, so it
  !     runs on a nuclei+AO EZFIO.
  !
  ! ``G1``
  !     C^T S C = 1. Orthonormality of the MO set in the AO metric; no group
  !     element enters it, and it runs whether or not the MO space is
  !     complete.
  !
  ! ``G2``
  !     X_g X_h = X_gh. Only meaningful when mo_num == ao_num.
  !
  ! delta_MO = max_g ``X_g^T X_g - 1`` is printed and NOT gated: it is G-iso
  ! transported into the MO metric (see the identity at the measurement).
  !
  ! Each gate prints the measured deviation and stops on failure. The 1e-7
  ! bound of G-iso is not machine precision: :c:data:`ao_overlap` itself is only
  ! group symmetric to ~3e-9 (see the comment at the measurement below).
  !
  ! ``delta_mo`` is returned so that the MO report can size its G5 window:
  ! it is the accuracy floor of every MO purity.
  END_DOC
  logical :: has_ao
  integer :: g, h, k, n_pair, ne, g_worst, h_worst, iw, jw, loc(2)
  integer :: i, j, it
  double precision, intent(out) :: delta_mo
  double precision :: dev, dev_g1, dev_g2, s, cnorm
  double precision, allocatable :: m(:,:), u(:,:), x(:,:), y(:,:), csc(:,:)
  double precision, allocatable :: cv(:), cw(:)

  delta_mo = 0.d0
  write(6,'(A)') 'AO basis'
  call ezfio_has_ao_basis_ao_num(has_ao)
  if (.not. has_ao) then
    write(6,'(A)') '  no AO basis in this EZFIO (ao_basis/ao_num missing)'
    write(6,'(A)') '  gates G-iso, G1, G2 skipped'
    return
  endif
  write(6,'(A,I0,A,I0,A,I0)') '  ao_num = ', ao_num, '   shells = ', &
       shell_num, '   mo_num = ', mo_num
  if (pg_basis_ierr /= 0) then
    write(6,'(A,A)') '  basis check: FAILED - ', trim(pg_basis_msg)
    write(6,'(A)') '  gates G-iso, G1, G2 skipped'
    return
  endif
  write(6,'(A)') '  basis check: every shell has a unique image shell'

  ! G-iso: with exact integrals and exact AOs, M_g^T S M_g = S to machine
  ! precision. This is the authoritative test of Sym^l(R_g) and of the
  ! B_p/B_q rescaling, and it needs no MOs.
  !
  ! Its floor is NOT machine precision but the symmetry-breaking noise of
  ! QP2's own overlap integrals: :math:`|S|` is only group symmetric to ~3e-9
  ! (formaldehyde: the C 1s / O 2p_x element that must vanish is 3.2e-9),
  ! and ``M_g^T S M_g - S`` inherits that one for one. The gate is 1e-7 on
  ! purpose: a wrong monomial representation, a wrong normalization factor
  ! or a wrong power ordering is O(1) here, three orders of magnitude above
  ! the noise, while the noise itself is not a property of the transform.
  allocate(m(ao_num,ao_num), u(ao_num,ao_num))
  dev = 0.d0
  g_worst = 0
  iw = 0
  jw = 0
  do g = 1, pg_n_elem
    call pg_ao_rep(g, m)
    u(1:ao_num,1:ao_num) = matmul(matmul(transpose(m), ao_overlap), m) &
         - ao_overlap
    loc = maxloc(abs(u(1:ao_num,1:ao_num)))
    if (abs(u(loc(1),loc(2))) > dev) then
      dev = abs(u(loc(1),loc(2)))
      g_worst = g
      iw = loc(1)
      jw = loc(2)
    endif
  enddo
  deallocate(m, u)
  write(6,'(A,ES12.4)') '  G-iso max |M_g^T S M_g - S| : ', dev
  if (dev > 1.d-7) then
    write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') '  worst element: g = ', &
         g_worst, ', AO (', iw, ',', jw, ') shells (', ao_shell(iw), ',', &
         ao_shell(jw), ')'
    write(6,'(A)') '  G-iso FAILED: the AO transform is not orthogonal'
    stop 1
  endif

  if (mo_num == 0) then
    write(6,'(A)') '  no MOs: gates G1, G2 skipped'
    return
  endif

  !---------------------------------------------------------------------
  ! G1: the MO set is orthonormal in the AO metric. No group
  ! element enters this test, it is well conditioned, and it must run whether
  ! or not the MO space is complete, so it comes first in the MO part.
  !---------------------------------------------------------------------
  allocate(csc(mo_num,mo_num))
  csc(1:mo_num,1:mo_num) = matmul(transpose(mo_coef), &
       matmul(ao_overlap, mo_coef))
  do k = 1, mo_num
    csc(k,k) = csc(k,k) - 1.d0
  enddo
  dev_g1 = maxval(abs(csc(1:mo_num,1:mo_num)))
  write(6,'(A,ES12.4)') '  G1    max |C^T S C - 1|      : ', dev_g1
  if (dev_g1 > 1.d-8) then
    write(6,'(A)') '  G1 FAILED: mo_coef is not orthonormal in the AO metric'
    stop 1
  endif

  !---------------------------------------------------------------------
  ! :math:`||C||`_2^2 = sigma_max(C)^2 = largest eigenvalue of C^T C (mo_num x
  ! mo_num), by power iteration from v = 1: deterministic, and no LAPACK.
  ! It is the amplification factor of the identity below.
  !---------------------------------------------------------------------
  csc(1:mo_num,1:mo_num) = matmul(transpose(mo_coef), mo_coef)
  allocate(cv(mo_num), cw(mo_num))
  cv(1:mo_num) = 1.d0
  cnorm = 0.d0
  do it = 1, 50
    cw(1:mo_num) = matmul(csc(1:mo_num,1:mo_num), cv(1:mo_num))
    dev = maxval(abs(cw(1:mo_num)))
    if (dev == 0.d0) exit
    cnorm = dev
    cv(1:mo_num) = cw(1:mo_num)/dev
  enddo
  deallocate(cv, cw, csc)

  !---------------------------------------------------------------------
  ! delta_MO: DIAGNOSTIC, never a gate. For square C with
  ! C^T S C = 1 one has C C^T = S^{-1}, hence EXACTLY
  !     X_g^T X_g - 1 = C^T (M_g^T S M_g - S) C = C^T deltaS_g C ,
  ! i.e. G-iso transported into the MO metric and amplified by at most
  ! :math:`||C||`_2^2. It is therefore 1 for ANY orthonormal MO set, localized or
  ! not, and it is not a symmetry-adaptation test. An absolute threshold on
  ! it is meaningless: benzene 6-31G measures 2.7e-7 at G-iso = 4.0e-9 on a
  ! geometry symmetric to 2.2e-15 bohr. What it does measure is the
  ! representation floor of the MO purity gates, so the report prints it and
  ! uses it to size its G5 window.
  !---------------------------------------------------------------------
  allocate(x(mo_num,mo_num), y(mo_num,mo_num))
  delta_mo = 0.d0
  do g = 1, pg_n_elem
    call pg_get_mo_rep(g, x)
    y(1:mo_num,1:mo_num) = matmul(transpose(x), x)
    do k = 1, mo_num
      y(k,k) = y(k,k) - 1.d0
    enddo
    delta_mo = max(delta_mo, maxval(abs(y(1:mo_num,1:mo_num))))
  enddo
  write(6,'(A,ES12.4,A)') '  delta_MO  max_g |X_g^T X_g - 1| : ', delta_mo, &
       '   (diagnostic only: exactly max_g |C^T (M_g^T S M_g - S) C|, implied by G-iso)'
  write(6,'(A,F12.1)') '  ||C||_2^2                     : ', cnorm
  write(6,'(A,ES12.4)') '  geom residual (bohr)          : ', pg_geom_residual
  if (pg_thresh < delta_mo) then
    write(6,'(A,ES10.4,A)') '  WARNING: pg_thresh (', pg_thresh, &
         ') is below the measured representation floor delta_MO; '// &
         'labels are not reliable on this basis'
  endif

  if (mo_num /= ao_num) then
    write(6,'(A,I0,A,I0,A)') '  G2 skipped: MO space incomplete (mo_num = ', &
         mo_num, ' /= ao_num = ', ao_num, '), representation property untested'
    deallocate(x, y)
    return
  endif

  ! G2: X_g X_h = X_{gh}. All pairs when the group is small; otherwise the 8
  ! first elements (which generate every group of the data set) plus 100
  ! fixed pseudo-random pairs: :math:`|G|^2 dgemms at :math:`|G|` = 128 is not a test, it
  ! is a benchmark.
  !
  ! The bound is 1e-7, not 1e-12, for two independent reasons, both larger
  ! than machine precision and neither of them a property of X:
  !   * the same integral noise that floors G-iso and delta_MO (~3e-9);
  !   * the group itself is known only to :c:data:`pg_tol`. A stored element is a
  !     representative whose matrix is accurate to :c:data:`pg_tol` / lever arm
  !     (~1e-8; measured 4.9e-9 on the C3 axis of NH3), so R_g R_h differs
  !     from the stored R_gh by that much, and X = C^T S M C inherits it one
  !     for one. The *table* is exact (products come from the permutations,
  !     never from a matrix comparison); only the matrices are approximate.
  ! A wrong product, a wrong monomial representation or a missing element is
  ! O(1) here.
  ne = pg_n_elem
  n_pair = 0
  dev_g2 = 0.d0
  g_worst = 0
  h_worst = 0
  if (ne <= 8) then
    do g = 1, ne
      do h = 1, ne
        call pg_g2_pair(g, h, dev)
        if (dev > dev_g2) then
          dev_g2 = dev
          g_worst = g
          h_worst = h
        endif
        n_pair = n_pair + 1
      enddo
    enddo
  else
    do g = 1, min(8, ne)
      do h = 1, min(8, ne)
        call pg_g2_pair(g, h, dev)
        if (dev > dev_g2) then
          dev_g2 = dev
          g_worst = g
          h_worst = h
        endif
        n_pair = n_pair + 1
      enddo
    enddo
    do k = 1, 100
      g = 1 + mod(k*7919, ne)
      h = 1 + mod(k*k*104729 + 13, ne)
      call pg_g2_pair(g, h, dev)
      if (dev > dev_g2) then
        dev_g2 = dev
        g_worst = g
        h_worst = h
      endif
      n_pair = n_pair + 1
    enddo
  endif
  write(6,'(A,I0,A,ES12.4)') '  G2    max |X_g X_h - X_gh|   : ', &
       n_pair, ' pairs, ', dev_g2
  if (dev_g2 > 1.d-7) then
    write(6,'(A,I0,A,I0,A,I0)') '  worst pair: g = ', g_worst, ', h = ', &
         h_worst, ', gh = ', pg_elem_mul(g_worst,h_worst)
    write(6,'(A)') '  G2 FAILED: X is not a representation of the group'
    stop 1
  endif
  deallocate(x, y)
end subroutine pg_report_ao

!-----------------------------------------------------------------------
! ``X_g X_h - X_{gh}`` for one pair, computed from scratch.
!-----------------------------------------------------------------------
subroutine pg_g2_pair(g, h, dev)
  implicit none
  BEGIN_DOC
  ! ``dev`` = max ``X_g X_h - X_{gh}`` for one pair. The product element comes
  ! from :c:data:`pg_elem_mul`, which is built from the permutations: no matrix is
  ! ever matched.
  END_DOC
  integer, intent(in) :: g, h
  double precision, intent(out) :: dev
  double precision, allocatable :: xg(:,:), xh(:,:), xp(:,:), d(:,:)
  integer :: gh, i, j

  gh = pg_elem_mul(g,h)
  allocate(xg(mo_num,mo_num), xh(mo_num,mo_num), xp(mo_num,mo_num), &
       d(mo_num,mo_num))
  call pg_get_mo_rep(g, xg)
  call pg_get_mo_rep(h, xh)
  call pg_get_mo_rep(gh, xp)
  do j = 1, mo_num
    do i = 1, mo_num
      d(i,j) = abs(sum(xg(i,1:mo_num)*xh(1:mo_num,j)) - xp(i,j))
    enddo
  enddo
  dev = maxval(d)
  deallocate(xg, xh, xp, d)
end subroutine pg_g2_pair

!-----------------------------------------------------------------------
! MO irrep labels, and gates G3..G6.
!-----------------------------------------------------------------------
subroutine pg_report_mo(delta_mo, write_mo)
  implicit none
  BEGIN_DOC
  ! Prints the label of every MO and the last four MO purity gates:
  !
  ! ``G3``
  !     every tr P^a is an integer (1e-6) and a multiple of d_a, i.e. n_a
  !     copies of irrep a. Authoritative test of the column matching:
  !     a sv/sd swap in an Abelian group changes n_a and is invisible to
  !     every earlier test.
  !
  ! ``G4``
  !     sum_a tr P^a = mo_num (1e-6).
  !
  ! ``G5``
  !     0 <= p_a(i) <= 1 (w5), sum_a p_a(i) = 1 (1e-6) and
  !     max_a p_a(i) = 1 (w5), with w5 = max(1e-8, 10*delta_MO). The last
  !     one is the MO-set-is-symmetry-adapted test.
  !
  ! ``G6``
  !     P^a P^a v = P^a v and sum_a P^a v = v (1e-6) for a fixed
  !     pseudo-random v: one pass over the group to build P^a v, then a
  !     second pass for P^a(P^a v). One X_g per element, never stored.
  !     For an irrep with w_a = 0 the ratio is pure representation noise,
  !     so the deviation of that projector is tested absolutely against
  !     the same 1e-6 (stricter than the ratio).
  !
  ! Every gate prints its measured value; a failing gate prints FAILED and
  ! the subroutine keeps going, so that one run shows every number. It stops
  ! with status 1 only at the very end.
  !
  ! The 1e-6 bound is not machine precision and must not be lowered: every
  ! X_g is built from the stored *representative* of an element, which is
  ! only as good as the geometry tolerance (measured: 2.3e-8 on X_g X_gh for
  ! the C3 axis of NH3 at pg_tol = 1e-4 bohr), and P^a is a linear
  ! combination of X_g, so :math:`||P^2 - P||` grows as (d_a/mult_a)^2 times that
  ! floor. A wrong column matching, a wrong mult_a or a
  ! missing irrep is O(1) in these gates.
  !
  ! ``delta_mo`` is the representation floor measured by the AO report (the
  ! diagnostic max_g ``X_g^T X_g - 1``, implied by G-iso). It sizes the G5
  ! window, so that a well-conditioned but wide basis cannot fail G5 on
  ! integral noise alone.
  END_DOC
  logical :: has_ao, failed
  integer :: a, c, g, i, j, m, n_mixed, a1, a2
  integer :: b, a0, nb, n_imp, n_bad, sdn
  double precision, intent(in) :: delta_mo
  logical, intent(in) :: write_mo
  double precision :: dev, s, wsum, na, pvnorm, w5
  double precision :: sn
  logical :: pure
  double precision :: dev_g3, dev_g3d, dev_g4, dev_g5a, dev_g5b, dev_g5c
  double precision :: dev_g6a, dev_g6e, dev_g6b
  integer :: n_g6e
  double precision, allocatable :: wraw(:), x(:,:), v(:), y(:), pv(:,:), pw(:,:)

  write(6,'(A)') ''
  write(6,'(A)') 'MO symmetry'
  call ezfio_has_ao_basis_ao_num(has_ao)
  if (.not. has_ao) then
    write(6,'(A)') '  no AO basis in this EZFIO (ao_basis/ao_num missing)'
    write(6,'(A)') '  MO labels and gates G3..G6 skipped'
    return
  endif
  if (mo_num == 0) then
    write(6,'(A)') '  no MOs in this EZFIO (mo_basis/mo_num is absent or 0)'
    write(6,'(A)') '  MO labels and gates G3..G6 skipped'
    return
  endif
  if (pg_basis_ierr /= 0) then
    write(6,'(A,A)') '  basis check: FAILED - ', trim(pg_basis_msg)
    write(6,'(A)') '  MO labels and gates G3..G6 skipped: pg_get_mo_rep '// &
         'would stop in pg_ao_rep_ex'
    return
  endif

  failed = .false.
  write(6,'(A,I0,A,I0,A,I0)') '  mo_num = ', mo_num, '   blocks = ', &
       pg_n_mo_block, '   |G| = ', pg_n_elem

  !---------------------------------------------------------------------
  ! Per-irrep table. w_a is the SUM OF THE ROWS of pg_mo_weight, the very
  ! number pg_irrep_weight rounds: taking the integer provider here would
  ! make the integer test of G3 vacuous.
  !---------------------------------------------------------------------
  allocate(wraw(pg_n_irrep))
  do a = 1, pg_n_irrep
    wsum = 0.d0
    do i = 1, mo_num
      wsum = wsum + pg_mo_weight(i,a)
    enddo
    wraw(a) = wsum
    if (abs(wraw(a)) < 1.d-10) wraw(a) = 0.d0
  enddo

  write(6,'(A)') ''
  write(6,'(A)') '  d_a = chi_a(E); mult_a = 1 for a real row and 2 for a row '// &
       'that is the sum'
  write(6,'(A)') '  of a complex-conjugate pair; w_a = tr P^a = sum_i p_a(i); '// &
       'n_a = w_a/d_a'
  write(6,'(A)') '     a  irrep        d_a  mult_a            w_a              n_a'
  do a = 1, pg_n_irrep
    na = wraw(a)/dble(pg_irrep_dim(a))
    write(6,'(I6,3X,A8,I5,I7,2(2X,F14.6))') a, &
         trim(adjustl(pg_irrep_name(a))), pg_irrep_dim(a), pg_irrep_mult(a), &
         wraw(a), na
  enddo

  !---------------------------------------------------------------------
  ! G3: tr P^a is an integer and a multiple of d_a.
  !---------------------------------------------------------------------
  write(6,'(A)') ''
  dev_g3 = 0.d0
  dev_g3d = 0.d0
  do a = 1, pg_n_irrep
    dev_g3 = max(dev_g3, abs(wraw(a) - dble(nint(wraw(a)))))
    na = wraw(a)/dble(pg_irrep_dim(a))
    dev_g3d = max(dev_g3d, abs(na - dble(nint(na))))
  enddo
  write(6,'(A,ES12.4)') '  G3    max_a |w_a - nint(w_a)|        : ', dev_g3
  write(6,'(A,ES12.4)') '  G3    max_a |w_a/d_a - nint(w_a/d_a)|: ', dev_g3d
  if (dev_g3 > 1.d-6) then
    write(6,'(A)') '  G3 FAILED: tr P^a is not an integer'
    failed = .true.
  endif
  if (dev_g3d > 1.d-6) then
    write(6,'(A)') '  G3 FAILED: tr P^a is not a multiple of d_a'
    failed = .true.
  endif

  !---------------------------------------------------------------------
  ! G4: the projectors account for every MO.
  !---------------------------------------------------------------------
  wsum = 0.d0
  do a = 1, pg_n_irrep
    do i = 1, mo_num
      wsum = wsum + pg_mo_weight(i,a)
    enddo
  enddo
  dev_g4 = abs(wsum - dble(mo_num))
  write(6,'(A,ES12.4)') '  G4    |sum_a w_a - mo_num|           : ', dev_g4
  if (dev_g4 > 1.d-6) then
    write(6,'(A)') '  G4 FAILED: the projectors do not sum to the identity'
    failed = .true.
  endif

  !---------------------------------------------------------------------
  ! G5: bounds, trace and completeness of the diagonal of the projectors.
  ! The window of the two purity clauses is w5 = max(1e-8, 10*delta_MO):
  ! delta_MO is the measured floor of the representation, so a
  ! well-conditioned but wide basis cannot fail G5 on integral noise alone.
  ! The 1e-6 of the sum rule is an identity of the characters, not a
  ! representation quantity, and stays absolute.
  !---------------------------------------------------------------------
  w5 = max(1.d-8, 10.d0*delta_mo)
  write(6,'(A,ES12.4)') '  G5    w5 = max(1e-8, 10*delta_MO)    : ', w5
  dev_g5a = 0.d0
  dev_g5b = 0.d0
  dev_g5c = 0.d0
  do i = 1, mo_num
    s = 0.d0
    do a = 1, pg_n_irrep
      s = s + pg_mo_weight(i,a)
      dev_g5c = max(dev_g5c, max(-pg_mo_weight(i,a), pg_mo_weight(i,a)-1.d0))
    enddo
    dev_g5a = max(dev_g5a, abs(s - 1.d0))
    dev_g5b = max(dev_g5b, abs(pg_mo_purity(i) - 1.d0))
  enddo
  write(6,'(A,ES12.4)') '  G5    max_i |sum_a p_a(i) - 1|       : ', dev_g5a
  write(6,'(A,ES12.4)') '  G5    max_i |max_a p_a(i) - 1|       : ', dev_g5b
  write(6,'(A,ES12.4)') '  G5    max_i,a max(-p_a, p_a-1)       : ', dev_g5c
  if (dev_g5a > 1.d-6) then
    write(6,'(A)') '  G5 FAILED: sum_a p_a(i) is not 1'
    failed = .true.
  endif
  if (dev_g5b > w5) then
    write(6,'(A)') '  G5 FAILED: the MO set is not symmetry adapted '// &
         '(some MO is not an eigenvector of the projectors)'
    failed = .true.
  endif
  if (dev_g5c > w5) then
    write(6,'(A)') '  G5 FAILED: a projector diagonal element is outside [0,1]'
    failed = .true.
  endif

  !---------------------------------------------------------------------
  ! G6: idempotence and completeness of the projectors, applied to a fixed
  ! pseudo-random v. v is cos(i) normalised: deterministic, machine
  ! independent, and it overlaps every block.
  !---------------------------------------------------------------------
  allocate(x(mo_num,mo_num), v(mo_num), y(mo_num))
  allocate(pv(mo_num,pg_n_irrep), pw(mo_num,pg_n_irrep))
  pvnorm = 0.d0
  do i = 1, mo_num
    v(i) = cos(dble(i))
    pvnorm = pvnorm + v(i)*v(i)
  enddo
  pvnorm = sqrt(pvnorm)
  v(1:mo_num) = v(1:mo_num)/pvnorm

  ! First pass: P^a v = (d_a/(mult_a*:math:`|G|`)) sum_c chi_a(c) sum_{g in c} X_g v.
  ! One build of X_g per element; the n_irrep accumulation of X_g v is done
  ! on the single vector y = X_g v, which is the same sum for less work.
  pv(1:mo_num,1:pg_n_irrep) = 0.d0
  do c = 1, pg_n_col
    do m = pg_col_start(c), pg_col_start(c) + pg_col_len(c) - 1
      g = pg_col_member(m)
      call pg_get_mo_rep(g, x)
      do i = 1, mo_num
        s = 0.d0
        do j = 1, mo_num
          s = s + x(i,j)*v(j)
        enddo
        y(i) = s
      enddo
      do a = 1, pg_n_irrep
        if (pg_char(a,c) == 0.d0) cycle
        do i = 1, mo_num
          pv(i,a) = pv(i,a) + pg_char(a,c)*y(i)
        enddo
      enddo
    enddo
  enddo
  do a = 1, pg_n_irrep
    s = dble(pg_irrep_dim(a))/(dble(pg_irrep_mult(a))*dble(pg_n_elem))
    do i = 1, mo_num
      pv(i,a) = s*pv(i,a)
    enddo
  enddo

  ! Second pass: the same for P^a(P^a v), i.e. n_irrep matvecs per element.
  pw(1:mo_num,1:pg_n_irrep) = 0.d0
  do c = 1, pg_n_col
    do m = pg_col_start(c), pg_col_start(c) + pg_col_len(c) - 1
      g = pg_col_member(m)
      call pg_get_mo_rep(g, x)
      do a = 1, pg_n_irrep
        if (pg_char(a,c) == 0.d0) cycle
        do j = 1, mo_num
          s = 0.d0
          do i = 1, mo_num
            s = s + x(j,i)*pv(i,a)
          enddo
          pw(j,a) = pw(j,a) + pg_char(a,c)*s
        enddo
      enddo
    enddo
  enddo
  do a = 1, pg_n_irrep
    s = dble(pg_irrep_dim(a))/(dble(pg_irrep_mult(a))*dble(pg_n_elem))
    do i = 1, mo_num
      pw(i,a) = s*pw(i,a)
    enddo
  enddo

  dev_g6a = 0.d0
  dev_g6e = 0.d0
  n_g6e = 0
  do a = 1, pg_n_irrep
    pvnorm = 0.d0
    dev = 0.d0
    do i = 1, mo_num
      pvnorm = pvnorm + pv(i,a)*pv(i,a)
      dev = dev + (pw(i,a)-pv(i,a))**2
    enddo
    pvnorm = sqrt(pvnorm)
    dev = sqrt(dev)
    ! An irrep absent from the MO set (w_a = 0, e.g. A2 for water) has
    ! P^a = 0 in exact arithmetic. Numerically P^a is only the
    ! representation error (~1e-9 here), so P^a v and P^a(P^a v) are two
    ! independent noise vectors of that size and their ratio is 1 whatever
    ! the transform does: the relative clause of G6 would test the noise.
    ! For those irreps the SAME 1e-6 bound is applied to the absolute
    ! deviation, which is at least as strict because :math:`||P^a v||` <= :math:`||v||` = 1.
    ! The test is selected by w_a, the projected rank printed above, not by
    ! :math:`||P^a v||`: a populated irrep whose image happens to be small must not
    ! escape the relative test.
    if (abs(wraw(a)) <= 1.d-6) then
      dev_g6e = max(dev_g6e, dev)
      n_g6e = n_g6e + 1
    else
      dev = dev/pvnorm
      dev_g6a = max(dev_g6a, dev)
    endif
  enddo
  deallocate(wraw)
  ! :math:`||v||` = 1 by construction, so the denominator below is 1.
  dev = 0.d0
  do i = 1, mo_num
    s = 0.d0
    do a = 1, pg_n_irrep
      s = s + pv(i,a)
    enddo
    dev = dev + (s - v(i))**2
  enddo
  dev_g6b = sqrt(dev)
  write(6,'(A,ES12.4)') '  G6    max_a ||P^a(P^a v)-P^a v||/||P^a v||: ', &
       dev_g6a
  if (n_g6e > 0) then
    write(6,'(A,I0,A,ES12.4)') '  G6    ', n_g6e, ' empty irrep(s), max_a '// &
         '||P^a(P^a v)-P^a v||: ', dev_g6e
  endif
  write(6,'(A,ES12.4)') '  G6    ||sum_a P^a v - v|| / ||v||    : ', dev_g6b
  if (dev_g6a > 1.d-6) then
    write(6,'(A)') '  G6 FAILED: a projector is not idempotent'
    failed = .true.
  endif
  if (dev_g6e > 1.d-6) then
    write(6,'(A)') '  G6 FAILED: a projector with no MO in this set '// &
         'does not vanish'
    failed = .true.
  endif
  if (dev_g6b > 1.d-6) then
    write(6,'(A)') '  G6 FAILED: the projectors do not sum to the identity'
    failed = .true.
  endif
  deallocate(x, v, y, pv, pw)

  !---------------------------------------------------------------------
  ! The MO table.
  !---------------------------------------------------------------------
  write(6,'(A)') ''
  write(6,'(A)') '  MO labels: i, irrep index (0 = no definite irrep), irrep '// &
       'name,'
  write(6,'(A)') '  symmetry block, block dimension, purity max_a p_a(i); a '// &
       'mixed MO adds'
  write(6,'(A)') '  its two largest irrep weights p_a(i).'
  write(6,'(A)') '      i  a  irrep       block  dim     purity'
  n_mixed = pg_n_mo_mixed
  do i = 1, mo_num
    if (pg_mo_irrep(i) == 0) then
      call pg_mo_top2(i, a1, a2)
      if (a2 > 0) then
        write(6,'(I7,I3,2X,A8,I6,I5,F11.6,2X,A,1X,F8.6,1X,A,1X,F8.6)') i, &
             0, 'mixed', pg_mo_block(i), pg_mo_block_dim(i), &
             pg_mo_purity(i), trim(pg_irrep_name(a1)), pg_mo_weight(i,a1), &
             trim(pg_irrep_name(a2)), pg_mo_weight(i,a2)
      else
        write(6,'(I7,I3,2X,A8,I6,I5,F11.6)') i, 0, 'mixed', &
             pg_mo_block(i), pg_mo_block_dim(i), pg_mo_purity(i)
      endif
    else
      write(6,'(I7,I3,2X,A8,I6,I5,F11.6)') i, pg_mo_irrep(i), &
           trim(adjustl(pg_irrep_name(pg_mo_irrep(i)))), pg_mo_block(i), &
           pg_mo_block_dim(i), pg_mo_purity(i)
    endif
  enddo
  write(6,'(A)') ''
  write(6,'(A,I0,A,I0,A)') '  MOs with no definite irrep: ', n_mixed, &
       ' of ', mo_num, '   (pg_mo_irrep = 0: max_a p_a < 1 - pg_thresh)'
  write(6,'(A)') '  The members of a degenerate block share ONE label and '// &
       'have no individual'
  write(6,'(A)') '  label: inside a block P^a is the d_a x d_a identity, so '// &
       'every partner has'
  write(6,'(A)') '  p_a = 1 and the a of the block is the a of each of its '// &
       'members.'

  !---------------------------------------------------------------------
  ! The invariant that catches a wrong block partition. It is
  ! a property of the BLOCKS, not of the labels: the label test acts on the
  ! weights, so it stays correct (every purity 1.000000) while the block
  ! structure is wrong. A block is one irrep space, and its size is a
  ! positive multiple of d_a of that label.
  !
  ! sum_a n_a counts the copies of each irrep and is the number of blocks an
  ! exact partition must give; fewer blocks than that can only be two copies
  ! of the SAME irrep accidentally degenerate, which SCF mixes into one
  ! invariant subspace. That is legitimate, but only once every block is
  ! pure and a multiple of d_a.
  !---------------------------------------------------------------------
  sn = 0.d0
  sdn = 0
  do a = 1, pg_n_irrep
    sn = sn + dble(pg_irrep_weight(a))/dble(pg_irrep_dim(a))
    sdn = sdn + pg_irrep_weight(a)
  enddo
  n_imp = 0
  n_bad = 0
  do b = 1, pg_n_mo_block
    nb = 0
    a0 = 0
    pure = .true.
    do i = 1, mo_num
      if (pg_mo_block(i) /= b) cycle
      if (nb == 0) then
        a0 = pg_mo_irrep(i)
      elseif (pg_mo_irrep(i) /= a0) then
        pure = .false.
      endif
      nb = nb + 1
    enddo
    if (nb == 0 .or. a0 <= 0) pure = .false.
    if (.not. pure) then
      n_imp = n_imp + 1
    elseif (mod(nb, pg_irrep_dim(a0)) /= 0) then
      n_bad = n_bad + 1
    endif
  enddo

  write(6,'(A)') ''
  write(6,'(A,I0,A,F9.4,A,I0)') '  blocks: pg_n_mo_block = ', &
       pg_n_mo_block, '   sum_a n_a =', sn, '   sum_a d_a n_a = ', sdn
  if (n_imp == 0 .and. n_bad == 0 .and. pg_n_mo_block == nint(sn)) then
    write(6,'(A)') '  every block carries one irrep and has a size that '// &
         'is a multiple of d_a'
  elseif (n_imp == 0 .and. n_bad == 0) then
    write(6,'(A)') '  blocks pure, fewer blocks than copies: accidental '// &
         'degeneracy of two copies of the same irrep (not a failure)'
  else
    write(6,'(A,I0,A,I0,A)') '  BLOCKS MISMATCH: ', n_imp, &
         ' block(s) mix irreps, ', n_bad, ' block(s) have a size that '// &
         'is not a multiple of d_a'
    failed = .true.
  endif

  if (failed) then
    write(6,'(A)') ''
    write(6,'(A)') '  MO symmetry: at least one gate FAILED (see above)'
    stop 1
  endif

  !---------------------------------------------------------------------
  ! Persistence. Only a checked labelling is written.
  !---------------------------------------------------------------------
  if (write_mo) then
    if (n_mixed /= 0) then
      write(6,'(A,I0,A)') '  PG_WRITE_MO_SYMMETRY ignored: ', n_mixed, &
           ' MO(s) have no definite irrep, and mo_symmetry would group them under 0'
    else
      call pg_write_mo_symmetry()
      write(6,'(A,I0,A)') '  mo_basis/mo_symmetry written: ', mo_num, &
           ' integers (irrep index, same integer = same irrep)'
    endif
  endif
end subroutine pg_report_mo
