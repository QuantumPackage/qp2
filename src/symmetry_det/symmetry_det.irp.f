program symmetry_det
  implicit none
  BEGIN_DOC
  ! Determinant and state symmetry driver.
  !
  ! Reads the wave function stored in the EZFIO directory given on the
  ! command line, labels every determinant with the irrep of the
  ! MO projector, and projects every CI state on the irreps.
  !
  ! Prints the determinant table, the per-irrep counts, the per-state
  ! weight vectors, and runs the gates D1-D7. Exits 1 if a gate fails, 0
  ! otherwise, so that it can be chained in a test script.
  END_DOC

  ! Without this the determinants module takes N_det = 1 and seeds
  ! psi_det with the Hartree-Fock determinant (read_wf = False means
  ! "build the space", not "read it"). Same pattern as src/tools/print_wf.
  read_wf = .True.
  touch read_wf
  call pg_det_report
end program symmetry_det

!-----------------------------------------------------------------------
! Report.
!-----------------------------------------------------------------------
subroutine pg_det_report
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Prints the labels and the gates, and stops with status 1 if any gate
  ! fails. Everything measured is printed before the verdict: one run
  ! shows the whole story.
  END_DOC
  integer, parameter :: ntmax = 200
  character*(256) :: soa, sob
  character*(8) :: lab
  integer :: i, a, s, k, ipos, nlab, nunl, jmax, nmax, nshow, a1, a2, nw
  integer :: list(N_int*bit_kind_size,2), n_el(2)
  double precision :: c, cmax, q, cache_mb, wv
  double precision :: w(pg_n_irrep)
  logical :: failed
  logical, external :: pg_det_gates

  if (pg_detection_failed) then
    write(6,'(A)') 'symmetry_det: point-group detection failed'
    stop 1
  endif

  write(6,'(A)') '============================================================'
  write(6,'(A)') ' symmetry_det driver: determinant and state symmetry'
  write(6,'(A)') '============================================================'
  write(6,'(A,A)')        'EZFIO directory : ', trim(ezfio_filename)
  write(6,'(A,A)')        'point group     : ', trim(pg_group_name)
  write(6,'(A,I0,A,I0)')  '|G|             : ', pg_n_elem, &
                          '    irreps: ', pg_n_irrep
  if (pg_n_mo_mixed == 0) then
    write(6,'(A,I0,A,I0,A)') 'MO basis        : ', mo_num, ' MOs in ', &
                            pg_n_mo_block, ' symmetry blocks, all labelled'
  else
    write(6,'(A,I0,A,I0,A,I0,A)') 'MO basis        : ', mo_num, ' MOs in ', &
                            pg_n_mo_block, ' symmetry blocks, ', &
                            pg_n_mo_mixed, ' without a definite irrep'
    write(6,'(A)') '  A mixed MO is a superposition: its two largest '// &
                   'irrep weights are'
    nshow = 0
    do i = 1, mo_num
      if (pg_mo_irrep(i) /= 0) cycle
      if (nshow >= 10) then
        write(6,'(A,I0,A)') '    ... (', pg_n_mo_mixed - 10, ' more)'
        exit
      endif
      call pg_mo_top2(i, a1, a2)
      if (a2 == 0) then
        write(6,'(A,I0,A,A,1X,F9.6)') '    MO ', i, ': ', &
             trim(pg_irrep_name(a1)), pg_mo_weight(i,a1)
      else
        write(6,'(A,I0,A,A,1X,F9.6,2X,A,1X,F9.6)') '    MO ', i, ': ', &
             trim(pg_irrep_name(a1)), pg_mo_weight(i,a1), &
             trim(pg_irrep_name(a2)), pg_mo_weight(i,a2)
      endif
      nshow = nshow + 1
    enddo
    write(6,'(A)') '  Determinant and state irreps below are '// &
                   'probabilities wherever they are not 1.'
  endif
  cache_mb = 8.d0*dble(mo_num)*dble(mo_num)*dble(pg_n_elem)/1048576.d0
  write(6,'(A,I0,A,I0,A,I0,A,F8.2,A)') 'X_g cache       : ', mo_num, ' x ', &
                          mo_num, ' x ', pg_n_elem, ' = ', cache_mb, ' MB'
  write(6,'(A,I0,A,I0)')  'N_det           : ', N_det, '    N_states: ', N_states
  write(6,'(A,L1)')       'read_wf         : ', read_wf
  nmax = 0
  if (N_det > 0) nmax = maxval(pg_det_block(1:N_det))
  write(6,'(A,I0,A,I0,A,I0)') 'families        : ', pg_n_family, &
                          '    largest family: ', nmax, &
                          '    pg_family_cap: ', pg_family_cap

  if (pg_label_ierr /= 0) then
    write(6,'(A,A)') 'symmetry_det: labelling failed: ', trim(pg_label_msg)
    stop 1
  endif

  !---------------------------------------------------------------------
  ! Determinant table.
  !---------------------------------------------------------------------
  nlab = 0
  nunl = 0
  do i = 1, N_det
    if (pg_det_irrep(i) /= 0) then
      nlab = nlab + 1
    else
      nunl = nunl + 1
    endif
  enddo
  write(6,'(A)') ''
  write(6,'(A,I0,A,I0,A)') 'determinants: ', nlab, ' with a definite irrep, ', &
       nunl, ' without'
  if (N_det <= ntmax) then
    write(6,'(A)') ''
    write(6,'(A)') ' det  family  |F|  irrep     purity    max|c|  st  '// &
         'alpha occupation            beta occupation'
    do i = 1, N_det
      if (pg_det_irrep(i) /= 0) then
        lab = pg_irrep_name(pg_det_irrep(i))
      else
        lab = '-'
      endif
      cmax = 0.d0
      jmax = 1
      do s = 1, N_states
        c = abs(psi_coef(i,s))
        if (c > cmax) then
          cmax = c
          jmax = s
        endif
      enddo
      call bitstring_to_list_ab(psi_det(1,1,i), list, n_el, N_int)
      soa = ''
      sob = ''
      ipos = 1
      do k = 1, min(n_el(1),12)
        write(soa(ipos:),'(I4)') list(k,1)
        ipos = ipos + 4
      enddo
      if (n_el(1) > 12) soa(ipos:) = ' ...'
      ipos = 1
      do k = 1, min(n_el(2),12)
        write(sob(ipos:),'(I4)') list(k,2)
        ipos = ipos + 4
      enddo
      if (n_el(2) > 12) sob(ipos:) = ' ...'
      write(6,'(I5,I8,I5,2X,A8,F9.6,1X,ES10.3,I4,2X,A25,1X,A25)') &
           i, pg_det_family(i), pg_det_block(i), trim(lab), &
           pg_det_purity(i), cmax, jmax, trim(soa), trim(sob)
    enddo
  else
    write(6,'(A,I0,A,I0,A)') '  (table of the first ', ntmax, &
         ' determinants hidden: N_det = ', N_det, '; the summary below is complete)'
  endif
  if (nunl > 0) then
    write(6,'(A)') ' determinants with no definite irrep (probability of each irrep):'
    write(6,'(A)') '    det  family  |F|   purity  composition'
    nunl = 0
    do i = 1, N_det
      if (pg_det_irrep(i) /= 0) cycle
      if (nunl >= 50) then
        write(6,'(A)') '   ...'
        exit
      endif
      call pg_det_weights(i, w)
      ! Self-check: pg_det_purity comes from the full family projector,
      ! w from its diagonal.
      if (abs(maxval(w(1:pg_n_irrep)) - pg_det_purity(i)) > 1.d-8) then
        write(6,'(A,I0,A,ES12.4)') '   (bookkeeping) determinant ', i, &
             ': diagonal and full projector differ by ', &
             abs(maxval(w(1:pg_n_irrep)) - pg_det_purity(i))
      endif
      write(6,'(I7,I8,I5,1X,F9.6,3X)',advance='no') i, pg_det_family(i), &
           pg_det_block(i), pg_det_purity(i)
      nw = 0
      do a = 1, pg_n_irrep
        if (w(a) < 1.d-4) cycle
        if (nw >= 4) then
          write(6,'(A)',advance='no') ' ...'
          exit
        endif
        write(6,'(A,1X,F8.6,2X)',advance='no') trim(pg_irrep_name(a)), w(a)
        nw = nw + 1
      enddo
      write(6,*)
      nunl = nunl + 1
    enddo
    write(6,'(A)') '   (weights below 1e-4 omitted, at most 4 shown). A '// &
         'determinant that does not fill'
    write(6,'(A)') '   whole blocks is a superposition of irreps; the '// &
         'distribution is the honest description.'
  endif

  !---------------------------------------------------------------------
  ! Counts per irrep.
  !---------------------------------------------------------------------
  write(6,'(A)') ''
  write(6,'(A)') ' determinants and weight per irrep (state 1):'
  do a = 1, pg_n_irrep
    nunl = 0
    q = 0.d0
    if (N_states > 0) then
      do i = 1, N_det
        if (pg_det_irrep(i) /= a) cycle
        nunl = nunl + 1
        q = q + psi_coef(i,1)*psi_coef(i,1)
      enddo
    endif
    write(6,'(2X,A8,A,I5,A,F12.8)') trim(pg_irrep_name(a)), '  n_det =', &
         nunl, '   sum_D c_D^2 =', q
  enddo

  !---------------------------------------------------------------------
  ! States.
  !---------------------------------------------------------------------
  write(6,'(A)') ''
  write(6,'(A)') ' states:'
  do s = 1, N_states
    q = 0.d0
    do a = 1, pg_n_irrep
      q = q + pg_state_irrep_weight(a,s)
    enddo
    if (pg_state_irrep(s) /= 0) then
      lab = pg_irrep_name(pg_state_irrep(s))
    else if (pg_state_dom(s) > 0) then
      ! No definite irrep: name the most probable one; w_max below says how
      ! much of the state it carries.
      lab = pg_irrep_name(pg_state_dom(s))
    else
      lab = 'none'
    endif
    write(6,'(A,I0,A,ES14.6,A,ES14.6)') '  state ', s, &
         ':  sum_a w_a =', q, '   sum_D c_D^2 =', pg_state_norm(s)
    write(6,'(A,A,A,F10.6,A,F10.6)') '    irrep = ', trim(lab), &
         '   w_max =', pg_state_purity(s), &
         '   mixed determinant fraction =', pg_det_mixed_frac(s)
    ipos = 0
    do a = 1, pg_n_irrep
      ! A projector weight is non-negative; a value that would print as
      ! -0.000000 is roundoff below the resolution of the line. Clamp only
      ! that, so a real negative weight is still visible.
      wv = pg_state_irrep_weight(a,s)
      if (wv < 0.d0 .and. wv > -5.d-7) wv = 0.d0
      write(6,'(A,A,A,F9.6,A)',advance='no') '    ', trim(pg_irrep_name(a)), &
           ':', wv, '  '
      ipos = ipos + 1
      if (mod(ipos,4) == 0 .and. a < pg_n_irrep) write(6,'(A)') ''
    enddo
    write(6,'(A)') ''
  enddo

  !---------------------------------------------------------------------
  ! Gates.
  !---------------------------------------------------------------------
  failed = .not. pg_det_gates()

  write(6,'(A)') ''
  if (failed) then
    write(6,'(A)') '  determinant/state symmetry: at least one gate FAILED'
    stop 1
  endif
  write(6,'(A)') '  all gates passed'
end subroutine pg_det_report

!-----------------------------------------------------------------------
! Gates D1-D7, plus the state gates.
!-----------------------------------------------------------------------
logical function pg_det_gates()
  use bitmasks
  implicit none
  BEGIN_DOC
  ! ``D1``
  !     sum_a T^a = 1 on every family (validates chi, mult_a and the
  !     family structure at once).
  !
  ! ``D2``
  !     (T^a)^2 = T^a and T^a Hermitian, every family, every a.
  !
  ! D2b tests the production element kernel (:c:func:`pg_family_vg`, which
  ! computes the ordered pairs), while the diagonal used for labels and
  ! purity is independently cross-checked against :c:func:`pg_det_weights`
  ! -- which still runs on :c:func:`pg_elem_me` -- by the report's
  ! self-check.
  !
  ! ``D3``
  !     -1e-8 <= purity <= 1+1e-8 for every determinant.
  !
  ! ``D4``
  !     ``sum_a w_a(s) - sum_D c_D(s)^2`` <= 1e-8: the accumulation loop
  !     reproduces the sum rule of the state weights.
  !
  ! ``D5``
  !     every labelled state has w_max >= 1-1e-6, and the mixed
  !     determinant fraction is reported.
  !
  ! ``D6``
  !     the dominant determinant of state 1 has irrep 1 for a closed-shell
  !     Hartree-Fock wave function (the cross-check that the determinant
  !     and MO labellings agree). Gated only when the
  !     determinant's family has one member: a block that is partially
  !     filled makes even a closed-shell occupation a mixture, so such a
  !     configuration has no one-determinant irrep to compare with. That
  !     case is printed, and the state weights carry the answer.
  !
  ! Families of size 1 are skipped by D1/D2: there T^a = delta is exact
  ! by construction, and the label already comes from that diagonal.
  !
  ! ``D7``
  !     the tensor contraction of :c:func:`pg_family_wadd` against the
  !     reference double sum ``sum_{ip,jp} c(ip,s) V_g(ip,jp) c(jp,s)``
  !     built from :c:func:`pg_family_vg`, on a deterministic sample of
  !     families. This is the gate that makes the W3 rewrite honest: the
  !     contraction is a reassociation of that sum, so it may differ in
  !     the last bits but nowhere else. At most ``d7_nfam`` families with
  !     ``np > 1`` are sampled (the first ones, so the sample is a fixed
  !     function of the wave function), and each family is checked on
  !     ``max(2, min(|G|, d7_pairs/np**2))`` evenly strided elements --
  !     g = 1 and at least one non-identity element are always included.
  !     The relative measure is ``max|q_contract - q_ref| / max|q_ref|``
  !     over the sample, gated at 1e-13.
  END_DOC
  integer :: ifam, i0, np, ip, jp, a, i, s, idom, jmax
  integer :: nsel(2,pg_n_mo_block)
  integer, allocatable :: mask(:,:,:)
  double precision, allocatable :: t(:,:,:), tt(:,:)
  double precision :: d1, d2a, d2b, d3min, d3max, d4, v, q, wmin
  double precision :: cmax, c
  logical :: closed, ok
  character*(32) :: nm
  ! D7: the production contraction against the reference kernel.
  ! At most this many families are sampled, and each family gets enough
  ! group elements that its reference costs about ``d7_pairs`` pairs
  ! (at least 2, so a non-identity element is always included, at most
  ! the whole group). Measured at 1 thread with callgrind:
  !   cisd_base      44.8M Ir =  6.7% of the 668.0M driver (200 families,
  !                  52,480 reference pairs; 82.5% of the cost is q_new)
  !   n2_rot_small   1.55G Ir =  4.3% of the 36.38G driver (101 families,
  !                  3,962,500 reference pairs; 1.27G of it is q_new,
  !                  because the contraction always runs all |G| elements
  !                  while only ng = 2 of them are compared)
  ! The plan's "first 200 families, all elements" would cost 245.7M
  ! reference pairs on n2_rot_small, about 16.5G Ir and +1.6 s there.
  ! D7 is reached only from :c:func:`pg_det_gates`, so ``cisd`` never
  ! pays for it.
  integer, parameter :: d7_nfam = 200, d7_pairs = 256
  integer :: nsam, ng, g, k, imin, npart, nfull
  integer :: fsec(2,2*pg_n_mo_block), psec(2,2*pg_n_mo_block)
  integer :: pn(2*pg_n_mo_block)
  integer :: mrad(2*pg_n_mo_block), strd(2*pg_n_mo_block)
  integer :: aoff(2*pg_n_mo_block)
  integer, allocatable :: pos(:,:,:), xrow(:,:,:), xcol(:,:,:), idx(:)
  integer, allocatable :: rep(:,:), si(:,:)
  double precision, allocatable :: vv(:,:), amat(:), qtot(:,:), cc(:,:)
  logical, allocatable :: gsel(:)
  double precision :: d7, d7ref, d7rel, qref
  integer, external :: pg_fam_amat_size

  pg_det_gates = .true.
  d1 = 0.d0
  d2a = 0.d0
  d2b = 0.d0
  do ifam = 1, pg_n_family
    i0 = pg_family_first(ifam)
    np = pg_det_block(i0)
    if (np <= 1) cycle
    call pg_det_counts(i0, nsel)
    allocate(mask(2,pg_n_mo_block,np), t(np,np,pg_n_irrep))
    allocate(tt(np,np))
    call pg_family_fill(nsel, np, mask)
    call pg_family_projector(np, nsel, mask, t)
    do jp = 1, np
      do ip = 1, np
        v = 0.d0
        do a = 1, pg_n_irrep
          v = v + t(ip,jp,a)
        enddo
        if (ip == jp) v = v - 1.d0
        d1 = max(d1, abs(v))
      enddo
    enddo
    do a = 1, pg_n_irrep
      do jp = 1, np
        do ip = 1, np
          d2b = max(d2b, abs(t(ip,jp,a) - t(jp,ip,a)))
        enddo
      enddo
      tt(1:np,1:np) = matmul(t(1:np,1:np,a), t(1:np,1:np,a))
      do jp = 1, np
        do ip = 1, np
          d2a = max(d2a, abs(tt(ip,jp) - t(ip,jp,a)))
        enddo
      enddo
    enddo
    deallocate(mask, t, tt)
  enddo

  d3min = 1.d0
  d3max = 0.d0
  do i = 1, N_det
    d3min = min(d3min, pg_det_purity(i))
    d3max = max(d3max, pg_det_purity(i))
  enddo

  d4 = 0.d0
  wmin = 2.d0
  do s = 1, N_states
    q = 0.d0
    do a = 1, pg_n_irrep
      q = q + pg_state_irrep_weight(a,s)
    enddo
    d4 = max(d4, abs(q - pg_state_norm(s)))
    if (pg_state_irrep(s) /= 0) wmin = min(wmin, pg_state_purity(s))
  enddo

  idom = 0
  jmax = 1
  if (N_states > 0 .and. N_det > 0) then
    cmax = -1.d0
    do i = 1, N_det
      c = abs(psi_coef(i,jmax))
      if (c > cmax) then
        cmax = c
        idom = i
      endif
    enddo
  endif

  !---------------------------------------------------------------------
  ! D7: the production contraction against the reference element kernel.
  !---------------------------------------------------------------------
  d7 = 0.d0
  d7ref = 0.d0
  nsam = 0
  if (N_states >= 1) then
    do ifam = 1, pg_n_family
      if (nsam >= d7_nfam) exit
      i0 = pg_family_first(ifam)
      np = pg_det_block(i0)
      if (np <= 1) cycle
      call pg_det_counts(i0, nsel)
      allocate(mask(2,pg_n_mo_block,np), idx(np), rep(2*pg_n_mo_block,np), &
               si(np,2*pg_n_mo_block), &
               pos(pg_det_nmax,np,2*pg_n_mo_block), &
               xrow(pg_det_nmax,np,2*pg_n_mo_block), &
               xcol(pg_det_nmax,np,2*pg_n_mo_block))
      call pg_family_fill(nsel, np, mask)
      if (pg_label_ierr /= 0) then
        deallocate(mask, idx, rep, si, pos, xrow, xcol)
        exit
      endif
      call pg_family_index(np, mask, idx, imin)
      call pg_fam_setup(np, nsel, mask, npart, nfull, fsec, psec, pn, pos, &
                        xrow, xcol, mrad, strd, rep, aoff, si)
      if (pg_label_ierr /= 0) then
        deallocate(mask, idx, rep, si, pos, xrow, xcol)
        exit
      endif
      allocate(cc(np,N_states), qtot(pg_n_elem,N_states), vv(np,np), &
               gsel(pg_n_elem), &
               amat(max(1, pg_fam_amat_size(npart, mrad))))
      cc(1:np,1:N_states) = 0.d0
      do s = 1, N_states
        do ip = 1, np
          if (idx(ip) > 0) cc(ip,s) = psi_coef(idx(ip),s)
        enddo
      enddo
      qtot(1:pg_n_elem,1:N_states) = 0.d0
      ! Reference on an evenly strided subset of the group: k = 1 gives
      ! g = 1, k = 2 already gives a non-identity element because ng <= |G|.
      ! The mask is filled with the same g expression as the reference loop
      ! below, so the elements the contraction computes and the elements
      ! the kernel checks cannot drift apart; everything else is skipped.
      ng = min(pg_n_elem, max(2, d7_pairs/(np*np)))
      gsel(1:pg_n_elem) = .false.
      do k = 1, ng
        g = 1 + (k-1)*pg_n_elem/ng
        gsel(g) = .true.
      enddo
      call pg_family_wadd(np, nsel, mask, cc, np, N_states, pg_xg(1,1,1), &
                          gsel, qtot)
      do k = 1, ng
        g = 1 + (k-1)*pg_n_elem/ng
        call pg_family_vg(g, np, nfull, fsec, npart, psec, pn, pos, xrow, &
                          xcol, mrad, rep, aoff, si, amat, pg_xg(1,1,1), vv)
        do s = 1, N_states
          qref = 0.d0
          do jp = 1, np
            do ip = 1, np
              qref = qref + cc(ip,s)*vv(ip,jp)*cc(jp,s)
            enddo
          enddo
          d7 = max(d7, abs(qtot(g,s) - qref))
          d7ref = max(d7ref, abs(qref))
        enddo
      enddo
      deallocate(mask, idx, rep, si, pos, xrow, xcol, cc, qtot, vv, amat, gsel)
      nsam = nsam + 1
    enddo
  endif

  write(6,'(A)') ''
  write(6,'(A)') ' gates:'
  ok = d1 <= 1.d-8
  write(6,'(A,ES12.4,A,L1)') '  D1  max_{F,D,D''} |sum_a T^a - 1|          = ', &
       d1, '   <= 1e-8   ', ok
  pg_det_gates = pg_det_gates .and. ok
  ok = d2a <= 1.d-6
  write(6,'(A,ES12.4,A,L1)') '  D2  max ||(T^a)^2 - T^a||                 = ', &
       d2a, '   <= 1e-6   ', ok
  pg_det_gates = pg_det_gates .and. ok
  ok = d2b <= 1.d-6
  write(6,'(A,ES12.4,A,L1)') '  D2  max ||T^a - (T^a)^T||                 = ', &
       d2b, '   <= 1e-6   ', ok
  pg_det_gates = pg_det_gates .and. ok
  ok = d3min >= -1.d-8 .and. d3max <= 1.d0 + 1.d-8
  write(6,'(A,2F12.8,A,L1)') '  D3  determinant purity range            = [', &
       d3min, d3max, ']  ', ok
  pg_det_gates = pg_det_gates .and. ok
  ok = d4 <= 1.d-8
  write(6,'(A,ES12.4,A,L1)') '  D4  max_s |sum_a w_a - sum_D c_D^2|       = ', &
       d4, '   <= 1e-8   ', ok
  pg_det_gates = pg_det_gates .and. ok
  if (wmin > 1.5d0) then
    write(6,'(A)') '  D5  no state has a definite irrep: nothing to check'
  else
    ok = wmin >= 1.d0 - 1.d-6
    write(6,'(A,F10.6,A,L1)') '  D5  min_a w_a of the labelled states      = ', &
         wmin, '   >= 1-1e-6  ', ok
    pg_det_gates = pg_det_gates .and. ok
  endif

  closed = .false.
  if (idom > 0 .and. N_det == 1) then
    closed = (elec_alpha_num == elec_beta_num)
    if (closed) closed = all(psi_det(1:N_int,1,idom) == psi_det(1:N_int,2,idom))
  endif
  if (idom == 0) then
    write(6,'(A)') '  D6  no dominant determinant (empty wave function)'
    pg_det_gates = .false.
  elseif (closed) then
    if (pg_det_block(idom) == 1) then
      if (pg_det_irrep(idom) > 0) then
        nm = trim(pg_irrep_name(pg_det_irrep(idom)))
      else
        nm = 'none'
      endif
      ok = (pg_det_irrep(idom) == 1)
      write(6,'(A,I3,A,A,A,L1)') '  D6  closed-shell HF determinant irrep     = ', &
           pg_det_irrep(idom), '  (', trim(nm), ')   ', ok
      pg_det_gates = pg_det_gates .and. ok
    else
      write(6,'(A,I0,A)') '  D6  closed-shell HF determinant: family of ', &
           pg_det_block(idom), ' members, no one-determinant irrep'
      write(6,'(A)') '      one block is partially filled: a closed-shell '// &
           'occupation of a degenerate shell is not a symmetry eigen-'// &
           'state, so the SCF solution itself breaks the point group '// &
           '(reported, not gated)'
    endif
  else
    write(6,'(A,I0,A,A)') '  D6  dominant determinant of state 1: irrep = ', &
         pg_det_irrep(idom), '  (not a closed-shell N_det=1 wave function: '// &
         'reported, not gated)'
  endif

  if (nsam == 0) then
    write(6,'(A)') '  D7  no family with more than one member: nothing to check'
  else
    ! d7ref is the largest |c^T V c| in the sample: the natural scale of the
    ! quantity. A vanishing reference reports the absolute deviation, which
    ! is the honest measure when there is no scale to divide by.
    if (d7ref > 0.d0) then
      d7rel = d7/d7ref
    else
      d7rel = d7
    endif
    ok = d7rel <= 1.d-13
    write(6,'(A,ES12.4,A,L1)') '  D7  max |q_contract - q_reference|/max|q_ref| = ', &
         d7rel, '   <= 1e-13   ', ok
    pg_det_gates = pg_det_gates .and. ok
  endif

  if (N_det > 1) then
    write(6,'(A)') '  note: a single-determinant irrep is a basis-state '// &
         'label; the state symmetry is the w_a vector above.'
  endif
end function pg_det_gates
