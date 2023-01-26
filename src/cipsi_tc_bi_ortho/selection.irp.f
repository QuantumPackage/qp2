use bitmasks

! ---

subroutine select_connected(i_generator, E0, pt2_data, b, subset, csubset)

  use bitmasks
  use selection_types

  implicit none
  integer,                intent(in)    :: i_generator, subset, csubset
  double precision,       intent(in)    :: E0(N_states)
  type(selection_buffer), intent(inout) :: b
  type(pt2_type),         intent(inout) :: pt2_data

  integer                               :: k, l
  integer(bit_kind)                     :: hole_mask(N_int,2), particle_mask(N_int,2)
  double precision, allocatable         :: fock_diag_tmp(:,:)

  allocate(fock_diag_tmp(2,mo_num+1))

  call build_fock_tmp(fock_diag_tmp, psi_det_generators(1,1,i_generator), N_int)

  do k = 1, N_int
      hole_mask(k,1)     = iand(generators_bitmask(k,1,s_hole), psi_det_generators(k,1,i_generator))
      hole_mask(k,2)     = iand(generators_bitmask(k,2,s_hole), psi_det_generators(k,2,i_generator))
      particle_mask(k,1) = iand(generators_bitmask(k,1,s_part), not(psi_det_generators(k,1,i_generator)) )
      particle_mask(k,2) = iand(generators_bitmask(k,2,s_part), not(psi_det_generators(k,2,i_generator)) )
  enddo
  call select_singles_and_doubles(i_generator, hole_mask, particle_mask, fock_diag_tmp, E0, pt2_data, b, subset, csubset)

  deallocate(fock_diag_tmp)

end subroutine select_connected

! ---

subroutine select_singles_and_doubles(i_generator, hole_mask,particle_mask, fock_diag_tmp, E0, pt2_data, buf, subset, csubset)

  BEGIN_DOC
  !  WARNING /!\ : It is assumed that the generators and selectors are psi_det_sorted_tc
  END_DOC

  use bitmasks
  use selection_types
  implicit none

  integer,           intent(in)         :: i_generator, subset, csubset
  integer(bit_kind), intent(in)         :: hole_mask(N_int,2), particle_mask(N_int,2)
  double precision,  intent(in)         :: fock_diag_tmp(mo_num)
  double precision,  intent(in)         :: E0(N_states)
  type(pt2_type),         intent(inout) :: pt2_data
  type(selection_buffer), intent(inout) :: buf

  double precision, parameter           :: norm_thr = 1.d-16

  integer                               :: h1, h2, s1, s2, s3, i1, i2, ib, sp, k, i, j, nt, ii, sze
  integer                               :: maskInd
  integer                               :: N_holes(2), N_particles(2)
  integer                               :: hole_list(N_int*bit_kind_size,2)
  integer                               :: particle_list(N_int*bit_kind_size,2)
  integer                               :: l_a, nmax, idx
  integer                               :: nb_count, maskInd_save
  integer(bit_kind)                     :: hole(N_int,2), particle(N_int,2), mask(N_int, 2), pmask(N_int, 2)
  integer(bit_kind)                     :: mobMask(N_int, 2), negMask(N_int, 2)
  logical                               :: fullMatch, ok
  logical                               :: monoAdo, monoBdo
  logical                               :: monoBdo_save
  logical                               :: found

  integer, allocatable                  :: preinteresting(:), prefullinteresting(:)
  integer, allocatable                  :: interesting(:), fullinteresting(:)
  integer, allocatable                  :: tmp_array(:)
  integer, allocatable                  :: indices(:), exc_degree(:), iorder(:)
  integer(bit_kind), allocatable        :: minilist(:, :, :), fullminilist(:, :, :)
  logical, allocatable                  :: banned(:,:,:), bannedOrb(:,:)
  double precision, allocatable         :: coef_fullminilist_rev(:,:)
  double precision, allocatable         :: mat(:,:,:), mat_p(:,:,:), mat_m(:,:,:)


  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_tc_order psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order psi_selectors_coef_transp_tc
  PROVIDE psi_selectors_rcoef_bi_orth_transp psi_selectors_lcoef_bi_orth_transp

  PROVIDE banned_excitation

  monoAdo = .true.
  monoBdo = .true.

  do k = 1, N_int
    hole    (k,1) = iand(psi_det_generators(k,1,i_generator), hole_mask(k,1))
    hole    (k,2) = iand(psi_det_generators(k,2,i_generator), hole_mask(k,2))
    particle(k,1) = iand(not(psi_det_generators(k,1,i_generator)), particle_mask(k,1))
    particle(k,2) = iand(not(psi_det_generators(k,2,i_generator)), particle_mask(k,2))
  enddo

  call bitstring_to_list_ab(hole    , hole_list    , N_holes    , N_int)
  call bitstring_to_list_ab(particle, particle_list, N_particles, N_int)

  allocate( indices(N_det), exc_degree( max(N_det_alpha_unique, N_det_beta_unique) ) )

  ! Pre-compute excitation degrees wrt alpha determinants
  k = 1
  do i = 1, N_det_alpha_unique
    call get_excitation_degree_spin(psi_det_alpha_unique(1,i), psi_det_generators(1,1,i_generator), exc_degree(i), N_int)
  enddo

  ! Iterate on 0SD beta, and find alphas 0SDTQ such that exc_degree <= 4
  do j = 1, N_det_beta_unique
    call get_excitation_degree_spin(psi_det_beta_unique(1,j), psi_det_generators(1,2,i_generator), nt, N_int)
    if (nt > 2) cycle
    do l_a = psi_bilinear_matrix_columns_loc(j), psi_bilinear_matrix_columns_loc(j+1)-1
      i = psi_bilinear_matrix_rows(l_a)
      if(nt + exc_degree(i) <= 4) then
        idx = psi_det_sorted_tc_order(psi_bilinear_matrix_order(l_a))
        if (psi_average_norm_contrib_sorted_tc(idx) > norm_thr) then
          indices(k) = idx
          k = k + 1
        endif
      endif
    enddo
  enddo

  ! Pre-compute excitation degrees wrt beta determinants
  do i = 1, N_det_beta_unique
    call get_excitation_degree_spin(psi_det_beta_unique(1,i), psi_det_generators(1,2,i_generator), exc_degree(i), N_int)
  enddo

  ! Iterate on 0S alpha, and find betas TQ such that exc_degree <= 4
  ! Remove also contributions < 1.d-20)
  do j = 1, N_det_alpha_unique
    call get_excitation_degree_spin(psi_det_alpha_unique(1,j), psi_det_generators(1,1,i_generator), nt, N_int)
    if (nt > 1) cycle
    do l_a = psi_bilinear_matrix_transp_rows_loc(j), psi_bilinear_matrix_transp_rows_loc(j+1)-1
      i = psi_bilinear_matrix_transp_columns(l_a)
      if(exc_degree(i) < 3) cycle
      if(nt + exc_degree(i) <= 4) then
        idx = psi_det_sorted_tc_order(                                  &
            psi_bilinear_matrix_order(                               &
            psi_bilinear_matrix_transp_order(l_a)))
        if(psi_average_norm_contrib_sorted_tc(idx) > norm_thr) then
          indices(k) = idx
          k = k + 1
        endif
      endif
    enddo
  enddo

  deallocate(exc_degree)
  nmax = k - 1

  call isort_noidx(indices,nmax)

  ! Start with 32 elements. Size will double along with the filtering.
  allocate(preinteresting(0:32), prefullinteresting(0:32), interesting(0:32), fullinteresting(0:32))
  preinteresting(:) = 0
  prefullinteresting(:) = 0

  do i = 1, N_int
    negMask(i,1) = not(psi_det_generators(i,1,i_generator))
    negMask(i,2) = not(psi_det_generators(i,2,i_generator))
  enddo

  do k = 1, nmax

    i = indices(k)
    mobMask(1,1) = iand(negMask(1,1), psi_det_sorted_tc(1,1,i))
    mobMask(1,2) = iand(negMask(1,2), psi_det_sorted_tc(1,2,i))
    nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
    do j = 2, N_int
      mobMask(j,1) = iand(negMask(j,1), psi_det_sorted_tc(j,1,i))
      mobMask(j,2) = iand(negMask(j,2), psi_det_sorted_tc(j,2,i))
      nt = nt + popcnt(mobMask(j, 1)) + popcnt(mobMask(j, 2))
    enddo

    if(nt <= 4) then
      if(i <= N_det_selectors) then
        sze = preinteresting(0)
        if(sze+1 == size(preinteresting)) then
          allocate(tmp_array(0:sze))
          tmp_array(0:sze) = preinteresting(0:sze)
          deallocate(preinteresting)
          allocate(preinteresting(0:2*sze))
          preinteresting(0:sze) = tmp_array(0:sze)
          deallocate(tmp_array)
        endif
        preinteresting(0) = sze+1
        preinteresting(sze+1) = i
      elseif(nt <= 2) then
        sze = prefullinteresting(0)
        if(sze+1 == size(prefullinteresting)) then
          allocate (tmp_array(0:sze))
          tmp_array(0:sze) = prefullinteresting(0:sze)
          deallocate(prefullinteresting)
          allocate(prefullinteresting(0:2*sze))
          prefullinteresting(0:sze) = tmp_array(0:sze)
          deallocate(tmp_array)
        endif
        prefullinteresting(0) = sze+1
        prefullinteresting(sze+1) = i
      endif
    endif

  enddo
  deallocate(indices)

  allocate( banned(mo_num, mo_num,2), bannedOrb(mo_num, 2) )
  allocate( mat(N_states, mo_num, mo_num) )
  allocate( mat_p(N_states, mo_num, mo_num), mat_m(N_states, mo_num, mo_num) )
  maskInd = -1

  do s1 = 1, 2
    do i1 = N_holes(s1), 1, -1   ! Generate low excitations first

      found = .False.
      monoBdo_save = monoBdo
      maskInd_save = maskInd
      do s2 = s1, 2
        ib = 1
        if(s1 == s2) ib = i1+1
        do i2 = N_holes(s2), ib, -1
          maskInd = maskInd + 1
          if(mod(maskInd, csubset) == (subset-1)) then
            found = .True.
          end if
        enddo
        if(s1 /= s2) monoBdo = .false.
      enddo

      if (.not.found) cycle
      monoBdo = monoBdo_save
      maskInd = maskInd_save

      h1 = hole_list(i1,s1)
      call apply_hole(psi_det_generators(1,1,i_generator), s1, h1, pmask, ok, N_int)

      negMask = not(pmask)

      interesting(0) = 0
      fullinteresting(0) = 0

      do ii = 1, preinteresting(0)
        i = preinteresting(ii)
        select case(N_int)
          case(1)
            mobMask(1,1) = iand(negMask(1,1), psi_det_sorted_tc(1,1,i))
            mobMask(1,2) = iand(negMask(1,2), psi_det_sorted_tc(1,2,i))
            nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
          case(2)
            mobMask(1:2,1) = iand(negMask(1:2,1), psi_det_sorted_tc(1:2,1,i))
            mobMask(1:2,2) = iand(negMask(1:2,2), psi_det_sorted_tc(1:2,2,i))
            nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2)) +     &
                popcnt(mobMask(2, 1)) + popcnt(mobMask(2, 2))
          case(3)
            mobMask(1:3,1) = iand(negMask(1:3,1), psi_det_sorted_tc(1:3,1,i))
            mobMask(1:3,2) = iand(negMask(1:3,2), psi_det_sorted_tc(1:3,2,i))
            nt = 0
            do j = 3, 1, -1
              if (mobMask(j,1) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 1))
                if (nt > 4) exit
              endif
              if (mobMask(j,2) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 2))
                if (nt > 4) exit
              endif
            enddo
          case(4)
            mobMask(1:4,1) = iand(negMask(1:4,1), psi_det_sorted_tc(1:4,1,i))
            mobMask(1:4,2) = iand(negMask(1:4,2), psi_det_sorted_tc(1:4,2,i))
            nt = 0
            do j = 4, 1, -1
              if (mobMask(j,1) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 1))
                if (nt > 4) exit
              endif
              if (mobMask(j,2) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 2))
                if (nt > 4) exit
              endif
            enddo
          case default
            mobMask(1:N_int,1) = iand(negMask(1:N_int,1), psi_det_sorted_tc(1:N_int,1,i))
            mobMask(1:N_int,2) = iand(negMask(1:N_int,2), psi_det_sorted_tc(1:N_int,2,i))
            nt = 0
            do j = N_int, 1, -1
              if (mobMask(j,1) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 1))
                if (nt > 4) exit
              endif
              if (mobMask(j,2) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 2))
                if (nt > 4) exit
              endif
            enddo
        end select

        if(nt <= 4) then
          sze = interesting(0)
          if(sze+1 == size(interesting)) then
            allocate (tmp_array(0:sze))
            tmp_array(0:sze) = interesting(0:sze)
            deallocate(interesting)
            allocate(interesting(0:2*sze))
            interesting(0:sze) = tmp_array(0:sze)
            deallocate(tmp_array)
          endif
          interesting(0) = sze+1
          interesting(sze+1) = i
          if(nt <= 2) then
            sze = fullinteresting(0)
            if(sze+1 == size(fullinteresting)) then
              allocate (tmp_array(0:sze))
              tmp_array(0:sze) = fullinteresting(0:sze)
              deallocate(fullinteresting)
              allocate(fullinteresting(0:2*sze))
              fullinteresting(0:sze) = tmp_array(0:sze)
              deallocate(tmp_array)
            endif
            fullinteresting(0) = sze+1
            fullinteresting(sze+1) = i
          endif
        endif

      enddo

      do ii = 1, prefullinteresting(0)
        i = prefullinteresting(ii)
        nt = 0
        mobMask(1,1) = iand(negMask(1,1), psi_det_sorted_tc(1,1,i))
        mobMask(1,2) = iand(negMask(1,2), psi_det_sorted_tc(1,2,i))
        nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
        if (nt > 2) cycle
        do j=N_int,2,-1
          mobMask(j,1) = iand(negMask(j,1), psi_det_sorted_tc(j,1,i))
          mobMask(j,2) = iand(negMask(j,2), psi_det_sorted_tc(j,2,i))
          nt = nt+ popcnt(mobMask(j, 1)) + popcnt(mobMask(j, 2))
          if (nt > 2) exit
        end do

        if(nt <= 2) then
          sze = fullinteresting(0)
          if (sze+1 == size(fullinteresting)) then
            allocate (tmp_array(0:sze))
            tmp_array(0:sze) = fullinteresting(0:sze)
            deallocate(fullinteresting)
            allocate(fullinteresting(0:2*sze))
            fullinteresting(0:sze) = tmp_array(0:sze)
            deallocate(tmp_array)
          endif
          fullinteresting(0) = sze+1
          fullinteresting(sze+1) = i
        endif
      enddo

      allocate( fullminilist (N_int, 2, fullinteresting(0)), &
                    minilist (N_int, 2,     interesting(0))  )

      do i = 1, fullinteresting(0)
        do k = 1, N_int
          fullminilist(k,1,i) = psi_det_sorted_tc(k,1,fullinteresting(i))
          fullminilist(k,2,i) = psi_det_sorted_tc(k,2,fullinteresting(i))
        enddo
      enddo

      do i = 1, interesting(0)
        do k = 1, N_int
          minilist(k,1,i) = psi_det_sorted_tc(k,1,interesting(i))
          minilist(k,2,i) = psi_det_sorted_tc(k,2,interesting(i))
        enddo
      enddo

      do s2 = s1, 2
        sp = s1

        if(s1 /= s2) sp = 3

        ib = 1
        if(s1 == s2) ib = i1+1
        monoAdo = .true.
        do i2 = N_holes(s2), ib, -1   ! Generate low excitations first

          h2 = hole_list(i2,s2)
          call apply_hole(pmask, s2,h2, mask, ok, N_int)
          banned(:,:,1) = banned_excitation(:,:)
          banned(:,:,2) = banned_excitation(:,:)
          do j = 1, mo_num
            bannedOrb(j, 1) = .true.
            bannedOrb(j, 2) = .true.
          enddo
          do s3 = 1, 2
            do i = 1, N_particles(s3)
              bannedOrb(particle_list(i,s3), s3) = .false.
            enddo
          enddo
          if(s1 /= s2) then
            if(monoBdo) then
              bannedOrb(h1,s1) = .false.
            endif
            if(monoAdo) then
              bannedOrb(h2,s2) = .false.
              monoAdo = .false.
            endif
          endif

          maskInd = maskInd + 1
          if(mod(maskInd, csubset) == (subset-1)) then

            call spot_isinwf(mask, fullminilist, i_generator, fullinteresting(0), banned, fullMatch, fullinteresting)
            if(fullMatch) cycle

            call splash_pq(mask, sp, minilist, i_generator, interesting(0), bannedOrb, banned, mat, interesting, mat_p, mat_m)

            call fill_buffer_double(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2_data, mat, buf, mat_p, mat_m)
          endif

        enddo

        if(s1 /= s2) monoBdo = .false.
      enddo

      deallocate(fullminilist, minilist)

    enddo
  enddo

  deallocate(preinteresting, prefullinteresting, interesting, fullinteresting)
  deallocate(banned, bannedOrb,mat)
  deallocate(mat_p, mat_m)

end subroutine select_singles_and_doubles

! ---

subroutine spot_isinwf(mask, det, i_gen, N, banned, fullMatch, interesting)

  use bitmasks
  implicit none

  BEGIN_DOC
  ! Identify the determinants in det which are in the internal space. These are
  ! the determinants that can be produced by creating two particles on the mask.
  END_DOC

  integer, intent(in) :: i_gen, N
  integer, intent(in) :: interesting(0:N)
  integer(bit_kind),intent(in) :: mask(N_int, 2), det(N_int, 2, N)
  logical, intent(inout) :: banned(mo_num, mo_num)
  logical, intent(out) :: fullMatch

  integer :: i, j, na, nb, list(3)
  integer(bit_kind) :: myMask(N_int, 2), negMask(N_int, 2)

  fullMatch = .false.

  do i=1,N_int
    negMask(i,1) = not(mask(i,1))
    negMask(i,2) = not(mask(i,2))
  end do

  genl : do i=1, N
    ! If det(i) can't be generated by the mask, cycle
    do j=1, N_int
      if(iand(det(j,1,i), mask(j,1)) /= mask(j, 1)) cycle genl
      if(iand(det(j,2,i), mask(j,2)) /= mask(j, 2)) cycle genl
    end do

    ! If det(i) < det(i_gen), it hs already been considered
    if(interesting(i) < i_gen) then
      fullMatch = .true.
      return
    end if

    ! Identify the particles
    do j=1, N_int
      myMask(j, 1) = iand(det(j, 1, i), negMask(j, 1))
      myMask(j, 2) = iand(det(j, 2, i), negMask(j, 2))
    end do

    call bitstring_to_list_in_selection(myMask(1,1), list(1), na, N_int)
    call bitstring_to_list_in_selection(myMask(1,2), list(na+1), nb, N_int)
    banned(list(1), list(2)) = .true.
  end do genl

end subroutine spot_isinwf

! ---

subroutine splash_pq(mask, sp, det, i_gen, N_sel, bannedOrb, banned, mat, interesting, mat_p, mat_m)

  BEGIN_DOC
  ! Computes the contributions A(r,s) by
  ! comparing the external determinant to all the internal determinants det(i).
  ! an applying two particles (r,s) to the mask.
  END_DOC

  use bitmasks
  implicit none

  integer, intent(in)             :: sp, i_gen, N_sel
  integer, intent(in)             :: interesting(0:N_sel)
  integer(bit_kind),intent(in)    :: mask(N_int, 2), det(N_int, 2, N_sel)
  logical, intent(inout)          :: bannedOrb(mo_num, 2), banned(mo_num, mo_num, 2)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  double precision, intent(inout) :: mat_p(N_states, mo_num, mo_num), mat_m(N_states, mo_num, mo_num)

  integer                         :: i, ii, j, k, l, h(0:2,2), p(0:4,2), nt
  integer(bit_kind)               :: perMask(N_int, 2), mobMask(N_int, 2), negMask(N_int, 2)
  integer(bit_kind)               :: phasemask(N_int,2)


  PROVIDE psi_selectors_coef_transp_tc psi_det_sorted_tc
  PROVIDE psi_selectors_rcoef_bi_orth_transp psi_selectors_lcoef_bi_orth_transp


  mat   = 0d0
  mat_p = 0d0
  mat_m = 0d0

  do i = 1, N_int
    negMask(i,1) = not(mask(i,1))
    negMask(i,2) = not(mask(i,2))
  end do

  do i = 1, N_sel
    if(interesting(i) < 0) then
      stop 'prefetch interesting(i) and det(i)'
    endif

    mobMask(1,1) = iand(negMask(1,1), det(1,1,i))
    mobMask(1,2) = iand(negMask(1,2), det(1,2,i))
    nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))

    if(nt > 4) cycle

    do j = 2, N_int
      mobMask(j,1) = iand(negMask(j,1), det(j,1,i))
      mobMask(j,2) = iand(negMask(j,2), det(j,2,i))
      nt = nt + popcnt(mobMask(j, 1)) + popcnt(mobMask(j, 2))
    enddo

    if(nt > 4) cycle

    if (interesting(i) == i_gen) then
      if(sp == 3) then
        do k = 1, mo_num
          do j = 1, mo_num
            banned(j,k,2) = banned(k,j,1)
          enddo
        enddo
      else
        do k = 1, mo_num
          do l = k+1, mo_num
            banned(l,k,1) = banned(k,l,1)
          enddo
        enddo
      endif
    endif

    if (interesting(i) >= i_gen) then

        call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
        call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)

        perMask(1,1) = iand(mask(1,1), not(det(1,1,i)))
        perMask(1,2) = iand(mask(1,2), not(det(1,2,i)))
        do j=2,N_int
          perMask(j,1) = iand(mask(j,1), not(det(j,1,i)))
          perMask(j,2) = iand(mask(j,2), not(det(j,2,i)))
        end do
!        call get_d3_h  ( det(1,1,i), bannedOrb, banned, mat         , mask, p, sp, psi_selectors_coef_transp_tc (1, interesting(i)) )
!        call get_d3_htc( det(1,1,i), bannedOrb, banned, mat_m, mat_p, mask, p, sp, psi_selectors_rcoef_bi_orth_transp(1, interesting(i)) &
!                       , psi_selectors_lcoef_bi_orth_transp(1, interesting(i)) )

        call bitstring_to_list_in_selection(perMask(1,1), h(1,1), h(0,1), N_int)
        call bitstring_to_list_in_selection(perMask(1,2), h(1,2), h(0,2), N_int)

        call get_mask_phase(psi_det_sorted_tc(1,1,interesting(i)), phasemask,N_int)
        if(nt == 4) then
          call get_d2  (det(1,1,i), phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, psi_selectors_coef_transp_tc(1, 1, interesting(i)))
!          call get_pm2(det(1,1,i), phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, psi_selectors_coef_transp_tc(1, interesting(i)))
        elseif(nt == 3) then
          call get_d1 (det(1,1,i), phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, psi_selectors_coef_transp_tc(1, 1, interesting(i)))
!          call get_pm1(det(1,1,i), phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, psi_selectors_coef_transp_tc(1, interesting(i)))
        else
          call get_d0 (det(1,1,i), phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, psi_selectors_coef_transp_tc(1, 1, interesting(i)))
!          call get_pm0(det(1,1,i), phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, psi_selectors_coef_transp_tc(1, interesting(i)))
        endif
    elseif(nt == 4) then
        call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
        call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)
        call past_d2(banned, p, sp)
    elseif(nt == 3) then
        call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
        call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)
        call past_d1(bannedOrb, p)
    endif
  enddo

end subroutine splash_pq

! ---

subroutine fill_buffer_double(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2_data, mat, buf, mat_p, mat_m)

  use bitmasks
  use selection_types
  implicit none

  integer,          intent(in)          :: i_generator, sp, h1, h2
  double precision, intent(in)          :: mat(N_states, mo_num, mo_num)
  double precision, intent(in)          :: mat_p(N_states, mo_num, mo_num), mat_m(N_states, mo_num, mo_num)
  logical,          intent(in)          :: bannedOrb(mo_num, 2), banned(mo_num, mo_num)
  double precision, intent(in)          :: fock_diag_tmp(mo_num)
  double precision, intent(in)          :: E0(N_states)
  type(pt2_type),         intent(inout) :: pt2_data
  type(selection_buffer), intent(inout) :: buf

  integer                               :: iii, s, degree
  integer                               :: s1, s2, p1, p2, ib, j, istate, jstate
  integer                               :: info, k , iwork(N_states+1)
  integer(bit_kind)                     :: occ(N_int,2), n
  integer(bit_kind)                     :: mask(N_int, 2), det(N_int, 2)
  logical                               :: do_cycle, ok, do_diag
  double precision                      :: delta_E, val, Hii, w, tmp, alpha_h_psi
  double precision                      :: E_shift
  double precision                      :: i_h_alpha, alpha_h_i, psi_h_alpha
  double precision                      :: e_pert(N_states), coef(N_states)
  double precision                      :: s_weight(N_states,N_states)
  double precision                      :: eigvalues(N_states+1)
  double precision                      :: work(1+6*(N_states+1)+2*(N_states+1)**2)

  integer,          external            :: number_of_holes, number_of_particles
  logical,          external            :: is_a_two_holes_two_particles
  logical,          external            :: is_a_1h1p
  double precision, external            :: diag_H_mat_elem_fock


  PROVIDE dominant_dets_of_cfgs N_dominant_dets_of_cfgs

  do jstate = 1, N_states
    do istate = 1, N_states
      s_weight(istate,jstate) = dsqrt(selection_weight(istate)*selection_weight(jstate))
    enddo
  enddo

  if(sp == 3) then
    s1 = 1
    s2 = 2
  else
    s1 = sp
    s2 = sp
  end if
  call apply_holes(psi_det_generators(1,1,i_generator), s1, h1, s2, h2, mask, ok, N_int)
  E_shift = 0.d0

  if (h0_type == 'CFG') then
    j = det_to_configuration(i_generator)
    E_shift = psi_det_Hii(i_generator) - psi_configuration_Hii(j)
  endif

  do p1 = 1, mo_num

    if(bannedOrb(p1, s1)) cycle
    ib = 1
    if(sp /= 3) ib = p1+1

    do p2 = ib, mo_num

      if(bannedOrb(p2, s2)) cycle
      if(banned(p1,p2)) cycle

      ! TODO ??
      !if(pseudo_sym)then
      !  if(dabs(mat(1, p1, p2)).lt.thresh_sym)then
      !    w = 0.d0
      !  endif
      !endif

      ! MANU: ERREUR dans les calculs puisque < I | H | J > = 0 
      ! n'implique pas < I | H_TC | J > = 0 ??
      !val = maxval(abs(mat(1:N_states, p1, p2)))
      !if( val == 0d0) cycle

      call apply_particles(mask, s1, p1, s2, p2, det, ok, N_int)

      if(do_only_cas) then
        if( number_of_particles(det) > 0 ) cycle
        if( number_of_holes(det) > 0 ) cycle
      endif

      if(do_ddci) then
        if(is_a_two_holes_two_particles(det)) cycle
      endif

      if(do_only_1h1p) then
        if(.not.is_a_1h1p(det)) cycle
      endif

      if(seniority_max >= 0) then
        s = 0
        do k = 1, N_int
          s = s + popcnt(ieor(det(k,1),det(k,2)))
        enddo
        if (s > seniority_max) cycle
      endif

      if(excitation_max >= 0) then
        do_cycle = .True.
        if(excitation_ref == 1) then
          call get_excitation_degree(HF_bitmask, det(1,1), degree, N_int)
          do_cycle = do_cycle .and. (degree > excitation_max)
        elseif(excitation_ref == 2) then
          do k = 1, N_dominant_dets_of_cfgs
            call get_excitation_degree(dominant_dets_of_cfgs(1,1,k), det(1,1), degree, N_int)
            do_cycle = do_cycle .and. (degree > excitation_max)
          enddo
        endif
        if(do_cycle) cycle
      endif

      if(excitation_alpha_max >= 0) then
        do_cycle = .True.
        if(excitation_ref == 1) then
          call get_excitation_degree_spin(HF_bitmask, det(1,1), degree, N_int)
          do_cycle = do_cycle .and. (degree > excitation_max)
        elseif (excitation_ref == 2) then
          do k = 1, N_dominant_dets_of_cfgs
            call get_excitation_degree_spin(dominant_dets_of_cfgs(1,1,k), det(1,1), degree, N_int)
            do_cycle = do_cycle .and. (degree > excitation_alpha_max)
          enddo
        endif
        if(do_cycle) cycle
      endif

      if(excitation_beta_max >= 0) then
        do_cycle = .True.
        if(excitation_ref == 1) then
          call get_excitation_degree_spin(HF_bitmask, det(1,2), degree, N_int)
          do_cycle = do_cycle .and. (degree > excitation_max)
        elseif(excitation_ref == 2) then
          do k = 1, N_dominant_dets_of_cfgs
            call get_excitation_degree(dominant_dets_of_cfgs(1,2,k), det(1,2), degree, N_int)
            do_cycle = do_cycle .and. (degree > excitation_beta_max)
          enddo
        endif
        if(do_cycle) cycle
      endif


      w = 0.d0

      e_pert = 0.d0
      coef = 0.d0
      do_diag = .False.

      ! psi_det_generators  --> |i> of psi_0
      ! psi_coef_generators --> c_i of psi_0
      !
      ! <alpha|H|psi_0> = \sum_i c_i <alpha|H|i>

        ! -------------------------------------------
        ! Non hermitian 
        ! c_alpha = <alpha|H(j)|psi_0>/delta_E(alpha)
        ! e_alpha = c_alpha * <psi_0|H(j)|alpha>
        ! <alpha|H|psi_0> and <psi_0|H|alpha>
        ! <det|H(j)|psi_0> and transpose 
        ! -------------------------------------------

        istate = 1
        call htilde_mu_mat_bi_ortho_tot(det, det, N_int, Hii)
        delta_E = E0(istate) - Hii + E_shift
        !delta_E = 1.d0

!        call get_excitation_degree( HF_bitmask, det, degree, N_int)

        double precision :: alpha_h_psi_tmp, psi_h_alpha_tmp
        psi_h_alpha_tmp = mat_m(istate, p1, p2)
        alpha_h_psi_tmp = mat_p(istate, p1, p2)
!
        psi_h_alpha = 0.d0
        alpha_h_psi = 0.d0
        do iii = 1, N_det
          call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,iii), det, N_int, i_h_alpha)
          call htilde_mu_mat_bi_ortho_tot(det, psi_det(1,1,iii), N_int, alpha_h_i)
!          psi_h_alpha += i_h_alpha * leigvec_tc_bi_orth(iii,1)
!          alpha_h_psi += alpha_h_i * reigvec_tc_bi_orth(iii,1) 
          psi_h_alpha += i_h_alpha * 1.d0
          alpha_h_psi += alpha_h_i * 1.d0
        enddo
!          print*,'---',p1,p2
!          call debug_det(det,N_int)
!          print*,psi_h_alpha    *alpha_h_psi,    psi_h_alpha,    alpha_h_psi  
!          print*,psi_h_alpha_tmp*alpha_h_psi_tmp,psi_h_alpha_tmp,alpha_h_psi_tmp  
!         if(dabs(psi_h_alpha - psi_h_alpha_tmp).gt.1.d-10 .or. dabs(alpha_h_psi - alpha_h_psi_tmp).gt.1.d-10)then
!        if(dabs(psi_h_alpha_tmp*alpha_h_psi_tmp).gt.1.d+10)then
        if(dabs(psi_h_alpha*alpha_h_psi - psi_h_alpha_tmp*alpha_h_psi_tmp).gt.1.d-10)then
!          print*,'---'
!          print*,psi_h_alpha    *alpha_h_psi,    psi_h_alpha,    alpha_h_psi  
!          print*,psi_h_alpha_tmp*alpha_h_psi_tmp,psi_h_alpha_tmp,alpha_h_psi_tmp  
         call debug_det(det,N_int)
          print*,dabs(psi_h_alpha*alpha_h_psi - psi_h_alpha_tmp*alpha_h_psi_tmp),psi_h_alpha    *alpha_h_psi,psi_h_alpha_tmp*alpha_h_psi_tmp
          print*,'-- Good '
          print*,   psi_h_alpha,    alpha_h_psi  
          print*,'-- bad '
          print*,psi_h_alpha_tmp,alpha_h_psi_tmp  
          print*,'-- details good'
        double precision :: accu_1, accu_2
        accu_1 = 0.d0
        accu_2 = 0.d0
        do iii = 1, N_det
          call get_excitation_degree( psi_det(1,1,iii), det, degree, N_int)
          call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,iii), det, N_int, i_h_alpha)
          call htilde_mu_mat_bi_ortho_tot(det, psi_det(1,1,iii), N_int, alpha_h_i)
          print*,iii,degree,i_h_alpha,alpha_h_i
          accu_1 += i_h_alpha
          accu_2 += alpha_h_i
          print*,accu_1,accu_2
          
        enddo
!          if(dabs(psi_h_alpha*alpha_h_psi).gt.1.d-10)then
!          print*,p1,p2
!          print*,det(1,1), det(1,2)
!          call debug_det(det,N_int)
!          print*,psi_h_alpha    *alpha_h_psi,    psi_h_alpha,    alpha_h_psi  
!          print*,psi_h_alpha_tmp*alpha_h_psi_tmp,psi_h_alpha_tmp,alpha_h_psi_tmp  
!          print*, dabs(psi_h_alpha*alpha_h_psi - psi_h_alpha_tmp*alpha_h_psi_tmp),& 
!                   psi_h_alpha    *alpha_h_psi,psi_h_alpha_tmp*alpha_h_psi_tmp
        stop
          endif
!          endif
!          stop
!        endif

        !if(alpha_h_psi*psi_h_alpha/delta_E.gt.1.d-10)then
        !  print*, 'E0,Hii,E_shift'
        !  print*, E0(istate), Hii, E_shift
        !  print*, psi_h_alpha, alpha_h_psi, delta_E
        !  print*, psi_h_alpha * alpha_h_psi / delta_E
        !  !if(Hii .lt. E0(istate)) then
        !  !  call debug_det(det, N_int)
        !  !  print*, ' |E0| < |Hii| !!!'
        !  !  print*, ' E0  = ', E0(istate)
        !  !  print*, ' Hii = ', Hii
        !  !endif
        !endif

        coef(istate)   = alpha_h_psi / delta_E 
        e_pert(istate) = coef(istate) * psi_h_alpha
        if(selection_tc     ==  1 )then
         if(e_pert(istate).lt.0.d0)then
          e_pert(istate) = 0.d0
         endif
        else if(selection_tc == -1)then
         if(e_pert(istate).gt.0.d0)then
          e_pert(istate) = 0.d0
         endif
        endif
         

        !if(e_pert(istate) .gt. 1.d-15) then
        !  print*, 'E0,Hii,E_shift'
        !  print*, E0(istate), Hii, E_shift
        !  print*, psi_h_alpha, alpha_h_psi, delta_E
        !  print*, psi_h_alpha*alpha_h_psi/delta_E
        !endif

!      elseif(cipsi_tc == "h_tc_2x2") then


      do_diag = sum(dabs(coef)) > 0.001d0 .and. N_states > 1

      do istate = 1, N_states

        alpha_h_psi = mat(istate, p1, p2)

        pt2_data % overlap(:,istate) = pt2_data % overlap(:,istate) + coef(:) * coef(istate)
        pt2_data % variance(istate)  = pt2_data % variance(istate) + dabs(e_pert(istate))
        pt2_data % pt2(istate)       = pt2_data % pt2(istate)      + e_pert(istate)

        select case (weight_selection)
          case(5)
            ! Variance selection
            if (h0_type == 'CFG') then
              w = min(w, - alpha_h_psi * alpha_h_psi * s_weight(istate,istate)) & 
                / c0_weight(istate)
            else
              w = min(w, - alpha_h_psi * alpha_h_psi * s_weight(istate,istate))
            endif
          case(6)
            if (h0_type == 'CFG') then
              w = min(w,- coef(istate) * coef(istate) * s_weight(istate,istate)) &
                / c0_weight(istate)
            else
              w = min(w,- coef(istate) * coef(istate) * s_weight(istate,istate))
            endif
          case default
            ! Energy selection
            if (h0_type == 'CFG') then
              !w = min(w, e_pert(istate) * s_weight(istate,istate)) / c0_weight(istate)
              w = min(w, -dabs(e_pert(istate)) * s_weight(istate,istate)) / c0_weight(istate)
            else
              !w = min(w, e_pert(istate) * s_weight(istate,istate))
              w = min(w, -dabs( e_pert(istate) ) * s_weight(istate,istate))
            endif
        endselect
      enddo

      if(h0_type == 'CFG') then
        do k = 1, N_int
          occ(k,1) = ieor(det(k,1), det(k,2))
          occ(k,2) = iand(det(k,1), det(k,2))
        enddo
        call configuration_to_dets_size(occ, n, elec_alpha_num, N_int)
        n = max(n,1)
        w *= dsqrt(dble(n))
      endif

      if(w <= buf%mini) then
        call add_to_selection_buffer(buf, det, w)
      endif

    enddo ! end do p2
  enddo   ! end do p1

end subroutine fill_buffer_double

! ---

subroutine get_mask_phase(det1, pm, Nint)

  use bitmasks
  implicit none

  integer, intent(in) :: Nint
  integer(bit_kind), intent(in) :: det1(Nint,2)
  integer(bit_kind), intent(out) :: pm(Nint,2)
  integer(bit_kind) :: tmp1, tmp2
  integer :: i
  tmp1 = 0_8
  tmp2 = 0_8
  select case (Nint)

BEGIN_TEMPLATE
    case ($Nint)
      do i=1,$Nint
        pm(i,1) = ieor(det1(i,1), shiftl(det1(i,1), 1))
        pm(i,2) = ieor(det1(i,2), shiftl(det1(i,2), 1))
        pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 2))
        pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 2))
        pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 4))
        pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 4))
        pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 8))
        pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 8))
        pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 16))
        pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 16))
        pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 32))
        pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 32))
        pm(i,1) = ieor(pm(i,1), tmp1)
        pm(i,2) = ieor(pm(i,2), tmp2)
        if(iand(popcnt(det1(i,1)), 1) == 1) tmp1 = not(tmp1)
        if(iand(popcnt(det1(i,2)), 1) == 1) tmp2 = not(tmp2)
      end do
SUBST [ Nint ]
1;;
2;;
3;;
4;;
END_TEMPLATE
    case default
      do i=1,Nint
        pm(i,1) = ieor(det1(i,1), shiftl(det1(i,1), 1))
        pm(i,2) = ieor(det1(i,2), shiftl(det1(i,2), 1))
        pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 2))
        pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 2))
        pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 4))
        pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 4))
        pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 8))
        pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 8))
        pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 16))
        pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 16))
        pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 32))
        pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 32))
        pm(i,1) = ieor(pm(i,1), tmp1)
        pm(i,2) = ieor(pm(i,2), tmp2)
        if(iand(popcnt(det1(i,1)), 1) == 1) tmp1 = not(tmp1)
        if(iand(popcnt(det1(i,2)), 1) == 1) tmp2 = not(tmp2)
      end do
  end select

end subroutine get_mask_phase

! ---

subroutine past_d1(bannedOrb, p)

  use bitmasks
  implicit none

  logical, intent(inout) :: bannedOrb(mo_num, 2)
  integer, intent(in) :: p(0:4, 2)
  integer :: i,s

  do s = 1, 2
    do i = 1, p(0, s)
      bannedOrb(p(i, s), s) = .true.
    end do
  end do

end subroutine past_d1

! ---

subroutine past_d2(banned, p, sp)

  use bitmasks
  implicit none

  logical, intent(inout) :: banned(mo_num, mo_num)
  integer, intent(in) :: p(0:4, 2), sp
  integer :: i,j

  if(sp == 3) then
    do j=1,p(0,2)
      do i=1,p(0,1)
        banned(p(i,1), p(j,2)) = .true.
      end do
    end do
  else
    do i=1,p(0, sp)
      do j=1,i-1
        banned(p(j,sp), p(i,sp)) = .true.
        banned(p(i,sp), p(j,sp)) = .true.
      end do
    end do
  end if

end subroutine past_d2

! ---

subroutine bitstring_to_list_in_selection( string, list, n_elements, Nint)

  BEGIN_DOC
  ! Gives the inidices(+1) of the bits set to 1 in the bit string
  END_DOC

  use bitmasks
  implicit none

  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint)
  integer, intent(out)           :: list(Nint*bit_kind_size)
  integer, intent(out)           :: n_elements

  integer                        :: i, ishift
  integer(bit_kind)              :: l

  n_elements = 0
  ishift = 2
  do i=1,Nint
    l = string(i)
    do while (l /= 0_bit_kind)
      n_elements = n_elements+1
      list(n_elements) = ishift+popcnt(l-1_bit_kind) - popcnt(l)
      l = iand(l,l-1_bit_kind)
    enddo
    ishift = ishift + bit_kind_size
  enddo

end subroutine bitstring_to_list_in_selection

! ---

