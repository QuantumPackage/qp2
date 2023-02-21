use bitmasks

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

end subroutine


subroutine select_connected(i_generator,E0,pt2_data,b,subset,csubset)
  use bitmasks
  use selection_types
  implicit none
  integer, intent(in)            :: i_generator, subset, csubset
  type(selection_buffer), intent(inout) :: b
  type(pt2_type), intent(inout)   :: pt2_data
  integer :: k,l
  double precision, intent(in)   :: E0(N_states)

  integer(bit_kind)              :: hole_mask(N_int,2), particle_mask(N_int,2)

  double precision, allocatable  :: fock_diag_tmp(:,:)

  allocate(fock_diag_tmp(2,mo_num+1))

  call build_fock_tmp(fock_diag_tmp,psi_det_generators(1,1,i_generator),N_int)

  do k=1,N_int
      hole_mask(k,1) = iand(generators_bitmask(k,1,s_hole), psi_det_generators(k,1,i_generator))
      hole_mask(k,2) = iand(generators_bitmask(k,2,s_hole), psi_det_generators(k,2,i_generator))
      particle_mask(k,1) = iand(generators_bitmask(k,1,s_part), not(psi_det_generators(k,1,i_generator)) )
      particle_mask(k,2) = iand(generators_bitmask(k,2,s_part), not(psi_det_generators(k,2,i_generator)) )
  enddo
  call select_singles_and_doubles(i_generator,hole_mask,particle_mask,fock_diag_tmp,E0,pt2_data,b,subset,csubset)
  deallocate(fock_diag_tmp)
end subroutine


double precision function get_phase_bi(phasemask, s1, s2, h1, p1, h2, p2, Nint)
  use bitmasks
  implicit none

  integer, intent(in) :: Nint
  integer(bit_kind), intent(in) :: phasemask(Nint,2)
  integer, intent(in) :: s1, s2, h1, h2, p1, p2
  logical :: change
  integer :: np
  double precision, save :: res(0:1) = (/1d0, -1d0/)

  integer :: h1_int, h2_int
  integer :: p1_int, p2_int
  integer :: h1_bit, h2_bit
  integer :: p1_bit, p2_bit
  h1_int = shiftr(h1-1,bit_kind_shift)+1
  h1_bit = h1 - shiftl(h1_int-1,bit_kind_shift)-1

  h2_int = shiftr(h2-1,bit_kind_shift)+1
  h2_bit = h2 - shiftl(h2_int-1,bit_kind_shift)-1

  p1_int = shiftr(p1-1,bit_kind_shift)+1
  p1_bit = p1 - shiftl(p1_int-1,bit_kind_shift)-1

  p2_int = shiftr(p2-1,bit_kind_shift)+1
  p2_bit = p2 - shiftl(p2_int-1,bit_kind_shift)-1


  ! Put the phasemask bits at position 0, and add them all
  h1_bit = int(shiftr(phasemask(h1_int,s1),h1_bit))
  p1_bit = int(shiftr(phasemask(p1_int,s1),p1_bit))
  h2_bit = int(shiftr(phasemask(h2_int,s2),h2_bit))
  p2_bit = int(shiftr(phasemask(p2_int,s2),p2_bit))

  np = h1_bit + p1_bit + h2_bit + p2_bit

  if(p1 < h1) np = np + 1
  if(p2 < h2) np = np + 1

  if(s1 == s2 .and. max(h1, p1) > min(h2, p2)) np = np + 1
  get_phase_bi = res(iand(np,1))
end


subroutine select_singles_and_doubles(i_generator,hole_mask,particle_mask,fock_diag_tmp,E0,pt2_data,buf,subset,csubset)
  use bitmasks
  use selection_types
  implicit none
  BEGIN_DOC
!            WARNING /!\ : It is assumed that the generators and selectors are psi_det_sorted
  END_DOC

  integer, intent(in)            :: i_generator, subset, csubset
  integer(bit_kind), intent(in)  :: hole_mask(N_int,2), particle_mask(N_int,2)
  double precision, intent(in)   :: fock_diag_tmp(mo_num)
  double precision, intent(in)   :: E0(N_states)
  type(pt2_type), intent(inout)   :: pt2_data
  type(selection_buffer), intent(inout) :: buf

  integer                         :: h1,h2,s1,s2,s3,i1,i2,ib,sp,k,i,j,nt,ii,sze
  integer(bit_kind)               :: hole(N_int,2), particle(N_int,2), mask(N_int, 2), pmask(N_int, 2)
  logical                         :: fullMatch, ok

  integer(bit_kind) :: mobMask(N_int, 2), negMask(N_int, 2)
  integer,allocatable               :: preinteresting(:), prefullinteresting(:)
  integer,allocatable               :: interesting(:), fullinteresting(:)
  integer,allocatable               :: tmp_array(:)
  integer(bit_kind), allocatable :: minilist(:, :, :), fullminilist(:, :, :)
  logical, allocatable           :: banned(:,:,:), bannedOrb(:,:)
  double precision, allocatable  :: coef_fullminilist_rev(:,:)


  double precision, allocatable   :: mat(:,:,:)

  logical :: monoAdo, monoBdo
  integer :: maskInd

  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_order psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order psi_selectors_coef_transp
  PROVIDE banned_excitation

  monoAdo = .true.
  monoBdo = .true.


  do k=1,N_int
    hole    (k,1) = iand(psi_det_generators(k,1,i_generator), hole_mask(k,1))
    hole    (k,2) = iand(psi_det_generators(k,2,i_generator), hole_mask(k,2))
    particle(k,1) = iand(not(psi_det_generators(k,1,i_generator)), particle_mask(k,1))
    particle(k,2) = iand(not(psi_det_generators(k,2,i_generator)), particle_mask(k,2))
  enddo


  integer                        :: N_holes(2), N_particles(2)
  integer                        :: hole_list(N_int*bit_kind_size,2)
  integer                        :: particle_list(N_int*bit_kind_size,2)

  call bitstring_to_list_ab(hole    , hole_list    , N_holes    , N_int)
  call bitstring_to_list_ab(particle, particle_list, N_particles, N_int)

  integer                        :: l_a, nmax, idx
  integer, allocatable           :: indices(:), exc_degree(:), iorder(:)

  ! Removed to avoid introducing determinants already presents in the wf
  !double precision, parameter :: norm_thr = 1.d-16

  allocate (indices(N_det),                                          &
      exc_degree(max(N_det_alpha_unique,N_det_beta_unique)))

  ! Pre-compute excitation degrees wrt alpha determinants
  k=1
  do i=1,N_det_alpha_unique
    call get_excitation_degree_spin(psi_det_alpha_unique(1,i),       &
        psi_det_generators(1,1,i_generator), exc_degree(i), N_int)
  enddo

  ! Iterate on 0SD beta, and find alphas 0SDTQ such that exc_degree <= 4
  do j=1,N_det_beta_unique
    call get_excitation_degree_spin(psi_det_beta_unique(1,j),        &
        psi_det_generators(1,2,i_generator), nt, N_int)
    if (nt > 2) cycle
    do l_a=psi_bilinear_matrix_columns_loc(j), psi_bilinear_matrix_columns_loc(j+1)-1
      i = psi_bilinear_matrix_rows(l_a)
      if (nt + exc_degree(i) <= 4) then
        idx = psi_det_sorted_order(psi_bilinear_matrix_order(l_a))
        ! Removed to avoid introducing determinants already presents in the wf
        !if (psi_average_norm_contrib_sorted(idx) > norm_thr) then
          indices(k) = idx
          k=k+1
        !endif
      endif
    enddo
  enddo

  ! Pre-compute excitation degrees wrt beta determinants
  do i=1,N_det_beta_unique
    call get_excitation_degree_spin(psi_det_beta_unique(1,i),        &
        psi_det_generators(1,2,i_generator), exc_degree(i), N_int)
  enddo

  ! Iterate on 0S alpha, and find betas TQ such that exc_degree <= 4
  ! Remove also contributions < 1.d-20)
  do j=1,N_det_alpha_unique
    call get_excitation_degree_spin(psi_det_alpha_unique(1,j),       &
        psi_det_generators(1,1,i_generator), nt, N_int)
    if (nt > 1) cycle
    do l_a=psi_bilinear_matrix_transp_rows_loc(j), psi_bilinear_matrix_transp_rows_loc(j+1)-1
      i = psi_bilinear_matrix_transp_columns(l_a)
      if (exc_degree(i) < 3) cycle
      if (nt + exc_degree(i) <= 4) then
        idx = psi_det_sorted_order(                                  &
            psi_bilinear_matrix_order(                               &
            psi_bilinear_matrix_transp_order(l_a)))
        ! Removed to avoid introducing determinants already presents in the wf
        !if (psi_average_norm_contrib_sorted(idx) > norm_thr) then
          indices(k) = idx
          k=k+1
        !endif
      endif
    enddo
  enddo

  deallocate(exc_degree)
  nmax=k-1

  call isort_noidx(indices,nmax)

  ! Start with 32 elements. Size will double along with the filtering.
  allocate(preinteresting(0:32), prefullinteresting(0:32),     &
      interesting(0:32), fullinteresting(0:32))
  preinteresting(:) = 0
  prefullinteresting(:) = 0

  do i=1,N_int
    negMask(i,1) = not(psi_det_generators(i,1,i_generator))
    negMask(i,2) = not(psi_det_generators(i,2,i_generator))
  end do

  do k=1,nmax
    i = indices(k)
    mobMask(1,1) = iand(negMask(1,1), psi_det_sorted(1,1,i))
    mobMask(1,2) = iand(negMask(1,2), psi_det_sorted(1,2,i))
    nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
    do j=2,N_int
      mobMask(j,1) = iand(negMask(j,1), psi_det_sorted(j,1,i))
      mobMask(j,2) = iand(negMask(j,2), psi_det_sorted(j,2,i))
      nt = nt + popcnt(mobMask(j, 1)) + popcnt(mobMask(j, 2))
    end do

    if(nt <= 4) then
      if(i <= N_det_selectors) then
        sze = preinteresting(0)
        if (sze+1 == size(preinteresting)) then
          allocate (tmp_array(0:sze))
          tmp_array(0:sze) = preinteresting(0:sze)
          deallocate(preinteresting)
          allocate(preinteresting(0:2*sze))
          preinteresting(0:sze) = tmp_array(0:sze)
          deallocate(tmp_array)
        endif
        preinteresting(0) = sze+1
        preinteresting(sze+1) = i
      else if(nt <= 2) then
        sze = prefullinteresting(0)
        if (sze+1 == size(prefullinteresting)) then
          allocate (tmp_array(0:sze))
          tmp_array(0:sze) = prefullinteresting(0:sze)
          deallocate(prefullinteresting)
          allocate(prefullinteresting(0:2*sze))
          prefullinteresting(0:sze) = tmp_array(0:sze)
          deallocate(tmp_array)
        endif
        prefullinteresting(0) = sze+1
        prefullinteresting(sze+1) = i
      end if
    end if
  end do
  deallocate(indices)

!  !$OMP CRITICAL
!  print *,  'Step1: ', i_generator, preinteresting(0)
!  !$OMP END CRITICAL

  allocate(banned(mo_num, mo_num,2), bannedOrb(mo_num, 2))
  allocate (mat(N_states, mo_num, mo_num))
  maskInd = -1

  integer                        :: nb_count, maskInd_save
  logical                        :: monoBdo_save
  logical                        :: found
  do s1=1,2
    do i1=N_holes(s1),1,-1   ! Generate low excitations first

      found = .False.
      monoBdo_save = monoBdo
      maskInd_save = maskInd
      do s2=s1,2
        ib = 1
        if(s1 == s2) ib = i1+1
        do i2=N_holes(s2),ib,-1
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
      call apply_hole(psi_det_generators(1,1,i_generator), s1,h1, pmask, ok, N_int)

      negMask = not(pmask)

      interesting(0) = 0
      fullinteresting(0) = 0

      do ii=1,preinteresting(0)
        i = preinteresting(ii)
        select case (N_int)
          case (1)
            mobMask(1,1) = iand(negMask(1,1), psi_det_sorted(1,1,i))
            mobMask(1,2) = iand(negMask(1,2), psi_det_sorted(1,2,i))
            nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
          case (2)
            mobMask(1:2,1) = iand(negMask(1:2,1), psi_det_sorted(1:2,1,i))
            mobMask(1:2,2) = iand(negMask(1:2,2), psi_det_sorted(1:2,2,i))
            nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2)) +     &
                popcnt(mobMask(2, 1)) + popcnt(mobMask(2, 2))
          case (3)
            mobMask(1:3,1) = iand(negMask(1:3,1), psi_det_sorted(1:3,1,i))
            mobMask(1:3,2) = iand(negMask(1:3,2), psi_det_sorted(1:3,2,i))
            nt = 0
            do j=3,1,-1
              if (mobMask(j,1) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 1))
                if (nt > 4) exit
              endif
              if (mobMask(j,2) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 2))
                if (nt > 4) exit
              endif
            end do
          case (4)
            mobMask(1:4,1) = iand(negMask(1:4,1), psi_det_sorted(1:4,1,i))
            mobMask(1:4,2) = iand(negMask(1:4,2), psi_det_sorted(1:4,2,i))
            nt = 0
            do j=4,1,-1
              if (mobMask(j,1) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 1))
                if (nt > 4) exit
              endif
              if (mobMask(j,2) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 2))
                if (nt > 4) exit
              endif
            end do
          case default
            mobMask(1:N_int,1) = iand(negMask(1:N_int,1), psi_det_sorted(1:N_int,1,i))
            mobMask(1:N_int,2) = iand(negMask(1:N_int,2), psi_det_sorted(1:N_int,2,i))
            nt = 0
            do j=N_int,1,-1
              if (mobMask(j,1) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 1))
                if (nt > 4) exit
              endif
              if (mobMask(j,2) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 2))
                if (nt > 4) exit
              endif
            end do
        end select

        if(nt <= 4) then
          sze = interesting(0)
          if (sze+1 == size(interesting)) then
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
          end if
        end if

      end do

      do ii=1,prefullinteresting(0)
        i = prefullinteresting(ii)
        nt = 0
        mobMask(1,1) = iand(negMask(1,1), psi_det_sorted(1,1,i))
        mobMask(1,2) = iand(negMask(1,2), psi_det_sorted(1,2,i))
        nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
        if (nt > 2) cycle
        do j=N_int,2,-1
          mobMask(j,1) = iand(negMask(j,1), psi_det_sorted(j,1,i))
          mobMask(j,2) = iand(negMask(j,2), psi_det_sorted(j,2,i))
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
        end if
      end do

      allocate (fullminilist (N_int, 2, fullinteresting(0)), &
                    minilist (N_int, 2,     interesting(0)) )
!      if(pert_2rdm)then
!        allocate(coef_fullminilist_rev(N_states,fullinteresting(0)))
!        do i=1,fullinteresting(0)
!          do j = 1, N_states
!            coef_fullminilist_rev(j,i) = psi_coef_sorted(fullinteresting(i),j)
!          enddo
!        enddo
!      endif

      do i=1,fullinteresting(0)
        fullminilist(:,:,i) = psi_det_sorted(:,:,fullinteresting(i))
      enddo

      do i=1,interesting(0)
        minilist(:,:,i) = psi_det_sorted(:,:,interesting(i))
      enddo

      do s2=s1,2
        sp = s1

        if(s1 /= s2) sp = 3

        ib = 1
        if(s1 == s2) ib = i1+1
        monoAdo = .true.
        do i2=N_holes(s2),ib,-1   ! Generate low excitations first

          h2 = hole_list(i2,s2)
          call apply_hole(pmask, s2,h2, mask, ok, N_int)
          banned(:,:,1) = banned_excitation(:,:)
          banned(:,:,2) = banned_excitation(:,:)
          do j=1,mo_num
            bannedOrb(j, 1) = .true.
            bannedOrb(j, 2) = .true.
          enddo
          do s3=1,2
            do i=1,N_particles(s3)
              bannedOrb(particle_list(i,s3), s3) = .false.
            enddo
          enddo
          if(s1 /= s2) then
            if(monoBdo) then
              bannedOrb(h1,s1) = .false.
            end if
            if(monoAdo) then
              bannedOrb(h2,s2) = .false.
              monoAdo = .false.
            end if
          end if

          maskInd = maskInd + 1
          if(mod(maskInd, csubset) == (subset-1)) then

            call spot_isinwf(mask, fullminilist, i_generator, fullinteresting(0), banned, fullMatch, fullinteresting)
            if(fullMatch) cycle
! !$OMP CRITICAL
!  print *,  'Step3: ', i_generator, h1, interesting(0)
! !$OMP END CRITICAL

            call splash_pq(mask, sp, minilist, i_generator, interesting(0), bannedOrb, banned, mat, interesting)

!            if(.not.pert_2rdm)then
             call fill_buffer_double(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2_data, mat, buf)
!            else
!             call fill_buffer_double_rdm(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2_data, mat, buf,fullminilist, coef_fullminilist_rev, fullinteresting(0))
!            endif
          end if
        enddo
        if(s1 /= s2) monoBdo = .false.
      enddo
      deallocate(fullminilist,minilist)
!      if(pert_2rdm)then
!        deallocate(coef_fullminilist_rev)
!      endif
    enddo
  enddo
  deallocate(preinteresting, prefullinteresting, interesting, fullinteresting)
  deallocate(banned, bannedOrb,mat)
end subroutine



subroutine fill_buffer_double(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2_data, mat, buf)
  use bitmasks
  use selection_types
  implicit none

  integer, intent(in) :: i_generator, sp, h1, h2
  double precision, intent(in) :: mat(N_states, mo_num, mo_num)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num)
  double precision, intent(in)    :: fock_diag_tmp(mo_num)
  double precision, intent(in)    :: E0(N_states)
  type(pt2_type), intent(inout)   :: pt2_data
  type(selection_buffer), intent(inout) :: buf
  logical :: ok
  integer :: s1, s2, p1, p2, ib, j, istate, jstate
  integer(bit_kind) :: mask(N_int, 2), det(N_int, 2)
  double precision :: e_pert(N_states), coef(N_states)
  double precision :: delta_E, val, Hii, w, tmp, alpha_h_psi
  double precision, external :: diag_H_mat_elem_fock
  double precision :: E_shift
  double precision :: s_weight(N_states,N_states)
  PROVIDE dominant_dets_of_cfgs N_dominant_dets_of_cfgs
  do jstate=1,N_states
    do istate=1,N_states
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

  do p1=1,mo_num
    if(bannedOrb(p1, s1)) cycle
    ib = 1
    if(sp /= 3) ib = p1+1

    do p2=ib,mo_num

! -----
! /!\ Generating only single excited determinants doesn't work because a
! determinant generated by a single excitation may be doubly excited wrt
! to a determinant of the future. In that case, the determinant will be
! detected as already generated when generating in the future with a
! double excitation.
!
!      if (.not.do_singles) then
!        if ((h1 == p1) .or. (h2 == p2)) then
!          cycle
!        endif
!      endif
!
!      if (.not.do_doubles) then
!        if ((h1 /= p1).and.(h2 /= p2)) then
!          cycle
!        endif
!      endif
! -----

      if(bannedOrb(p2, s2)) cycle
      if(banned(p1,p2)) cycle

      if(pseudo_sym)then
        if(dabs(mat(1, p1, p2)).lt.thresh_sym)then
          w = 0.d0
        endif
      endif

      val = maxval(abs(mat(1:N_states, p1, p2)))
      if( val == 0d0) cycle
      call apply_particles(mask, s1, p1, s2, p2, det, ok, N_int)

      if (do_only_cas) then
        integer, external :: number_of_holes, number_of_particles
        if (number_of_particles(det)>0) then
          cycle
        endif
        if (number_of_holes(det)>0) then
          cycle
        endif
      endif

      if (do_ddci) then
        logical, external  :: is_a_two_holes_two_particles
        if (is_a_two_holes_two_particles(det)) then
          cycle
        endif
      endif

      if (do_only_1h1p) then
        logical, external :: is_a_1h1p
        if (.not.is_a_1h1p(det)) cycle
      endif

      if (seniority_max >= 0) then
        integer :: s
        s = 0
        do k=1,N_int
          s = s + popcnt(ieor(det(k,1),det(k,2)))
        enddo

        if (s > seniority_max) cycle
      endif


      integer :: degree
      logical :: do_cycle
      if (excitation_max >= 0) then
        do_cycle = .True.
        if (excitation_ref == 1) then
          call get_excitation_degree(HF_bitmask,det(1,1),degree,N_int)
          do_cycle = do_cycle .and. (degree > excitation_max)
        else if (excitation_ref == 2) then
          do k=1,N_dominant_dets_of_cfgs
            call get_excitation_degree(dominant_dets_of_cfgs(1,1,k),det(1,1),degree,N_int)
            do_cycle = do_cycle .and. (degree > excitation_max)
          enddo
        endif
        if (do_cycle) cycle
      endif


      if (excitation_alpha_max >= 0) then
        do_cycle = .True.
        if (excitation_ref == 1) then
          call get_excitation_degree_spin(HF_bitmask,det(1,1),degree,N_int)
          do_cycle = do_cycle .and. (degree > excitation_max)
        else if (excitation_ref == 2) then
          do k=1,N_dominant_dets_of_cfgs
            call get_excitation_degree_spin(dominant_dets_of_cfgs(1,1,k),det(1,1),degree,N_int)
            do_cycle = do_cycle .and. (degree > excitation_alpha_max)
          enddo
        endif
        if (do_cycle) cycle
      endif


      if (excitation_beta_max >= 0) then
        do_cycle = .True.
        if (excitation_ref == 1) then
          call get_excitation_degree_spin(HF_bitmask,det(1,2),degree,N_int)
          do_cycle = do_cycle .and. (degree > excitation_max)
        else if (excitation_ref == 2) then
          do k=1,N_dominant_dets_of_cfgs
            call get_excitation_degree(dominant_dets_of_cfgs(1,2,k),det(1,2),degree,N_int)
            do_cycle = do_cycle .and. (degree > excitation_beta_max)
          enddo
        endif
        if (do_cycle) cycle
      endif

      if (twice_hierarchy_max >= 0) then
        s = 0
        do k=1,N_int
          s = s + popcnt(ieor(det(k,1),det(k,2)))
        enddo
        if ( mod(s,2)>0 ) stop 'For now, hierarchy CI is defined only for an even number of electrons'
        if (excitation_ref == 1) then
          call get_excitation_degree(HF_bitmask,det(1,1),degree,N_int)
        else if (excitation_ref == 2) then
          stop 'For now, hierarchy CI is defined only for a single reference determinant'
!         do k=1,N_dominant_dets_of_cfgs
!           call get_excitation_degree(dominant_dets_of_cfgs(1,1,k),det(1,1),degree,N_int)
!         enddo
        endif
        integer :: twice_hierarchy
        twice_hierarchy = degree + s/2
        if (twice_hierarchy > twice_hierarchy_max) cycle
      endif

      Hii = diag_H_mat_elem_fock(psi_det_generators(1,1,i_generator),det,fock_diag_tmp,N_int)

      w = 0d0

      e_pert = 0.d0
      coef = 0.d0
      logical :: do_diag
      do_diag = .False.

      do istate=1,N_states
        delta_E = E0(istate) - Hii + E_shift
        alpha_h_psi = mat(istate, p1, p2)
        if (alpha_h_psi == 0.d0) cycle

        val = alpha_h_psi + alpha_h_psi
        tmp = dsqrt(delta_E * delta_E + val * val)
        if (delta_E < 0.d0) then
            tmp = -tmp
        endif

        !e_pert(istate) = alpha_h_psi * alpha_h_psi / (E0(istate) - Hii)
        e_pert(istate) = 0.5d0 * (tmp - delta_E)

        if (dabs(alpha_h_psi) > 1.d-4) then
          coef(istate) = e_pert(istate) / alpha_h_psi
        else
          coef(istate) = alpha_h_psi / delta_E
        endif
      enddo

      do_diag = sum(dabs(coef)) > 0.001d0 .and. N_states > 1

      double precision :: eigvalues(N_states+1)
      double precision :: work(1+6*(N_states+1)+2*(N_states+1)**2)
      integer :: info, k , iwork(N_states+1)

      if (do_diag) then
        double precision :: pt2_matrix(N_states+1,N_states+1)
        pt2_matrix(N_states+1,N_states+1) = Hii+E_shift
        do istate=1,N_states
          pt2_matrix(:,istate) = 0.d0
          pt2_matrix(istate,istate) = E0(istate)
          pt2_matrix(istate,N_states+1) = mat(istate,p1,p2)
          pt2_matrix(N_states+1,istate) = mat(istate,p1,p2)
        enddo

        call DSYEV( 'V', 'U', N_states+1, pt2_matrix, N_states+1, eigvalues, &
                     work, size(work), info )
        if (info /= 0) then
          print *, 'error in '//irp_here
          stop -1
        endif
        pt2_matrix = dabs(pt2_matrix)
        iwork(1:N_states+1) = maxloc(pt2_matrix,DIM=1)
        do k=1,N_states
          e_pert(k) = eigvalues(iwork(k)) - E0(k)
        enddo
      endif



!      ! Gram-Schmidt using input overlap matrix
!      do istate=1,N_states
!        do jstate=1,istate-1
!          if ( (pt2_overlap(jstate,istate) == 0.d0).or.(pt2_overlap(jstate,jstate) == 0.d0) ) cycle
!          coef(istate) = coef(istate) - pt2_overlap(jstate,istate)/pt2_overlap(jstate,jstate) * coef(jstate)
!        enddo
!      enddo

      do istate=1, N_states

        alpha_h_psi = mat(istate, p1, p2)

        pt2_data % overlap(:,istate) = pt2_data % overlap(:,istate) + coef(:) * coef(istate)
        pt2_data % variance(istate)  = pt2_data % variance(istate) + alpha_h_psi * alpha_h_psi
        pt2_data % pt2(istate)       = pt2_data % pt2(istate)      + e_pert(istate)

!!!DEBUG
!        delta_E = E0(istate) - Hii + E_shift
!        pt2_data % pt2(istate) = pt2_data % pt2(istate) + alpha_h_psi**2/delta_E
!
!        integer :: k
!        double precision :: alpha_h_psi_2,hij
!        alpha_h_psi_2 = 0.d0
!        do k = 1,N_det_selectors
!         call i_H_j(det,psi_selectors(1,1,k),N_int,hij)
!         alpha_h_psi_2 = alpha_h_psi_2 + psi_selectors_coef(k,istate) * hij
!        enddo
!        if(dabs(alpha_h_psi_2 - alpha_h_psi).gt.1.d-12)then
!         call debug_det(psi_det_generators(1,1,i_generator),N_int)
!         call debug_det(det,N_int)
!         print*,'alpha_h_psi,alpha_h_psi_2 = ',alpha_h_psi,alpha_h_psi_2
!         stop
!        endif
!!!DEBUG

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
              w = min(w, e_pert(istate) * s_weight(istate,istate)) / c0_weight(istate)
            else
              w = min(w, e_pert(istate) * s_weight(istate,istate))
            endif

        end select

        ! To force the inclusion of determinants with a positive pt2 contribution
        if (e_pert(istate) > 1d-8) then
          w = -huge(1.0)
        endif

      end do

!!!BEGIN_DEBUG
!      ! To check if the pt2 is taking determinants already in the wf
!      if (is_in_wavefunction(det(N_int,1),N_int)) then
!        logical, external :: is_in_wavefunction
!        print*, 'A determinant contributing to the pt2 is already in'
!        print*, 'the wave function:'
!        call  print_det(det(N_int,1),N_int)
!        print*,'contribution to the pt2 for the states:', e_pert(:)
!        print*,'error in the filtering in'
!        print*, 'cipsi/selection.irp.f sub:  selecte_singles_and_doubles'
!        print*, 'abort'
!        call abort
!      endif
!!!END_DEBUG

      integer(bit_kind) :: occ(N_int,2), n
      if (h0_type == 'CFG') then
        do k=1,N_int
          occ(k,1) = ieor(det(k,1),det(k,2))
          occ(k,2) = iand(det(k,1),det(k,2))
        enddo
        call configuration_to_dets_size(occ,n,elec_alpha_num,N_int)
        n = max(n,1)
        w *= dsqrt(dble(n))
      endif

      if(w <= buf%mini) then
        call add_to_selection_buffer(buf, det, w)
      end if
    end do
  end do
end

subroutine splash_pq(mask, sp, det, i_gen, N_sel, bannedOrb, banned, mat, interesting)
  use bitmasks
  implicit none
  BEGIN_DOC
! Computes the contributions A(r,s) by
! comparing the external determinant to all the internal determinants det(i).
! an applying two particles (r,s) to the mask.
  END_DOC

  integer, intent(in)            :: sp, i_gen, N_sel
  integer, intent(in)            :: interesting(0:N_sel)
  integer(bit_kind),intent(in)   :: mask(N_int, 2), det(N_int, 2, N_sel)
  logical, intent(inout)         :: bannedOrb(mo_num, 2), banned(mo_num, mo_num, 2)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)

  integer                        :: i, ii, j, k, l, h(0:2,2), p(0:4,2), nt
  integer(bit_kind)              :: perMask(N_int, 2), mobMask(N_int, 2), negMask(N_int, 2)
  integer(bit_kind)             :: phasemask(N_int,2)

  PROVIDE psi_selectors_coef_transp psi_det_sorted
  mat = 0d0

  do i=1,N_int
    negMask(i,1) = not(mask(i,1))
    negMask(i,2) = not(mask(i,2))
  end do

  do i=1, N_sel
    if (interesting(i) < 0) then
      stop 'prefetch interesting(i) and det(i)'
    endif

    mobMask(1,1) = iand(negMask(1,1), det(1,1,i))
    mobMask(1,2) = iand(negMask(1,2), det(1,2,i))
    nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))

    if(nt > 4) cycle

    do j=2,N_int
      mobMask(j,1) = iand(negMask(j,1), det(j,1,i))
      mobMask(j,2) = iand(negMask(j,2), det(j,2,i))
      nt = nt + popcnt(mobMask(j, 1)) + popcnt(mobMask(j, 2))
    end do

    if(nt > 4) cycle

    if (interesting(i) == i_gen) then
        if(sp == 3) then
          do k=1,mo_num
            do j=1,mo_num
              banned(j,k,2) = banned(k,j,1)
            enddo
          enddo
        else
          do k=1,mo_num
          do l=k+1,mo_num
            banned(l,k,1) = banned(k,l,1)
          end do
          end do
        end if
    end if

    if (interesting(i) >= i_gen) then
        call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
        call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)

        perMask(1,1) = iand(mask(1,1), not(det(1,1,i)))
        perMask(1,2) = iand(mask(1,2), not(det(1,2,i)))
        do j=2,N_int
          perMask(j,1) = iand(mask(j,1), not(det(j,1,i)))
          perMask(j,2) = iand(mask(j,2), not(det(j,2,i)))
        end do

        call bitstring_to_list_in_selection(perMask(1,1), h(1,1), h(0,1), N_int)
        call bitstring_to_list_in_selection(perMask(1,2), h(1,2), h(0,2), N_int)

        call get_mask_phase(psi_det_sorted(1,1,interesting(i)), phasemask,N_int)
        if(nt == 4) then
!          call get_d2_reference(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp(1, interesting(i)))
          call get_d2(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp(1, interesting(i)))
        else if(nt == 3) then
!          call get_d1_reference(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp(1, interesting(i)))
          call get_d1(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp(1, interesting(i)))
        else
!          call get_d0_reference(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp(1, interesting(i)))
          call get_d0(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp(1, interesting(i)))
        end if
    else if(nt == 4) then
        call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
        call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)
        call past_d2(banned, p, sp)
    else if(nt == 3) then
        call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
        call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)
        call past_d1(bannedOrb, p)
    end if
  end do

end


subroutine get_d2(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  double precision, intent(in) :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  double precision, external :: get_phase_bi, mo_two_e_integral

  integer :: i, j, k, tip, ma, mi, puti, putj
  integer :: h1, h2, p1, p2, i1, i2
  double precision :: hij, phase

  integer, parameter:: turn2d(2,3,4) = reshape((/0,0, 0,0, 0,0,  3,4, 0,0, 0,0,  2,4, 1,4, 0,0,  2,3, 1,3, 1,2 /), (/2,3,4/))
  integer, parameter :: turn2(2) = (/2, 1/)
  integer, parameter :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer :: bant
  bant = 1

  tip = p(0,1) * p(0,2)

  ma = sp
  if(p(0,1) > p(0,2)) ma = 1
  if(p(0,1) < p(0,2)) ma = 2
  mi = mod(ma, 2) + 1

  if(sp == 3) then
    if(ma == 2) bant = 2

    if(tip == 3) then
      puti = p(1, mi)
      if(bannedOrb(puti, mi)) return
      h1 = h(1, ma)
      h2 = h(2, ma)

      do i = 1, 3
        putj = p(i, ma)
        if(banned(putj,puti,bant)) cycle
        i1 = turn3(1,i)
        i2 = turn3(2,i)
        p1 = p(i1, ma)
        p2 = p(i2, ma)

        hij = mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2, p1, h1, h2)
        if (hij == 0.d0) cycle

        hij = hij * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)

        if(ma == 1) then
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, putj, puti) = mat(k, putj, puti) + coefs(k) * hij
          enddo
        else
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
          enddo
        end if
      end do
    else
      h1 = h(1,1)
      h2 = h(1,2)
      do j = 1,2
        putj = p(j, 2)
        if(bannedOrb(putj, 2)) cycle
        p2 = p(turn2(j), 2)
        do i = 1,2
          puti = p(i, 1)

          if(banned(puti,putj,bant) .or. bannedOrb(puti,1)) cycle
          p1 = p(turn2(i), 1)

          hij = mo_two_e_integral(p1, p2, h1, h2)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
            enddo
          endif
        end do
      end do
    end if

  else
    if(tip == 0) then
      h1 = h(1, ma)
      h2 = h(2, ma)
      do i=1,3
        puti = p(i, ma)
        if(bannedOrb(puti,ma)) cycle
        do j=i+1,4
          putj = p(j, ma)
          if(bannedOrb(putj,ma)) cycle
          if(banned(puti,putj,1)) cycle

          i1 = turn2d(1, i, j)
          i2 = turn2d(2, i, j)
          p1 = p(i1, ma)
          p2 = p(i2, ma)
          hij = mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2,p1, h1, h2)
          if (hij == 0.d0) cycle

          hij = hij * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, putj) = mat(k, puti, putj) +coefs(k) * hij
          enddo
        end do
      end do
    else if(tip == 3) then
      h1 = h(1, mi)
      h2 = h(1, ma)
      p1 = p(1, mi)
      do i=1,3
        puti = p(turn3(1,i), ma)
        if(bannedOrb(puti,ma)) cycle
        putj = p(turn3(2,i), ma)
        if(bannedOrb(putj,ma)) cycle
        if(banned(puti,putj,1)) cycle
        p2 = p(i, ma)

        hij = mo_two_e_integral(p1, p2, h1, h2)
        if (hij == 0.d0) cycle

        hij = hij * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2, N_int)
        if (puti < putj) then
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
          enddo
        else
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, putj, puti) = mat(k, putj, puti) + coefs(k) * hij
          enddo
        endif
      end do
    else ! tip == 4
      puti = p(1, sp)
      putj = p(2, sp)
      if(.not. banned(puti,putj,1)) then
        p1 = p(1, mi)
        p2 = p(2, mi)
        h1 = h(1, mi)
        h2 = h(2, mi)
        hij = (mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2,p1, h1, h2))
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
          enddo
        end if
      end if
    end if
  end if
end


subroutine get_d1(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  use bitmasks
  implicit none

  integer(bit_kind), intent(in)  :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in)  :: phasemask(N_int,2)
  logical, intent(in)            :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind)              :: det(N_int, 2)
  double precision, intent(in)   :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in)            :: h(0:2,2), p(0:4,2), sp
  double precision, external     :: get_phase_bi, mo_two_e_integral
  logical                        :: ok

  logical, allocatable           :: lbanned(:,:)
  integer                        :: puti, putj, ma, mi, s1, s2, i, i1, i2, j
  integer                        :: hfix, pfix, h1, h2, p1, p2, ib, k, l

  integer, parameter             :: turn2(2) = (/2,1/)
  integer, parameter             :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer                        :: bant
  double precision, allocatable :: hij_cache(:,:)
  double precision               :: hij, tmp_row(N_states, mo_num), tmp_row2(N_states, mo_num)
  PROVIDE mo_integrals_map N_int

  allocate (lbanned(mo_num, 2))
  allocate (hij_cache(mo_num,2))
  lbanned = bannedOrb

  do i=1, p(0,1)
    lbanned(p(i,1), 1) = .true.
  end do
  do i=1, p(0,2)
    lbanned(p(i,2), 2) = .true.
  end do

  ma = 1
  if(p(0,2) >= 2) ma = 2
  mi = turn2(ma)

  bant = 1

  if(sp == 3) then
    !move MA
    if(ma == 2) bant = 2
    puti = p(1,mi)
    hfix = h(1,ma)
    p1 = p(1,ma)
    p2 = p(2,ma)
    if(.not. bannedOrb(puti, mi)) then
      call get_mo_two_e_integrals(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map)
      call get_mo_two_e_integrals(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map)
      tmp_row = 0d0
      do putj=1, hfix-1
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hij = hij_cache(putj,1) - hij_cache(putj,2)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row(k,putj) = tmp_row(k,putj) + hij * coefs(k)
          enddo
        endif
      end do
      do putj=hfix+1, mo_num
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hij = hij_cache(putj,2) - hij_cache(putj,1)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row(k,putj) = tmp_row(k,putj) + hij * coefs(k)
          enddo
        endif
      end do

      if(ma == 1) then
        mat(1:N_states,1:mo_num,puti) = mat(1:N_states,1:mo_num,puti) + tmp_row(1:N_states,1:mo_num)
      else
        do l=1,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k,puti,l) = mat(k,puti,l) + tmp_row(k,l)
          enddo
        enddo
      end if
    end if

    !MOVE MI
    pfix = p(1,mi)
    tmp_row = 0d0
    tmp_row2 = 0d0
    call get_mo_two_e_integrals(hfix,pfix,p1,mo_num,hij_cache(1,1),mo_integrals_map)
    call get_mo_two_e_integrals(hfix,pfix,p2,mo_num,hij_cache(1,2),mo_integrals_map)
    putj = p1
    do puti=1,mo_num !HOT
      if(lbanned(puti,mi)) cycle
      !p1 fixed
      putj = p1
      if(.not. banned(putj,puti,bant)) then
        hij = hij_cache(puti,2)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row(k,puti) = tmp_row(k,puti) + hij * coefs(k)
          enddo
        endif
      end if
!    enddo
!
      putj = p2
!    do puti=1,mo_num !HOT
      if(.not. banned(putj,puti,bant)) then
        hij = hij_cache(puti,1)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
          do k=1,N_states
            tmp_row2(k,puti) = tmp_row2(k,puti) + hij * coefs(k)
          enddo
        endif
      end if
    end do

    if(mi == 1) then
      mat(:,:,p1) = mat(:,:,p1) + tmp_row(:,:)
      mat(:,:,p2) = mat(:,:,p2) + tmp_row2(:,:)
    else
      do l=1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k,p1,l) = mat(k,p1,l) + tmp_row(k,l)
          mat(k,p2,l) = mat(k,p2,l) + tmp_row2(k,l)
        enddo
      enddo
    end if

  else  ! sp /= 3

    if(p(0,ma) == 3) then
      do i=1,3
        hfix = h(1,ma)
        puti = p(i, ma)
        p1 = p(turn3(1,i), ma)
        p2 = p(turn3(2,i), ma)
        call get_mo_two_e_integrals(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map)
        call get_mo_two_e_integrals(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map)
        tmp_row = 0d0
        do putj=1,hfix-1
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache(putj,1) - hij_cache(putj,2)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
            tmp_row(:,putj) = tmp_row(:,putj) + hij * coefs(:)
          endif
        end do
        do putj=hfix+1,mo_num
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache(putj,2) - hij_cache(putj,1)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
            tmp_row(:,putj) = tmp_row(:,putj) + hij * coefs(:)
          endif
        end do

        mat(:, :puti-1, puti) = mat(:, :puti-1, puti) + tmp_row(:,:puti-1)
        do l=puti,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, l) = mat(k, puti,l) + tmp_row(k,l)
          enddo
        enddo
      end do
    else
      hfix = h(1,mi)
      pfix = p(1,mi)
      p1 = p(1,ma)
      p2 = p(2,ma)
      tmp_row = 0d0
      tmp_row2 = 0d0
      call get_mo_two_e_integrals(hfix,p1,pfix,mo_num,hij_cache(1,1),mo_integrals_map)
      call get_mo_two_e_integrals(hfix,p2,pfix,mo_num,hij_cache(1,2),mo_integrals_map)
      putj = p2
      do puti=1,mo_num
        if(lbanned(puti,ma)) cycle
        putj = p2
        if(.not. banned(puti,putj,1)) then
          hij = hij_cache(puti,1)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              tmp_row(k,puti) = tmp_row(k,puti) + hij * coefs(k)
            enddo
          endif
        end if

        putj = p1
        if(.not. banned(puti,putj,1)) then
          hij = hij_cache(puti,2)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
            do k=1,N_states
              tmp_row2(k,puti) = tmp_row2(k,puti) + hij * coefs(k)
            enddo
          endif
        end if
      end do
      mat(:,:p2-1,p2) = mat(:,:p2-1,p2) + tmp_row(:,:p2-1)
      do l=p2,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k,p2,l) = mat(k,p2,l) + tmp_row(k,l)
        enddo
      enddo
      mat(:,:p1-1,p1) = mat(:,:p1-1,p1) + tmp_row2(:,:p1-1)
      do l=p1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k,p1,l) = mat(k,p1,l) + tmp_row2(k,l)
        enddo
      enddo
    end if
  end if
  deallocate(lbanned,hij_cache)

 !! MONO
    if(sp == 3) then
      s1 = 1
      s2 = 2
    else
      s1 = sp
      s2 = sp
    end if

    do i1=1,p(0,s1)
      ib = 1
      if(s1 == s2) ib = i1+1
      do i2=ib,p(0,s2)
        p1 = p(i1,s1)
        p2 = p(i2,s2)
        if(bannedOrb(p1, s1) .or. bannedOrb(p2, s2) .or. banned(p1, p2, 1)) cycle
        call apply_particles(mask, s1, p1, s2, p2, det, ok, N_int)
        call i_h_j(gen, det, N_int, hij)
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k, p1, p2) = mat(k, p1, p2) + coefs(k) * hij
        enddo
      end do
    end do
end




subroutine get_d0(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: gen(N_int, 2), mask(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind) :: det(N_int, 2)
  double precision, intent(in) :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  integer :: i, j, k, s, h1, h2, p1, p2, puti, putj
  double precision :: hij, phase
  double precision, external :: get_phase_bi, mo_two_e_integral
  logical :: ok

  integer, parameter :: bant=1
  double precision, allocatable :: hij_cache1(:), hij_cache2(:)
  allocate (hij_cache1(mo_num),hij_cache2(mo_num))


  if(sp == 3) then ! AB
    h1 = p(1,1)
    h2 = p(1,2)
    do p1=1, mo_num
      if(bannedOrb(p1, 1)) cycle
      call get_mo_two_e_integrals(p1,h2,h1,mo_num,hij_cache1,mo_integrals_map)
      do p2=1, mo_num
        if(bannedOrb(p2,2)) cycle
        if(banned(p1, p2, bant)) cycle ! rentable?
        if(p1 == h1 .or. p2 == h2) then
          call apply_particles(mask, 1,p1,2,p2, det, ok, N_int)
          call i_h_j(gen, det, N_int, hij)
        else
          phase = get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
          hij = hij_cache1(p2) * phase
        end if
        if (hij == 0.d0) cycle
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k, p1, p2) = mat(k, p1, p2) + coefs(k) * hij  ! HOTSPOT
        enddo
      end do
    end do

  else ! AA BB
    p1 = p(1,sp)
    p2 = p(2,sp)
    do puti=1, mo_num
      if(bannedOrb(puti, sp)) cycle
      call get_mo_two_e_integrals(puti,p2,p1,mo_num,hij_cache1,mo_integrals_map)
      call get_mo_two_e_integrals(puti,p1,p2,mo_num,hij_cache2,mo_integrals_map)
      do putj=puti+1, mo_num
        if(bannedOrb(putj, sp)) cycle
        if(banned(puti, putj, bant)) cycle ! rentable?
        if(puti == p1 .or. putj == p2 .or. puti == p2 .or. putj == p1) then
          call apply_particles(mask, sp,puti,sp,putj, det, ok, N_int)
          call i_h_j(gen, det, N_int, hij)
          if (hij == 0.d0) cycle
        else
          hij = (mo_two_e_integral(p1, p2, puti, putj) -  mo_two_e_integral(p2, p1, puti, putj))
          if (hij == 0.d0) cycle
          hij = hij * get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
        end if
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
        enddo
      end do
    end do
  end if

  deallocate(hij_cache1,hij_cache2)
end


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
end


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
end



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
end


subroutine bitstring_to_list_in_selection( string, list, n_elements, Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Gives the indices(+1) of the bits set to 1 in the bit string
  END_DOC
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

end
!




! OLD unoptimized routines for debugging
! ======================================

subroutine get_d0_reference(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: gen(N_int, 2), mask(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind) :: det(N_int, 2)
  double precision, intent(in) :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  integer :: i, j, s, h1, h2, p1, p2, puti, putj
  double precision :: hij, phase
  double precision, external :: get_phase_bi, mo_two_e_integral
  logical :: ok

  integer :: bant
  bant = 1


  if(sp == 3) then ! AB
    h1 = p(1,1)
    h2 = p(1,2)
    do p1=1, mo_num
      if(bannedOrb(p1, 1)) cycle
      do p2=1, mo_num
        if(bannedOrb(p2,2)) cycle
        if(banned(p1, p2, bant)) cycle ! rentable?
        if(p1 == h1 .or. p2 == h2) then
          call apply_particles(mask, 1,p1,2,p2, det, ok, N_int)
          call i_h_j(gen, det, N_int, hij)
        else
          phase = get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
          hij = mo_two_e_integral(p1, p2, h1, h2) * phase
        end if
        mat(:, p1, p2) = mat(:, p1, p2) + coefs(:) * hij
      end do
    end do
  else ! AA BB
    p1 = p(1,sp)
    p2 = p(2,sp)
    do puti=1, mo_num
      if(bannedOrb(puti, sp)) cycle
      do putj=puti+1, mo_num
        if(bannedOrb(putj, sp)) cycle
        if(banned(puti, putj, bant)) cycle ! rentable?
        if(puti == p1 .or. putj == p2 .or. puti == p2 .or. putj == p1) then
          call apply_particles(mask, sp,puti,sp,putj, det, ok, N_int)
          call i_h_j(gen, det, N_int, hij)
        else
          hij = (mo_two_e_integral(p1, p2, puti, putj) -  mo_two_e_integral(p2, p1, puti, putj))* get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
        end if
        mat(:, puti, putj) = mat(:, puti, putj) + coefs(:) * hij
      end do
    end do
  end if
end

subroutine get_d1_reference(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  use bitmasks
  implicit none

  integer(bit_kind), intent(in)  :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in)  :: phasemask(N_int,2)
  logical, intent(in)            :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind)              :: det(N_int, 2)
  double precision, intent(in)   :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in)            :: h(0:2,2), p(0:4,2), sp
  double precision               :: hij, tmp_row(N_states, mo_num), tmp_row2(N_states, mo_num)
  double precision, external     :: get_phase_bi, mo_two_e_integral
  logical                        :: ok

  logical, allocatable           :: lbanned(:,:)
  integer                        :: puti, putj, ma, mi, s1, s2, i, i1, i2, j
  integer                        :: hfix, pfix, h1, h2, p1, p2, ib

  integer, parameter             :: turn2(2) = (/2,1/)
  integer, parameter             :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer                        :: bant


  allocate (lbanned(mo_num, 2))
  lbanned = bannedOrb

  do i=1, p(0,1)
    lbanned(p(i,1), 1) = .true.
  end do
  do i=1, p(0,2)
    lbanned(p(i,2), 2) = .true.
  end do

  ma = 1
  if(p(0,2) >= 2) ma = 2
  mi = turn2(ma)

  bant = 1

  if(sp == 3) then
    !move MA
    if(ma == 2) bant = 2
    puti = p(1,mi)
    hfix = h(1,ma)
    p1 = p(1,ma)
    p2 = p(2,ma)
    if(.not. bannedOrb(puti, mi)) then
      tmp_row = 0d0
      do putj=1, hfix-1
        if(lbanned(putj, ma) .or. banned(putj, puti,bant)) cycle
        hij = (mo_two_e_integral(p1, p2, putj, hfix)-mo_two_e_integral(p2,p1,putj,hfix)) * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
        tmp_row(1:N_states,putj) = tmp_row(1:N_states,putj) + hij * coefs(1:N_states)
      end do
      do putj=hfix+1, mo_num
        if(lbanned(putj, ma) .or. banned(putj, puti,bant)) cycle
        hij = (mo_two_e_integral(p1, p2, hfix, putj)-mo_two_e_integral(p2,p1,hfix,putj)) * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
        tmp_row(1:N_states,putj) = tmp_row(1:N_states,putj) + hij * coefs(1:N_states)
      end do

      if(ma == 1) then
        mat(1:N_states,1:mo_num,puti) = mat(1:N_states,1:mo_num,puti) + tmp_row(1:N_states,1:mo_num)
      else
        mat(1:N_states,puti,1:mo_num) = mat(1:N_states,puti,1:mo_num) + tmp_row(1:N_states,1:mo_num)
      end if
    end if

    !MOVE MI
    pfix = p(1,mi)
    tmp_row = 0d0
    tmp_row2 = 0d0
    do puti=1,mo_num
      if(lbanned(puti,mi)) cycle
      !p1 fixed
      putj = p1
      if(.not. banned(putj,puti,bant)) then
        hij = mo_two_e_integral(p2,pfix,hfix,puti) * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
        tmp_row(:,puti) = tmp_row(:,puti) + hij * coefs(:)
      end if

      putj = p2
      if(.not. banned(putj,puti,bant)) then
        hij = mo_two_e_integral(p1,pfix,hfix,puti) * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
        tmp_row2(:,puti) = tmp_row2(:,puti) + hij * coefs(:)
      end if
    end do

    if(mi == 1) then
      mat(:,:,p1) = mat(:,:,p1) + tmp_row(:,:)
      mat(:,:,p2) = mat(:,:,p2) + tmp_row2(:,:)
    else
      mat(:,p1,:) = mat(:,p1,:) + tmp_row(:,:)
      mat(:,p2,:) = mat(:,p2,:) + tmp_row2(:,:)
    end if
  else
    if(p(0,ma) == 3) then
      do i=1,3
        hfix = h(1,ma)
        puti = p(i, ma)
        p1 = p(turn3(1,i), ma)
        p2 = p(turn3(2,i), ma)
        tmp_row = 0d0
        do putj=1,hfix-1
          if(lbanned(putj,ma) .or. banned(puti,putj,1)) cycle
          hij = (mo_two_e_integral(p1, p2, putj, hfix)-mo_two_e_integral(p2,p1,putj,hfix)) * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
          tmp_row(:,putj) = tmp_row(:,putj) + hij * coefs(:)
        end do
        do putj=hfix+1,mo_num
          if(lbanned(putj,ma) .or. banned(puti,putj,1)) cycle
          hij = (mo_two_e_integral(p1, p2, hfix, putj)-mo_two_e_integral(p2,p1,hfix,putj)) * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
          tmp_row(:,putj) = tmp_row(:,putj) + hij * coefs(:)
        end do

        mat(:, :puti-1, puti) = mat(:, :puti-1, puti) + tmp_row(:,:puti-1)
        mat(:, puti, puti:) = mat(:, puti, puti:) + tmp_row(:,puti:)
      end do
    else
      hfix = h(1,mi)
      pfix = p(1,mi)
      p1 = p(1,ma)
      p2 = p(2,ma)
      tmp_row = 0d0
      tmp_row2 = 0d0
      do puti=1,mo_num
        if(lbanned(puti,ma)) cycle
        putj = p2
        if(.not. banned(puti,putj,1)) then
          hij = mo_two_e_integral(pfix, p1, hfix, puti) * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
          tmp_row(:,puti) = tmp_row(:,puti) + hij * coefs(:)
        end if

        putj = p1
        if(.not. banned(puti,putj,1)) then
          hij = mo_two_e_integral(pfix, p2, hfix, puti) * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
          tmp_row2(:,puti) = tmp_row2(:,puti) + hij * coefs(:)
        end if
      end do
      mat(:,:p2-1,p2) = mat(:,:p2-1,p2) + tmp_row(:,:p2-1)
      mat(:,p2,p2:) = mat(:,p2,p2:) + tmp_row(:,p2:)
      mat(:,:p1-1,p1) = mat(:,:p1-1,p1) + tmp_row2(:,:p1-1)
      mat(:,p1,p1:) = mat(:,p1,p1:) + tmp_row2(:,p1:)
    end if
  end if
  deallocate(lbanned)

 !! MONO
    if(sp == 3) then
      s1 = 1
      s2 = 2
    else
      s1 = sp
      s2 = sp
    end if

    do i1=1,p(0,s1)
      ib = 1
      if(s1 == s2) ib = i1+1
      do i2=ib,p(0,s2)
        p1 = p(i1,s1)
        p2 = p(i2,s2)
        if(bannedOrb(p1, s1) .or. bannedOrb(p2, s2) .or. banned(p1, p2, 1)) cycle
        call apply_particles(mask, s1, p1, s2, p2, det, ok, N_int)
        call i_h_j(gen, det, N_int, hij)
        mat(:, p1, p2) = mat(:, p1, p2) + coefs(:) * hij
      end do
    end do
end

subroutine get_d2_reference(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(2,N_int)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  double precision, intent(in) :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  double precision, external :: get_phase_bi, mo_two_e_integral

  integer :: i, j, tip, ma, mi, puti, putj
  integer :: h1, h2, p1, p2, i1, i2
  double precision :: hij, phase

  integer, parameter:: turn2d(2,3,4) = reshape((/0,0, 0,0, 0,0,  3,4, 0,0, 0,0,  2,4, 1,4, 0,0,  2,3, 1,3, 1,2 /), (/2,3,4/))
  integer, parameter :: turn2(2) = (/2, 1/)
  integer, parameter :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer :: bant
  bant = 1

  tip = p(0,1) * p(0,2)

  ma = sp
  if(p(0,1) > p(0,2)) ma = 1
  if(p(0,1) < p(0,2)) ma = 2
  mi = mod(ma, 2) + 1

  if(sp == 3) then
    if(ma == 2) bant = 2

    if(tip == 3) then
      puti = p(1, mi)
      do i = 1, 3
        putj = p(i, ma)
        if(banned(putj,puti,bant)) cycle
        i1 = turn3(1,i)
        i2 = turn3(2,i)
        p1 = p(i1, ma)
        p2 = p(i2, ma)
        h1 = h(1, ma)
        h2 = h(2, ma)

        hij = (mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2,p1, h1, h2)) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
        if(ma == 1) then
          mat(:, putj, puti) = mat(:, putj, puti) + coefs(:) * hij
        else
          mat(:, puti, putj) = mat(:, puti, putj) + coefs(:) * hij
        end if
      end do
    else
      h1 = h(1,1)
      h2 = h(1,2)
      do j = 1,2
        putj = p(j, 2)
        p2 = p(turn2(j), 2)
        do i = 1,2
          puti = p(i, 1)

          if(banned(puti,putj,bant)) cycle
          p1 = p(turn2(i), 1)

          hij = mo_two_e_integral(p1, p2, h1, h2) * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2,N_int)
          mat(:, puti, putj) = mat(:, puti, putj) + coefs(:) * hij
        end do
      end do
    end if

  else
    if(tip == 0) then
      h1 = h(1, ma)
      h2 = h(2, ma)
      do i=1,3
      puti = p(i, ma)
      do j=i+1,4
        putj = p(j, ma)
        if(banned(puti,putj,1)) cycle

        i1 = turn2d(1, i, j)
        i2 = turn2d(2, i, j)
        p1 = p(i1, ma)
        p2 = p(i2, ma)
        hij = (mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2,p1, h1, h2)) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2,N_int)
        mat(:, puti, putj) = mat(:, puti, putj) + coefs(:) * hij
      end do
      end do
    else if(tip == 3) then
      h1 = h(1, mi)
      h2 = h(1, ma)
      p1 = p(1, mi)
      do i=1,3
        puti = p(turn3(1,i), ma)
        putj = p(turn3(2,i), ma)
        if(banned(puti,putj,1)) cycle
        p2 = p(i, ma)

        hij = mo_two_e_integral(p1, p2, h1, h2) * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2,N_int)
        mat(:, min(puti, putj), max(puti, putj)) = mat(:, min(puti, putj), max(puti, putj)) + coefs(:) * hij
      end do
    else ! tip == 4
      puti = p(1, sp)
      putj = p(2, sp)
      if(.not. banned(puti,putj,1)) then
        p1 = p(1, mi)
        p2 = p(2, mi)
        h1 = h(1, mi)
        h2 = h(2, mi)
        hij = (mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2,p1, h1, h2)) * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2,N_int)
        mat(:, puti, putj) = mat(:, puti, putj) + coefs(:) * hij
      end if
    end if
  end if
end


