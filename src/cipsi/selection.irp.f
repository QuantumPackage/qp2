use bitmasks

BEGIN_PROVIDER [logical , pert_2rdm ]
 implicit none
 pert_2rdm = .False.
END_PROVIDER 

BEGIN_PROVIDER [ double precision, pt2_match_weight, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Weights adjusted along the selection to make the PT2 contributions
 ! of each state coincide.
 END_DOC
 pt2_match_weight(:) = 1.d0
END_PROVIDER

BEGIN_PROVIDER [ double precision, variance_match_weight, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Weights adjusted along the selection to make the variances 
 ! of each state coincide.
 END_DOC
 variance_match_weight(:) = 1.d0
END_PROVIDER

subroutine update_pt2_and_variance_weights(pt2, variance, norm, N_st)
  implicit none
  BEGIN_DOC
! Updates the rPT2- and Variance- matching weights.
  END_DOC
  integer, intent(in)          :: N_st
  double precision, intent(in) :: pt2(N_st)
  double precision, intent(in) :: variance(N_st)
  double precision, intent(in) :: norm(N_st)

  double precision :: avg, rpt2(N_st), element, dt, x
  integer          :: k
  integer, save    :: i_iter=0
  integer, parameter :: i_itermax = 3
  double precision, allocatable, save :: memo_variance(:,:), memo_pt2(:,:)

  if (i_iter == 0) then
    allocate(memo_variance(N_st,i_itermax), memo_pt2(N_st,i_itermax))
    memo_pt2(:,:) = 1.d0
    memo_variance(:,:) = 1.d0
  endif

  i_iter = i_iter+1
  if (i_iter > i_itermax) then
    i_iter = 1
  endif

  dt = 4.d0 

  do k=1,N_st
    rpt2(k) = pt2(k)/(1.d0 + norm(k))                                                     
  enddo

  avg = sum(rpt2(1:N_st)) / dble(N_st)
  do k=1,N_st
    element = exp(dt*(rpt2(k)/avg -1.d0))
    element = min(1.5d0 , element)
    element = max(0.5d0 , element)
    memo_pt2(k,i_iter) = element
    pt2_match_weight(k) = product(memo_pt2(k,:))
  enddo

  avg = sum(variance(1:N_st)) / dble(N_st)
  do k=1,N_st
    element = exp(dt*(variance(k)/avg -1.d0))
    element = min(1.5d0 , element)
    element = max(0.5d0 , element)
    memo_variance(k,i_iter) = element
    variance_match_weight(k) = product(memo_variance(k,:))
  enddo

  print *, '# PT2 weight ', real(pt2_match_weight(:),4)
  print *, '# var weight ', real(variance_match_weight(:),4)
  SOFT_TOUCH pt2_match_weight variance_match_weight
end


BEGIN_PROVIDER [ double precision, selection_weight, (N_states) ]
   implicit none
   BEGIN_DOC
   ! Weights used in the selection criterion
   END_DOC
   select case (weight_selection)

     case (0)
      print *,  'Using input weights in selection'
      selection_weight(1:N_states) = state_average_weight(1:N_states)

     case (1)
      print *,  'Using 1/c_max^2 weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) 

     case (2)
      print *,  'Using pt2-matching weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * pt2_match_weight(1:N_states)

     case (3)
      print *,  'Using variance-matching weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * variance_match_weight(1:N_states)

     case (4)
      print *,  'Using variance- and pt2-matching weights in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * variance_match_weight(1:N_states) * pt2_match_weight(1:N_states)

     case (5)
      print *,  'Using variance-matching weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * variance_match_weight(1:N_states)

    end select

END_PROVIDER


subroutine get_mask_phase(det1, pm, Nint)
  use bitmasks
  implicit none
  integer, intent(in) :: Nint
  integer(bit_kind), intent(in) :: det1(Nint,2)
  integer(bit_kind), intent(out) :: pm(Nint,2)
  integer(bit_kind) :: tmp1, tmp2
  integer :: i
  pm(1:Nint,1:2) = det1(1:Nint,1:2)
  tmp1 = 0_8
  tmp2 = 0_8
  do i=1,Nint
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 1))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 1))
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

end subroutine


subroutine select_connected(i_generator,E0,pt2,variance,norm,b,subset,csubset)
  use bitmasks
  use selection_types
  implicit none
  integer, intent(in)            :: i_generator, subset, csubset
  type(selection_buffer), intent(inout) :: b
  double precision, intent(inout)  :: pt2(N_states)
  double precision, intent(inout)  :: variance(N_states)
  double precision, intent(inout)  :: norm(N_states)
  integer :: k,l
  double precision, intent(in)   :: E0(N_states)

  integer(bit_kind)              :: hole_mask(N_int,2), particle_mask(N_int,2)

  double precision, allocatable  :: fock_diag_tmp(:,:)

  allocate(fock_diag_tmp(2,mo_num+1))

  call build_fock_tmp(fock_diag_tmp,psi_det_generators(1,1,i_generator),N_int)

  do l=1,N_generators_bitmask
    do k=1,N_int
      hole_mask(k,1) = iand(generators_bitmask(k,1,s_hole,l), psi_det_generators(k,1,i_generator))
      hole_mask(k,2) = iand(generators_bitmask(k,2,s_hole,l), psi_det_generators(k,2,i_generator))
      particle_mask(k,1) = iand(generators_bitmask(k,1,s_part,l), not(psi_det_generators(k,1,i_generator)) )
      particle_mask(k,2) = iand(generators_bitmask(k,2,s_part,l), not(psi_det_generators(k,2,i_generator)) )
    enddo
    call select_singles_and_doubles(i_generator,hole_mask,particle_mask,fock_diag_tmp,E0,pt2,variance,norm,b,subset,csubset)
  enddo
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


subroutine select_singles_and_doubles(i_generator,hole_mask,particle_mask,fock_diag_tmp,E0,pt2,variance,norm,buf,subset,csubset)
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
  double precision, intent(inout) :: pt2(N_states)
  double precision, intent(inout)  :: variance(N_states)
  double precision, intent(inout)  :: norm(N_states)
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
  double precision, allocatable  :: coef_fullminilist(:,:)


  double precision, allocatable   :: mat(:,:,:)

  logical :: monoAdo, monoBdo
  integer :: maskInd

  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_order psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order psi_selectors_coef_transp

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
  allocate (indices(N_det),                                          &
      exc_degree(max(N_det_alpha_unique,N_det_beta_unique)))
  
  k=1
  do i=1,N_det_alpha_unique
    call get_excitation_degree_spin(psi_det_alpha_unique(1,i),       &
        psi_det_generators(1,1,i_generator), exc_degree(i), N_int)
  enddo
  
  do j=1,N_det_beta_unique
    call get_excitation_degree_spin(psi_det_beta_unique(1,j),        &
        psi_det_generators(1,2,i_generator), nt, N_int)
    if (nt > 2) cycle
    do l_a=psi_bilinear_matrix_columns_loc(j), psi_bilinear_matrix_columns_loc(j+1)-1
      i = psi_bilinear_matrix_rows(l_a)
      if (nt + exc_degree(i) <= 4) then
        idx = psi_det_sorted_order(psi_bilinear_matrix_order(l_a))
        if (psi_average_norm_contrib_sorted(idx) < 1.d-12) cycle
        indices(k) = idx
        k=k+1
      endif
    enddo
  enddo
  
  do i=1,N_det_beta_unique
    call get_excitation_degree_spin(psi_det_beta_unique(1,i),        &
        psi_det_generators(1,2,i_generator), exc_degree(i), N_int)
  enddo
  
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
        if (psi_average_norm_contrib_sorted(idx) < 1.d-12) cycle
        indices(k) = idx
        k=k+1
      endif
    enddo
  enddo
  
  deallocate(exc_degree)
  nmax=k-1

  allocate(iorder(nmax))
  do i=1,nmax
    iorder(i) = i
  enddo
  call isort(indices,iorder,nmax)
  deallocate(iorder)
  
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
        select case (N_int)
          case (1)
            mobMask(1,1) = iand(negMask(1,1), psi_det_sorted(1,1,preinteresting(ii)))
            mobMask(1,2) = iand(negMask(1,2), psi_det_sorted(1,2,preinteresting(ii)))
            nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
          case (2)
            mobMask(1:2,1) = iand(negMask(1:2,1), psi_det_sorted(1:2,1,preinteresting(ii)))
            mobMask(1:2,2) = iand(negMask(1:2,2), psi_det_sorted(1:2,2,preinteresting(ii)))
            nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2)) +     &
                popcnt(mobMask(2, 1)) + popcnt(mobMask(2, 2))
          case (3)
            mobMask(1:3,1) = iand(negMask(1:3,1), psi_det_sorted(1:3,1,preinteresting(ii)))
            mobMask(1:3,2) = iand(negMask(1:3,2), psi_det_sorted(1:3,2,preinteresting(ii)))
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
            mobMask(1:4,1) = iand(negMask(1:4,1), psi_det_sorted(1:4,1,preinteresting(ii)))
            mobMask(1:4,2) = iand(negMask(1:4,2), psi_det_sorted(1:4,2,preinteresting(ii)))
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
            mobMask(1:N_int,1) = iand(negMask(1:N_int,1), psi_det_sorted(1:N_int,1,preinteresting(ii)))
            mobMask(1:N_int,2) = iand(negMask(1:N_int,2), psi_det_sorted(1:N_int,2,preinteresting(ii)))
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
          i = preinteresting(ii)
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
      if(pert_2rdm)then
       allocate(coef_fullminilist(fullinteresting(0),N_states)) 
       do i=1,fullinteresting(0)
         coef_fullminilist(i,:) = psi_coef_sorted(fullinteresting(i),:)
       enddo
      endif
      do i=1,fullinteresting(0)
        fullminilist(1:N_int,1:2,i) = psi_det_sorted(1:N_int,1:2,fullinteresting(i))
      enddo
      
      do i=1,interesting(0)
        minilist(1:N_int,1:2,i) = psi_det_sorted(1:N_int,1:2,interesting(i))
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
          banned = .false.
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
            
            if(.not.pert_2rdm)then
             call fill_buffer_double(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2, variance, norm, mat, buf)
            else 
             call fill_buffer_double_rdm(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2, variance, norm, mat, buf,fullminilist, coef_fullminilist, fullinteresting(0))
            endif
          end if
        enddo
        if(s1 /= s2) monoBdo = .false.
      enddo
      deallocate(fullminilist,minilist)
      if(pert_2rdm)then
       deallocate(coef_fullminilist)
      endif
    enddo
  enddo
  deallocate(preinteresting, prefullinteresting, interesting, fullinteresting)
  deallocate(banned, bannedOrb,mat)
end subroutine



subroutine fill_buffer_double(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2, variance, norm, mat, buf)
  use bitmasks
  use selection_types
  implicit none

  integer, intent(in) :: i_generator, sp, h1, h2
  double precision, intent(in) :: mat(N_states, mo_num, mo_num)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num)
  double precision, intent(in)           :: fock_diag_tmp(mo_num)
  double precision, intent(in)    :: E0(N_states)
  double precision, intent(inout) :: pt2(N_states)
  double precision, intent(inout) :: variance(N_states)
  double precision, intent(inout) :: norm(N_states)
  type(selection_buffer), intent(inout) :: buf
  logical :: ok
  integer :: s1, s2, p1, p2, ib, j, istate
  integer(bit_kind) :: mask(N_int, 2), det(N_int, 2)
  double precision :: e_pert, delta_E, val, Hii, sum_e_pert, tmp, alpha_h_psi, coef
  double precision, external :: diag_H_mat_elem_fock
  double precision :: E_shift

  logical, external :: detEq
  double precision, allocatable :: values(:)
  integer, allocatable          :: keys(:,:)
  integer                       :: nkeys
  

  if(sp == 3) then
    s1 = 1
    s2 = 2
  else
    s1 = sp
    s2 = sp
  end if
  call apply_holes(psi_det_generators(1,1,i_generator), s1, h1, s2, h2, mask, ok, N_int)
  E_shift = 0.d0

  if (h0_type == 'SOP') then
    j = det_to_occ_pattern(i_generator)
    E_shift = psi_det_Hii(i_generator) - psi_occ_pattern_Hii(j)
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


      if( sum(abs(mat(1:N_states, p1, p2))) == 0d0) cycle
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


      Hii = diag_H_mat_elem_fock(psi_det_generators(1,1,i_generator),det,fock_diag_tmp,N_int)

      sum_e_pert = 0d0

      do istate=1,N_states
        delta_E = E0(istate) - Hii + E_shift
        alpha_h_psi = mat(istate, p1, p2)
        val = alpha_h_psi + alpha_h_psi
        tmp = dsqrt(delta_E * delta_E + val * val)
        if (delta_E < 0.d0) then
            tmp = -tmp
        endif
        e_pert = 0.5d0 * (tmp - delta_E)
        coef = e_pert / alpha_h_psi
        pt2(istate) = pt2(istate) + e_pert
        variance(istate) = variance(istate) + alpha_h_psi * alpha_h_psi
        norm(istate) = norm(istate) + coef * coef

        if (weight_selection /= 5) then
          ! Energy selection
          sum_e_pert = sum_e_pert + e_pert * selection_weight(istate)
        else
          ! Variance selection
          sum_e_pert = sum_e_pert - alpha_h_psi * alpha_h_psi * selection_weight(istate)
        endif
      end do
      if(pseudo_sym)then
       if(dabs(mat(1, p1, p2)).lt.thresh_sym)then 
        sum_e_pert = 10.d0
       endif
      endif

      if(sum_e_pert <= buf%mini) then
        call add_to_selection_buffer(buf, det, sum_e_pert)
      end if
    end do
  end do
end


subroutine fill_buffer_double_rdm(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2, variance, norm, mat, buf, psi_det_connection, psi_coef_connection, n_det_connection)
  use bitmasks
  use selection_types
  implicit none
  
  integer, intent(in)           :: n_det_connection
  double precision, intent(in)  :: psi_coef_connection(n_det_connection,N_states)
  integer(bit_kind), intent(in) :: psi_det_connection(N_int,2,n_det_connection)
  integer, intent(in) :: i_generator, sp, h1, h2
  double precision, intent(in) :: mat(N_states, mo_num, mo_num)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num)
  double precision, intent(in)           :: fock_diag_tmp(mo_num)
  double precision, intent(in)    :: E0(N_states)
  double precision, intent(inout) :: pt2(N_states)
  double precision, intent(inout) :: variance(N_states)
  double precision, intent(inout) :: norm(N_states)
  type(selection_buffer), intent(inout) :: buf
  logical :: ok
  integer :: s1, s2, p1, p2, ib, j, istate
  integer(bit_kind) :: mask(N_int, 2), det(N_int, 2)
  double precision :: e_pert, delta_E, val, Hii, sum_e_pert, tmp, alpha_h_psi, coef(N_states)
  double precision, external :: diag_H_mat_elem_fock
  double precision :: E_shift

  logical, external :: detEq
  double precision, allocatable :: values(:)
  integer, allocatable          :: keys(:,:)
  integer                       :: nkeys
  integer :: sze_buff
  sze_buff = 5 * mo_num ** 2 
  allocate(keys(4,sze_buff),values(sze_buff))
  nkeys = 0
  if(sp == 3) then
    s1 = 1
    s2 = 2
  else
    s1 = sp
    s2 = sp
  end if
  call apply_holes(psi_det_generators(1,1,i_generator), s1, h1, s2, h2, mask, ok, N_int)
  E_shift = 0.d0

  if (h0_type == 'SOP') then
    j = det_to_occ_pattern(i_generator)
    E_shift = psi_det_Hii(i_generator) - psi_occ_pattern_Hii(j)
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


      if( sum(abs(mat(1:N_states, p1, p2))) == 0d0) cycle
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


      Hii = diag_H_mat_elem_fock(psi_det_generators(1,1,i_generator),det,fock_diag_tmp,N_int)

      sum_e_pert = 0d0

      do istate=1,N_states
        delta_E = E0(istate) - Hii + E_shift
        alpha_h_psi = mat(istate, p1, p2)
        val = alpha_h_psi + alpha_h_psi
        tmp = dsqrt(delta_E * delta_E + val * val)
        if (delta_E < 0.d0) then
            tmp = -tmp
        endif
        e_pert = 0.5d0 * (tmp - delta_E)
        coef(istate) = e_pert / alpha_h_psi
        pt2(istate) = pt2(istate) + e_pert
        variance(istate) = variance(istate) + alpha_h_psi * alpha_h_psi
        norm(istate) = norm(istate) + coef(istate) * coef(istate)

        if (weight_selection /= 5) then
          ! Energy selection
          sum_e_pert = sum_e_pert + e_pert * selection_weight(istate)
        else
          ! Variance selection
          sum_e_pert = sum_e_pert - alpha_h_psi * alpha_h_psi * selection_weight(istate)
        endif
      end do
      
      call give_2rdm_pert_contrib(det,coef,psi_det_connection,psi_coef_connection,n_det_connection,nkeys,keys,values,sze_buff)

      if(sum_e_pert <= buf%mini) then
        call add_to_selection_buffer(buf, det, sum_e_pert)
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
          call get_d2(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp(1, interesting(i)))
        else if(nt == 3) then
          call get_d1(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp(1, interesting(i)))
        else
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
          do k=1,N_states
            mat(k, putj, puti) = mat(k, putj, puti) +coefs(k) * hij
          enddo
        else
          do k=1,N_states
            mat(k, puti, putj) = mat(k, puti, putj) +coefs(k) * hij
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

          hij = mo_two_e_integral(p1, p2, h1, h2) * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
          do k=1,N_states
            mat(k, puti, putj) = mat(k, puti, putj) +coefs(k) * hij
          enddo
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
        hij = (mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2,p1, h1, h2)) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
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

        hij = mo_two_e_integral(p1, p2, h1, h2) * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2, N_int)
        do k=1,N_states
          mat(k, min(puti, putj), max(puti, putj)) = mat(k, min(puti, putj), max(puti, putj)) + coefs(k) * hij
        enddo
      end do
    else ! tip == 4
      puti = p(1, sp)
      putj = p(2, sp)
      if(.not. banned(puti,putj,1)) then
        p1 = p(1, mi)
        p2 = p(2, mi)
        h1 = h(1, mi)
        h2 = h(2, mi)
        hij = (mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2,p1, h1, h2)) * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2, N_int)
        do k=1,N_states
          mat(k, puti, putj) = mat(k, puti, putj) +coefs(k) * hij
        enddo
      end if
    end if
  end if
end


subroutine get_d1(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  use bitmasks
  implicit none

  integer(bit_kind), intent(in)  :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in)            :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind)              :: det(N_int, 2)
  double precision, intent(in)   :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in)            :: h(0:2,2), p(0:4,2), sp
  double precision, external     :: get_phase_bi, mo_two_e_integral
  logical                        :: ok

  logical, allocatable           :: lbanned(:,:)
  integer                        :: puti, putj, ma, mi, s1, s2, i, i1, i2, j
  integer                        :: hfix, pfix, h1, h2, p1, p2, ib, k

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
          tmp_row(1:N_states,putj) = tmp_row(1:N_states,putj) + hij * coefs(1:N_states)
        endif
      end do
      do putj=hfix+1, mo_num
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hij = hij_cache(putj,2) - hij_cache(putj,1)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
          tmp_row(1:N_states,putj) = tmp_row(1:N_states,putj) + hij * coefs(1:N_states)
        endif
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
    call get_mo_two_e_integrals(hfix,pfix,p1,mo_num,hij_cache(1,1),mo_integrals_map)
    call get_mo_two_e_integrals(hfix,pfix,p2,mo_num,hij_cache(1,2),mo_integrals_map)
    putj = p1
    do puti=1,mo_num
      !p1 fixed
      if(.not.(banned(putj,puti,bant).or.lbanned(puti,mi))) then
        hij = hij_cache(puti,2)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
          do k=1,N_states
            tmp_row(k,puti) = tmp_row(k,puti) + hij * coefs(k)
          enddo
        endif
      end if
    enddo

    putj = p2
    do puti=1,mo_num
      if(.not.(banned(putj,puti,bant)).or.(lbanned(puti,mi))) then
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
      mat(:,p1,:) = mat(:,p1,:) + tmp_row(:,:)
      mat(:,p2,:) = mat(:,p2,:) + tmp_row2(:,:)
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
        mat(:, puti, puti:) = mat(:, puti,puti:) + tmp_row(:,puti:)
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
        if(.not. banned(puti,putj,1)) then
          if(lbanned(puti,ma)) cycle
          hij = hij_cache(puti,1)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
            do k=1,N_states
              tmp_row(k,puti) = tmp_row(k,puti) + hij * coefs(k)
            enddo
          endif
        end if
      enddo

      putj = p1
      do puti=1,mo_num
        if(.not. banned(puti,putj,1)) then
          if(lbanned(puti,ma)) cycle
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
      mat(:,p2,p2:) = mat(:,p2,p2:) + tmp_row(:,p2:)
      mat(:,:p1-1,p1) = mat(:,:p1-1,p1) + tmp_row2(:,:p1-1)
      mat(:,p1,p1:) = mat(:,p1,p1:) + tmp_row2(:,p1:)
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
      p1 = p(i1,s1)
      if(s1 == s2) ib = i1+1
      if(bannedOrb(p1, s1)) cycle
      do i2=ib,p(0,s2)
        p2 = p(i2,s2)
        if(bannedOrb(p2, s2) .or. banned(p1, p2, 1)) cycle
        call apply_particles(mask, s1, p1, s2, p2, det, ok, N_int)
        call i_h_j(gen, det, N_int, hij)
        mat(:, p1, p2) = mat(:, p1, p2) + coefs(:) * hij
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
    do p2=1, mo_num
      if(bannedOrb(p2,2)) cycle
      call get_mo_two_e_integrals(p2,h1,h2,mo_num,hij_cache1,mo_integrals_map)
      do p1=1, mo_num
        if(bannedOrb(p1, 1) .or. banned(p1, p2, bant)) cycle
        if(p1 /= h1 .and. p2 /= h2) then
          if (hij_cache1(p1) == 0.d0) cycle
          phase = get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
          hij = hij_cache1(p1) * phase
        else
          call apply_particles(mask, 1,p1,2,p2, det, ok, N_int)
          call i_h_j(gen, det, N_int, hij)
          if (hij == 0.d0) cycle
        end if
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
        if(bannedOrb(putj, sp) .or. banned(putj, sp, bant)) cycle
        if(puti /= p1 .and. putj /= p2 .and. puti /= p2 .and. putj /= p1) then
          hij = hij_cache1(putj) -  hij_cache2(putj)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
            do k=1,N_states
              mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
            enddo
          endif
        else
          call apply_particles(mask, sp,puti,sp,putj, det, ok, N_int)
          call i_h_j(gen, det, N_int, hij)
          if (hij /= 0.d0) then
            do k=1,N_states
              mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
            enddo
          endif
        end if
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
    do i=1,p(0,1)
      do j=1,p(0,2)
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
  ! Gives the inidices(+1) of the bits set to 1 in the bit string
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
