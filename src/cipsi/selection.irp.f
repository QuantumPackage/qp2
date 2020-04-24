
use bitmasks

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

  avg = sum(rpt2(1:N_st)) / dble(N_st) - 1.d-32 ! Avoid future division by zero
  do k=1,N_st
    element = exp(dt*(rpt2(k)/avg -1.d0))
    element = min(1.5d0 , element)
    element = max(0.5d0 , element)
    memo_pt2(k,i_iter) = element
    pt2_match_weight(k) = product(memo_pt2(k,:))
  enddo

  avg = sum(variance(1:N_st)) / dble(N_st) + 1.d-32 ! Avoid future division by zero
  do k=1,N_st
    element = exp(dt*(variance(k)/avg -1.d0))
    element = min(1.5d0 , element)
    element = max(0.5d0 , element)
    memo_variance(k,i_iter) = element
    variance_match_weight(k) = product(memo_variance(k,:))
  enddo

  threshold_davidson_pt2 = min(1.d-6, &
     max(threshold_davidson, 1.e-1 * PT2_relative_error * minval(abs(rpt2(1:N_states)))) )

  SOFT_TOUCH pt2_match_weight variance_match_weight threshold_davidson_pt2
end


BEGIN_PROVIDER [ double precision, selection_weight, (N_states) ]
   implicit none
   BEGIN_DOC
   ! Weights used in the selection criterion
   END_DOC
   select case (weight_selection)

     case (0)
      print *,  'Using input weights in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * state_average_weight(1:N_states)

     case (1)
      print *,  'Using 1/c_max^2 weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) 

     case (2)
      print *,  'Using pt2-matching weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * pt2_match_weight(1:N_states)
      print *, '# PT2 weight ', real(pt2_match_weight(:),4)

     case (3)
      print *,  'Using variance-matching weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * variance_match_weight(1:N_states)
      print *, '# var weight ', real(variance_match_weight(:),4)

     case (4)
      print *,  'Using variance- and pt2-matching weights in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * sqrt(variance_match_weight(1:N_states) * pt2_match_weight(1:N_states))
      print *, '# PT2 weight ', real(pt2_match_weight(:),4)
      print *, '# var weight ', real(variance_match_weight(:),4)

     case (5)
      print *,  'Using variance-matching weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * variance_match_weight(1:N_states)
      print *, '# var weight ', real(variance_match_weight(:),4)

     case (6)
      print *,  'Using CI coefficient-based selection' 
      selection_weight(1:N_states) = c0_weight(1:N_states)

     case (7)
      print *,  'Input weights multiplied by variance- and pt2-matching'
      selection_weight(1:N_states) = c0_weight(1:N_states) * sqrt(variance_match_weight(1:N_states) * pt2_match_weight(1:N_states)) * state_average_weight(1:N_states)
      print *, '# PT2 weight ', real(pt2_match_weight(:),4)
      print *, '# var weight ', real(variance_match_weight(:),4)

     case (8)
      print *,  'Input weights multiplied by pt2-matching'
      selection_weight(1:N_states) = c0_weight(1:N_states) * pt2_match_weight(1:N_states) * state_average_weight(1:N_states)
      print *, '# PT2 weight ', real(pt2_match_weight(:),4)

     case (9)
      print *,  'Input weights multiplied by variance-matching'
      selection_weight(1:N_states) = c0_weight(1:N_states) * variance_match_weight(1:N_states) * state_average_weight(1:N_states)
      print *, '# var weight ', real(variance_match_weight(:),4)

    end select
     print *, '# Total weight ', real(selection_weight(:),4)

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
  !todo: simplify for kpts
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

  ! possible holes and particles for this generator
  ! hole_mask: occupied in this generator   .AND. occupied in generators_bitmask_hole
  ! part_mask: unoccupied in this generator .AND. occupied in generators_bitmask_part
  do k=1,N_int
      hole_mask(k,1) = iand(generators_bitmask(k,1,s_hole), psi_det_generators(k,1,i_generator))
      hole_mask(k,2) = iand(generators_bitmask(k,2,s_hole), psi_det_generators(k,2,i_generator))
      particle_mask(k,1) = iand(generators_bitmask(k,1,s_part), not(psi_det_generators(k,1,i_generator)) )
      particle_mask(k,2) = iand(generators_bitmask(k,2,s_part), not(psi_det_generators(k,2,i_generator)) )
  enddo
  call select_singles_and_doubles(i_generator,hole_mask,particle_mask,fock_diag_tmp,E0,pt2,variance,norm,b,subset,csubset)
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
  integer   :: kh1,kh2,kpt12,kk1,kk2,ik01,ik02,ik1,ik2
  integer(bit_kind)               :: hole(N_int,2), particle(N_int,2), mask(N_int, 2), pmask(N_int, 2)
  logical                         :: fullMatch, ok

  integer(bit_kind) :: mobMask(N_int, 2), negMask(N_int, 2)
  integer,allocatable               :: preinteresting(:), prefullinteresting(:)
  integer,allocatable               :: interesting(:), fullinteresting(:)
  integer,allocatable               :: tmp_array(:)
  integer(bit_kind), allocatable :: minilist(:, :, :), fullminilist(:, :, :)
  logical, allocatable           :: banned(:,:,:), bannedOrb(:,:)
  double precision, allocatable  :: coef_fullminilist_rev(:,:)
  complex*16, allocatable  :: coef_fullminilist_rev_complex(:,:)


  double precision, allocatable   :: mat(:,:,:)
  complex*16, allocatable   :: mat_complex(:,:,:)

  logical :: monoAdo, monoBdo
  integer :: maskInd

  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_order psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order
  if (is_complex) then
    PROVIDE psi_selectors_coef_transp_complex
  else
    PROVIDE psi_selectors_coef_transp
  endif

  monoAdo = .true.
  monoBdo = .true.
  
  !todo: this is already done in select_connected? why repeat?
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
  
  
  ! S_s = selectors
  ! S_0 = {|D_G>}                                       (i_generator determinant)
  ! S_j = {|D_k> : |D_k> \in T_j|D_G> }                 (i.e. S_2 is all dets connected to |D_G> by a double excitation)
  ! S_2b = S_2 \intersection {|D_k> : a_{h1}|D_k> != 0} (in S_2 and h1 is occupied)
  ! S_2' = S_2 \ {|D_k> : a_{h1}|D_k> != 0}             (in S_2 and h1 is not occupied)
  ! S_4b = S_4 \intersection {|D_k> : a_{h1}|D_k> != 0} (in S_4 and h1 is occupied)
  ! S_4' = S_4 \ {|D_k> : a_{h1}|D_k> != 0}             (in S_4 and h1 is not occupied)

  ! construct the following sets of determinants:
  !   preinteresting: S_pi = (U_{j=0..4} S_j) \intersection S_s
  !   prefullinteresting: S_pfi = (U_{j=0..2} S_j) \ S_s
  !   interesting: S_i = S_pi \ S_4b = ( (U_{j=0..3} S_j) U S_4' ) \intersection S_s
  !   fullinteresting: S_fi = S_i U (S_pfi \ S_2b) = (S_0 U S_1 U S_2') 
  !     (in order, first elements are in S_s, later elements are not in S_s) 


  ! get indices of all unique dets for which total excitation degree (relative to i_generator) is <= 4
  k=1
  ! get exc_degree(i) for each unique alpha det(i) from i_generator(alpha)
  do i=1,N_det_alpha_unique
    call get_excitation_degree_spin(psi_det_alpha_unique(1,i),       &
        psi_det_generators(1,1,i_generator), exc_degree(i), N_int)
  enddo
  
  ! get exc_degree (= nt) for each unique beta det(j) from i_generator(beta)
  do j=1,N_det_beta_unique
    call get_excitation_degree_spin(psi_det_beta_unique(1,j),        &
        psi_det_generators(1,2,i_generator), nt, N_int)
    if (nt > 2) cycle ! don't keep anything more than double beta exc
    do l_a=psi_bilinear_matrix_columns_loc(j), psi_bilinear_matrix_columns_loc(j+1)-1
      i = psi_bilinear_matrix_rows(l_a)
      if (nt + exc_degree(i) <= 4) then ! don't keep anything more than 4-fold total exc
        idx = psi_det_sorted_order(psi_bilinear_matrix_order(l_a))
        if (psi_average_norm_contrib_sorted(idx) > 0.d0) then
          indices(k) = idx
          k=k+1
        endif
      endif
    enddo
  enddo
  

  !  indices now contains det indices (in psi_det_sorted) of dets which differ from generator by:
  !  (exc_alpha,exc_beta) in
  !  (4,0)
  !  (3,0), (3,1)
  !  (2,0), (2,1), (2,2)
  !  (1,0), (1,1), (1,2)
  !  (0,0), (0,1), (0,2)
  !
  !                              (4,0)
  !                       (3,0), (3,1)
  !                (2,0), (2,1), (2,2)
  !         (1,0), (1,1), (1,2)
  !  (0,0), (0,1), (0,2)
  !
  !  below, add (0,3), (0,4), (1,3)

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
        if (psi_average_norm_contrib_sorted(idx) > 0.d0) then
          indices(k) = idx
          k=k+1
        endif
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
  ! sort indices by location in psi_det_sorted
  
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
    ! mobMask in psi_det(i) but not in i_generator
    ! nt = popcnt(mobMask)
    mobMask(1,1) = iand(negMask(1,1), psi_det_sorted(1,1,i))
    mobMask(1,2) = iand(negMask(1,2), psi_det_sorted(1,2,i))
    nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
    do j=2,N_int
      mobMask(j,1) = iand(negMask(j,1), psi_det_sorted(j,1,i))
      mobMask(j,2) = iand(negMask(j,2), psi_det_sorted(j,2,i))
      nt = nt + popcnt(mobMask(j, 1)) + popcnt(mobMask(j, 2))
    end do
    
    ! preinteresting:     within a 4-fold excitation from i_generator; in selectors
    ! prefullinteresting: within a double excitation from i_generator; not in selectors

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
  if (is_complex) then
  allocate (mat_complex(N_states, mo_num, mo_num))
  else
  allocate (mat(N_states, mo_num, mo_num))
  endif
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
!todo: kpts
      kh1 = (h1-1)/mo_num_per_kpt + 1
      ! pmask is i_generator det with bit at h1 set to zero
      call apply_hole(psi_det_generators(1,1,i_generator), s1,h1, pmask, ok, N_int)
      
      negMask = not(pmask)
!      
      ! see set definitions above 
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
        
        ! nt = ( orbs occupied in preinteresting(ii) and not occupied in i_gen(after removing elec from h1) )
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
      if(pert_2rdm)then
        if (is_complex) then
          print*,irp_here,' not implemented for complex: pert_2rdm'
          stop -1
        else
          allocate(coef_fullminilist_rev(N_states,fullinteresting(0))) 
          do i=1,fullinteresting(0)
            do j = 1, N_states
              coef_fullminilist_rev(j,i) = psi_coef_sorted(fullinteresting(i),j)
            enddo
          enddo
        endif
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
          if (is_complex) then
!=============================================================
!!todo use this once kpts are implemented
            kh2 = (h2-1)/mo_num_per_kpt + 1
            kpt12 = kconserv(kh1,kh2,1)
            ! mask is gen_i with (h1,s1),(h2,s2) removed
            call apply_hole(pmask, s2,h2, mask, ok, N_int)
            banned = .true.
            ! only allow excitations that conserve momentum
            do kk1=1,kpt_num
              ! equivalent to kk2 = kconserv(kh1,kh2,kk1)
              kk2 = kconserv(kpt12,1,kk1)
              ik01 = (kk1-1) * mo_num_per_kpt + 1 !first mo in kk1
              ik02 = (kk2-1) * mo_num_per_kpt + 1 !first mo in kk2
              do ik1 = ik01, ik01 + mo_num_per_kpt - 1 !loop over mos in kk1
                do ik2 = ik02, ik02 + mo_num_per_kpt - 1 !loop over mos in kk2
                  ! depending on sp, might not need both of these?
                  ! sp=1 (a,a) or sp=2 (b,b): only use banned(:,:,1)
                  ! sp=3 (a,b): banned(alpha,beta,1) is transpose of banned(beta,alpha,2)
                  banned(ik1,ik2,1) = .false.
                  banned(ik1,ik2,2) = .false.
                enddo
              enddo
            enddo
!=============================================================
!            ! mask is gen_i with (h1,s1),(h2,s2) removed
!            call apply_hole(pmask, s2,h2, mask, ok, N_int)
!            banned = .false.
!=============================================================
          else
            call apply_hole(pmask, s2,h2, mask, ok, N_int)
            banned = .false.
          endif
          do j=1,mo_num
            bannedOrb(j, 1) = .true.
            bannedOrb(j, 2) = .true.
          enddo
          do s3=1,2
            do i=1,N_particles(s3)
              bannedOrb(particle_list(i,s3), s3) = .false. ! allow excitation into orbitals in particle_list
            enddo
          enddo
          if(s1 /= s2) then
            if(monoBdo) then
              bannedOrb(h1,s1) = .false. ! allow alpha elec to go back into alpha hole
            end if
            if(monoAdo) then
              bannedOrb(h2,s2) = .false. ! allow beta elec to go back into beta hole
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
            if (is_complex) then
              call splash_pq_complex(mask, sp, minilist, i_generator, interesting(0), bannedOrb, banned, mat_complex, interesting)
              
              if(.not.pert_2rdm)then
               call fill_buffer_double_complex(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2, variance, norm, mat_complex, buf)
              else
                print*,irp_here,' not implemented for complex (fill_buffer_double_rdm_complex)'
                stop -1
               !call fill_buffer_double_rdm_complex(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2, variance, norm, mat_complex, buf,fullminilist, coef_fullminilist_rev_complex, fullinteresting(0))
              endif
            else
            call splash_pq(mask, sp, minilist, i_generator, interesting(0), bannedOrb, banned, mat, interesting)
            
            if(.not.pert_2rdm)then
             call fill_buffer_double(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2, variance, norm, mat, buf)
            else 
             call fill_buffer_double_rdm(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2, variance, norm, mat, buf,fullminilist, coef_fullminilist_rev, fullinteresting(0))
            endif
            endif!complex
          end if
        enddo !i2
        if(s1 /= s2) monoBdo = .false.
      enddo !s2
      deallocate(fullminilist,minilist)
      if(pert_2rdm)then
        if (is_complex) then
          print*,irp_here,' not implemented for complex: pert_2rdm'
          stop -1
        else
          deallocate(coef_fullminilist_rev)
        endif
      endif
    enddo ! i1
  enddo ! s1
  deallocate(preinteresting, prefullinteresting, interesting, fullinteresting)
  deallocate(banned, bannedOrb)
  if (is_complex) then
    deallocate(mat_complex)
  else
    deallocate(mat)
  endif
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
  double precision :: e_pert, delta_E, val, Hii, w, tmp, alpha_h_psi, coef
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

      Hii = diag_h_mat_elem_fock(psi_det_generators(1,1,i_generator),det,fock_diag_tmp,N_int)

      w = 0d0

!      integer(bit_kind) :: occ(N_int,2), n
!      call occ_pattern_of_det(det,occ,N_int)
!      call occ_pattern_to_dets_size(occ,n,elec_alpha_num,N_int)


      do istate=1,N_states
        delta_E = E0(istate) - Hii + E_shift
        alpha_h_psi = mat(istate, p1, p2)
        val = alpha_h_psi + alpha_h_psi
        tmp = dsqrt(delta_E * delta_E + val * val)
        if (delta_E < 0.d0) then
            tmp = -tmp
        endif
        e_pert = 0.5d0 * (tmp - delta_E)
        if (dabs(alpha_h_psi) > 1.d-4) then
          coef = e_pert / alpha_h_psi
        else
          coef = alpha_h_psi / delta_E
        endif
        pt2(istate) = pt2(istate) + e_pert
        variance(istate) = variance(istate) + alpha_h_psi * alpha_h_psi
        norm(istate) = norm(istate) + coef * coef

!!!DEBUG
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
            w = w - alpha_h_psi * alpha_h_psi * selection_weight(istate)

          case(6)
            w = w - coef * coef * selection_weight(istate)

          case default
            ! Energy selection
            w = w + e_pert * selection_weight(istate)

        end select
      end do


      if(pseudo_sym)then
        if(dabs(mat(1, p1, p2)).lt.thresh_sym)then 
          w = 0.d0
        endif
      endif

!      w = dble(n) * w

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
    do j=1, N_int  ! if all occupied orbs in mask are not also occupied in det(i), go to next det
      if(iand(det(j,1,i), mask(j,1)) /= mask(j, 1)) cycle genl
      if(iand(det(j,2,i), mask(j,2)) /= mask(j, 2)) cycle genl
    end do

    ! If det(i) < det(i_gen), it hs already been considered
    if(interesting(i) < i_gen) then
      fullMatch = .true.
      return
    end if

    ! Identify the particles
    do j=1, N_int ! if electrons are excited into the orbs given by myMask, resulting determinant will be det(i)
      myMask(j, 1) = iand(det(j, 1, i), negMask(j, 1))
      myMask(j, 2) = iand(det(j, 2, i), negMask(j, 2))
    end do

    ! don't allow excitations into this pair of orbitals?
    ! should 'banned' have dimensions (mo_num,mo_num,2)?
    ! is it always true that popcnt(myMask) = 2 ? (sum over N_int and alpha/beta spins)
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
        mat(:, p1, p2) += coefs(:) * hij
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
        mat(:, puti, putj) += coefs(:) * hij
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
        tmp_row(1:N_states,putj) += hij * coefs(1:N_states)
      end do
      do putj=hfix+1, mo_num
        if(lbanned(putj, ma) .or. banned(putj, puti,bant)) cycle
        hij = (mo_two_e_integral(p1, p2, hfix, putj)-mo_two_e_integral(p2,p1,hfix,putj)) * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
        tmp_row(1:N_states,putj) += hij * coefs(1:N_states)
      end do

      if(ma == 1) then           
        mat(1:N_states,1:mo_num,puti) += tmp_row(1:N_states,1:mo_num)
      else
        mat(1:N_states,puti,1:mo_num) += tmp_row(1:N_states,1:mo_num)
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
        tmp_row(:,puti) += hij * coefs(:)
      end if
      
      putj = p2
      if(.not. banned(putj,puti,bant)) then
        hij = mo_two_e_integral(p1,pfix,hfix,puti) * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
        tmp_row2(:,puti) += hij * coefs(:)
      end if
    end do
    
    if(mi == 1) then
      mat(:,:,p1) += tmp_row(:,:)
      mat(:,:,p2) += tmp_row2(:,:)
    else
      mat(:,p1,:) += tmp_row(:,:)
      mat(:,p2,:) += tmp_row2(:,:)
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
          tmp_row(:,putj) += hij * coefs(:)
        end do
        do putj=hfix+1,mo_num
          if(lbanned(putj,ma) .or. banned(puti,putj,1)) cycle
          hij = (mo_two_e_integral(p1, p2, hfix, putj)-mo_two_e_integral(p2,p1,hfix,putj)) * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
          tmp_row(:,putj) += hij * coefs(:)
        end do

        mat(:, :puti-1, puti) += tmp_row(:,:puti-1)
        mat(:, puti, puti:) += tmp_row(:,puti:)
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
          tmp_row(:,puti) += hij * coefs(:)
        end if
        
        putj = p1
        if(.not. banned(puti,putj,1)) then
          hij = mo_two_e_integral(pfix, p2, hfix, puti) * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
          tmp_row2(:,puti) += hij * coefs(:)
        end if
      end do
      mat(:,:p2-1,p2) += tmp_row(:,:p2-1)
      mat(:,p2,p2:) += tmp_row(:,p2:)
      mat(:,:p1-1,p1) += tmp_row2(:,:p1-1)
      mat(:,p1,p1:) += tmp_row2(:,p1:)
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
        mat(:, p1, p2) += coefs(:) * hij
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
          mat(:, putj, puti) += coefs(:) * hij
        else
          mat(:, puti, putj) += coefs(:) * hij
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
          mat(:, puti, putj) += coefs(:) * hij
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
        mat(:, puti, putj) += coefs(:) * hij
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
        mat(:, min(puti, putj), max(puti, putj)) += coefs(:) * hij
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
        mat(:, puti, putj) += coefs(:) * hij
      end if
    end if
  end if
end 


!==============================================================================!
!                                                                              !
!                                    Complex                                   !
!                                                                              !
!==============================================================================!

subroutine fill_buffer_double_complex(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2, variance, norm, mat, buf)
  !todo: should be okay for complex
  use bitmasks
  use selection_types
  implicit none

  integer, intent(in) :: i_generator, sp, h1, h2
  complex*16, intent(in) :: mat(N_states, mo_num, mo_num)
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
  double precision :: e_pert, delta_E, val, Hii, w, tmp
  complex*16 ::  alpha_h_psi, coef, val_c
  double precision, external :: diag_H_mat_elem_fock
  double precision :: E_shift

!  logical, external :: detEq
!  double precision, allocatable :: values(:)
!  integer, allocatable          :: keys(:,:)
!  integer                       :: nkeys
  

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

      val = maxval(cdabs(mat(1:N_states, p1, p2)))
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

      Hii = diag_h_mat_elem_fock(psi_det_generators(1,1,i_generator),det,fock_diag_tmp,N_int)

      w = 0d0

!      integer(bit_kind) :: occ(N_int,2), n
!      call occ_pattern_of_det(det,occ,N_int)
!      call occ_pattern_to_dets_size(occ,n,elec_alpha_num,N_int)


      do istate=1,N_states
        delta_E = E0(istate) - Hii + E_shift
        alpha_h_psi = mat(istate, p1, p2)
        val_c = alpha_h_psi + alpha_h_psi
        tmp = dsqrt(delta_E * delta_E + cdabs(val_c * val_c))
        if (delta_E < 0.d0) then
            tmp = -tmp
        endif
        e_pert = 0.5d0 * (tmp - delta_E)
        if (cdabs(alpha_h_psi) > 1.d-4) then
          coef = e_pert / alpha_h_psi
        else
          coef = alpha_h_psi / delta_E
        endif
        pt2(istate) = pt2(istate) + e_pert
        variance(istate) = variance(istate) + cdabs(alpha_h_psi * alpha_h_psi)
        norm(istate) = norm(istate) + cdabs(coef * coef)

!!!DEBUG
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
            w = w - cdabs(alpha_h_psi * alpha_h_psi) * selection_weight(istate)

          case(6)
            w = w - cdabs(coef * coef) * selection_weight(istate)

          case default
            ! Energy selection
            w = w + e_pert * selection_weight(istate)

        end select
      end do


      if(pseudo_sym)then
        if(cdabs(mat(1, p1, p2)).lt.thresh_sym)then 
          w = 0.d0
        endif
      endif

!      w = dble(n) * w

      if(w <= buf%mini) then
        call add_to_selection_buffer(buf, det, w)
      end if
    end do
  end do
end

subroutine splash_pq_complex(mask, sp, det, i_gen, N_sel, bannedOrb, banned, mat, interesting)
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
  ! mat should be out, not inout? (if only called from select_singles_and_doubles)
  complex*16, intent(inout) :: mat(N_states, mo_num, mo_num)

  integer                        :: i, ii, j, k, l, h(0:2,2), p(0:4,2), nt
  integer(bit_kind)              :: perMask(N_int, 2), mobMask(N_int, 2), negMask(N_int, 2)
  integer(bit_kind)             :: phasemask(N_int,2)

  PROVIDE psi_selectors_coef_transp_complex psi_det_sorted
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

    ! p contains orbs in det that are not in the doubly ionized generator
    if (interesting(i) >= i_gen) then ! det past i_gen
        call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
        call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)

        perMask(1,1) = iand(mask(1,1), not(det(1,1,i)))
        perMask(1,2) = iand(mask(1,2), not(det(1,2,i)))
        do j=2,N_int
          perMask(j,1) = iand(mask(j,1), not(det(j,1,i)))
          perMask(j,2) = iand(mask(j,2), not(det(j,2,i)))
        end do

      ! h contains orbs in the doubly ionized generator that are not in det
        call bitstring_to_list_in_selection(perMask(1,1), h(1,1), h(0,1), N_int)
        call bitstring_to_list_in_selection(perMask(1,2), h(1,2), h(0,2), N_int)

        call get_mask_phase(psi_det_sorted(1,1,interesting(i)), phasemask,N_int)
        if(nt == 4) then ! differ by 6 (2,4)
          call get_d2_complex(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp_complex(1, interesting(i)))
        else if(nt == 3) then ! differ by 4 (1,3)
          call get_d1_complex(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp_complex(1, interesting(i)))
        else ! differ by 2 (0,2) 
          call get_d0_complex(det(1,1,i), phasemask, bannedOrb, banned, mat, mask, h, p, sp, psi_selectors_coef_transp_complex(1, interesting(i)))
        end if
    else if(nt == 4) then ! differ by 6 (2,4); i_gen past det
        call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
        call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)
        call past_d2(banned, p, sp)
    else if(nt == 3) then ! differ by 4 (1,3); i_gen past det
        call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
        call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)
        call past_d1(bannedOrb, p)
    end if
  end do

end


subroutine get_d2_complex(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  !todo: indices/conjg should be correct for complex
  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  complex*16, intent(in) :: coefs(N_states)
  complex*16, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  double precision, external :: get_phase_bi
  complex*16, external :: mo_two_e_integral_complex

  integer :: i, j, k, tip, ma, mi, puti, putj
  integer :: h1, h2, p1, p2, i1, i2
  double precision :: phase
  complex*16 :: hij

  integer, parameter:: turn2d(2,3,4) = reshape((/0,0, 0,0, 0,0,  3,4, 0,0, 0,0,  2,4, 1,4, 0,0,  2,3, 1,3, 1,2 /), (/2,3,4/))
  integer, parameter :: turn2(2) = (/2, 1/)
  integer, parameter :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer :: bant
  bant = 1

  tip = p(0,1) * p(0,2) ! number of alpha particles times number of beta particles

  ma = sp !1:(alpha,alpha); 2:(b,b); 3:(a,b)
  if(p(0,1) > p(0,2)) ma = 1 ! more alpha particles than beta particles
  if(p(0,1) < p(0,2)) ma = 2 ! fewer alpha particles than beta particles
  mi = mod(ma, 2) + 1

  if(sp == 3) then ! if one alpha and one beta xhole 
    !(where xholes refer to the ionizations from the generator, not the holes occupied in the ionized generator)
    if(ma == 2) bant = 2 ! if more beta particles than alpha particles

    if(tip == 3) then ! if 3 of one particle spin and 1 of the other particle spin
      puti = p(1, mi)
      if(bannedOrb(puti, mi)) return
      h1 = h(1, ma)
      h2 = h(2, ma)

      do i = 1, 3    ! loop over all 3 combinations of 2 particles with spin ma
        putj = p(i, ma)
        if(banned(putj,puti,bant)) cycle
        i1 = turn3(1,i)
        i2 = turn3(2,i)
        p1 = p(i1, ma)
        p2 = p(i2, ma)
        
     ! |G> = |psi_{gen,i}>
     ! |G'> = a_{x1} a_{x2} |G>
     ! |alpha> = a_{puti}^{\dagger} a_{putj}^{\dagger} |G'>
     ! |alpha> = t_{x1,x2}^{puti,putj} |G>
     ! hij = <psi_{selectors,i}|H|alpha>
     ! |alpha> = t_{p1,p2}^{h1,h2}|psi_{selectors,i}>
        !todo: <i|H|j>  =  (<h1,h2|p1,p2> - <h1,h2|p2,p1>) * phase
        !    <psi|H|j> +=  dconjg(c_i) * <i|H|j>
        !      <j|H|i>  =  (<p1,p2|h1,h2> - <p2,p1|h1,h2>) * phase
        !    <j|H|psi> +=  <j|H|i> * c_i
        hij = mo_two_e_integral_complex(p1, p2, h1, h2) - mo_two_e_integral_complex(p2, p1, h1, h2)
        if (hij == (0.d0,0.d0)) cycle

        ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
        hij = dconjg(hij) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)

        if(ma == 1) then ! if particle spins are (alpha,alpha,alpha,beta), then puti is beta and putj is alpha
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, putj, puti) = mat(k, putj, puti) + coefs(k) * hij
          enddo
        else            ! if particle spins are (beta,beta,beta,alpha), then puti is alpha and putj is beta
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
          enddo
        end if
      end do
    else ! if 2 alpha and 2 beta particles
      h1 = h(1,1)
      h2 = h(1,2)
      do j = 1,2 ! loop over all 4 combinations of one alpha and one beta particle
        putj = p(j, 2)
        if(bannedOrb(putj, 2)) cycle
        p2 = p(turn2(j), 2)
        do i = 1,2
          puti = p(i, 1)

          if(banned(puti,putj,bant) .or. bannedOrb(puti,1)) cycle
          p1 = p(turn2(i), 1)

    ! hij = <psi_{selectors,i}|H|alpha> 
          hij = mo_two_e_integral_complex(p1, p2, h1, h2)
          if (hij /= (0.d0,0.d0)) then
            ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
            hij = dconjg(hij) * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
            enddo
          endif
        end do
      end do
    end if

  else ! if holes are (a,a) or (b,b)
    if(tip == 0) then ! if particles are (a,a,a,a) or (b,b,b,b)
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
          hij = mo_two_e_integral_complex(p1, p2, h1, h2) - mo_two_e_integral_complex(p2,p1, h1, h2)
          if (hij == (0.d0,0.d0)) cycle

          ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
          hij = dconjg(hij) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, putj) = mat(k, puti, putj) +coefs(k) * hij
          enddo
        end do
      end do
    else if(tip == 3) then ! if particles are (a,a,a,b) (ma=1,mi=2) or (a,b,b,b) (ma=2,mi=1)
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

        hij = mo_two_e_integral_complex(p1, p2, h1, h2)
        if (hij == (0.d0,0.d0)) cycle

        ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
        hij = dconjg(hij) * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2, N_int)
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
    else ! tip == 4  (a,a,b,b)
      puti = p(1, sp)
      putj = p(2, sp)
      if(.not. banned(puti,putj,1)) then
        p1 = p(1, mi)
        p2 = p(2, mi)
        h1 = h(1, mi)
        h2 = h(2, mi)
        hij = (mo_two_e_integral_complex(p1, p2, h1, h2) - mo_two_e_integral_complex(p2,p1, h1, h2))
        if (hij /= (0.d0,0.d0)) then
          ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
          hij = dconjg(hij) * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
          enddo
        end if
      end if
    end if
  end if
end


subroutine get_d1_complex(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  !todo: indices should be okay for complex?
  use bitmasks
  implicit none

  integer(bit_kind), intent(in)  :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in)  :: phasemask(N_int,2)
  logical, intent(in)            :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind)              :: det(N_int, 2)
  complex*16, intent(in)   :: coefs(N_states)
  complex*16, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in)            :: h(0:2,2), p(0:4,2), sp
  double precision, external     :: get_phase_bi
  complex*16, external     :: mo_two_e_integral_complex
  logical                        :: ok

  logical, allocatable           :: lbanned(:,:)
  integer                        :: puti, putj, ma, mi, s1, s2, i, i1, i2, j
  integer                        :: hfix, pfix, h1, h2, p1, p2, ib, k, l

  integer :: kp1,ip1, kp2,ip2, khfix,ihfix, kputi,iputi, kputj,iputj, putj0
  integer :: kpfix, ipfix, puti0
  integer :: kputi1,kputi2,puti01,puti02
  integer :: ii0

  integer, parameter             :: turn2(2) = (/2,1/)
  integer, parameter             :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer                        :: bant
  complex*16, allocatable :: hij_cache(:,:),hij_cache2(:,:)
  complex*16               :: hij, tmp_row(N_states, mo_num), tmp_row2(N_states, mo_num)
  complex*16               :: tmp_row_kpts(N_states, mo_num), tmp_row2_kpts(N_states, mo_num)
  complex*16               :: tmp_row_kpts2(N_states, mo_num_per_kpt), tmp_row2_kpts2(N_states,mo_num_per_kpt)
  complex*16 :: tmp_mat1(N_states,mo_num,mo_num), tmp_mat2(N_states,mo_num,mo_num)
  PROVIDE mo_integrals_map N_int

  allocate (lbanned(mo_num, 2))
  allocate (hij_cache(mo_num,2),hij_cache2(mo_num_per_kpt,2))
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
    kputi = (puti-1)/mo_num_per_kpt + 1
    khfix = (hfix-1)/mo_num_per_kpt + 1
    kp1   =   (p1-1)/mo_num_per_kpt + 1
    kp2   =   (p2-1)/mo_num_per_kpt + 1
    iputi = mod(puti-1,mo_num_per_kpt) + 1
    ihfix = mod(hfix-1,mo_num_per_kpt) + 1
    ip1   = mod(p1-1,  mo_num_per_kpt) + 1
    ip2   = mod(p2-1,  mo_num_per_kpt) + 1

    if(.not. bannedOrb(puti, mi)) then
      !call get_mo_two_e_integrals_complex(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
      !call get_mo_two_e_integrals_complex(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
      call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p1,ip1,kp1,p2,ip2,kp2,mo_num_per_kpt,hij_cache2(1,1),mo_integrals_map,mo_integrals_map_2)
      call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p2,ip2,kp2,p1,ip1,kp1,mo_num_per_kpt,hij_cache2(1,2),mo_integrals_map,mo_integrals_map_2)
      tmp_row = (0.d0,0.d0)
      tmp_row_kpts2 = (0.d0,0.d0)
      kputj = kconserv(kp1,kp2,khfix)
      putj0 = (kputj-1)*mo_num_per_kpt
      !do putj=1, hfix-1
      !  if(lbanned(putj, ma)) cycle
      !  if(banned(putj, puti,bant)) cycle
      !  hij = hij_cache(putj,1) - hij_cache(putj,2)
      !  if (hij /= (0.d0,0.d0)) then
      !    hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
      !    !DIR$ LOOP COUNT AVG(4)
      !    do k=1,N_states
      !      tmp_row(k,putj) = tmp_row(k,putj) + hij * coefs(k)
      !    enddo
      !  endif
      !end do
      !do putj=hfix+1, mo_num
      !  if(lbanned(putj, ma)) cycle
      !  if(banned(putj, puti,bant)) cycle
      !  hij = hij_cache(putj,2) - hij_cache(putj,1)
      !  if (hij /= (0.d0,0.d0)) then
      !    hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
      !    !DIR$ LOOP COUNT AVG(4)
      !    do k=1,N_states
      !      tmp_row(k,putj) = tmp_row(k,putj) + hij * coefs(k)
      !    enddo
      !  endif
      !end do
      !===========================
      ! begin kpts testing
      do putj = putj0+1, hfix-1
        iputj = putj-putj0
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hij = hij_cache2(iputj,1) - hij_cache2(iputj,2)
        if (hij /= (0.d0,0.d0)) then
          hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            !tmp_row_kpts(k,putj) = tmp_row_kpts(k,putj) + hij * coefs(k)
            tmp_row_kpts2(k,iputj) = tmp_row_kpts2(k,iputj) + hij * coefs(k)
          enddo
        endif
      end do
      do putj = hfix+1,putj0+mo_num_per_kpt
        iputj = putj - putj0
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hij = hij_cache2(iputj,2) - hij_cache2(iputj,1)
        if (hij /= (0.d0,0.d0)) then
          hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            !tmp_row_kpts(k,putj) = tmp_row_kpts(k,putj) + hij * coefs(k)
            tmp_row_kpts2(k,iputj) = tmp_row_kpts2(k,iputj) + hij * coefs(k)
          enddo
        endif
      end do
      ! end kpts testing
      !===========================================================
      !print*,'tmp_row_k,tmp_row'
      !do ii0=1,mo_num
      !  if (cdabs(tmp_row_kpts(1,ii0)-tmp_row(1,ii0)).gt.1.d-12) then
      !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG, ',ii0,hfix,p1,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
      !  endif
      !enddo
      !===========================================================
      if(ma == 1) then
        !mat(1:N_states,1:mo_num,puti) = mat(1:N_states,1:mo_num,puti) + tmp_row(1:N_states,1:mo_num)
        mat(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) = mat(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) + &
                                         tmp_row_kpts2(1:N_states,1:mo_num_per_kpt)
      else
        do l=1,mo_num_per_kpt
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k,puti,l+putj0) = mat(k,puti,l+putj0) + tmp_row_kpts2(k,l)
          enddo
        enddo
      end if
    end if

    !MOVE MI
    pfix = p(1,mi)
    kpfix = (pfix-1)/mo_num_per_kpt + 1
    ipfix = mod(pfix-1,mo_num_per_kpt) + 1
    !tmp_row = (0.d0,0.d0)
    !tmp_row2 = (0.d0,0.d0)
    !tmp_row_kpts = (0.d0,0.d0)
    !tmp_row2_kpts = (0.d0,0.d0)
    tmp_row_kpts2 = (0.d0,0.d0)
    tmp_row2_kpts2 = (0.d0,0.d0)
    !call get_mo_two_e_integrals_complex(hfix,pfix,p1,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
    !call get_mo_two_e_integrals_complex(hfix,pfix,p2,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
    call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,pfix,ipfix,kpfix,p1,ip1,kp1,mo_num_per_kpt,hij_cache2(1,1),mo_integrals_map,mo_integrals_map_2)
    call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,pfix,ipfix,kpfix,p2,ip2,kp2,mo_num_per_kpt,hij_cache2(1,2),mo_integrals_map,mo_integrals_map_2)
    putj = p1
    !============
    !begin ref
    !do puti=1,mo_num !HOT
    !  if(lbanned(puti,mi)) cycle
    !  !p1 fixed
    !  putj = p1
    !  if(.not. banned(putj,puti,bant)) then
    !    hij = hij_cache(puti,2)
    !    if (hij /= (0.d0,0.d0)) then
    !      hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
    !      !DIR$ LOOP COUNT AVG(4)
    !      do k=1,N_states
    !        tmp_row(k,puti) = tmp_row(k,puti) + hij * coefs(k)
    !      enddo
    !    endif
    !  end if
!   ! enddo
!   !   
    !  putj = p2
!   ! do puti=1,mo_num !HOT
    !  if(.not. banned(putj,puti,bant)) then
    !    hij = hij_cache(puti,1)
    !    if (hij /= (0.d0,0.d0)) then
    !      hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
    !      do k=1,N_states
    !        tmp_row2(k,puti) = tmp_row2(k,puti) + hij * coefs(k)
    !      enddo
    !    endif
    !  end if
    !end do
    !end ref
    !===================
    !begin kpts
    if (kp1.eq.kp2) then
    !if (.False.) then
      kputi1 = kconserv(kpfix,kp1,khfix)
      kputi2 = kputi1
      puti01 = (kputi1-1)*mo_num_per_kpt
      puti02 = puti01
      do iputi=1,mo_num_per_kpt !HOT
        puti = puti01 + iputi
        if(lbanned(puti,mi)) cycle
        !p1 fixed
        putj = p1
        if(.not. banned(putj,puti,bant)) then
          hij = hij_cache2(iputi,2)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              tmp_row_kpts2(k,iputi) = tmp_row_kpts2(k,iputi) + hij * coefs(k)
              !tmp_row_kpts(k,puti) = tmp_row_kpts(k,puti) + hij * coefs(k)
            enddo
          endif
        end if
!      enddo
!        
        putj = p2
!      do puti=1,mo_num !HOT
        if(.not. banned(putj,puti,bant)) then
          hij = hij_cache2(iputi,1)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
            do k=1,N_states
              tmp_row2_kpts2(k,iputi) = tmp_row2_kpts2(k,iputi) + hij * coefs(k)
              !tmp_row2_kpts(k,puti) = tmp_row2_kpts(k,puti) + hij * coefs(k)
            enddo
          endif
        end if
      end do
    else !kp1.ne.kp2
      kputi2 = kconserv(kpfix,kp2,khfix)
      puti02 = (kputi2-1)*mo_num_per_kpt
      putj = p1
      do iputi=1,mo_num_per_kpt !HOT
        puti = puti02 + iputi
        if(lbanned(puti,mi)) cycle
        !p1 fixed
        if(.not. banned(putj,puti,bant)) then
          hij = hij_cache2(iputi,2)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              tmp_row_kpts2(k,iputi) = tmp_row_kpts2(k,iputi) + hij * coefs(k)
              !tmp_row_kpts(k,puti) = tmp_row_kpts(k,puti) + hij * coefs(k)
            enddo
          endif
        end if
      enddo
!        
      putj = p2
      kputi1 = kconserv(kpfix,kp1,khfix)
      puti01 = (kputi1-1)*mo_num_per_kpt
      do iputi=1,mo_num_per_kpt !HOT
        puti = puti01 + iputi
        if(lbanned(puti,mi)) cycle
        if(.not. banned(putj,puti,bant)) then
          hij = hij_cache2(iputi,1)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
            do k=1,N_states
              tmp_row2_kpts2(k,iputi) = tmp_row2_kpts2(k,iputi) + hij * coefs(k)
              !tmp_row2_kpts(k,puti) = tmp_row2_kpts(k,puti) + hij * coefs(k)
            enddo
          endif
        end if
      end do
    endif
    !end kpts
    !===================
    !test printing
    !print'((A),5(I5))','kpt info1: ',kconserv(kpfix,kp2,khfix),khfix,kpfix,kp2,kputi2
    !print'((A),5(I5))','kpt info2: ',kconserv(kpfix,kp1,khfix),khfix,kpfix,kp1,kputi1
    !do ii0=1,mo_num
    !  if (cdabs(tmp_row_kpts(1,ii0)-tmp_row(1,ii0)).gt.1.d-12) then
    !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1a, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
    !!  else if ((cdabs(tmp_row_kpts(1,ii0))+cdabs(tmp_row(1,ii0))).gt.1.d-12) then
    !!    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1b, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
    !  endif
    !  if (cdabs(tmp_row2_kpts(1,ii0)-tmp_row2(1,ii0)).gt.1.d-12) then
    !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 2a, ',ii0,hfix,pfix,p1,tmp_row2_kpts(1,ii0),tmp_row2(1,ii0)
    !!  else if ((cdabs(tmp_row2_kpts(1,ii0))+cdabs(tmp_row2(1,ii0))).gt.1.d-12) then
    !!    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 2b, ',ii0,hfix,pfix,p1,tmp_row2_kpts(1,ii0),tmp_row2(1,ii0)
    !  endif
    !enddo
    !===================

    if(mi == 1) then
      !mat(:,:,p1) = mat(:,:,p1) + tmp_row(:,:)
      !mat(:,:,p2) = mat(:,:,p2) + tmp_row2(:,:)
      mat(:,puti02+1:puti02+mo_num_per_kpt,p1) = mat(:,puti02+1:puti02+mo_num_per_kpt,p1) + tmp_row_kpts2(:,:)
      mat(:,puti01+1:puti01+mo_num_per_kpt,p2) = mat(:,puti01+1:puti01+mo_num_per_kpt,p2) + tmp_row2_kpts2(:,:)
    else
      do l=1,mo_num_per_kpt
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k,p1,l+puti02) = mat(k,p1,l+puti02) + tmp_row_kpts2(k,l)
          mat(k,p2,l+puti01) = mat(k,p2,l+puti01) + tmp_row2_kpts2(k,l)
        enddo
      enddo
    end if
    !todo: kpts okay up to this point in get_d1_complex

  else  ! sp /= 3

    if(p(0,ma) == 3) then
      do i=1,3
        hfix = h(1,ma)
        puti = p(i, ma)
        p1 = p(turn3(1,i), ma)
        p2 = p(turn3(2,i), ma)
        kputi = (puti-1)/mo_num_per_kpt + 1
        khfix = (hfix-1)/mo_num_per_kpt + 1
        kp1   =   (p1-1)/mo_num_per_kpt + 1
        kp2   =   (p2-1)/mo_num_per_kpt + 1
        iputi = mod(puti-1,mo_num_per_kpt) + 1
        ihfix = mod(hfix-1,mo_num_per_kpt) + 1
        ip1   = mod(p1-1,  mo_num_per_kpt) + 1
        ip2   = mod(p2-1,  mo_num_per_kpt) + 1
        call get_mo_two_e_integrals_complex(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
        call get_mo_two_e_integrals_complex(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
        call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p1,ip1,kp1,p2,ip2,kp2,mo_num_per_kpt,hij_cache2(1,1),mo_integrals_map,mo_integrals_map_2)
        call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p2,ip2,kp2,p1,ip1,kp1,mo_num_per_kpt,hij_cache2(1,2),mo_integrals_map,mo_integrals_map_2)
        tmp_row = (0.d0,0.d0)
        !tmp_row_kpts = (0.d0,0.d0)
        tmp_row_kpts2 = (0.d0,0.d0)
        !===================
        !begin ref
        do putj=1,hfix-1
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache(putj,1) - hij_cache(putj,2)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
            tmp_row(:,putj) = tmp_row(:,putj) + hij * coefs(:)
          endif
        end do
        do putj=hfix+1,mo_num
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache(putj,2) - hij_cache(putj,1)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
            tmp_row(:,putj) = tmp_row(:,putj) + hij * coefs(:)
          endif
        end do
        !end ref
        !=================
        !begin kpts
        kputj = kconserv(kp1,kp2,khfix)
        putj0 = (kputj-1)*mo_num_per_kpt
        do putj = putj0+1,hfix-1
          iputj = putj - putj0
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache2(iputj,1) - hij_cache2(iputj,2)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
            !tmp_row_kpts(:,putj) = tmp_row_kpts(:,putj) + hij * coefs(:)
            tmp_row_kpts2(:,iputj) = tmp_row_kpts2(:,iputj) + hij * coefs(:)
          endif
        end do
        do putj=hfix+1,putj0+mo_num_per_kpt
          iputj = putj - putj0
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache2(iputj,2) - hij_cache2(iputj,1)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
            !tmp_row_kpts(:,putj) = tmp_row_kpts(:,putj) + hij * coefs(:)
            tmp_row_kpts2(:,iputj) = tmp_row_kpts2(:,iputj) + hij * coefs(:)
          endif
        end do

        !end kpts
    !do ii0=1,mo_num
    !  if (cdabs(tmp_row_kpts(1,ii0)-tmp_row(1,ii0)).gt.1.d-12) then
    !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1a, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
    !!  else if ((cdabs(tmp_row_kpts(1,ii0))+cdabs(tmp_row(1,ii0))).gt.1.d-12) then
    !!    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1b, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
    !  endif
    !enddo
        !=================
        tmp_mat1 = (0.d0,0.d0)
        tmp_mat2 = (0.d0,0.d0)
        tmp_mat1(:, :puti-1, puti) = tmp_mat1(:, :puti-1, puti) + tmp_row(:,:puti-1)
        do l=puti,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_mat1(k, puti, l) = tmp_mat1(k, puti,l) + tmp_row(k,l)
          enddo
        enddo
        !=================
        if (kputj.lt.kputi) then
          tmp_mat2(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) =  &
                  tmp_mat2(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) + &
                  tmp_row_kpts2(1:N_states,1:mo_num_per_kpt)
        else if (kputj.gt.kputi) then
          do l=1,mo_num_per_kpt
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              tmp_mat2(k, puti, l+putj0) = tmp_mat2(k, puti,l+putj0) + tmp_row_kpts2(k,l)
            enddo
          enddo
        else !kputj == kputi
          tmp_mat2(1:N_states,putj0+1:puti-1,puti) =  &
                  tmp_mat2(1:N_states,putj0+1:puti-1,puti) + &
                  tmp_row_kpts2(1:N_states,1:iputi-1)
          do l=iputi,mo_num_per_kpt
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              tmp_mat2(k, puti, l+putj0) = tmp_mat2(k, puti,l+putj0) + tmp_row_kpts2(k,l)
            enddo
          enddo
        endif
        !=================
        do k=1,N_states
          do l=1,mo_num
            do ii0=1,mo_num
              if (cdabs(tmp_mat2(k,l,ii0)-tmp_mat1(k,l,ii0)).gt.1.d-12) then
                print'((A),6(I5),2(2(E25.15),2X))','WarNInG 3a, ',k,l,ii0,hfix,p1,p2,tmp_mat2(k,l,ii0),tmp_mat1(k,l,ii0)
            !  else if ((cdabs(tmp_row_kpts(1,ii0))+cdabs(tmp_row(1,ii0))).gt.1.d-12) then
            !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1b, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
              endif
            enddo
          enddo
        enddo

        !=================
        mat(:, :puti-1, puti) = mat(:, :puti-1, puti) + tmp_row(:,:puti-1)
        do l=puti,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, l) = mat(k, puti,l) + tmp_row(k,l)
          enddo
        enddo
        !!=================
        !!todo: check for iputi=1,2
        !if (kputj.lt.kputi) then
        !  mat(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) =  &
        !          mat(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) + &
        !          tmp_row_kpts2(1:N_states,1:mo_num_per_kpt)
        !else if (kputj.gt.kputi) then
        !  do l=1,mo_num_per_kpt
        !    !DIR$ LOOP COUNT AVG(4)
        !    do k=1,N_states
        !      mat(k, puti, l+putj0) = mat(k, puti,l+putj0) + tmp_row_kpts2(k,l)
        !    enddo
        !  enddo
        !else !kputj == kputi
        !  mat(1:N_states,putj0+1:puti-1,puti) =  &
        !          mat(1:N_states,putj0+1:puti-1,puti) + &
        !          tmp_row_kpts2(1:N_states,1:iputi-1)
        !  do l=iputi,mo_num_per_kpt
        !    !DIR$ LOOP COUNT AVG(4)
        !    do k=1,N_states
        !      mat(k, puti, l+putj0) = mat(k, puti,l+putj0) + tmp_row_kpts2(k,l)
        !    enddo
        !  enddo
        !endif
      end do
    else
      hfix = h(1,mi)
      pfix = p(1,mi)
      p1 = p(1,ma)
      p2 = p(2,ma)
      kpfix = (pfix-1)/mo_num_per_kpt + 1
      khfix = (hfix-1)/mo_num_per_kpt + 1
      kp1   =   (p1-1)/mo_num_per_kpt + 1
      kp2   =   (p2-1)/mo_num_per_kpt + 1
      ipfix = mod(pfix-1,mo_num_per_kpt) + 1
      ihfix = mod(hfix-1,mo_num_per_kpt) + 1
      ip1   = mod(p1-1,  mo_num_per_kpt) + 1
      ip2   = mod(p2-1,  mo_num_per_kpt) + 1
      tmp_row = (0.d0,0.d0)
      tmp_row2 = (0.d0,0.d0)
      !tmp_row_kpts = (0.d0,0.d0)
      !tmp_row2_kpts = (0.d0,0.d0)
      call get_mo_two_e_integrals_complex(hfix,p1,pfix,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
      call get_mo_two_e_integrals_complex(hfix,p2,pfix,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
      !call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p1,ip1,kp1,pfix,ipfix,kpfix,mo_num_per_kpt,hij_cache2(1,1),mo_integrals_map,mo_integrals_map_2)
      !call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p2,ip2,kp2,pfix,ipfix,kpfix,mo_num_per_kpt,hij_cache2(1,2),mo_integrals_map,mo_integrals_map_2)
      !===============
      !begin ref
      putj = p2
      do puti=1,mo_num
        if(lbanned(puti,ma)) cycle
        putj = p2
        if(.not. banned(puti,putj,1)) then
          hij = hij_cache(puti,1)
          if (hij /= (0.d0,0.d0)) then
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
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
            do k=1,N_states
              tmp_row2(k,puti) = tmp_row2(k,puti) + hij * coefs(k)
            enddo
          endif
        end if
      end do
      !end ref
      !===============
      !begin kpts
      !todo: combine if kp1==kp2
 !     putj = p2
 !     kputi1 = kconserv(kp1,kpfix,khfix)
 !     puti01 = (kputi1-1)*mo_num_per_kpt
 !     do iputi=1,mo_num_per_kpt
 !       puti = puti01 + iputi
 !       if(lbanned(puti,ma)) cycle
 !       if(.not. banned(puti,putj,1)) then
 !         hij = hij_cache2(iputi,1)
 !         if (hij /= (0.d0,0.d0)) then
 !           hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
 !           !DIR$ LOOP COUNT AVG(4)
 !           do k=1,N_states
 !             tmp_row_kpts(k,puti) = tmp_row_kpts(k,puti) + hij * coefs(k)
 !           enddo
 !         endif
 !       end if
 !     enddo
 !     putj = p1
 !     kputi2 = kconserv(kp2,kpfix,khfix)
 !     puti02 = (kputi2-1)*mo_num_per_kpt
 !     do iputi=1,mo_num_per_kpt
 !       puti = puti02 + iputi
 !       if(lbanned(puti,ma)) cycle
 !       if(.not. banned(puti,putj,1)) then
 !         hij = hij_cache2(iputi,2)
 !         if (hij /= (0.d0,0.d0)) then
 !           hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
 !           do k=1,N_states
 !             tmp_row2_kpts(k,puti) = tmp_row2_kpts(k,puti) + hij * coefs(k)
 !           enddo
 !         endif
 !       end if
 !     end do
 !     !end kpts
 !     !===============
 !   !test printing
 !   !print'((A),5(I5))','kpt info1: ',kconserv(kpfix,kp2,khfix),khfix,kpfix,kp2,kputi2
 !   !print'((A),5(I5))','kpt info2: ',kconserv(kpfix,kp1,khfix),khfix,kpfix,kp1,kputi1
 !   do ii0=1,mo_num
 !     if (cdabs(tmp_row_kpts(1,ii0)-tmp_row(1,ii0)).gt.1.d-12) then
 !       print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1a, ',ii0,hfix,p1,pfix,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
 !   !  else if ((cdabs(tmp_row_kpts(1,ii0))+cdabs(tmp_row(1,ii0))).gt.1.d-12) then
 !   !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1b, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
 !     endif
 !     if (cdabs(tmp_row2_kpts(1,ii0)-tmp_row2(1,ii0)).gt.1.d-12) then
 !       print'((A),4(I5),2(2(E25.15),2X))','WarNInG 2a, ',ii0,hfix,p2,pfix,tmp_row2_kpts(1,ii0),tmp_row2(1,ii0)
 !   !  else if ((cdabs(tmp_row2_kpts(1,ii0))+cdabs(tmp_row2(1,ii0))).gt.1.d-12) then
 !   !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 2b, ',ii0,hfix,pfix,p1,tmp_row2_kpts(1,ii0),tmp_row2(1,ii0)
 !     endif
 !   enddo
    !===================
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
        ! gen is a selector; mask is ionized generator; det is alpha
        ! hij is contribution to <psi|H|alpha>
        call i_h_j_complex(gen, det, N_int, hij)
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
          mat(k, p1, p2) = mat(k, p1, p2) + coefs(k) * dconjg(hij)
        enddo
      end do
    end do
end




subroutine get_d0_complex(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  !todo: indices/conjg should be okay for complex
  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: gen(N_int, 2), mask(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind) :: det(N_int, 2)
  complex*16, intent(in) :: coefs(N_states)
  complex*16, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  integer :: i, j, k, s, h1, h2, p1, p2, puti, putj
  double precision :: phase
  complex*16 :: hij
  double precision, external :: get_phase_bi
  complex*16, external :: mo_two_e_integral_complex
  logical :: ok

  integer, parameter :: bant=1
  complex*16, allocatable :: hij_cache1(:), hij_cache2(:)
  allocate (hij_cache1(mo_num),hij_cache2(mo_num))


  if(sp == 3) then ! AB
    h1 = p(1,1)
    h2 = p(1,2)
    do p1=1, mo_num
      if(bannedOrb(p1, 1)) cycle
      call get_mo_two_e_integrals_complex(p1,h2,h1,mo_num,hij_cache1,mo_integrals_map,mo_integrals_map_2)
      do p2=1, mo_num
        if(bannedOrb(p2,2)) cycle
        if(banned(p1, p2, bant)) cycle ! rentable?
        if(p1 == h1 .or. p2 == h2) then
          call apply_particles(mask, 1,p1,2,p2, det, ok, N_int)
          ! call i_h_j_complex(gen, det, N_int, hij) ! need to take conjugate of this
          call i_h_j_complex(det, gen, N_int, hij)
        else
          phase = get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
          hij = hij_cache1(p2) * phase
        end if
        if (hij == (0.d0,0.d0)) cycle
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
      call get_mo_two_e_integrals_complex(puti,p2,p1,mo_num,hij_cache1,mo_integrals_map,mo_integrals_map_2)
      call get_mo_two_e_integrals_complex(puti,p1,p2,mo_num,hij_cache2,mo_integrals_map,mo_integrals_map_2)
      do putj=puti+1, mo_num
        if(bannedOrb(putj, sp)) cycle
        if(banned(puti, putj, bant)) cycle ! rentable?
        if(puti == p1 .or. putj == p2 .or. puti == p2 .or. putj == p1) then
          call apply_particles(mask, sp,puti,sp,putj, det, ok, N_int)
          !call i_h_j_complex(gen, det, N_int, hij) ! need to take conjugate of this
          call i_h_j_complex(det, gen, N_int, hij)
          if (hij == (0.d0,0.d0)) cycle
        else
          hij = (mo_two_e_integral_complex(p1, p2, puti, putj) -  mo_two_e_integral_complex(p2, p1, puti, putj))
          if (hij == (0.d0,0.d0)) cycle
          hij = dconjg(hij) * get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
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

