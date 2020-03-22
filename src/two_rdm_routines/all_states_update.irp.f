 subroutine orb_range_diagonal_contrib_to_all_2_rdm_dm_all_states(det_1,c_1,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the DIAGONAL PART of the two body rdms in a specific range of orbitals for a given determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1,N_st) is the two-body rdm to be updated in physicist notation 
!
! orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
!
! ispin determines which spin-spin component of the two-rdm you will update 
!
! ispin   == 1 :: alpha/ alpha
! ispin   == 2 :: beta / beta
! ispin   == 3 :: alpha/ beta
! ispin   == 4 :: spin traced <=> total two-rdm 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,N_st,ispin
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1,N_st)
 integer(bit_kind), intent(in)  :: det_1(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 double precision, intent(in)   :: c_1(N_st)

 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate
 integer(bit_kind) :: det_1_act(N_int,2)
 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
 do i = 1, N_int
  det_1_act(i,1) =   iand(det_1(i,1),orb_bitmask(i)) 
  det_1_act(i,2) =   iand(det_1(i,2),orb_bitmask(i)) 
 enddo

 alpha_alpha = .False.   
 beta_beta   = .False.
 alpha_beta  = .False.
 spin_trace  = .False.
 if(     ispin == 1)then
  alpha_alpha = .True.
 else if(ispin == 2)then
  beta_beta   = .True.
 else if(ispin == 3)then
  alpha_beta  = .True.
 else if(ispin == 4)then
  spin_trace  = .True.
 endif
 call bitstring_to_list_ab(det_1_act, occ, n_occ_ab, N_int)
 logical :: is_integer_in_string
 integer :: i1,i2
 if(alpha_beta)then
  do istate = 1, N_st
   do i = 1, n_occ_ab(1)
    i1 = occ(i,1)
    do j = 1, n_occ_ab(2)
     i2 = occ(j,2)
     h1 = list_orb_reverse(i1)
     h2 = list_orb_reverse(i2)
     ! If alpha/beta, electron 1 is alpha, electron 2 is beta
     ! Therefore you don't necessayr have symmetry between electron 1 and 2 
     big_array(h1,h2,h1,h2,istate) += 1.0d0 * c_1(istate) 
    enddo 
   enddo
  enddo
 else if (alpha_alpha)then
  do istate = 1, N_st
   do i = 1, n_occ_ab(1)
    i1 = occ(i,1)
    do j = 1, n_occ_ab(1)
     i2 = occ(j,1)
     h1 = list_orb_reverse(i1)
     h2 = list_orb_reverse(i2)
     big_array(h1,h2,h1,h2,istate) += 0.5d0 * c_1(istate) 
     big_array(h1,h2,h2,h1,istate) -= 0.5d0 * c_1(istate) 
    enddo
   enddo
  enddo
!  pause
 else if (beta_beta)then
  do istate = 1, N_st
   do i = 1, n_occ_ab(2)
    i1 = occ(i,2)
    do j = 1, n_occ_ab(2)
     i2 = occ(j,2)
     h1 = list_orb_reverse(i1)
     h2 = list_orb_reverse(i2)
     big_array(h1,h2,h1,h2,istate) += 0.5d0 * c_1(istate) 
     big_array(h1,h2,h2,h1,istate) -= 0.5d0 * c_1(istate) 
    enddo
   enddo
  enddo
 else if(spin_trace)then
  ! 0.5 * (alpha beta + beta alpha)
  do istate = 1, N_st
   do i = 1, n_occ_ab(1)
    i1 = occ(i,1)
    do j = 1, n_occ_ab(2)
     i2 = occ(j,2)
     h1 = list_orb_reverse(i1)
     h2 = list_orb_reverse(i2)
     big_array(h1,h2,h1,h2,istate) += 0.5d0 * c_1(istate)
     big_array(h2,h1,h2,h1,istate) += 0.5d0 * c_1(istate)
    enddo 
   enddo
   do i = 1, n_occ_ab(1)
    i1 = occ(i,1)
    do j = 1, n_occ_ab(1)
     i2 = occ(j,1)
     h1 = list_orb_reverse(i1)
     h2 = list_orb_reverse(i2)
     big_array(h1,h2,h1,h2,istate) += 0.5d0 * c_1(istate) 
     big_array(h1,h2,h2,h1,istate) -= 0.5d0 * c_1(istate) 
    enddo
   enddo
   do i = 1, n_occ_ab(2)
    i1 = occ(i,2)
    do j = 1, n_occ_ab(2)
     i2 = occ(j,2)
     h1 = list_orb_reverse(i1)
     h2 = list_orb_reverse(i2)
     big_array(h1,h2,h1,h2,istate) += 0.5d0 * c_1(istate) 
     big_array(h1,h2,h2,h1,istate) -= 0.5d0 * c_1(istate) 
    enddo
   enddo
  enddo
 endif
 end


 subroutine orb_range_off_diagonal_double_to_2_rdm_ab_dm_all_states(det_1,det_2,c_1,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a alpha/beta DOUBLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1,N_st) is the two-body rdm to be updated in physicist notation 
!
! orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
!
! ispin determines which spin-spin component of the two-rdm you will update 
!
! ispin   == 1 :: alpha/ alpha
! ispin   == 2 :: beta / beta
! ispin   == 3 :: alpha/ beta
! ispin   == 4 :: spin traced <=> total two-rdm 
!
! here, only ispin == 3 or 4 will do something
 END_DOC
 implicit none
 integer, intent(in) :: dim1,N_st,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1,N_st)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1(N_st)
 integer :: i,j,h1,h2,p1,p2,istate
 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
 logical :: is_integer_in_string
 alpha_alpha = .False.   
 beta_beta   = .False.
 alpha_beta  = .False.
 spin_trace  = .False.
 if(     ispin == 1)then
  alpha_alpha = .True.
 else if(ispin == 2)then
  beta_beta   = .True.
 else if(ispin == 3)then
  alpha_beta  = .True.
 else if(ispin == 4)then
  spin_trace  = .True.
 endif
 call get_double_excitation(det_1,det_2,exc,phase,N_int)
 h1 = exc(1,1,1) 
 if(.not.is_integer_in_string(h1,orb_bitmask,N_int))return
 h1 = list_orb_reverse(h1)
 h2 = exc(1,1,2) 
 if(.not.is_integer_in_string(h2,orb_bitmask,N_int))return
 h2 = list_orb_reverse(h2)
 p1 = exc(1,2,1)
 if(.not.is_integer_in_string(p1,orb_bitmask,N_int))return
 p1 = list_orb_reverse(p1)
 p2 = exc(1,2,2)
 if(.not.is_integer_in_string(p2,orb_bitmask,N_int))return
 p2 = list_orb_reverse(p2)
 do istate = 1, N_st
  if(alpha_beta)then
    big_array(h1,h2,p1,p2,istate) += c_1(istate) * phase 
  else if(spin_trace)then
    big_array(h1,h2,p1,p2,istate) += 0.5d0 * c_1(istate) * phase 
    big_array(p1,p2,h1,h2,istate) += 0.5d0 * c_1(istate) * phase 
  endif
 enddo
 end

 subroutine orb_range_off_diagonal_single_to_2_rdm_ab_dm_all_states(det_1,det_2,c_1,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a SINGLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1,N_st) is the two-body rdm to be updated in physicist notation 
!
! orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
!
! ispin determines which spin-spin component of the two-rdm you will update 
!
! ispin   == 1 :: alpha/ alpha
! ispin   == 2 :: beta / beta
! ispin   == 3 :: alpha/ beta
! ispin   == 4 :: spin traced <=> total two-rdm 
!
! here, only ispin == 3 or 4 will do something
 END_DOC
 implicit none
 integer, intent(in) :: dim1,N_st,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1,N_st)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1(N_st)

 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate,p1
 integer                        :: exc(0:2,2,2)
 double precision               :: phase

 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
 logical :: is_integer_in_string
 alpha_alpha = .False.   
 beta_beta   = .False.
 alpha_beta  = .False.
 spin_trace  = .False.
 if(     ispin == 1)then
  alpha_alpha = .True.
 else if(ispin == 2)then
  beta_beta   = .True.
 else if(ispin == 3)then
  alpha_beta  = .True.
 else if(ispin == 4)then
  spin_trace  = .True.
 endif

 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 call get_single_excitation(det_1,det_2,exc,phase,N_int)
 if(alpha_beta)then
  do istate = 1, N_st
   if (exc(0,1,1) == 1) then
    ! Mono alpha
    h1 = exc(1,1,1)
    if(.not.is_integer_in_string(h1,orb_bitmask,N_int))return
    h1 = list_orb_reverse(h1)
    p1 = exc(1,2,1)
    if(.not.is_integer_in_string(p1,orb_bitmask,N_int))return
    p1 = list_orb_reverse(p1)
     do i = 1, n_occ_ab(2)
      h2 = occ(i,2)
      if(.not.is_integer_in_string(h2,orb_bitmask,N_int))cycle
      h2 = list_orb_reverse(h2)
      big_array(h1,h2,p1,h2,istate) += c_1(istate) * phase
     enddo 
   else 
    ! Mono beta
    h1 = exc(1,1,2)
    if(.not.is_integer_in_string(h1,orb_bitmask,N_int))return
    h1 = list_orb_reverse(h1)
    p1 = exc(1,2,2)
    if(.not.is_integer_in_string(p1,orb_bitmask,N_int))return
    p1 = list_orb_reverse(p1)
     do i = 1, n_occ_ab(1)
      h2 = occ(i,1)
      if(.not.is_integer_in_string(h2,orb_bitmask,N_int))cycle
      h2 = list_orb_reverse(h2)
      big_array(h2,h1,h2,p1,istate) += c_1(istate) * phase
     enddo 
   endif
  enddo
 else if(spin_trace)then
  if (exc(0,1,1) == 1) then
   ! Mono alpha
   h1 = exc(1,1,1)
   if(.not.is_integer_in_string(h1,orb_bitmask,N_int))return
   h1 = list_orb_reverse(h1)
   p1 = exc(1,2,1)
   if(.not.is_integer_in_string(p1,orb_bitmask,N_int))return
   p1 = list_orb_reverse(p1)
   do istate = 1, N_st
    do i = 1, n_occ_ab(2)
     h2 = occ(i,2)
     if(.not.is_integer_in_string(h2,orb_bitmask,N_int))cycle
     h2 = list_orb_reverse(h2)
     big_array(h1,h2,p1,h2,istate) +=  0.5d0 * c_1(istate) * phase
     big_array(h2,h1,h2,p1,istate) +=  0.5d0 * c_1(istate) * phase
    enddo 
   enddo
  else 
   ! Mono beta
   h1 = exc(1,1,2)
   if(.not.is_integer_in_string(h1,orb_bitmask,N_int))return
   h1 = list_orb_reverse(h1)
   p1 = exc(1,2,2)
   if(.not.is_integer_in_string(p1,orb_bitmask,N_int))return
   p1 = list_orb_reverse(p1)
   do istate = 1, N_st
    do i = 1, n_occ_ab(1)
     h2 = occ(i,1)
     if(.not.is_integer_in_string(h2,orb_bitmask,N_int))cycle
     h2 = list_orb_reverse(h2)
     big_array(h1,h2,p1,h2,istate) +=  0.5d0 * c_1(istate) * phase
     big_array(h2,h1,h2,p1,istate) +=  0.5d0 * c_1(istate) * phase
    enddo 
   enddo
  endif
 endif
 end

 subroutine orb_range_off_diagonal_single_to_2_rdm_aa_dm_all_states(det_1,det_2,c_1,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a ALPHA SINGLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1,N_st) is the two-body rdm to be updated in physicist notation 
!
! orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
!
! ispin determines which spin-spin component of the two-rdm you will update 
!
! ispin   == 1 :: alpha/ alpha
! ispin   == 2 :: beta / beta
! ispin   == 3 :: alpha/ beta
! ispin   == 4 :: spin traced <=> total two-rdm 
!
! here, only ispin == 1 or 4 will do something
 END_DOC
 use bitmasks
 implicit none
 integer, intent(in) :: dim1,N_st,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1,N_st)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1(N_st)

 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate,p1
 integer                        :: exc(0:2,2,2)
 double precision               :: phase

 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
 logical :: is_integer_in_string
 alpha_alpha = .False.   
 beta_beta   = .False.
 alpha_beta  = .False.
 spin_trace  = .False.
 if(     ispin == 1)then
  alpha_alpha = .True.
 else if(ispin == 2)then
  beta_beta   = .True.
 else if(ispin == 3)then
  alpha_beta  = .True.
 else if(ispin == 4)then
  spin_trace  = .True.
 endif

 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 call get_single_excitation(det_1,det_2,exc,phase,N_int)
 if(alpha_alpha.or.spin_trace)then
  if (exc(0,1,1) == 1) then
   ! Mono alpha
   h1 = exc(1,1,1)
   if(.not.is_integer_in_string(h1,orb_bitmask,N_int))return
   h1 = list_orb_reverse(h1)
   p1 = exc(1,2,1)
   if(.not.is_integer_in_string(p1,orb_bitmask,N_int))return
   p1 = list_orb_reverse(p1)
   do istate = 1, N_st
    do i = 1, n_occ_ab(1)
     h2 = occ(i,1)
     if(.not.is_integer_in_string(h2,orb_bitmask,N_int))cycle
     h2 = list_orb_reverse(h2)
     big_array(h1,h2,p1,h2,istate) += 0.5d0 * c_1(istate) * phase
     big_array(h1,h2,h2,p1,istate) -= 0.5d0 * c_1(istate) * phase
 
     big_array(h2,h1,h2,p1,istate) += 0.5d0 * c_1(istate) * phase
     big_array(h2,h1,p1,h2,istate) -= 0.5d0 * c_1(istate) * phase
    enddo 
   enddo
  else 
   return
  endif
 endif
 end

 subroutine orb_range_off_diagonal_single_to_2_rdm_bb_dm_all_states(det_1,det_2,c_1,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a BETA  SINGLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1,N_st) is the two-body rdm to be updated in physicist notation 
!
! orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
!
! ispin determines which spin-spin component of the two-rdm you will update 
!
! ispin   == 1 :: alpha/ alpha
! ispin   == 2 :: beta / beta
! ispin   == 3 :: alpha/ beta
! ispin   == 4 :: spin traced <=> total two-rdm 
!
! here, only ispin == 2 or 4 will do something
 END_DOC
 implicit none
 integer, intent(in) :: dim1,N_st,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1,N_st)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1(N_st)


 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate,p1
 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
 logical :: is_integer_in_string
 alpha_alpha = .False.   
 beta_beta   = .False.
 alpha_beta  = .False.
 spin_trace  = .False.
 if(     ispin == 1)then
  alpha_alpha = .True.
 else if(ispin == 2)then
  beta_beta   = .True.
 else if(ispin == 3)then
  alpha_beta  = .True.
 else if(ispin == 4)then
  spin_trace  = .True.
 endif


 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 call get_single_excitation(det_1,det_2,exc,phase,N_int)
 if(beta_beta.or.spin_trace)then
  if (exc(0,1,1) == 1) then
   return
  else
   ! Mono beta
   h1 =  exc(1,1,2)
   if(.not.is_integer_in_string(h1,orb_bitmask,N_int))return
   h1 = list_orb_reverse(h1)
   p1 =  exc(1,2,2)
   if(.not.is_integer_in_string(p1,orb_bitmask,N_int))return
   p1 = list_orb_reverse(p1)
   do istate = 1, N_st
    do i = 1, n_occ_ab(2)
     h2 = occ(i,2)
     if(.not.is_integer_in_string(h2,orb_bitmask,N_int))cycle
     h2 = list_orb_reverse(h2)
     big_array(h1,h2,p1,h2,istate) += 0.5d0 * c_1(istate) * phase
     big_array(h1,h2,h2,p1,istate) -= 0.5d0 * c_1(istate) * phase
 
     big_array(h2,h1,h2,p1,istate) += 0.5d0 * c_1(istate) * phase
     big_array(h2,h1,p1,h2,istate) -= 0.5d0 * c_1(istate) * phase
     enddo 
   enddo
  endif
 endif
 end


 subroutine orb_range_off_diagonal_double_to_2_rdm_aa_dm_all_states(det_1,det_2,c_1,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a ALPHA/ALPHA DOUBLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1,N_st) is the two-body rdm to be updated in physicist notation 
!
! orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
!
! ispin determines which spin-spin component of the two-rdm you will update 
!
! ispin   == 1 :: alpha/ alpha
! ispin   == 2 :: beta / beta
! ispin   == 3 :: alpha/ beta
! ispin   == 4 :: spin traced <=> total two-rdm 
!
! here, only ispin == 1 or 4 will do something
 END_DOC
 implicit none
 integer, intent(in) :: dim1,N_st,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1,N_st)
 integer(bit_kind), intent(in)  :: det_1(N_int),det_2(N_int)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1(N_st)

 integer :: i,j,h1,h2,p1,p2,istate
 integer                        :: exc(0:2,2)
 double precision               :: phase

 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
 logical :: is_integer_in_string
 alpha_alpha = .False.   
 beta_beta   = .False.
 alpha_beta  = .False.
 spin_trace  = .False.
 if(     ispin == 1)then
  alpha_alpha = .True.
 else if(ispin == 2)then
  beta_beta   = .True.
 else if(ispin == 3)then
  alpha_beta  = .True.
 else if(ispin == 4)then
  spin_trace  = .True.
 endif
 call get_double_excitation_spin(det_1,det_2,exc,phase,N_int)
 h1 =exc(1,1)
 if(.not.is_integer_in_string(h1,orb_bitmask,N_int))return
 h1 = list_orb_reverse(h1)
 h2 =exc(2,1)
 if(.not.is_integer_in_string(h2,orb_bitmask,N_int))return
 h2 = list_orb_reverse(h2)
 p1 =exc(1,2)
 if(.not.is_integer_in_string(p1,orb_bitmask,N_int))return
 p1 = list_orb_reverse(p1)
 p2 =exc(2,2)
 if(.not.is_integer_in_string(p2,orb_bitmask,N_int))return
 p2 = list_orb_reverse(p2)
 if(alpha_alpha.or.spin_trace)then
  do istate = 1, N_st
   big_array(h1,h2,p1,p2,istate) += 0.5d0 * c_1(istate) * phase
   big_array(h1,h2,p2,p1,istate) -= 0.5d0 * c_1(istate) * phase
                                         
   big_array(h2,h1,p2,p1,istate) += 0.5d0 * c_1(istate) * phase
   big_array(h2,h1,p1,p2,istate) -= 0.5d0 * c_1(istate) * phase
  enddo
 endif
 end

 subroutine orb_range_off_diagonal_double_to_2_rdm_bb_dm_all_states(det_1,det_2,c_1,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a BETA /BETA  DOUBLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1,N_st) is the two-body rdm to be updated in physicist notation 
!
! orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
!
! ispin determines which spin-spin component of the two-rdm you will update 
!
! ispin   == 1 :: alpha/ alpha
! ispin   == 2 :: beta / beta
! ispin   == 3 :: alpha/ beta
! ispin   == 4 :: spin traced <=> total two-rdm 
!
! here, only ispin == 2 or 4 will do something
 END_DOC
 implicit none

 integer, intent(in) :: dim1,N_st,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1,N_st)
 integer(bit_kind), intent(in)  :: det_1(N_int),det_2(N_int)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1(N_st)

 integer :: i,j,h1,h2,p1,p2,istate
 integer                        :: exc(0:2,2)
 double precision               :: phase
 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
 logical :: is_integer_in_string
 alpha_alpha = .False.   
 beta_beta   = .False.
 alpha_beta  = .False.
 spin_trace  = .False.
 if(     ispin == 1)then
  alpha_alpha = .True.
 else if(ispin == 2)then
  beta_beta   = .True.
 else if(ispin == 3)then
  alpha_beta  = .True.
 else if(ispin == 4)then
  spin_trace  = .True.
 endif

 call get_double_excitation_spin(det_1,det_2,exc,phase,N_int)
 h1 =exc(1,1)
 if(.not.is_integer_in_string(h1,orb_bitmask,N_int))return
 h1 = list_orb_reverse(h1)
 h2 =exc(2,1)
 if(.not.is_integer_in_string(h2,orb_bitmask,N_int))return
 h2 = list_orb_reverse(h2)
 p1 =exc(1,2)
 if(.not.is_integer_in_string(p1,orb_bitmask,N_int))return
 p1 = list_orb_reverse(p1)
 p2 =exc(2,2)
 if(.not.is_integer_in_string(p2,orb_bitmask,N_int))return
 p2 = list_orb_reverse(p2)
 do istate = 1, N_st
  if(beta_beta.or.spin_trace)then
   big_array(h1,h2,p1,p2,istate) += 0.5d0 * c_1(istate)* phase 
   big_array(h1,h2,p2,p1,istate) -= 0.5d0 * c_1(istate)* phase 
 
   big_array(h2,h1,p2,p1,istate) += 0.5d0 * c_1(istate)* phase 
   big_array(h2,h1,p1,p2,istate) -= 0.5d0 * c_1(istate)* phase 
  endif
 enddo
 end

