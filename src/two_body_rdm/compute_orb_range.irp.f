
 subroutine orb_range_diagonal_contrib_to_two_rdm_ab_dm(det_1,c_1,big_array,dim1,orb_bitmask)
 use bitmasks
 BEGIN_DOC
! routine that update the DIAGONAL PART of the alpha/beta two body rdm in a specific range of orbitals 
! c_1 is supposed to be a scalar quantity, such as state averaged coef 
 END_DOC
 implicit none
 integer, intent(in) :: dim1
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 double precision, intent(in)   :: c_1
 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 do i = 1, n_occ_ab(1)
  h1 = occ(i,1)
  do j = 1, n_occ_ab(2)
   h2 = occ(j,2)
   big_array(h1,h2,h1,h2) += c_1
  enddo 
 enddo
 end


 subroutine orb_range_diagonal_contrib_to_all_two_rdm_dm(det_1,c_1,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the DIAGONAL PART of the two body rdms in a specific range of orbitals for a given determinant det_1
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation 
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
 integer, intent(in) :: dim1,ispin
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 double precision, intent(in)   :: c_1

 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2
 integer(bit_kind) :: det_1_act(N_int,2)
 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
 do i = 1, N_int
  det_1_act(i,1) =   iand(det_1(i,1),orb_bitmask(i)) 
  det_1_act(i,2) =   iand(det_1(i,2),orb_bitmask(i)) 
 enddo
 
!print*,'ahah'
!call debug_det(det_1_act,N_int)
!pause
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
 BEGIN_DOC
! no factor 1/2 have to be taken into account as the permutations are already taken into account
 END_DOC
 call bitstring_to_list_ab(det_1_act, occ, n_occ_ab, N_int)
 logical :: is_integer_in_string
 integer :: i1,i2
 if(alpha_beta)then
  do i = 1, n_occ_ab(1)
   i1 = occ(i,1)
!  if(.not.is_integer_in_string(i1,orb_bitmask,N_int))cycle
   do j = 1, n_occ_ab(2)
!   if(.not.is_integer_in_string(i2,orb_bitmask,N_int))cycle
    i2 = occ(j,2)
    h1 = list_orb_reverse(i1)
    h2 = list_orb_reverse(i2)
    big_array(h1,h2,h1,h2) += c_1 
   enddo 
  enddo
 else if (alpha_alpha)then
  do i = 1, n_occ_ab(1)
   i1 = occ(i,1)
!  if(.not.is_integer_in_string(i1,orb_bitmask,N_int))cycle
   do j = 1, n_occ_ab(1)
    i2 = occ(j,1)
!   if(.not.is_integer_in_string(i2,orb_bitmask,N_int))cycle
    h1 = list_orb_reverse(i1)
    h2 = list_orb_reverse(i2)
    big_array(h1,h2,h1,h2) += 0.5d0 * c_1 
    big_array(h1,h2,h2,h1) -= 0.5d0 * c_1 
   enddo
  enddo
 else if (beta_beta)then
  do i = 1, n_occ_ab(2)
   i1 = occ(i,2)
!  if(.not.is_integer_in_string(i1,orb_bitmask,N_int))cycle
   do j = 1, n_occ_ab(2)
    i2 = occ(j,2)
!   if(.not.is_integer_in_string(i2,orb_bitmask,N_int))cycle
    h1 = list_orb_reverse(i1)
    h2 = list_orb_reverse(i2)
    big_array(h1,h2,h1,h2) += 0.5d0 * c_1 
    big_array(h1,h2,h2,h1) -= 0.5d0 * c_1 
   enddo
  enddo
 else if(spin_trace)then
  ! 0.5 * (alpha beta + beta alpha)
  do i = 1, n_occ_ab(1)
   i1 = occ(i,1)
!  if(.not.is_integer_in_string(i1,orb_bitmask,N_int))cycle
   do j = 1, n_occ_ab(2)
    i2 = occ(j,2)
!   if(.not.is_integer_in_string(i2,orb_bitmask,N_int))cycle
    h1 = list_orb_reverse(i1)
    h2 = list_orb_reverse(i2)
    big_array(h1,h2,h1,h2) += 0.5d0 * (c_1 )
    big_array(h2,h1,h2,h1) += 0.5d0 * (c_1 )
   enddo 
  enddo
 !stop
  do i = 1, n_occ_ab(1)
   i1 = occ(i,1)
!  if(.not.is_integer_in_string(i1,orb_bitmask,N_int))cycle
   do j = 1, n_occ_ab(1)
    i2 = occ(j,1)
!   if(.not.is_integer_in_string(i2,orb_bitmask,N_int))cycle
    h1 = list_orb_reverse(i1)
    h2 = list_orb_reverse(i2)
    big_array(h1,h2,h1,h2) += 0.5d0 * c_1 
    big_array(h1,h2,h2,h1) -= 0.5d0 * c_1 
   enddo
  enddo
  do i = 1, n_occ_ab(2)
   i1 = occ(i,2)
!  if(.not.is_integer_in_string(i1,orb_bitmask,N_int))cycle
   do j = 1, n_occ_ab(2)
    i2 = occ(j,2)
!   if(.not.is_integer_in_string(i2,orb_bitmask,N_int))cycle
    h1 = list_orb_reverse(i1)
    h2 = list_orb_reverse(i2)
    big_array(h1,h2,h1,h2) += 0.5d0 * c_1 
    big_array(h1,h2,h2,h1) -= 0.5d0 * c_1 
   enddo
  enddo
 endif
 end


 subroutine orb_range_off_diagonal_double_to_two_rdm_ab_dm(det_1,det_2,c_1,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a alpha/beta DOUBLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation 
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
 integer, intent(in) :: dim1,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1
 integer :: i,j,h1,h2,p1,p2
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
!print*,''
!do i = 1, mo_num
! print*,'list_orb',i,list_orb_reverse(i) 
!enddo
 call get_double_excitation(det_1,det_2,exc,phase,N_int)
 h1 = exc(1,1,1) 
!print*,'h1',h1
 if(.not.is_integer_in_string(h1,orb_bitmask,N_int))return
 h1 = list_orb_reverse(h1)
!print*,'passed h1 = ',h1
 h2 = exc(1,1,2) 
!print*,'h2',h2
 if(.not.is_integer_in_string(h2,orb_bitmask,N_int))return
 h2 = list_orb_reverse(h2)
!print*,'passed h2 = ',h2
 p1 = exc(1,2,1)
!print*,'p1',p1
 if(.not.is_integer_in_string(p1,orb_bitmask,N_int))return
 p1 = list_orb_reverse(p1)
!print*,'passed p1 = ',p1
 p2 = exc(1,2,2)
!print*,'p2',p2
 if(.not.is_integer_in_string(p2,orb_bitmask,N_int))return
 p2 = list_orb_reverse(p2)
!print*,'passed p2 = ',p2
 if(alpha_beta)then
   big_array(h1,h2,p1,p2) += c_1 * phase 
 else if(spin_trace)then
   big_array(h1,h2,p1,p2) += 0.5d0 * c_1 * phase 
   big_array(p1,p2,h1,h2) += 0.5d0 * c_1 * phase 
  !print*,'h1,h2,p1,p2',h1,h2,p1,p2
  !print*,'',big_array(h1,h2,p1,p2)
 endif
 end

 subroutine orb_range_off_diagonal_single_to_two_rdm_ab_dm(det_1,det_2,c_1,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a SINGLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation 
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
 integer, intent(in) :: dim1,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1

 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,p1
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
     big_array(h1,h2,p1,h2) += c_1 * phase
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
     big_array(h2,h1,h2,p1) += c_1 * phase
    enddo 
  endif
 else if(spin_trace)then
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
     big_array(h1,h2,p1,h2) +=  0.5d0 * c_1 * phase
     big_array(h2,h1,h2,p1) +=  0.5d0 * c_1 * phase
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
     big_array(h1,h2,p1,h2) +=  0.5d0 * c_1 * phase
     big_array(h2,h1,h2,p1) +=  0.5d0 * c_1 * phase
    enddo 
  endif
 endif
 end

 subroutine orb_range_off_diagonal_single_to_two_rdm_aa_dm(det_1,det_2,c_1,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a ALPHA SINGLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation 
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
 integer, intent(in) :: dim1,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1

 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,p1
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
    do i = 1, n_occ_ab(1)
     h2 = occ(i,1)
     if(.not.is_integer_in_string(h2,orb_bitmask,N_int))cycle
     h2 = list_orb_reverse(h2)
     big_array(h1,h2,p1,h2) += 0.5d0 * c_1 * phase
     big_array(h1,h2,h2,p1) -= 0.5d0 * c_1 * phase
 
     big_array(h2,h1,h2,p1) += 0.5d0 * c_1 * phase
     big_array(h2,h1,p1,h2) -= 0.5d0 * c_1 * phase
    enddo 
  else 
   return
  endif
 endif
 end

 subroutine orb_range_off_diagonal_single_to_two_rdm_bb_dm(det_1,det_2,c_1,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a BETA  SINGLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation 
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
 integer, intent(in) :: dim1,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1


 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,p1
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
   do i = 1, n_occ_ab(2)
    h2 = occ(i,2)
    if(.not.is_integer_in_string(h2,orb_bitmask,N_int))cycle
    h2 = list_orb_reverse(h2)
    big_array(h1,h2,p1,h2) += 0.5d0 * c_1 * phase
    big_array(h1,h2,h2,p1) -= 0.5d0 * c_1 * phase
 
    big_array(h2,h1,h2,p1) += 0.5d0 * c_1 * phase
    big_array(h2,h1,p1,h2) -= 0.5d0 * c_1 * phase
   enddo 
  endif
 endif
 end


 subroutine orb_range_off_diagonal_double_to_two_rdm_aa_dm(det_1,det_2,c_1,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a ALPHA/ALPHA DOUBLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation 
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
 integer, intent(in) :: dim1,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int),det_2(N_int)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1

 integer :: i,j,h1,h2,p1,p2
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
   big_array(h1,h2,p1,p2) += 0.5d0 * c_1 * phase
   big_array(h1,h2,p2,p1) -= 0.5d0 * c_1 * phase
                                         
   big_array(h2,h1,p2,p1) += 0.5d0 * c_1 * phase
   big_array(h2,h1,p1,p2) -= 0.5d0 * c_1 * phase
 endif
 end

 subroutine orb_range_off_diagonal_double_to_two_rdm_bb_dm(det_1,det_2,c_1,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for 
!
! a given couple of determinant det_1, det_2 being a BETA /BETA  DOUBLE excitation with respect to one another
! 
! c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
! 
! big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation 
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

 integer, intent(in) :: dim1,ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int),det_2(N_int)
 integer(bit_kind), intent(in)  :: orb_bitmask(N_int)
 integer, intent(in) :: list_orb_reverse(mo_num)
 double precision, intent(in)   :: c_1

 integer :: i,j,h1,h2,p1,p2
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
 if(beta_beta.or.spin_trace)then
  big_array(h1,h2,p1,p2) += 0.5d0 * c_1* phase 
  big_array(h1,h2,p2,p1) -= 0.5d0 * c_1* phase 

  big_array(h2,h1,p2,p1) += 0.5d0 * c_1* phase 
  big_array(h2,h1,p1,p2) -= 0.5d0 * c_1* phase 
 endif
 end

