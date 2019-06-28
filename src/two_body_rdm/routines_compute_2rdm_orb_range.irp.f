
 subroutine orb_range_diagonal_contrib_to_two_rdm_ab_dm(det_1,c_1,big_array,dim1,norb,list_orb)
 use bitmasks
 BEGIN_DOC
! routine that update the DIAGONAL PART of the alpha/beta two body rdm in a specific range of orbitals 
! c_1 is supposed to be a scalar quantity, such as state averaged coef 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,norb,list_orb(norb)
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2)
 double precision, intent(in)   :: c_1
 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 do i = 1, n_occ_ab(1)
  h1 = occ(i,1)
  do j = 1, n_occ_ab(2)
   h2 = occ(j,2)
   big_array(h1,h1,h2,h2) += c_1
  enddo 
 enddo
 end


 subroutine orb_range_diagonal_contrib_to_all_two_rdm_dm(det_1,c_1,big_array,dim1,norb,list_orb,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the DIAGONAL PART of ALL THREE two body rdm 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,norb,list_orb(norb),ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2)
 double precision, intent(in)   :: c_1

 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate
 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
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
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 if(alpha_beta)then
  do i = 1, n_occ_ab(1)
   h1 = occ(i,1)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    big_array(h1,h1,h2,h2) += c_1 
   enddo 
  enddo
 else if (alpha_alpha)then
  do i = 1, n_occ_ab(1)
   h1 = occ(i,1)
   do j = 1, n_occ_ab(1)
    h2 = occ(j,1)
    big_array(h1,h1,h2,h2) += 0.5d0 * c_1 
    big_array(h1,h2,h2,h1) -= 0.5d0 * c_1 
   enddo
  enddo
 else if (beta_beta)then
  do i = 1, n_occ_ab(2)
   h1 = occ(i,2)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    big_array(h1,h1,h2,h2) += 0.5d0 * c_1 
    big_array(h1,h2,h2,h1) -= 0.5d0 * c_1 
   enddo
  enddo
 else if(spin_trace)then
  ! 0.5 * (alpha beta + beta alpha)
  do i = 1, n_occ_ab(1)
   h1 = occ(i,1)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    big_array(h1,h1,h2,h2) += 0.5d0 * (c_1 )
    big_array(h2,h2,h1,h1) += 0.5d0 * (c_1 )
   enddo 
  enddo
  ! alpha alpha
  do i = 1, n_occ_ab(1)
   h1 = occ(i,1)
   do j = 1, n_occ_ab(1)
    h2 = occ(j,1)
    big_array(h1,h1,h2,h2) += 0.5d0 * c_1 
    big_array(h1,h2,h2,h1) -= 0.5d0 * c_1 
   enddo
  enddo
  ! beta beta
  do i = 1, n_occ_ab(2)
   h1 = occ(i,2)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    big_array(h1,h1,h2,h2) += 0.5d0 * c_1 
    big_array(h1,h2,h2,h1) -= 0.5d0 * c_1 
   enddo
  enddo
 endif
 end


 subroutine orb_range_off_diagonal_double_to_two_rdm_ab_dm(det_1,det_2,c_1,big_array,dim1,norb,list_orb,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the alpha/beta 2RDM only for DOUBLE EXCITATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,norb,list_orb(norb),ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 double precision, intent(in)   :: c_1
 integer :: i,j,h1,h2,p1,p2,istate
 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
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
 h2 = exc(1,1,2) 
 p1 = exc(1,2,1)
 p2 = exc(1,2,2)
 if(alpha_beta)then
   big_array(h1,p1,h2,p2) += c_1 * phase 
 else if(spin_trace)then
   big_array(h1,p1,h2,p2) += 0.5d0 * c_1 * phase 
   big_array(h2,p2,h1,p1) += 0.5d0 * c_1 * phase 
 endif
 end

 subroutine orb_range_off_diagonal_single_to_two_rdm_ab_dm(det_1,det_2,c_1,big_array,dim1,norb,list_orb,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the alpha/beta 2RDM only for SINGLE EXCITATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,norb,list_orb(norb),ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 double precision, intent(in)   :: c_1

 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate,p1
 integer                        :: exc(0:2,2,2)
 double precision               :: phase

 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
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
   p1 = exc(1,2,1)
    do i = 1, n_occ_ab(2)
     h2 = occ(i,2)
     big_array(h1,p1,h2,h2) += c_1 * phase
    enddo 
  else 
   ! Mono beta
   h1 = exc(1,1,2)
   p1 = exc(1,2,2)
    do i = 1, n_occ_ab(1)
     h2 = occ(i,1)
     big_array(h2,h2,h1,p1) += c_1 * phase
    enddo 
  endif
 else if(spin_trace)then
  if (exc(0,1,1) == 1) then
   ! Mono alpha
   h1 = exc(1,1,1)
   p1 = exc(1,2,1)
    do i = 1, n_occ_ab(2)
     h2 = occ(i,2)
     big_array(h1,p1,h2,h2) +=  0.5d0 * c_1 * phase
     big_array(h2,h2,h1,p1) +=  0.5d0 * c_1 * phase
    enddo 
  else 
   ! Mono beta
   h1 = exc(1,1,2)
   p1 = exc(1,2,2)
    do i = 1, n_occ_ab(1)
     h2 = occ(i,1)
     big_array(h1,p1,h2,h2) +=  0.5d0 * c_1 * phase
     big_array(h2,h2,h1,p1) +=  0.5d0 * c_1 * phase
    enddo 
  endif
 endif
 end

 subroutine orb_range_off_diagonal_single_to_two_rdm_aa_dm(det_1,det_2,c_1,big_array,dim1,norb,list_orb,ispin)
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the alpha/alpha 2RDM only for SINGLE EXCITATIONS 
 END_DOC
 use bitmasks
 implicit none
 integer, intent(in) :: dim1,norb,list_orb(norb),ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 double precision, intent(in)   :: c_1

 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate,p1
 integer                        :: exc(0:2,2,2)
 double precision               :: phase

 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
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
   p1 = exc(1,2,1)
    do i = 1, n_occ_ab(1)
     h2 = occ(i,1)
     big_array(h1,p1,h2,h2) += 0.5d0 * c_1 * phase
     big_array(h1,h2,h2,p1) -= 0.5d0 * c_1 * phase
 
     big_array(h2,h2,h1,p1) += 0.5d0 * c_1 * phase
     big_array(h2,p1,h1,h2) -= 0.5d0 * c_1 * phase
    enddo 
  else 
   return
  endif
 endif
 end

 subroutine orb_range_off_diagonal_single_to_two_rdm_bb_dm(det_1,det_2,c_1,big_array,dim1,norb,list_orb,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the beta /beta  2RDM only for SINGLE EXCITATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,norb,list_orb(norb),ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 double precision, intent(in)   :: c_1


 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate,p1
 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
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
   p1 =  exc(1,2,2)
   do istate = 1, N_states
    do i = 1, n_occ_ab(2)
     h2 = occ(i,2)
     big_array(h1,p1,h2,h2) += 0.5d0 * c_1 * phase
     big_array(h1,h2,h2,p1) -= 0.5d0 * c_1 * phase
 
     big_array(h2,h2,h1,p1) += 0.5d0 * c_1 * phase
     big_array(h2,p1,h1,h2) -= 0.5d0 * c_1 * phase
     enddo 
   enddo
  endif
 endif
 end


 subroutine orb_range_off_diagonal_double_to_two_rdm_aa_dm(det_1,det_2,c_1,big_array,dim1,norb,list_orb,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the alpha/alpha 2RDM only for DOUBLE EXCITATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,norb,list_orb(norb),ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int),det_2(N_int)
 double precision, intent(in)   :: c_1

 integer :: i,j,h1,h2,p1,p2,istate
 integer                        :: exc(0:2,2)
 double precision               :: phase

 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
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
 h2 =exc(2,1)
 p1 =exc(1,2)
 p2 =exc(2,2)
 if(alpha_alpha.or.spin_trace)then
  do istate = 1, N_states
   big_array(h1,p1,h2,p2) += 0.5d0 * c_1 * phase
   big_array(h1,p2,h2,p1) -= 0.5d0 * c_1 * phase
                                         
   big_array(h2,p2,h1,p1) += 0.5d0 * c_1 * phase
   big_array(h2,p1,h1,p2) -= 0.5d0 * c_1 * phase
  enddo
 endif
 end

 subroutine orb_range_off_diagonal_double_to_two_rdm_bb_dm(det_1,det_2,c_1,big_array,dim1,norb,list_orb,ispin)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the beta /beta  2RDM only for DOUBLE EXCITATIONS 
 END_DOC
 implicit none

 integer, intent(in) :: dim1,norb,list_orb(norb),ispin
 double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1)
 integer(bit_kind), intent(in)  :: det_1(N_int),det_2(N_int)
 double precision, intent(in)   :: c_1

 integer :: i,j,h1,h2,p1,p2,istate
 integer                        :: exc(0:2,2)
 double precision               :: phase
 logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
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
 h2 =exc(2,1)
 p1 =exc(1,2)
 p2 =exc(2,2)
 if(beta_beta.or.spin_trace)then
  big_array(h1,p1,h2,p2) += 0.5d0 * c_1* phase 
  big_array(h1,p2,h2,p1) -= 0.5d0 * c_1* phase 

  big_array(h2,p2,h1,p1) += 0.5d0 * c_1* phase 
  big_array(h2,p1,h1,p2) -= 0.5d0 * c_1* phase 
 endif
 end

