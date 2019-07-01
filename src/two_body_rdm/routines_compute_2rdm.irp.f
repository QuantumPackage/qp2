

 subroutine diagonal_contrib_to_two_rdm_ab_dm(det_1,c_1,big_array,dim1,dim2,dim3,dim4)
 use bitmasks
 BEGIN_DOC
! routine that update the DIAGONAL PART of the alpha/beta two body rdm IN CHEMIST NOTATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int,2)
 double precision, intent(in)   :: c_1(N_states)
 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate
 double precision               :: c_1_bis
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 do istate = 1, N_states
  c_1_bis = c_1(istate) * c_1(istate)
  do i = 1, n_occ_ab(1)
   h1 = occ(i,1)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    big_array(h1,h1,h2,h2,istate) += c_1_bis 
   enddo 
  enddo
 enddo
 end


 subroutine diagonal_contrib_to_all_two_rdm_dm(det_1,c_1,big_array_aa,big_array_bb,big_array_ab,dim1,dim2,dim3,dim4)
 use bitmasks
 BEGIN_DOC
! routine that update the DIAGONAL PART of ALL THREE two body rdm IN CHEMIST NOTATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array_ab(dim1,dim2,dim3,dim4,N_states)
 double precision, intent(inout) :: big_array_aa(dim1,dim2,dim3,dim4,N_states)
 double precision, intent(inout) :: big_array_bb(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int,2)
 double precision, intent(in)   :: c_1(N_states)
 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate
 double precision               :: c_1_bis
 BEGIN_DOC
! no factor 1/2 have to be taken into account as the permutations are already taken into account
 END_DOC
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 do istate = 1, N_states
  c_1_bis = c_1(istate) * c_1(istate)
  do i = 1, n_occ_ab(1)
   h1 = occ(i,1)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    big_array_ab(h1,h1,h2,h2,istate) += c_1_bis 
   enddo 
   do j = 1, n_occ_ab(1)
    h2 = occ(j,1)
    big_array_aa(h1,h1,h2,h2,istate) += 0.5d0 * c_1_bis 
    big_array_aa(h1,h2,h2,h1,istate) -= 0.5d0 * c_1_bis 
   enddo
  enddo
  do i = 1, n_occ_ab(2)
   h1 = occ(i,2)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    big_array_bb(h1,h1,h2,h2,istate) += 0.5d0 * c_1_bis 
    big_array_bb(h1,h2,h2,h1,istate) -= 0.5d0 * c_1_bis 
   enddo
  enddo
 enddo
 end


 subroutine off_diagonal_double_to_two_rdm_ab_dm(det_1,det_2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the alpha/beta 2RDM only for DOUBLE EXCITATIONS  IN CHEMIST NOTATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 double precision, intent(in)   :: c_1(N_states),c_2(N_states)
 integer :: i,j,h1,h2,p1,p2,istate
 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 call get_double_excitation(det_1,det_2,exc,phase,N_int)
 h1 = exc(1,1,1) 
 h2 = exc(1,1,2) 
 p1 = exc(1,2,1)
 p2 = exc(1,2,2)
 do istate = 1, N_states
  big_array(h1,p1,h2,p2,istate) += c_1(istate) * phase * c_2(istate)
! big_array(p1,h1,p2,h2,istate) += c_1(istate) * phase * c_2(istate)
 enddo
 end

 subroutine off_diagonal_single_to_two_rdm_ab_dm(det_1,det_2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the alpha/beta 2RDM only for SINGLE EXCITATIONS IN CHEMIST NOTATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 double precision, intent(in)   :: c_1(N_states),c_2(N_states)
 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate,p1
 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 call get_single_excitation(det_1,det_2,exc,phase,N_int)
 if (exc(0,1,1) == 1) then
  ! Mono alpha
  h1 = exc(1,1,1)
  p1 = exc(1,2,1)
  do istate = 1, N_states
   do i = 1, n_occ_ab(2)
    h2 = occ(i,2)
    big_array(h1,p1,h2,h2,istate) += 1.d0 * c_1(istate) * c_2(istate) * phase
    enddo 
  enddo
 else 
  ! Mono beta
  h1 = exc(1,1,2)
  p1 = exc(1,2,2)
  do istate = 1, N_states
   do i = 1, n_occ_ab(1)
    h2 = occ(i,1)
    big_array(h2,h2,h1,p1,istate) += 1.d0 * c_1(istate) * c_2(istate) * phase
    enddo 
  enddo
 endif
 end

 subroutine off_diagonal_single_to_two_rdm_aa_dm(det_1,det_2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the alpha/alpha 2RDM only for SINGLE EXCITATIONS IN CHEMIST NOTATIONS 
 END_DOC
 use bitmasks
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 double precision, intent(in)   :: c_1(N_states),c_2(N_states)
 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate,p1
 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 call get_single_excitation(det_1,det_2,exc,phase,N_int)
 if (exc(0,1,1) == 1) then
  ! Mono alpha
  h1 = exc(1,1,1)
  p1 = exc(1,2,1)
  do istate = 1, N_states
   do i = 1, n_occ_ab(1)
    h2 = occ(i,1)
    big_array(h1,p1,h2,h2,istate) += 0.5d0 * c_1(istate) * c_2(istate) * phase
    big_array(h1,h2,h2,p1,istate) -= 0.5d0 * c_1(istate) * c_2(istate) * phase

    big_array(h2,h2,h1,p1,istate) += 0.5d0 * c_1(istate) * c_2(istate) * phase
    big_array(h2,p1,h1,h2,istate) -= 0.5d0 * c_1(istate) * c_2(istate) * phase
    enddo 
  enddo
 else 
  return
 endif
 end

 subroutine off_diagonal_single_to_two_rdm_bb_dm(det_1,det_2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the beta /beta  2RDM only for SINGLE EXCITATIONS IN CHEMIST NOTATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 double precision, intent(in)   :: c_1(N_states),c_2(N_states)
 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate,p1
 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 call get_single_excitation(det_1,det_2,exc,phase,N_int)
 if (exc(0,1,1) == 1) then
  return
 else
  ! Mono beta
  h1 =  exc(1,1,2)
  p1 =  exc(1,2,2)
  do istate = 1, N_states
   do i = 1, n_occ_ab(2)
    h2 = occ(i,2)
    big_array(h1,p1,h2,h2,istate) += 0.5d0 * c_1(istate) * c_2(istate) * phase
    big_array(h1,h2,h2,p1,istate) -= 0.5d0 * c_1(istate) * c_2(istate) * phase

    big_array(h2,h2,h1,p1,istate) += 0.5d0 * c_1(istate) * c_2(istate) * phase
    big_array(h2,p1,h1,h2,istate) -= 0.5d0 * c_1(istate) * c_2(istate) * phase
    enddo 
  enddo
 endif
 end


 subroutine off_diagonal_double_to_two_rdm_aa_dm(det_1,det_2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the alpha/alpha 2RDM only for DOUBLE EXCITATIONS IN CHEMIST NOTATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int),det_2(N_int)
 double precision, intent(in)   :: c_1(N_states),c_2(N_states)
 integer :: i,j,h1,h2,p1,p2,istate
 integer                        :: exc(0:2,2)
 double precision               :: phase
 call get_double_excitation_spin(det_1,det_2,exc,phase,N_int)
 h1 =exc(1,1)
 h2 =exc(2,1)
 p1 =exc(1,2)
 p2 =exc(2,2)
!print*,'h1,p1,h2,p2',h1,p1,h2,p2,c_1(istate) * phase * c_2(istate)
 do istate = 1, N_states
  big_array(h1,p1,h2,p2,istate) += 0.5d0 * c_1(istate) * phase * c_2(istate)
  big_array(h1,p2,h2,p1,istate) -= 0.5d0 * c_1(istate) * phase * c_2(istate)

  big_array(h2,p2,h1,p1,istate) += 0.5d0 * c_1(istate) * phase * c_2(istate)
  big_array(h2,p1,h1,p2,istate) -= 0.5d0 * c_1(istate) * phase * c_2(istate)
 enddo
 end

 subroutine off_diagonal_double_to_two_rdm_bb_dm(det_1,det_2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)
 use bitmasks
 BEGIN_DOC
! routine that update the OFF DIAGONAL PART of the beta /beta  2RDM only for DOUBLE EXCITATIONS 
 END_DOC
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int),det_2(N_int)
 double precision, intent(in)   :: c_1(N_states),c_2(N_states)
 integer :: i,j,h1,h2,p1,p2,istate
 integer                        :: exc(0:2,2)
 double precision               :: phase
 call get_double_excitation_spin(det_1,det_2,exc,phase,N_int)
 h1 =exc(1,1)
 h2 =exc(2,1)
 p1 =exc(1,2)
 p2 =exc(2,2)
!print*,'h1,p1,h2,p2',h1,p1,h2,p2,c_1(istate) * phase * c_2(istate)
 do istate = 1, N_states
  big_array(h1,p1,h2,p2,istate) += 0.5d0 * c_1(istate) * phase * c_2(istate)
  big_array(h1,p2,h2,p1,istate) -= 0.5d0 * c_1(istate) * phase * c_2(istate)

  big_array(h2,p2,h1,p1,istate) += 0.5d0 * c_1(istate) * phase * c_2(istate)
  big_array(h2,p1,h1,p2,istate) -= 0.5d0 * c_1(istate) * phase * c_2(istate)
 enddo
 end

