
! This file is an example of the kind of manipulations that you can do with providers 
!

!!!!!!!!!!!!!!!!!!!!!!!!!! Main providers useful for the program !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!               type             name 
BEGIN_PROVIDER [ double precision, trace_mo_one_e_ints]
 implicit none
 BEGIN_DOC
! trace_mo_one_e_ints = Trace of the one-electron integrals on the MO basis
!
!                     = sum_i mo_one_e_integrals(i,i)
 END_DOC
 integer :: i
 trace_mo_one_e_ints = 0.d0
 do i = 1, mo_num
  trace_mo_one_e_ints += mo_one_e_integrals(i,i)
 enddo
END_PROVIDER  

BEGIN_PROVIDER [ double precision, trace_ao_one_e_ints]
 implicit none
 BEGIN_DOC
! trace_ao_one_e_ints = Trace of the one-electron integrals on the AO basis taking into account the non orthogonality
!
! Be aware that the trace of an operator in a non orthonormal basis is Tr(A S^{-1}) = \sum_{m,n}(A_mn S^{-1}_mn) 
! 
! WARNING: it is equal to the trace on the MO basis if and only if the AO basis and MO basis 
!          have the same number of functions
 END_DOC
 integer :: i,j
 double precision :: accu
 double precision, allocatable :: inv_overlap_times_integrals(:,:) ! = h S^{-1}
 allocate(inv_overlap_times_integrals(ao_num,ao_num))
 ! routine that computes the product of two matrices, you can check it with 
 ! irpman get_AB_prod
 call get_AB_prod(ao_one_e_integrals,ao_num,ao_num,s_inv,ao_num,inv_overlap_times_integrals)
 ! Tr(inv_overlap_times_integrals) = Tr(h S^{-1})
 trace_ao_one_e_ints = 0.d0
 do i = 1, ao_num
  trace_ao_one_e_ints += inv_overlap_times_integrals(i,i)
 enddo
 !
 ! testing the formula Tr(A S^{-1}) = \sum_{m,n}(A_mn S^{-1}_mn)
 double precision :: test
 test = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   test += ao_one_e_integrals(j,i) * s_inv(i,j)
  enddo
 enddo
 if(dabs(accu - trace_ao_one_e_ints).gt.1.d-12)then
  print*,'Warning ! '
  print*,'Something is wrong because Tr(AB) \ne sum_{mn}A_mn B_nm'
 endif
END_PROVIDER

BEGIN_PROVIDER [ double precision, trace_ao_one_e_ints_from_mo]
 implicit none
 BEGIN_DOC
! trace_ao_one_e_ints_from_mo = Trace of the one-electron integrals on the AO basis after projection on the MO basis 
!
!                             = Tr([SC h {SC}^+] S^{-1})
!
!                             = Be aware that the trace of an operator in a non orthonormal basis is = Tr(A S^{-1}) where S is the metric
! Must be equal to the trace_mo_one_e_ints 
 END_DOC
 integer :: i
 double precision, allocatable :: inv_overlap_times_integrals(:,:)
 allocate(inv_overlap_times_integrals(ao_num,ao_num))
 ! Using the provider ao_one_e_integrals_from_mo = [SC h {SC}^+] 
 call get_AB_prod(ao_one_e_integrals_from_mo,ao_num,ao_num,s_inv,ao_num,inv_overlap_times_integrals)
 ! inv_overlap_times_integrals = [SC h {SC}^+] S^{-1}
 trace_ao_one_e_ints_from_mo = 0.d0
 ! Computing the trace
 do i = 1, ao_num
  trace_ao_one_e_ints_from_mo += inv_overlap_times_integrals(i,i)
 enddo
END_PROVIDER  

!!!!!!!!!!!!!!!!!!!!!!!!!!! Additional providers to check some stuffs !!!!!!!!!!!!!!!!!!!!!!!!!

BEGIN_PROVIDER [ double precision, ao_one_e_int_no_ov_from_mo, (ao_num, ao_num) ]
 BEGIN_DOC
 ! ao_one_e_int_no_ov_from_mo = C mo_one_e_integrals C^T 
 !
 ! WARNING : NON EQUAL TO ao_one_e_integrals due to the non orthogonality 
 END_DOC
 call mo_to_ao_no_overlap(mo_one_e_integrals,mo_num,ao_one_e_int_no_ov_from_mo,ao_num) 
END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_one_e_int_no_ov_from_mo_ov_ov, (ao_num, ao_num)]
 BEGIN_DOC
 ! ao_one_e_int_no_ov_from_mo_ov_ov = S ao_one_e_int_no_ov_from_mo S = SC mo_one_e_integrals (SC)^T
 !
 ! EQUAL TO ao_one_e_integrals ONLY IF ao_num = mo_num
 END_DOC
 double precision, allocatable :: tmp(:,:)
 allocate(tmp(ao_num, ao_num))
 call get_AB_prod(ao_overlap,ao_num,ao_num,ao_one_e_int_no_ov_from_mo,ao_num,tmp)
 call get_AB_prod(tmp,ao_num,ao_num,ao_overlap,ao_num,ao_one_e_int_no_ov_from_mo_ov_ov)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, c_t_s_c, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! C^T S C = should be the identity 
 END_DOC 
 call get_AB_prod(mo_coef_transp,mo_num,ao_num,S_mo_coef,mo_num,c_t_s_c)
END_PROVIDER 

