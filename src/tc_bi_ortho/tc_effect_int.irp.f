
  
BEGIN_PROVIDER [double precision, mo_tc_effec2e_int, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC 
  ! 
  ! mo_tc_effec2e_int(p,q,s,t) = < p q| V(12) | s t > + \sum_i < p q i | L(123)| s t i >
  !
  ! the potential V(12) contains ALL TWO-E CONTRIBUTION OF THE TC-HAMILTONIAN
  !
  END_DOC

  implicit none
  integer          :: i, j, k, l, ii
  double precision :: integral

  PROVIDE mo_bi_ortho_tc_two_e_chemist

  do j = 1, mo_num
    do i = 1, mo_num
      do l = 1, mo_num 
        do k = 1, mo_num
          mo_tc_effec2e_int(k,l,i,j) = mo_bi_ortho_tc_two_e_chemist(k,i,l,j)
 
          do ii = 1, elec_alpha_num
            call give_integrals_3_body_bi_ort(k, l, ii, i, j, ii, integral)
            mo_tc_effec2e_int(k,l,i,j) -= 2.d0 * integral 
          enddo
        enddo
      enddo
    enddo
  enddo
        
  FREE mo_bi_ortho_tc_two_e_chemist
        
END_PROVIDER 

! ---

