BEGIN_PROVIDER [ double precision, fock_a_abb_3e_bi_orth_old, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! fock_a_abb_3e_bi_orth_old(a,i) = bi-ortho 3-e Fock matrix for alpha electrons from alpha,beta,beta contribution
 END_DOC
 fock_a_abb_3e_bi_orth_old = 0.d0
 integer :: i,a,j,k
 double precision :: direct_int, exch_23_int
 do i = 1, mo_num
  do a = 1, mo_num
   
   do j = 1, elec_beta_num
    do k = j+1, elec_beta_num
      ! see contrib_3e_soo
      call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int) ! < a k j | i k j >
      call  give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)! < a k j | i j k > : E_23
      fock_a_abb_3e_bi_orth_old(a,i) += direct_int - exch_23_int 
    enddo
   enddo

  enddo
 enddo
 fock_a_abb_3e_bi_orth_old = - fock_a_abb_3e_bi_orth_old
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_a_aba_3e_bi_orth_old, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! fock_a_aba_3e_bi_orth_old(a,i) = bi-ortho 3-e Fock matrix for alpha electrons from alpha,alpha,beta contribution
 END_DOC
 fock_a_aba_3e_bi_orth_old = 0.d0
 integer :: i,a,j,k
 double precision :: direct_int, exch_13_int
 do i = 1, mo_num
  do a = 1, mo_num
   
   do j = 1, elec_alpha_num ! a
    do k = 1, elec_beta_num ! b
                                                                       !   a b a   a b a
      call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )! < a k j | i k j >
      call  give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)! < a k j | j k i > : E_13 
      fock_a_aba_3e_bi_orth_old(a,i) += direct_int - exch_13_int 
    enddo
   enddo

  enddo
 enddo
 fock_a_aba_3e_bi_orth_old = - fock_a_aba_3e_bi_orth_old
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_a_aaa_3e_bi_orth_old, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! fock_a_aaa_3e_bi_orth_old(a,i) = bi-ortho 3-e Fock matrix for alpha electrons from alpha,alpha,alpha contribution
 END_DOC
 fock_a_aaa_3e_bi_orth_old = 0.d0
 integer :: i,a,j,k
 double precision :: direct_int, exch_13_int, exch_23_int, exch_12_int, c_3_int, c_minus_3_int
 do i = 1, mo_num
  do a = 1, mo_num
   
   do j = 1, elec_alpha_num
    do k = j+1, elec_alpha_num
      ! positive terms :: cycle contrib 
      call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )!!! < a k j | i k j >
      call  give_integrals_3_body_bi_ort(a, k, j, j, i, k, c_3_int)      ! < a k j | j i k >
      call  give_integrals_3_body_bi_ort(a, k, j, k, j, i, c_minus_3_int)! < a k j | k j i >
      fock_a_aaa_3e_bi_orth_old(a,i) += direct_int + c_3_int + c_minus_3_int 
      ! negative terms :: exchange contrib
      call  give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)!!! < a k j | j k i > : E_13 
      call  give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)!!! < a k j | i j k > : E_23
      call  give_integrals_3_body_bi_ort(a, k, j, k, i, j, exch_12_int)!!! < a k j | k i j > : E_12
      fock_a_aaa_3e_bi_orth_old(a,i) += - exch_13_int - exch_23_int  - exch_12_int 
    enddo
   enddo

  enddo
 enddo
 fock_a_aaa_3e_bi_orth_old = - fock_a_aaa_3e_bi_orth_old
END_PROVIDER 

BEGIN_PROVIDER [double precision, fock_a_tot_3e_bi_orth_old, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! fock_a_tot_3e_bi_orth_old = bi-ortho 3-e Fock matrix for alpha electrons from all possible spin contributions 
 END_DOC
 fock_a_tot_3e_bi_orth_old = fock_a_abb_3e_bi_orth_old + fock_a_aba_3e_bi_orth_old + fock_a_aaa_3e_bi_orth_old 

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_b_baa_3e_bi_orth_old, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! fock_b_baa_3e_bi_orth_old(a,i) = bi-ortho 3-e Fock matrix for beta electrons from beta,alpha,alpha contribution
 END_DOC
 fock_b_baa_3e_bi_orth_old = 0.d0
 integer :: i,a,j,k
 double precision :: direct_int, exch_23_int
 do i = 1, mo_num
  do a = 1, mo_num
   
   do j = 1, elec_alpha_num
    do k = j+1, elec_alpha_num
      call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int) ! < a k j | i k j >
      call  give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)! < a k j | i j k > : E_23
      fock_b_baa_3e_bi_orth_old(a,i) += direct_int - exch_23_int 
    enddo
   enddo

  enddo
 enddo
 fock_b_baa_3e_bi_orth_old = - fock_b_baa_3e_bi_orth_old
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_b_bab_3e_bi_orth_old, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! fock_b_bab_3e_bi_orth_old(a,i) = bi-ortho 3-e Fock matrix for beta electrons from beta,alpha,beta contribution
 END_DOC
 fock_b_bab_3e_bi_orth_old = 0.d0
 integer :: i,a,j,k
 double precision :: direct_int, exch_13_int
 do i = 1, mo_num
  do a = 1, mo_num
   
   do j = 1, elec_beta_num
    do k = 1, elec_alpha_num
      !                                                                    b a b   b a b
      call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int) ! < a k j | i k j >
      call  give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)! < a k j | j k i > : E_13
      fock_b_bab_3e_bi_orth_old(a,i) += direct_int - exch_13_int 
    enddo
   enddo

  enddo
 enddo
 fock_b_bab_3e_bi_orth_old = - fock_b_bab_3e_bi_orth_old
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_b_bbb_3e_bi_orth_old, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! fock_b_bbb_3e_bi_orth_old(a,i) = bi-ortho 3-e Fock matrix for alpha electrons from alpha,alpha,alpha contribution
 END_DOC
 fock_b_bbb_3e_bi_orth_old = 0.d0
 integer :: i,a,j,k
 double precision :: direct_int, exch_13_int, exch_23_int, exch_12_int, c_3_int, c_minus_3_int
 do i = 1, mo_num
  do a = 1, mo_num
   
   do j = 1, elec_beta_num
    do k = j+1, elec_beta_num
      ! positive terms :: cycle contrib 
      call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )!!! < a k j | i k j >
      call  give_integrals_3_body_bi_ort(a, k, j, j, i, k, c_3_int)      ! < a k j | j i k >
      call  give_integrals_3_body_bi_ort(a, k, j, k, j, i, c_minus_3_int)! < a k j | k j i >
      fock_b_bbb_3e_bi_orth_old(a,i) += direct_int + c_3_int + c_minus_3_int 
      ! negative terms :: exchange contrib
      call  give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)!!! < a k j | j k i > : E_13 
      call  give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)!!! < a k j | i j k > : E_23
      call  give_integrals_3_body_bi_ort(a, k, j, k, i, j, exch_12_int)!!! < a k j | k i j > : E_12
      fock_b_bbb_3e_bi_orth_old(a,i) += - exch_13_int - exch_23_int  - exch_12_int 
    enddo
   enddo

  enddo
 enddo
 fock_b_bbb_3e_bi_orth_old = - fock_b_bbb_3e_bi_orth_old
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_b_tot_3e_bi_orth_old, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! fock_b_tot_3e_bi_orth_old = bi-ortho 3-e Fock matrix for alpha electrons from all possible spin contributions 
 END_DOC
 fock_b_tot_3e_bi_orth_old = fock_b_bbb_3e_bi_orth_old + fock_b_bab_3e_bi_orth_old + fock_b_baa_3e_bi_orth_old

END_PROVIDER 
