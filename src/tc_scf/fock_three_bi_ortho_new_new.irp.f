subroutine contrib_3e_sss(a,i,j,k,integral)
 integer, intent(in) :: a,i,j,k
 BEGIN_DOC
 ! returns the pure same spin contribution to F(a,i) from two orbitals j,k
 END_DOC
 double precision, intent(out) :: integral
 double precision :: direct_int, exch_13_int, exch_23_int, exch_12_int, c_3_int, c_minus_3_int
 call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )!!! < a k j | i k j >
 call  give_integrals_3_body_bi_ort(a, k, j, j, i, k, c_3_int)      ! < a k j | j i k >
 call  give_integrals_3_body_bi_ort(a, k, j, k, j, i, c_minus_3_int)! < a k j | k j i >
 integral = direct_int + c_3_int + c_minus_3_int 
 ! negative terms :: exchange contrib
 call  give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)!!! < a k j | j k i > : E_13 
 call  give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)!!! < a k j | i j k > : E_23
 call  give_integrals_3_body_bi_ort(a, k, j, k, i, j, exch_12_int)!!! < a k j | k i j > : E_12
 integral += - exch_13_int - exch_23_int  - exch_12_int 
 integral = -integral
end

subroutine contrib_3e_soo(a,i,j,k,integral)
 integer, intent(in) :: a,i,j,k
 BEGIN_DOC
 ! returns the same spin / opposite spin / opposite spin contribution to F(a,i) from two orbitals j,k
 END_DOC
 double precision, intent(out) :: integral
 double precision :: direct_int, exch_23_int
 call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int) ! < a k j | i k j >
 call  give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)! < a k j | i j k > : E_23
 integral = direct_int - exch_23_int 
 integral = -integral
end

subroutine contrib_3e_sos(a,i,j,k,integral)
 implicit none
 integer, intent(in) :: a,i,j,k
 BEGIN_DOC
 ! returns the same spin / opposite spin / same spin contribution to F(a,i) from two orbitals j,k
 END_DOC
 double precision, intent(out) :: integral
 double precision :: direct_int, exch_13_int
 call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )! < a k j | i k j >
 call  give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)! < a k j | j k i > : E_13 
 integral = direct_int - exch_13_int 
 integral = -integral
end

BEGIN_PROVIDER [double precision, fock_a_tot_3e_bi_orth_new, (mo_num, mo_num)]
 implicit none
 integer :: i,a,j,k
 double precision :: contrib_sss, contrib_sos, contrib_soo
 fock_a_tot_3e_bi_orth_new = 0.d0
 do i = 1, mo_num
  do a = 1, mo_num
   
!   do j = 1, elec_beta_num
!    do k = 1, elec_beta_num
!      call contrib_3e_sss(a,i,j,k,contrib_sss)
!      call contrib_3e_soo(a,i,j,k,contrib_soo)
!      call contrib_3e_sos(a,i,j,k,contrib_sos)
!      fock_a_tot_3e_bi_orth_new(a,i) += 0.5d0 * (contrib_sss + contrib_soo) + contrib_sos
!    enddo
!   enddo
   
   fock_a_tot_3e_bi_orth_new(a,i) += fock_cs_3e_bi_orth(a,i)

   do j = elec_beta_num + 1, elec_alpha_num 
    do k = 1, elec_beta_num
     call contrib_3e_sss(a,i,j,k,contrib_sss)
     call contrib_3e_sos(a,i,j,k,contrib_sos)
     fock_a_tot_3e_bi_orth_new(a,i) += 0.5d0 * contrib_sss + contrib_sos
    enddo
   enddo

   do j = 1, elec_alpha_num
    do k = elec_beta_num+1, elec_alpha_num
     call contrib_3e_sss(a,i,j,k,contrib_sss)
     fock_a_tot_3e_bi_orth_new(a,i) += 0.5d0 * contrib_sss 
    enddo
   enddo
!   do j = 1, elec_beta_num
!    do k = elec_beta_num+1, elec_alpha_num
!     call contrib_3e_sss(a,i,j,k,contrib_sss)
!     fock_a_tot_3e_bi_orth_new(a,i) += 0.5d0 * contrib_sss 
!    enddo
!   enddo
!  
!   do j = elec_beta_num+1, elec_alpha_num
!    do k = elec_beta_num+1, elec_alpha_num
!     call contrib_3e_sss(a,i,j,k,contrib_sss)
!     fock_a_tot_3e_bi_orth_new(a,i) += 0.5d0 * contrib_sss 
!    enddo
!   enddo

  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, fock_cs_3e_bi_orth, (mo_num, mo_num)]
 implicit none
 integer :: i,a,j,k
 double precision :: contrib_sss, contrib_sos, contrib_soo
 fock_cs_3e_bi_orth = 0.d0
 do i = 1, mo_num
  do a = 1, mo_num
   
   do j = 1, elec_beta_num
    do k = 1, elec_beta_num
      call contrib_3e_sss(a,i,j,k,contrib_sss)
      call contrib_3e_soo(a,i,j,k,contrib_soo)
      call contrib_3e_sos(a,i,j,k,contrib_sos)
      fock_cs_3e_bi_orth(a,i) += 0.5d0 * (contrib_sss + contrib_soo) + contrib_sos
    enddo
   enddo
  
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, fock_b_tot_3e_bi_orth_new, (mo_num, mo_num)]
 implicit none
 integer :: i,a,j,k
 double precision :: contrib_sss, contrib_sos, contrib_soo
 fock_b_tot_3e_bi_orth_new = 0.d0
 do i = 1, mo_num
  do a = 1, mo_num
   
!   do j = 1, elec_beta_num
!    do k = 1, elec_beta_num
!      call contrib_3e_sss(a,i,j,k,contrib_sss)
!      call contrib_3e_soo(a,i,j,k,contrib_soo)
!      call contrib_3e_sos(a,i,j,k,contrib_sos)
!      fock_b_tot_3e_bi_orth_new(a,i) += 0.5d0 * (contrib_sss + contrib_soo) + contrib_sos
!    enddo
!   enddo
   fock_b_tot_3e_bi_orth_new(a,i) += fock_cs_3e_bi_orth(a,i)

   do j = elec_beta_num + 1, elec_alpha_num 
    do k = 1, elec_alpha_num
      call contrib_3e_soo(a,i,j,k,contrib_soo)
      fock_b_tot_3e_bi_orth_new(a,i) += 0.5d0 * contrib_soo 
    enddo
   enddo
!   do j = elec_beta_num + 1, elec_alpha_num 
!    do k = 1, elec_beta_num
!      call contrib_3e_soo(a,i,j,k,contrib_soo)
!      fock_b_tot_3e_bi_orth_new(a,i) += 0.5d0 * contrib_soo 
!    enddo
!   enddo

   do j = 1, elec_beta_num
    do k = elec_beta_num+1, elec_alpha_num
     call contrib_3e_soo(a,i,j,k,contrib_soo)
     call contrib_3e_sos(a,i,j,k,contrib_sos)
     fock_b_tot_3e_bi_orth_new(a,i) += 0.5d0 * contrib_soo + contrib_sos
    enddo
   enddo

!   do j = elec_beta_num+1, elec_alpha_num
!    do k = elec_beta_num+1, elec_alpha_num
!     call contrib_3e_soo(a,i,j,k,contrib_soo)
!     fock_b_tot_3e_bi_orth_new(a,i) += 0.5d0 * contrib_soo
!    enddo
!   enddo

  enddo
 enddo

END_PROVIDER 
