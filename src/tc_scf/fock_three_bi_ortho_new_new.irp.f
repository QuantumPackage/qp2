
! ---

BEGIN_PROVIDER [double precision, fock_a_tot_3e_bi_orth, (mo_num, mo_num)]

  implicit none
  integer :: i, a

  PROVIDE mo_l_coef mo_r_coef

  fock_a_tot_3e_bi_orth = 0.d0

  do i = 1, mo_num
    do a = 1, mo_num
      fock_a_tot_3e_bi_orth(a,i) += fock_cs_3e_bi_orth  (a,i)
      fock_a_tot_3e_bi_orth(a,i) += fock_a_tmp1_bi_ortho(a,i)
      fock_a_tot_3e_bi_orth(a,i) += fock_a_tmp2_bi_ortho(a,i)
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_b_tot_3e_bi_orth, (mo_num, mo_num)]

  implicit none
  integer :: i, a

  PROVIDE mo_l_coef mo_r_coef

  fock_b_tot_3e_bi_orth = 0.d0

  do i = 1, mo_num
    do a = 1, mo_num
      fock_b_tot_3e_bi_orth(a,i) += fock_cs_3e_bi_orth  (a,i)
      fock_b_tot_3e_bi_orth(a,i) += fock_b_tmp2_bi_ortho(a,i)
      fock_b_tot_3e_bi_orth(a,i) += fock_b_tmp1_bi_ortho(a,i)
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_cs_3e_bi_orth, (mo_num, mo_num)]

  implicit none
  integer          :: i, a, j, k
  double precision :: contrib_sss, contrib_sos, contrib_soo, contrib
  double precision :: direct_int, exch_13_int, exch_23_int, exch_12_int, c_3_int, c_minus_3_int
  double precision :: new

  PROVIDE mo_l_coef mo_r_coef

  fock_cs_3e_bi_orth = 0.d0

  do i = 1, mo_num
    do a = 1, mo_num
    
      do j = 1, elec_beta_num
        do k = 1, elec_beta_num

          !!call contrib_3e_sss(a,i,j,k,contrib_sss)
          !!call contrib_3e_soo(a,i,j,k,contrib_soo)
          !!call contrib_3e_sos(a,i,j,k,contrib_sos)
          !!contrib = 0.5d0 * (contrib_sss + contrib_soo) + contrib_sos
 
          call give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )!!! < a k j | i k j >
          call give_integrals_3_body_bi_ort(a, k, j, j, i, k, c_3_int)      ! < a k j | j i k >
          call give_integrals_3_body_bi_ort(a, k, j, k, j, i, c_minus_3_int)! < a k j | k j i >

          ! negative terms :: exchange contrib
          call give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)!!! < a k j | j k i > : E_13 
          call give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)!!! < a k j | i j k > : E_23
          call give_integrals_3_body_bi_ort(a, k, j, k, i, j, exch_12_int)!!! < a k j | k i j > : E_12

          new = 2.d0 * direct_int + 0.5d0 * (c_3_int + c_minus_3_int - exch_12_int) -1.5d0 * exch_13_int - exch_23_int

          fock_cs_3e_bi_orth(a,i) += new
        enddo
      enddo
    enddo
  enddo
 
  fock_cs_3e_bi_orth = - fock_cs_3e_bi_orth

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_a_tmp1_bi_ortho, (mo_num, mo_num)]

  implicit none
  integer          :: i, a, j, k
  double precision :: contrib_sss, contrib_sos, contrib_soo, contrib
  double precision :: direct_int, exch_13_int, exch_23_int, exch_12_int, c_3_int, c_minus_3_int
  double precision :: new

  PROVIDE mo_l_coef mo_r_coef

  fock_a_tmp1_bi_ortho = 0.d0

  do i = 1, mo_num
    do a = 1, mo_num
    
      do j = elec_beta_num + 1, elec_alpha_num 
        do k = 1, elec_beta_num
          call give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )!!! < a k j | i k j >
          call give_integrals_3_body_bi_ort(a, k, j, j, i, k, c_3_int)      ! < a k j | j i k >
          call give_integrals_3_body_bi_ort(a, k, j, k, j, i, c_minus_3_int)! < a k j | k j i >
          call give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)!!! < a k j | j k i > : E_13 
          call give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)!!! < a k j | i j k > : E_23
          call give_integrals_3_body_bi_ort(a, k, j, k, i, j, exch_12_int)!!! < a k j | k i j > : E_12
          
          fock_a_tmp1_bi_ortho(a,i) += 1.5d0 * (direct_int - exch_13_int) + 0.5d0 * (c_3_int + c_minus_3_int - exch_23_int - exch_12_int)
        enddo
      enddo
    enddo
  enddo

  fock_a_tmp1_bi_ortho = - fock_a_tmp1_bi_ortho

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_a_tmp2_bi_ortho, (mo_num, mo_num)]

  implicit none
  integer          :: i, a, j, k
  double precision :: contrib_sss

  PROVIDE mo_l_coef mo_r_coef
 
  fock_a_tmp2_bi_ortho = 0.d0

  do i = 1, mo_num
    do a = 1, mo_num
      do j = 1, elec_alpha_num
        do k = elec_beta_num+1, elec_alpha_num
          call contrib_3e_sss(a, i, j, k, contrib_sss)

          fock_a_tmp2_bi_ortho(a,i) += 0.5d0 * contrib_sss
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_b_tmp1_bi_ortho, (mo_num, mo_num)]

  implicit none
  integer          :: i, a, j, k
  double precision :: direct_int, exch_13_int, exch_23_int, exch_12_int
  double precision :: new

  PROVIDE mo_l_coef mo_r_coef

  fock_b_tmp1_bi_ortho = 0.d0

  do i = 1, mo_num
    do a = 1, mo_num
      do j = 1, elec_beta_num
        do k = elec_beta_num+1, elec_alpha_num
          call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )!!! < a k j | i k j >
          call  give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)!!! < a k j | j k i > : E_13 
          call  give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)!!! < a k j | i j k > : E_23

          fock_b_tmp1_bi_ortho(a,i) += 1.5d0 * direct_int - 0.5d0 * exch_23_int - exch_13_int
        enddo
      enddo
    enddo
  enddo

  fock_b_tmp1_bi_ortho = - fock_b_tmp1_bi_ortho

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_b_tmp2_bi_ortho, (mo_num, mo_num)]

  implicit none
  integer          :: i, a, j, k
  double precision :: contrib_soo

  PROVIDE mo_l_coef mo_r_coef

  fock_b_tmp2_bi_ortho = 0.d0

  do i = 1, mo_num
    do a = 1, mo_num
      do j = elec_beta_num + 1, elec_alpha_num 
        do k = 1, elec_alpha_num
          call contrib_3e_soo(a, i, j, k, contrib_soo)

          fock_b_tmp2_bi_ortho(a,i) += 0.5d0 * contrib_soo
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

subroutine contrib_3e_sss(a, i, j, k, integral)

  BEGIN_DOC
  ! returns the pure same spin contribution to F(a,i) from two orbitals j,k
  END_DOC

  implicit none
  integer,          intent(in)  :: a, i, j, k
  double precision, intent(out) :: integral
  double precision              :: direct_int, exch_13_int, exch_23_int, exch_12_int, c_3_int, c_minus_3_int

  PROVIDE mo_l_coef mo_r_coef

  call give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )!!! < a k j | i k j >
  call give_integrals_3_body_bi_ort(a, k, j, j, i, k, c_3_int)      ! < a k j | j i k >
  call give_integrals_3_body_bi_ort(a, k, j, k, j, i, c_minus_3_int)! < a k j | k j i >
  integral = direct_int + c_3_int + c_minus_3_int 

  ! negative terms :: exchange contrib
  call give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)!!! < a k j | j k i > : E_13 
  call give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)!!! < a k j | i j k > : E_23
  call give_integrals_3_body_bi_ort(a, k, j, k, i, j, exch_12_int)!!! < a k j | k i j > : E_12
  integral += - exch_13_int - exch_23_int  - exch_12_int 

  integral = -integral

end

! ---

subroutine contrib_3e_soo(a,i,j,k,integral)

  BEGIN_DOC
  ! returns the same spin / opposite spin / opposite spin contribution to F(a,i) from two orbitals j,k
  END_DOC

  implicit none
  integer,          intent(in)  :: a, i, j, k
  double precision, intent(out) :: integral
  double precision              :: direct_int, exch_23_int

  PROVIDE mo_l_coef mo_r_coef

  call give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int) ! < a k j | i k j >
  call give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)! < a k j | i j k > : E_23
  integral = direct_int - exch_23_int 

  integral = -integral

end

! ---

subroutine contrib_3e_sos(a, i, j, k, integral)

  BEGIN_DOC
  ! returns the same spin / opposite spin / same spin contribution to F(a,i) from two orbitals j,k
  END_DOC

  PROVIDE mo_l_coef mo_r_coef

  implicit none
  integer,          intent(in)  :: a, i, j, k
  double precision, intent(out) :: integral
  double precision              :: direct_int, exch_13_int

  call give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )! < a k j | i k j >
  call give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)! < a k j | j k i > : E_13 
  integral = direct_int - exch_13_int 

  integral = -integral

end

! ---

