
! ---

BEGIN_PROVIDER [double precision, fock_a_tot_3e_bi_orth, (mo_num, mo_num)]

  BEGIN_DOC
  !
  ! Alpha part of the Fock matrix from three-electron terms 
  !
  ! WARNING :: non hermitian if bi-ortho MOS used 
  !
  ! This calculation becomes the dominant part one the integrals are provided
  ! 
  END_DOC

  implicit none
  integer          :: i, a
  double precision :: t0, t1

  !print*, ' Providing fock_a_tot_3e_bi_orth ...'
  !call wall_time(t0)

  PROVIDE mo_l_coef mo_r_coef
  PROVIDE fock_cs_3e_bi_orth fock_a_tmp1_bi_ortho fock_a_tmp2_bi_ortho

  fock_a_tot_3e_bi_orth = 0.d0

  do i = 1, mo_num
    do a = 1, mo_num
      fock_a_tot_3e_bi_orth(a,i) += fock_cs_3e_bi_orth  (a,i)
      fock_a_tot_3e_bi_orth(a,i) += fock_a_tmp1_bi_ortho(a,i)
      fock_a_tot_3e_bi_orth(a,i) += fock_a_tmp2_bi_ortho(a,i)
    enddo
  enddo

  !call wall_time(t1)
  !print*, ' Wall time for fock_a_tot_3e_bi_orth =', t1 - t0 

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_b_tot_3e_bi_orth, (mo_num, mo_num)]

  BEGIN_DOC
  !
  ! Beta part of the Fock matrix from three-electron terms 
  !
  ! WARNING :: non hermitian if bi-ortho MOS used 
  !
  ! This calculation becomes the dominant part one the integrals are provided
  !
  END_DOC

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
  integer                       :: i, a, j, k
  double precision              :: contrib_sss, contrib_sos, contrib_soo, contrib
  double precision              :: direct_int, exch_13_int, exch_23_int, exch_12_int, c_3_int, c_minus_3_int
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:,:)

  !print*, ' Providing fock_cs_3e_bi_orth ...'
  !call wall_time(t0)

  PROVIDE mo_l_coef mo_r_coef

  ! to PROVIDE stuffs
  call give_integrals_3_body_bi_ort(1, 1, 1, 1, 1, 1, contrib) 

  fock_cs_3e_bi_orth = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                                                                              &
  !$OMP PRIVATE (i, a, j, k, direct_int, c_3_int, c_minus_3_int, exch_13_int, exch_23_int, exch_12_int, tmp) &
  !$OMP SHARED  (mo_num, elec_beta_num, fock_cs_3e_bi_orth)

  allocate(tmp(mo_num,mo_num))
  tmp = 0.d0

  !$OMP DO
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

          tmp(a,i) += 2.d0 * direct_int + 0.5d0 * (c_3_int + c_minus_3_int - exch_12_int) -1.5d0 * exch_13_int - exch_23_int
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do i = 1, mo_num
    do a = 1, mo_num
      fock_cs_3e_bi_orth(a,i) += tmp(a,i)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmp)
  !$OMP END PARALLEL
 
  fock_cs_3e_bi_orth = - fock_cs_3e_bi_orth

  !call wall_time(t1)
  !print*, ' Wall time for fock_cs_3e_bi_orth =', t1-t0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_a_tmp1_bi_ortho, (mo_num, mo_num)]

  implicit none
  integer                       :: i, a, j, k, ee
  double precision              :: contrib_sss, contrib_sos, contrib_soo, contrib
  double precision              :: direct_int, exch_13_int, exch_23_int, exch_12_int, c_3_int, c_minus_3_int
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:,:)

  !print*, ' Providing fock_a_tmp1_bi_ortho ...'
  !call wall_time(t0)

  PROVIDE mo_l_coef mo_r_coef

  ! to PROVIDE stuffs
  call give_integrals_3_body_bi_ort(1, 1, 1, 1, 1, 1, contrib) 

  ee = elec_beta_num + 1
  fock_a_tmp1_bi_ortho = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                                                                              &
  !$OMP PRIVATE (i, a, j, k, direct_int, c_3_int, c_minus_3_int, exch_13_int, exch_23_int, exch_12_int, tmp) &
  !$OMP SHARED  (mo_num, elec_alpha_num, elec_beta_num, ee, fock_a_tmp1_bi_ortho)

  allocate(tmp(mo_num,mo_num))
  tmp = 0.d0

  !$OMP DO
  do i = 1, mo_num
    do a = 1, mo_num

      do j = ee, elec_alpha_num 
        do k = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )!!! < a k j | i k j >
          call give_integrals_3_body_bi_ort(a, k, j, j, i, k, c_3_int)      ! < a k j | j i k >
          call give_integrals_3_body_bi_ort(a, k, j, k, j, i, c_minus_3_int)! < a k j | k j i >
          call give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)!!! < a k j | j k i > : E_13 
          call give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)!!! < a k j | i j k > : E_23
          call give_integrals_3_body_bi_ort(a, k, j, k, i, j, exch_12_int)!!! < a k j | k i j > : E_12
          
          tmp(a,i) += 1.5d0 * (direct_int - exch_13_int) + 0.5d0 * (c_3_int + c_minus_3_int - exch_23_int - exch_12_int)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do i = 1, mo_num
    do a = 1, mo_num
      fock_a_tmp1_bi_ortho(a,i) += tmp(a,i)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmp)
  !$OMP END PARALLEL

  fock_a_tmp1_bi_ortho = - fock_a_tmp1_bi_ortho

  !call wall_time(t1)
  !print*, ' Wall time for fock_a_tmp1_bi_ortho =', t1-t0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_a_tmp2_bi_ortho, (mo_num, mo_num)]

  implicit none
  integer                       :: i, a, j, k, ee
  double precision              :: contrib_sss
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:,:)

  !print*, ' Providing fock_a_tmp2_bi_ortho ...'
  !call wall_time(t0)

  PROVIDE mo_l_coef mo_r_coef
 
  ! to PROVIDE stuffs
  call contrib_3e_sss(1, 1, 1, 1, contrib_sss)

  ee = elec_beta_num + 1
  fock_a_tmp2_bi_ortho = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                &
  !$OMP PRIVATE (i, a, j, k, contrib_sss, tmp) &
  !$OMP SHARED  (mo_num, elec_alpha_num, ee, fock_a_tmp2_bi_ortho)

  allocate(tmp(mo_num,mo_num))
  tmp = 0.d0

  !$OMP DO
  do i = 1, mo_num
    do a = 1, mo_num
      do j = 1, elec_alpha_num
        do k = ee, elec_alpha_num
          call contrib_3e_sss(a, i, j, k, contrib_sss)

          tmp(a,i) += 0.5d0 * contrib_sss
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do i = 1, mo_num
    do a = 1, mo_num
      fock_a_tmp2_bi_ortho(a,i) += tmp(a,i)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmp)
  !$OMP END PARALLEL

  !call wall_time(t1)
  !print*, ' Wall time for fock_a_tmp2_bi_ortho =', t1-t0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_b_tmp1_bi_ortho, (mo_num, mo_num)]

  implicit none
  integer                       :: i, a, j, k, ee
  double precision              :: direct_int, exch_13_int, exch_23_int, exch_12_int
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:,:)

  !print*, ' Providing fock_b_tmp1_bi_ortho ...'
  !call wall_time(t0)

  PROVIDE mo_l_coef mo_r_coef

  ! to PROVIDE stuffs
  call give_integrals_3_body_bi_ort(1, 1, 1, 1, 1, 1, direct_int) 

  ee = elec_beta_num + 1
  fock_b_tmp1_bi_ortho = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                                         &
  !$OMP PRIVATE (i, a, j, k, direct_int, exch_13_int, exch_23_int, tmp) &
  !$OMP SHARED  (mo_num, elec_beta_num, elec_alpha_num, ee, fock_b_tmp1_bi_ortho)

  allocate(tmp(mo_num,mo_num))
  tmp = 0.d0

  !$OMP DO
  do i = 1, mo_num
    do a = 1, mo_num
      do j = 1, elec_beta_num
        do k = ee, elec_alpha_num
          call  give_integrals_3_body_bi_ort(a, k, j, i, k, j, direct_int )!!! < a k j | i k j >
          call  give_integrals_3_body_bi_ort(a, k, j, j, k, i, exch_13_int)!!! < a k j | j k i > : E_13 
          call  give_integrals_3_body_bi_ort(a, k, j, i, j, k, exch_23_int)!!! < a k j | i j k > : E_23

          tmp(a,i) += 1.5d0 * direct_int - 0.5d0 * exch_23_int - exch_13_int
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do i = 1, mo_num
    do a = 1, mo_num
      fock_b_tmp1_bi_ortho(a,i) += tmp(a,i)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmp)
  !$OMP END PARALLEL

  fock_b_tmp1_bi_ortho = - fock_b_tmp1_bi_ortho

  !call wall_time(t1)
  !print*, ' Wall time for fock_b_tmp1_bi_ortho =', t1-t0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_b_tmp2_bi_ortho, (mo_num, mo_num)]

  implicit none
  integer                       :: i, a, j, k, ee
  double precision              :: contrib_soo
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:,:)

  !print*, ' Providing fock_b_tmp2_bi_ortho ...'
  !call wall_time(t0)

  PROVIDE mo_l_coef mo_r_coef

  ! to PROVIDE stuffs
  call contrib_3e_soo(1, 1, 1, 1, contrib_soo)

  ee = elec_beta_num + 1
  fock_b_tmp2_bi_ortho = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                &
  !$OMP PRIVATE (i, a, j, k, contrib_soo, tmp) &
  !$OMP SHARED  (mo_num, elec_alpha_num, ee, fock_b_tmp2_bi_ortho)

  allocate(tmp(mo_num,mo_num))
  tmp = 0.d0

  !$OMP DO
  do i = 1, mo_num
    do a = 1, mo_num
      do j = ee, elec_alpha_num 
        do k = 1, elec_alpha_num
          call contrib_3e_soo(a, i, j, k, contrib_soo)

          tmp(a,i) += 0.5d0 * contrib_soo
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do i = 1, mo_num
    do a = 1, mo_num
      fock_b_tmp2_bi_ortho(a,i) += tmp(a,i)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmp)
  !$OMP END PARALLEL

  !call wall_time(t1)
  !print*, ' Wall time for fock_b_tmp2_bi_ortho =', t1-t0

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

