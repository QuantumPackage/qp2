
! ---

BEGIN_PROVIDER [ double precision, three_e_4_idx_direct_bi_ort, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_4_idx_direct_bi_ort(m,j,k,i) = <mjk|-L|mji> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

 implicit none
 integer          :: i, j, k, m
 double precision :: integral, wall1, wall0

  three_e_4_idx_direct_bi_ort = 0.d0
  print *, ' Providing the three_e_4_idx_direct_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                   &
 !$OMP DEFAULT (NONE)             &
 !$OMP PRIVATE (i,j,k,m,integral) & 
 !$OMP SHARED (mo_num,three_e_4_idx_direct_bi_ort)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          call give_integrals_3_body_bi_ort(m, j, k, m, j, i, integral)
          three_e_4_idx_direct_bi_ort(m,j,k,i) = -1.d0 * integral 
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_4_idx_direct_bi_ort', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, three_e_4_idx_cycle_1_bi_ort, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_4_idx_cycle_1_bi_ort(m,j,k,i) = <mjk|-L|jim> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

  implicit none
  integer          :: i, j, k, m
  double precision :: integral, wall1, wall0

  three_e_4_idx_cycle_1_bi_ort = 0.d0
  print *, ' Providing the three_e_4_idx_cycle_1_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                   &
 !$OMP DEFAULT (NONE)             &
 !$OMP PRIVATE (i,j,k,m,integral) & 
 !$OMP SHARED (mo_num,three_e_4_idx_cycle_1_bi_ort)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          call give_integrals_3_body_bi_ort(m, j, k, j, i, m, integral)
          three_e_4_idx_cycle_1_bi_ort(m,j,k,i) = -1.d0 * integral 
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_4_idx_cycle_1_bi_ort', wall1 - wall0

END_PROVIDER 

! --

BEGIN_PROVIDER [ double precision, three_e_4_idx_cycle_2_bi_ort, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_4_idx_cycle_2_bi_ort(m,j,k,i) = <mjk|-L|imj> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

  implicit none
  integer          :: i, j, k, m
  double precision :: integral, wall1, wall0

  three_e_4_idx_cycle_2_bi_ort = 0.d0
  print *, ' Providing the three_e_4_idx_cycle_2_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                   &
 !$OMP DEFAULT (NONE)             &
 !$OMP PRIVATE (i,j,k,m,integral) & 
 !$OMP SHARED (mo_num,three_e_4_idx_cycle_2_bi_ort)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          call give_integrals_3_body_bi_ort(m, j, k, i, m, j, integral)
          three_e_4_idx_cycle_2_bi_ort(m,j,k,i) = -1.d0 * integral 
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_4_idx_cycle_2_bi_ort', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, three_e_4_idx_exch23_bi_ort, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_4_idx_exch23_bi_ort(m,j,k,i) = <mjk|-L|jmi> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

  implicit none
  integer          :: i, j, k, m
  double precision :: integral, wall1, wall0

  three_e_4_idx_exch23_bi_ort = 0.d0
  print *, ' Providing the three_e_4_idx_exch23_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                   &
 !$OMP DEFAULT (NONE)             &
 !$OMP PRIVATE (i,j,k,m,integral) & 
 !$OMP SHARED (mo_num,three_e_4_idx_exch23_bi_ort)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          call give_integrals_3_body_bi_ort(m, j, k, j, m, i, integral)
          three_e_4_idx_exch23_bi_ort(m,j,k,i) = -1.d0 * integral 
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_4_idx_exch23_bi_ort', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, three_e_4_idx_exch13_bi_ort, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_4_idx_exch13_bi_ort(m,j,k,i) = <mjk|-L|ijm> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  END_DOC

  implicit none
  integer          :: i, j, k, m
  double precision :: integral, wall1, wall0

  three_e_4_idx_exch13_bi_ort = 0.d0
  print *, ' Providing the three_e_4_idx_exch13_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                   &
 !$OMP DEFAULT (NONE)             &
 !$OMP PRIVATE (i,j,k,m,integral) & 
 !$OMP SHARED (mo_num,three_e_4_idx_exch13_bi_ort)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          call give_integrals_3_body_bi_ort(m, j, k, i, j, m, integral)
          three_e_4_idx_exch13_bi_ort(m,j,k,i) = -1.d0 * integral 
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_4_idx_exch13_bi_ort', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, three_e_4_idx_exch12_bi_ort, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  ! 
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_4_idx_exch12_bi_ort(m,j,k,i) = <mjk|-L|mij> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  ! 
  END_DOC

  implicit none
  integer          :: i, j, k, m
  double precision :: integral, wall1, wall0

  three_e_4_idx_exch12_bi_ort = 0.d0
  print *, ' Providing the three_e_4_idx_exch12_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                   &
 !$OMP DEFAULT (NONE)             &
 !$OMP PRIVATE (i,j,k,m,integral) & 
 !$OMP SHARED (mo_num,three_e_4_idx_exch12_bi_ort)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          call give_integrals_3_body_bi_ort(m, j, k, m, i, j, integral)
          three_e_4_idx_exch12_bi_ort(m,j,k,i) = -1.d0 * integral 
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_4_idx_exch12_bi_ort', wall1 - wall0

END_PROVIDER 

! ---

