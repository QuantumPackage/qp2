
BEGIN_PROVIDER [ double precision, three_e_diag_parrallel_spin_prov, (mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator ON A BI ORTHONORMAL BASIS 
  !
  ! three_e_diag_parrallel_spin_prov(m,j,i) = All combinations of the form <mji|-L|mji> for same spin matrix elements  
  ! 
  ! notice the -1 sign: in this way three_e_diag_parrallel_spin_prov can be directly used to compute Slater rules with a + sign
  !
  END_DOC

  implicit none
  integer          :: i, j, m
  double precision :: integral, wall1, wall0, three_e_diag_parrallel_spin

  three_e_diag_parrallel_spin_prov = 0.d0
  print *, ' Providing the three_e_diag_parrallel_spin_prov ...'

 integral = three_e_diag_parrallel_spin(1,1,1) ! to provide all stuffs
  call wall_time(wall0)
 !$OMP PARALLEL                 &
 !$OMP DEFAULT (NONE)           &
 !$OMP PRIVATE (i,j,m,integral) & 
 !$OMP SHARED (mo_num,three_e_diag_parrallel_spin_prov)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do j = 1, mo_num
      do m = j, mo_num
        three_e_diag_parrallel_spin_prov(m,j,i) =  three_e_diag_parrallel_spin(m,j,i)
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do i = 1, mo_num
    do j = 1, mo_num
      do m = 1, j
        three_e_diag_parrallel_spin_prov(m,j,i) = three_e_diag_parrallel_spin_prov(j,m,i)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print *, ' wall time for three_e_diag_parrallel_spin_prov', wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_e_single_parrallel_spin_prov, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_single_parrallel_spin_prov(m,j,k,i) = All combination of <mjk|-L|mji> for same spin matrix elements 
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

 implicit none
 integer          :: i, j, k, m
 double precision :: integral, wall1, wall0, three_e_single_parrallel_spin

  three_e_single_parrallel_spin_prov = 0.d0
  print *, ' Providing the three_e_single_parrallel_spin_prov ...'

  integral = three_e_single_parrallel_spin(1,1,1,1)
  call wall_time(wall0)
 !$OMP PARALLEL                   &
 !$OMP DEFAULT (NONE)             &
 !$OMP PRIVATE (i,j,k,m,integral) & 
 !$OMP SHARED (mo_num,three_e_single_parrallel_spin_prov)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          three_e_single_parrallel_spin_prov(m,j,k,i) = three_e_single_parrallel_spin(m,j,k,i)
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_single_parrallel_spin_prov', wall1 - wall0

END_PROVIDER 


! ---

BEGIN_PROVIDER [ double precision, three_e_double_parrallel_spin_prov, (mo_num, mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_double_parrallel_spin_prov(m,l,j,k,i) = <mlk|-L|mji> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  END_DOC

  implicit none
  integer          :: i, j, k, m, l
  double precision :: integral, wall1, wall0, three_e_double_parrallel_spin

  three_e_double_parrallel_spin_prov = 0.d0
  print *, ' Providing the three_e_double_parrallel_spin_prov ...'
  call wall_time(wall0)

 integral = three_e_double_parrallel_spin(1,1,1,1,1)
 !$OMP PARALLEL                     &
 !$OMP DEFAULT (NONE)               &
 !$OMP PRIVATE (i,j,k,m,l,integral) & 
 !$OMP SHARED (mo_num,three_e_double_parrallel_spin_prov)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
          do m = 1, mo_num
            three_e_double_parrallel_spin_prov(m,l,j,k,i) = three_e_double_parrallel_spin(m,l,j,k,i)
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_double_parrallel_spin_prov', wall1 - wall0

END_PROVIDER 

