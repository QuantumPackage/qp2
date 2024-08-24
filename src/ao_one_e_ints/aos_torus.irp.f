
! ---

 BEGIN_PROVIDER [ double precision, ao_overlap_torus  , (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_overlap_torus_x, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_overlap_torus_y, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_overlap_torus_z, (ao_num, ao_num) ]

  BEGIN_DOC
  ! Overlap between atomic basis functions:
  !
  ! TODO
  !
  ! :math:`\int \chi_i(r) \chi_j(r) dr`
  END_DOC

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta, c
  double precision :: A_center(3), B_center(3)

  PROVIDE torus_length
  print*, ' torus Lx = ', torus_length(1)
  print*, ' torus Ly = ', torus_length(2)
  print*, ' torus Lz = ', torus_length(3)

  ao_overlap_torus   = 0.d0
  ao_overlap_torus_x = 0.d0
  ao_overlap_torus_y = 0.d0
  ao_overlap_torus_z = 0.d0

!  dim1=100
!  !$OMP PARALLEL DO SCHEDULE(GUIDED) &
!  !$OMP DEFAULT(NONE) &
!  !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
!  !$OMP  overlap_x,overlap_y, overlap_z, overlap, &
!  !$OMP  alpha, beta,i,j,c) &
!  !$OMP SHARED(nucl_coord,ao_power,ao_prim_num, &
!  !$OMP  ao_overlap_x,ao_overlap_y,ao_overlap_z,ao_overlap,ao_num,ao_coef_normalized_ordered_transp,ao_nucl, &
!  !$OMP  ao_expo_ordered_transp,dim1)
!  do j=1,ao_num
!  A_center(1) = nucl_coord( ao_nucl(j), 1 )
!  A_center(2) = nucl_coord( ao_nucl(j), 2 )
!  A_center(3) = nucl_coord( ao_nucl(j), 3 )
!  power_A(1)  = ao_power( j, 1 )
!  power_A(2)  = ao_power( j, 2 )
!  power_A(3)  = ao_power( j, 3 )
!  do i= 1,ao_num
!    B_center(1) = nucl_coord( ao_nucl(i), 1 )
!    B_center(2) = nucl_coord( ao_nucl(i), 2 )
!    B_center(3) = nucl_coord( ao_nucl(i), 3 )
!    power_B(1)  = ao_power( i, 1 )
!    power_B(2)  = ao_power( i, 2 )
!    power_B(3)  = ao_power( i, 3 )
!    do n = 1,ao_prim_num(j)
!    alpha = ao_expo_ordered_transp(n,j)
!    do l = 1, ao_prim_num(i)
!      beta = ao_expo_ordered_transp(l,i)
!      call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
!      c = ao_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i)
!      ao_overlap(i,j) += c * overlap
!      if(isnan(ao_overlap(i,j)))then
!       print*,'i,j',i,j
!       print*,'l,n',l,n
!       print*,'c,overlap',c,overlap
!       print*,overlap_x,overlap_y,overlap_z
!       stop
!      endif
!      ao_overlap_x(i,j) += c * overlap_x
!      ao_overlap_y(i,j) += c * overlap_y
!      ao_overlap_z(i,j) += c * overlap_z
!    enddo
!    enddo
!  enddo
!  enddo
!  !$OMP END PARALLEL DO

END_PROVIDER

! ---

