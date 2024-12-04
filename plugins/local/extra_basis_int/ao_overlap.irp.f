
! ---

 BEGIN_PROVIDER [double precision, ao_extra_overlap  , (ao_extra_num, ao_extra_num)]

  BEGIN_DOC
  ! Overlap between atomic basis functions:
  !
  ! :math:`\int \chi_i(r) \chi_j(r) dr`
  END_DOC

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta, c
  double precision :: A_center(3), B_center(3)

  ao_extra_overlap   = 0.d0

   dim1=100
   !$OMP PARALLEL DO SCHEDULE(GUIDED) &
   !$OMP DEFAULT(NONE) &
   !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
   !$OMP  overlap_x,overlap_y, overlap_z, overlap, &
   !$OMP  alpha, beta,i,j,n,l,c) &
   !$OMP SHARED(extra_nucl_coord,ao_extra_power,ao_extra_prim_num, &
   !$OMP  ao_extra_overlap,ao_extra_num,ao_extra_coef_normalized_ordered_transp,ao_extra_nucl, &
   !$OMP  ao_extra_expo_ordered_transp,dim1)
   do j=1,ao_extra_num
   A_center(1) = extra_nucl_coord( ao_extra_nucl(j), 1 )
   A_center(2) = extra_nucl_coord( ao_extra_nucl(j), 2 )
   A_center(3) = extra_nucl_coord( ao_extra_nucl(j), 3 )
   power_A(1)  = ao_extra_power( j, 1 )
   power_A(2)  = ao_extra_power( j, 2 )
   power_A(3)  = ao_extra_power( j, 3 )
   do i= 1,ao_extra_num
     B_center(1) = extra_nucl_coord( ao_extra_nucl(i), 1 )
     B_center(2) = extra_nucl_coord( ao_extra_nucl(i), 2 )
     B_center(3) = extra_nucl_coord( ao_extra_nucl(i), 3 )
     power_B(1)  = ao_extra_power( i, 1 )
     power_B(2)  = ao_extra_power( i, 2 )
     power_B(3)  = ao_extra_power( i, 3 )
     do n = 1,ao_extra_prim_num(j)
     alpha = ao_extra_expo_ordered_transp(n,j)
     do l = 1, ao_extra_prim_num(i)
       beta = ao_extra_expo_ordered_transp(l,i)
       call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
       c = ao_extra_coef_normalized_ordered_transp(n,j) * ao_extra_coef_normalized_ordered_transp(l,i)
       ao_extra_overlap(i,j) += c * overlap
       if(isnan(ao_extra_overlap(i,j)))then
        print*,'i,j',i,j
        print*,'l,n',l,n
        print*,'c,overlap',c,overlap
        print*,overlap_x,overlap_y,overlap_z
        stop
       endif
     enddo
     enddo
   enddo
   enddo
   !$OMP END PARALLEL DO



END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, ao_extra_overlap_mixed  , (ao_num, ao_extra_num)]

  BEGIN_DOC
  ! Overlap between atomic basis functions:
  !
  ! <AO_i|AO_j extra basis>
  END_DOC

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta, c
  double precision :: A_center(3), B_center(3)

  ao_extra_overlap_mixed   = 0.d0

   dim1=100
   !$OMP PARALLEL DO SCHEDULE(GUIDED) &
   !$OMP DEFAULT(NONE) &
   !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
   !$OMP  overlap_x,overlap_y, overlap_z, overlap, &
   !$OMP  alpha, beta,i,j,n,l,c) &
   !$OMP SHARED(extra_nucl_coord,ao_extra_power,ao_extra_prim_num, &
   !$OMP  ao_extra_overlap_mixed,ao_extra_num,ao_extra_coef_normalized_ordered_transp,ao_extra_nucl, &
   !$OMP  ao_extra_expo_ordered_transp,dim1, &
   !$OMP  nucl_coord,ao_power,ao_prim_num, &
   !$OMP  ao_num,ao_coef_normalized_ordered_transp,ao_nucl, &
   !$OMP  ao_expo_ordered_transp)

   do j=1,ao_extra_num
   A_center(1) = extra_nucl_coord( ao_extra_nucl(j), 1 )
   A_center(2) = extra_nucl_coord( ao_extra_nucl(j), 2 )
   A_center(3) = extra_nucl_coord( ao_extra_nucl(j), 3 )
   power_A(1)  = ao_extra_power( j, 1 )
   power_A(2)  = ao_extra_power( j, 2 )
   power_A(3)  = ao_extra_power( j, 3 )
   do i= 1,ao_num
     B_center(1) = nucl_coord( ao_nucl(i), 1 )
     B_center(2) = nucl_coord( ao_nucl(i), 2 )
     B_center(3) = nucl_coord( ao_nucl(i), 3 )
     power_B(1)  = ao_power( i, 1 )
     power_B(2)  = ao_power( i, 2 )
     power_B(3)  = ao_power( i, 3 )
     do n = 1,ao_extra_prim_num(j)
      alpha = ao_extra_expo_ordered_transp(n,j)
      do l = 1, ao_prim_num(i)
        beta = ao_expo_ordered_transp(l,i)
        call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
        c = ao_extra_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i)
        ao_extra_overlap_mixed(i,j) += c * overlap
        if(isnan(ao_extra_overlap_mixed(i,j)))then
         print*,'i,j',i,j
         print*,'l,n',l,n
         print*,'c,overlap',c,overlap
         print*,overlap_x,overlap_y,overlap_z
         stop
        endif
      enddo
     enddo
    enddo
   enddo
   !$OMP END PARALLEL DO


END_PROVIDER

! ---

