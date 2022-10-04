subroutine overlap_gauss_xyz_r12_ao(D_center,delta,i,j,gauss_ints)
 implicit none
 BEGIN_DOC
! gauss_ints(m) = \int dr AO_i(r) AO_j(r) x/y/z e^{-delta |r-D_center|^2}
!
! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in)  :: D_center(3), delta
 double precision, intent(out) :: gauss_ints(3)

 integer :: num_a,num_b,power_A(3), power_B(3),l,k,m
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta,gauss_ints_tmp(3)
 gauss_ints = 0.d0
 if(ao_overlap_abs(j,i).lt.1.d-12)then
  return
 endif
 num_A = ao_nucl(i)
 power_A(1:3)= ao_power(i,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j)
 power_B(1:3)= ao_power(j,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 do l=1,ao_prim_num(i)
  alpha = ao_expo_ordered_transp(l,i)     
  do k=1,ao_prim_num(j)
   beta = ao_expo_ordered_transp(k,j)     
   call overlap_gauss_xyz_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,gauss_ints_tmp)
   do m = 1, 3
    gauss_ints(m) += gauss_ints_tmp(m) *  ao_coef_normalized_ordered_transp(l,i)             &
                                       *  ao_coef_normalized_ordered_transp(k,j) 
   enddo
  enddo 
 enddo
end



double precision function overlap_gauss_xyz_r12_ao_specific(D_center,delta,i,j,mx)
 implicit none
 BEGIN_DOC
! \int dr AO_i(r) AO_j(r) x/y/z e^{-delta |r-D_center|^2}
!
! with mx == 1 ==> x, mx == 2 ==> y, mx == 3 ==> z
 END_DOC
 integer, intent(in) :: i,j,mx
 double precision, intent(in)  :: D_center(3), delta

 integer :: num_a,num_b,power_A(3), power_B(3),l,k
 double precision :: gauss_int
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta
 double precision :: overlap_gauss_xyz_r12_specific
 overlap_gauss_xyz_r12_ao_specific = 0.d0
 if(ao_overlap_abs(j,i).lt.1.d-12)then
  return
 endif
 num_A = ao_nucl(i)
 power_A(1:3)= ao_power(i,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j)
 power_B(1:3)= ao_power(j,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 do l=1,ao_prim_num(i)
  alpha = ao_expo_ordered_transp(l,i)     
  do k=1,ao_prim_num(j)
   beta = ao_expo_ordered_transp(k,j)     
   gauss_int = overlap_gauss_xyz_r12_specific(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,mx)
   overlap_gauss_xyz_r12_ao_specific = gauss_int *  ao_coef_normalized_ordered_transp(l,i)             &
                                                 *  ao_coef_normalized_ordered_transp(k,j) 
  enddo 
 enddo
end


subroutine overlap_gauss_r12_all_ao(D_center,delta,aos_ints)
 implicit none
 double precision, intent(in) :: D_center(3), delta
 double precision, intent(out):: aos_ints(ao_num,ao_num)

 integer :: num_a,num_b,power_A(3), power_B(3),l,k,i,j
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta,analytical_j
 aos_ints = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   if(ao_overlap_abs(j,i).lt.1.d-12)cycle
   num_A = ao_nucl(i)
   power_A(1:3)= ao_power(i,1:3)
   A_center(1:3) = nucl_coord(num_A,1:3)
   num_B = ao_nucl(j)
   power_B(1:3)= ao_power(j,1:3)
   B_center(1:3) = nucl_coord(num_B,1:3)
   do l=1,ao_prim_num(i)
    alpha = ao_expo_ordered_transp(l,i)     
    do k=1,ao_prim_num(j)
     beta = ao_expo_ordered_transp(k,j)     
     analytical_j = overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
     aos_ints(j,i) += analytical_j *  ao_coef_normalized_ordered_transp(l,i)             &
                                   *  ao_coef_normalized_ordered_transp(k,j) 
    enddo 
   enddo
  enddo
 enddo
end

double precision function overlap_gauss_r12_ao(D_center,delta,i,j)
 implicit none
 BEGIN_DOC
! \int dr AO_i(r) AO_j(r) e^{-delta |r-D_center|^2}
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in) :: D_center(3), delta

 integer :: num_a,num_b,power_A(3), power_B(3),l,k
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta,analytical_j
 overlap_gauss_r12_ao = 0.d0
 if(ao_overlap_abs(j,i).lt.1.d-12)then
  return
 endif
 num_A = ao_nucl(i)
 power_A(1:3)= ao_power(i,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j)
 power_B(1:3)= ao_power(j,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 do l=1,ao_prim_num(i)
  alpha = ao_expo_ordered_transp(l,i)     
  do k=1,ao_prim_num(j)
   beta = ao_expo_ordered_transp(k,j)     
   analytical_j = overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
   overlap_gauss_r12_ao += analytical_j *  ao_coef_normalized_ordered_transp(l,i)             &
                                        *  ao_coef_normalized_ordered_transp(k,j) 
  enddo 
 enddo
end


