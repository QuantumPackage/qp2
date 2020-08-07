double precision function ao_value(i,r)
 implicit none
 BEGIN_DOC
! Returns the value of the i-th ao at point $\textbf{r}$
 END_DOC
 double precision, intent(in) :: r(3)
 integer, intent(in) :: i

 integer :: m,num_ao
 double precision :: center_ao(3)
 double precision :: beta
 integer :: power_ao(3)
 double precision :: accu,dx,dy,dz,r2
 num_ao = ao_nucl(i)
 power_ao(1:3)= ao_power(i,1:3)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 dx = (r(1) - center_ao(1))
 dy = (r(2) - center_ao(2))
 dz = (r(3) - center_ao(3))
 r2 = dx*dx + dy*dy + dz*dz
 dx = dx**power_ao(1)
 dy = dy**power_ao(2)
 dz = dz**power_ao(3)

 accu = 0.d0
 do m=1,ao_prim_num(i)
   beta = ao_expo_ordered_transp(m,i)
   accu += ao_coef_normalized_ordered_transp(m,i) * dexp(-beta*r2)
 enddo
 ao_value = accu * dx * dy * dz

end


double precision function primitive_value(i,j,r)
 implicit none
 BEGIN_DOC
! Returns the value of the j-th primitive of the i-th |AO| at point $\textbf{r}
! **without the coefficient**
 END_DOC
 double precision, intent(in) :: r(3)
 integer, intent(in) :: i,j

 integer :: m,num_ao
 double precision :: center_ao(3)
 double precision :: beta
 integer :: power_ao(3)
 double precision :: accu,dx,dy,dz,r2
 num_ao = ao_nucl(i)
 power_ao(1:3)= ao_power(i,1:3)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 dx = (r(1) - center_ao(1))
 dy = (r(2) - center_ao(2))
 dz = (r(3) - center_ao(3))
 r2 = dx*dx + dy*dy + dz*dz
 dx = dx**power_ao(1)
 dy = dy**power_ao(2)
 dz = dz**power_ao(3)

 accu = 0.d0
 m=j
 beta = ao_expo_ordered_transp(m,i)
 accu += dexp(-beta*r2)
 primitive_value = accu * dx * dy * dz

end


subroutine give_all_aos_at_r(r,aos_array)
 implicit none
 BEGIN_dOC
! input  : r == r(1) = x and so on
!
! output : aos_array(i) = aos(i) evaluated in $\textbf{r}$
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out):: aos_array(ao_num)

 integer :: power_ao(3)
 integer :: i,j,k,l,m
 double precision :: dx,dy,dz,r2
 double precision ::      dx2,dy2,dz2
 double precision :: center_ao(3)
 double precision :: beta
 do i = 1, nucl_num
  center_ao(1:3) = nucl_coord(i,1:3)
  dx = (r(1) - center_ao(1))
  dy = (r(2) - center_ao(2))
  dz = (r(3) - center_ao(3))
  r2 = dx*dx + dy*dy + dz*dz
  do j = 1,Nucl_N_Aos(i)
   k = Nucl_Aos_transposed(j,i) ! index of the ao in the ordered format
   aos_array(k) = 0.d0
   power_ao(1:3)= ao_power_ordered_transp_per_nucl(1:3,j,i)
   dx2 = dx**power_ao(1)
   dy2 = dy**power_ao(2)
   dz2 = dz**power_ao(3)
   do l = 1,ao_prim_num(k)
    beta = ao_expo_ordered_transp_per_nucl(l,j,i)
    if(dabs(beta*r2).gt.40.d0)cycle
    aos_array(k)+= ao_coef_normalized_ordered_transp_per_nucl(l,j,i) * dexp(-beta*r2)
   enddo
   aos_array(k) = aos_array(k) * dx2 * dy2 * dz2
  enddo
 enddo
end


subroutine give_all_aos_and_grad_at_r(r,aos_array,aos_grad_array)
 implicit none
 BEGIN_DOC
! input : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! output : 
!
! * aos_array(i) = ao(i) evaluated at ro
! * aos_grad_array(1,i) = gradient X of the ao(i) evaluated at $\textbf{r}$
!
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: aos_array(ao_num)
 double precision, intent(out) :: aos_grad_array(3,ao_num)

 integer :: power_ao(3)
 integer :: i,j,k,l,m
 double precision :: dx,dy,dz,r2
 double precision ::      dx2,dy2,dz2
 double precision ::      dx1,dy1,dz1
 double precision :: center_ao(3)
 double precision :: beta,accu_1,accu_2,contrib
 do i = 1, nucl_num
  center_ao(1:3) = nucl_coord(i,1:3)
  dx = (r(1) - center_ao(1))
  dy = (r(2) - center_ao(2))
  dz = (r(3) - center_ao(3))
  r2 = dx*dx + dy*dy + dz*dz
  do j = 1,Nucl_N_Aos(i)
   k = Nucl_Aos_transposed(j,i) ! index of the ao in the ordered format
   aos_array(k) = 0.d0
   aos_grad_array(1,k) = 0.d0
   aos_grad_array(2,k) = 0.d0
   aos_grad_array(3,k) = 0.d0
   power_ao(1:3)= ao_power_ordered_transp_per_nucl(1:3,j,i)
   dx2 = dx**power_ao(1)
   dy2 = dy**power_ao(2)
   dz2 = dz**power_ao(3)
   if(power_ao(1) .ne. 0)then
    dx1 = dble(power_ao(1)) * dx**(power_ao(1)-1)
   else
    dx1 = 0.d0
   endif
   if(power_ao(2) .ne. 0)then
    dy1 = dble(power_ao(2)) * dy**(power_ao(2)-1)
   else
    dy1 = 0.d0
   endif
   if(power_ao(3) .ne. 0)then
    dz1 = dble(power_ao(3)) * dz**(power_ao(3)-1)
   else
    dz1 = 0.d0
   endif
   accu_1 = 0.d0
   accu_2 = 0.d0
   do l = 1,ao_prim_num(k)
    beta = ao_expo_ordered_transp_per_nucl(l,j,i)
    contrib = 0.d0
    if(beta*r2.gt.50.d0)cycle
    contrib = ao_coef_normalized_ordered_transp_per_nucl(l,j,i) * dexp(-beta*r2)
    accu_1 += contrib
    accu_2 += contrib * beta
   enddo
   aos_array(k) = accu_1 * dx2 * dy2 * dz2
   aos_grad_array(1,k) = accu_1 * dx1  * dy2 * dz2- 2.d0 * dx2 * dx  * dy2 * dz2 * accu_2
   aos_grad_array(2,k) = accu_1 * dx2  * dy1 * dz2- 2.d0 * dx2 * dy2 * dy  * dz2 * accu_2
   aos_grad_array(3,k) = accu_1 * dx2  * dy2 * dz1- 2.d0 * dx2 * dy2 * dz2 * dz  * accu_2
  enddo
 enddo
end


subroutine give_all_aos_and_grad_and_lapl_at_r(r,aos_array,aos_grad_array,aos_lapl_array)
 implicit none
 BEGIN_DOC
! input  : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! output :
!
! * aos_array(i) = ao(i) evaluated at $\textbf{r}$
! * aos_grad_array(1,i) = $\nabla_x$ of the ao(i) evaluated at $\textbf{r}$
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: aos_array(ao_num)
 double precision, intent(out) :: aos_grad_array(3,ao_num)
 double precision, intent(out) :: aos_lapl_array(3,ao_num)

 integer :: power_ao(3)
 integer :: i,j,k,l,m
 double precision :: dx,dy,dz,r2
 double precision ::      dx2,dy2,dz2
 double precision ::      dx1,dy1,dz1
 double precision ::      dx3,dy3,dz3
 double precision ::      dx4,dy4,dz4
 double precision ::      dx5,dy5,dz5
 double precision :: center_ao(3)
 double precision :: beta,accu_1,accu_2,accu_3,contrib
 do i = 1, nucl_num
  center_ao(1:3) = nucl_coord(i,1:3)
  dx = (r(1) - center_ao(1))
  dy = (r(2) - center_ao(2))
  dz = (r(3) - center_ao(3))
  r2 = dx*dx + dy*dy + dz*dz
  do j = 1,Nucl_N_Aos(i)
   k = Nucl_Aos_transposed(j,i) ! index of the ao in the ordered format
   aos_array(k) = 0.d0
   aos_grad_array(1,k) = 0.d0
   aos_grad_array(2,k) = 0.d0
   aos_grad_array(3,k) = 0.d0

   aos_lapl_array(1,k) = 0.d0
   aos_lapl_array(2,k) = 0.d0
   aos_lapl_array(3,k) = 0.d0

   power_ao(1:3)= ao_power_ordered_transp_per_nucl(1:3,j,i)
   dx2 = dx**power_ao(1)
   dy2 = dy**power_ao(2)
   dz2 = dz**power_ao(3)
   if(power_ao(1) .ne. 0)then
    dx1 = dble(power_ao(1)) * dx**(power_ao(1)-1)
   else
    dx1 = 0.d0
   endif
   ! For the Laplacian
   if(power_ao(1) .ge. 2)then
    dx3 = dble(power_ao(1)) * dble((power_ao(1)-1))  * dx**(power_ao(1)-2)
   else
    dx3 = 0.d0
   endif
   if(power_ao(1) .ge. 1)then
    dx4 = dble((2 * power_ao(1) + 1))  * dx**(power_ao(1)) 
   else
    dx4 = dble((power_ao(1) + 1))  * dx**(power_ao(1)) 
   endif

   dx5 = dx**(power_ao(1)+2)

   if(power_ao(2) .ne. 0)then
    dy1 = dble(power_ao(2)) * dy**(power_ao(2)-1)
   else
    dy1 = 0.d0
   endif
   ! For the Laplacian
   if(power_ao(2) .ge. 2)then
    dy3 = dble(power_ao(2)) * dble((power_ao(2)-1))  * dy**(power_ao(2)-2)
   else
    dy3 = 0.d0
   endif

   if(power_ao(2) .ge. 1)then
    dy4 = dble((2 * power_ao(2) + 1))  * dy**(power_ao(2)) 
   else
    dy4 = dble((power_ao(2) + 1))  * dy**(power_ao(2)) 
   endif

   dy5 = dy**(power_ao(2)+2)


   if(power_ao(3) .ne. 0)then
    dz1 = dble(power_ao(3)) * dz**(power_ao(3)-1)
   else
    dz1 = 0.d0
   endif
   ! For the Laplacian
   if(power_ao(3) .ge. 2)then
    dz3 = dble(power_ao(3)) * dble((power_ao(3)-1))  * dz**(power_ao(3)-2)
   else
    dz3 = 0.d0
   endif

   if(power_ao(3) .ge. 1)then
    dz4 = dble((2 * power_ao(3) + 1))  * dz**(power_ao(3)) 
   else
    dz4 = dble((power_ao(3) + 1))  * dz**(power_ao(3)) 
   endif

   dz5 = dz**(power_ao(3)+2)


   accu_1 = 0.d0
   accu_2 = 0.d0
   accu_3 = 0.d0
   do l = 1,ao_prim_num(k)
    beta = ao_expo_ordered_transp_per_nucl(l,j,i)
    contrib = ao_coef_normalized_ordered_transp_per_nucl(l,j,i) * dexp(-beta*r2)
    accu_1 += contrib
    accu_2 += contrib * beta
    accu_3 += contrib * beta**2
   enddo
   aos_array(k) = accu_1 * dx2 * dy2 * dz2

   aos_grad_array(1,k) = accu_1 * dx1  * dy2 * dz2- 2.d0 * dx2 * dx  * dy2 * dz2 * accu_2
   aos_grad_array(2,k) = accu_1 * dx2  * dy1 * dz2- 2.d0 * dx2 * dy2 * dy  * dz2 * accu_2
   aos_grad_array(3,k) = accu_1 * dx2  * dy2 * dz1- 2.d0 * dx2 * dy2 * dz2 * dz  * accu_2

   aos_lapl_array(1,k) = accu_1 * dx3  * dy2 * dz2- 2.d0 * dx4 * dy2 * dz2* accu_2 +4.d0 * dx5 *dy2 * dz2* accu_3
   aos_lapl_array(2,k) = accu_1 * dx2  * dy3 * dz2- 2.d0 * dx2 * dy4 * dz2* accu_2 +4.d0 * dx2 *dy5 * dz2* accu_3
   aos_lapl_array(3,k) = accu_1 * dx2  * dy2 * dz3- 2.d0 * dx2 * dy2 * dz4* accu_2 +4.d0 * dx2 *dy2 * dz5* accu_3

  enddo
 enddo
end


