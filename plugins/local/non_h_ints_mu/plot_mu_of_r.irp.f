program plot_mu_of_r
 implicit none
 read_wf = .False.
 touch read_wf 
 call routine_print

end


subroutine routine_print
 implicit none
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.mu_of_r'
 i_unit_output = getUnitAndOpen(output,'w')
 integer :: ipoint,nx,i
 double precision :: xmax,xmin,r(3),dx,sigma
 double precision :: mu_val, mu_der(3),dm_a,dm_b,grad,grad_dm_a(3), grad_dm_b(3)
 xmax =  5.D0
 xmin = -5.D0
 nx = 10000
 dx = (xmax - xmin)/dble(nx)
 r = 0.d0
 r(1) = xmin
 do ipoint = 1, nx
  call mu_r_val_and_grad(r, r, mu_val, mu_der)
  call density_and_grad_alpha_beta(r,dm_a,dm_b, grad_dm_a, grad_dm_b)
  sigma = 0.d0
  do i = 1,3
   sigma += grad_dm_a(i)**2
  enddo
  sigma=dsqrt(sigma)
  grad = mu_der(1)**2 + mu_der(2)**2 + mu_der(3)**2 
  grad = dsqrt(grad)
  write(i_unit_output,'(100(F16.7,X))')r(1),mu_val,dm_a+dm_b,grad,sigma/dm_a
  r(1) += dx
 enddo
end
