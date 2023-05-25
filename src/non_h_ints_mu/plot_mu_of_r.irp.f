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
 integer :: ipoint,nx
 double precision :: xmax,xmin,r(3),dx
 double precision :: mu_val, mu_der(3),dm_a,dm_b,grad
 xmax =  5.D0
 xmin = -5.D0
 nx = 10000
 dx = (xmax - xmin)/dble(nx)
 r = 0.d0
 r(1) = xmin
 do ipoint = 1, nx
  call mu_r_val_and_grad(r, r, mu_val, mu_der)
  call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
  grad = mu_der(1)**2 + mu_der(2)**2 + mu_der(3)**2 
  grad = dsqrt(grad)
  write(i_unit_output,'(100(F16.7,X))')r(1),mu_val,dm_a+dm_b,grad
  r(1) += dx
 enddo
end
