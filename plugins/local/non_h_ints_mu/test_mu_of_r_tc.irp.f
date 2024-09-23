program test_mu_of_r_tc
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  ! You specify that you want to avoid any contribution from 
  ! orbitals coming from core 
 call test_grad_f_mean_field
 call test_grad_mu_mf
 call plot_mu_of_r_mf
end


subroutine test_grad_f_mean_field
 implicit none
 integer :: i_point,k
 double precision :: weight,r(3)
 double precision :: grad_f_mf_ab(3), grad_two_bod_dens(3)
 double precision :: grad_dm_a(3), grad_dm_b(3)
 double precision :: f_mf_ab,two_bod_dens, dm_a, dm_b

 double precision :: num_grad_f_mf_ab(3), num_grad_two_bod_dens(3)
 double precision :: num_grad_dm_a(3), num_grad_dm_b(3)
 double precision :: f_mf_ab_p,f_mf_ab_m
 double precision :: two_bod_dens_p, two_bod_dens_m
 double precision :: dm_a_p, dm_a_m
 double precision :: dm_b_p, dm_b_m
 double precision :: rbis(3), dr
 double precision :: accu_grad_f_mf_ab(3),accu_grad_two_bod_dens(3)
 double precision :: accu_grad_dm_a(3),accu_grad_dm_b(3)
 double precision :: accu_f_mf_ab, accu_two_bod_dens, accu_dm_a, accu_dm_b
 dr = 0.00001d0
 accu_f_mf_ab = 0.d0 
 accu_two_bod_dens = 0.d0 
 accu_dm_a = 0.d0 
 accu_dm_b = 0.d0

 accu_grad_f_mf_ab = 0.d0
 accu_grad_two_bod_dens = 0.d0
 accu_grad_dm_a = 0.d0
 accu_grad_dm_b = 0.d0
 do i_point = 1, n_points_final_grid
  r(1:3)   = final_grid_points(1:3,i_point)
  weight = final_weight_at_r_vector(i_point)
  call get_grad_f_mf_ab(r,grad_f_mf_ab, grad_two_bod_dens,f_mf_ab,two_bod_dens, dm_a, dm_b,grad_dm_a, grad_dm_b)
  call get_f_mf_ab(r,f_mf_ab_p,two_bod_dens_p, dm_a_p, dm_b_p)
  accu_f_mf_ab += weight * dabs(f_mf_ab - f_mf_ab_p)
  accu_two_bod_dens += weight * dabs(two_bod_dens - two_bod_dens_p)
  accu_dm_a += weight*dabs(dm_a - dm_a_p)
  accu_dm_b += weight*dabs(dm_b - dm_b_p)
  do k = 1, 3
   rbis = r
   rbis(k) += dr
   call get_f_mf_ab(rbis,f_mf_ab_p,two_bod_dens_p, dm_a_p, dm_b_p)
   rbis = r
   rbis(k) -= dr
   call get_f_mf_ab(rbis,f_mf_ab_m,two_bod_dens_m, dm_a_m, dm_b_m)
   num_grad_f_mf_ab(k) = (f_mf_ab_p - f_mf_ab_m)/(2.d0*dr)
   num_grad_two_bod_dens(k) = (two_bod_dens_p - two_bod_dens_m)/(2.d0*dr)
   num_grad_dm_a(k) = (dm_a_p - dm_a_m)/(2.d0*dr)
   num_grad_dm_b(k) = (dm_b_p - dm_b_m)/(2.d0*dr)
  enddo
  do k = 1, 3
   accu_grad_f_mf_ab(k) += weight * dabs(grad_f_mf_ab(k) - num_grad_f_mf_ab(k))
   accu_grad_two_bod_dens(k) += weight * dabs(grad_two_bod_dens(k) - num_grad_two_bod_dens(k))
   accu_grad_dm_a(k) += weight * dabs(grad_dm_a(k) - num_grad_dm_a(k))
   accu_grad_dm_b(k) += weight * dabs(grad_dm_b(k) - num_grad_dm_b(k))
  enddo
 enddo
 print*,'accu_f_mf_ab = ',accu_f_mf_ab
 print*,'accu_two_bod_dens = ',accu_two_bod_dens
 print*,'accu_dm_a = ',accu_dm_a
 print*,'accu_dm_b = ',accu_dm_b
 print*,'accu_grad_f_mf_ab = '
 print*,accu_grad_f_mf_ab
 print*,'accu_grad_two_bod_dens = '
 print*,accu_grad_two_bod_dens
 print*,'accu_dm_a = '
 print*,accu_grad_dm_a
 print*,'accu_dm_b = '
 print*,accu_grad_dm_b

end

subroutine test_grad_mu_mf
 implicit none
 integer :: i_point,k
 double precision :: weight,r(3),rbis(3)
 double precision :: mu_mf, dm,grad_mu_mf(3), grad_dm(3)
 double precision :: mu_mf_p, mu_mf_m, dm_m, dm_p, num_grad_mu_mf(3),dr, num_grad_dm(3)
 double precision :: accu_mu, accu_dm, accu_grad_dm(3), accu_grad_mu_mf(3)
 dr = 0.00001d0
 accu_grad_mu_mf = 0.d0
 accu_mu = 0.d0
 accu_grad_dm = 0.d0
 accu_dm = 0.d0
 do i_point = 1, n_points_final_grid
  r(1:3)   = final_grid_points(1:3,i_point)
  weight = final_weight_at_r_vector(i_point)
  call grad_mu_of_r_mean_field(r,mu_mf, dm, grad_mu_mf, grad_dm)
  call mu_of_r_mean_field(r,mu_mf_p, dm_p)
  accu_mu += weight*dabs(mu_mf_p - mu_mf)
  accu_dm += weight*dabs(dm_p - dm)
  do k = 1, 3
   rbis = r
   rbis(k) += dr
   call mu_of_r_mean_field(rbis,mu_mf_p, dm_p)
   rbis = r
   rbis(k) -= dr
   call mu_of_r_mean_field(rbis,mu_mf_m, dm_m)

   num_grad_mu_mf(k) = (mu_mf_p - mu_mf_m)/(2.d0*dr)
   num_grad_dm(k) = (dm_p - dm_m)/(2.d0*dr)
  enddo
  do k = 1, 3
   accu_grad_dm(k)+= weight *dabs(num_grad_dm(k) - grad_dm(k))
   accu_grad_mu_mf(k)+= weight *dabs(num_grad_mu_mf(k) - grad_mu_mf(k))
  enddo
 enddo
 print*,'accu_mu = ',accu_mu
 print*,'accu_dm = ',accu_dm
 print*,'accu_grad_dm = '
 print*, accu_grad_dm
 print*,'accu_grad_mu_mf = '
 print*, accu_grad_mu_mf

end

subroutine plot_mu_of_r_mf
 implicit none
  include 'constants.include.F'
 integer :: ipoint,npoint
 double precision :: dx,r(3),xmax,xmin
 double precision :: accu_mu,accu_nelec,mu_mf, dm,mu_mf_tc
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.mu_mf'
 i_unit_output = getUnitAndOpen(output,'w')
 xmax = 5.D0
 xmin = 0.d0
 npoint = 10000
 dx = (xmax - xmin)/dble(npoint)
 r = 0.d0
 r(1) = xmin
 accu_mu = 0.d0
 accu_nelec = 0.d0
 do ipoint = 1, npoint
  call mu_of_r_mean_field(r,mu_mf, dm)
  call mu_of_r_mean_field_tc(r,mu_mf_tc, dm)
  write(i_unit_output,'(100(F16.10,X))')r(1),mu_mf,mu_mf_tc,dm
  accu_mu    += mu_mf * dm * r(1)**2*dx*4.D0*pi
  accu_nelec +=         dm * r(1)**2*dx*4.D0*pi
  r(1) += dx
 enddo
 print*,'nelec      = ',accu_nelec
 print*,'mu average = ',accu_mu/accu_nelec
end
