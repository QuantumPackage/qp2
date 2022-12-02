
 BEGIN_PROVIDER [double precision, mu_of_r_extra_grid_test_func , (n_points_extra_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 double precision :: weight, rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs
 double precision :: elec_a,elec_b
 double precision :: r(3),dx,mu_test_func,mu_tmp,mu,damped_mu, mu_min
 dx = 1.d-5
 elec_a = 0.d0
 elec_b = 0.d0
 mu_min = mu_erf
 do ipoint = 1, n_points_extra_final_grid
  weight = final_weight_at_r_vector_extra(ipoint)
  r(:) = final_grid_points_extra(:,ipoint)
  mu = mu_test_func(r)
  mu_of_r_extra_grid_test_func(ipoint,1) = mu

  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight

 enddo 

END_PROVIDER 




 BEGIN_PROVIDER [double precision, average_mu_test_func      ]
&BEGIN_PROVIDER [double precision, mu_of_r_test_func , (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_test_func , (3,n_points_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 double precision :: weight, rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs
 double precision :: elec_a,elec_b,grad_mu(3),mu_test_func
 double precision :: mu, r(3),dx,damped_mu,mu_min,mu_tmp
 mu_min = mu_erf
 dx = 1.d-5
 average_mu_test_func = 0.d0
 elec_a = 0.d0
 elec_b = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  r(:) = final_grid_points(:,ipoint)
  mu_tmp = mu_test_func(r)

  mu_of_r_test_func(ipoint,1) = mu_tmp

  average_mu_test_func    +=  mu_of_r_test_func(ipoint,1) * weight * rho_hf
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight

  call get_grad_mu_test_func(r,dx,grad_mu)
  if(mu_tmp.lt.mu)then
   grad_mu = 0.d0
  endif
  grad_mu_of_r_test_func(:,ipoint,1) = grad_mu(:)

 enddo 
 average_mu_test_func    = average_mu_test_func   / dble(elec_a+ elec_b)

END_PROVIDER 


double precision function mu_test_func(r)
 implicit none
 double precision, intent(in) :: r(3)
 double precision :: x,y,z
 x = r(1) 
 y = r(2) 
 z = r(3)
 if(mu_test_choice == "cos")then
  mu_test_func = mu_erf + ampl_cos * dcos(omega_cos * x) * dcos(omega_cos * y) * dcos(omega_cos * z) !& 
!                        * dexp(-dexp_gauss * (x**2 + y**2 + z**2))
 else if(mu_test_choice== "gauss" )then
  mu_test_func = mu_erf + ampl_cos * dexp(-dexp_gauss * (x**2 + y**2 + z**2))
 endif
 mu_test_func = max(mu_of_r_min,mu_test_func)
end

subroutine get_grad_mu_test_func(r,dx,grad_mu)
 implicit none
 double precision, intent(in) :: r(3),dx
 double precision, intent(out):: grad_mu(3)
 double precision :: x,y,z
 if(mu_test_choice == "cos")then
  x = r(1) 
  y = r(2) 
  z = r(3)
  grad_mu(1) = -ampl_cos * omega_cos * dsin(omega_cos * x) * dcos(omega_cos * y) * dcos(omega_cos * z)
  grad_mu(2) = -ampl_cos * omega_cos * dcos(omega_cos * x) * dsin(omega_cos * y) * dcos(omega_cos * z)
  grad_mu(3) = -ampl_cos * omega_cos * dcos(omega_cos * x) * dcos(omega_cos * y) * dsin(omega_cos * z)
 else
  call get_grad_mu_test_num(r,dx,grad_mu) 
 endif
end

subroutine get_grad_mu_test_num(r,dx,grad_mu)
 implicit none
 double precision, intent(in) :: r(3), dx
 double precision, intent(out):: grad_mu(3)
 double precision :: r1(3),mu_plus,mu_minus,mu_test_func
 integer :: m
 do m = 1, 3 ! compute grad mu
  r1 = r
  r1(m) += dx 
  mu_plus = mu_test_func(r1)
  r1 = r 
  r1(m) -= dx 
  mu_minus = mu_test_func(r1)
  grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
 enddo

end
