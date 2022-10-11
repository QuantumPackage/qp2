
BEGIN_PROVIDER [double precision, mu_of_r_extra_grid_grad_n , (n_points_extra_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 double precision :: weight, rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs
 double precision :: elec_a,elec_b
 double precision :: r(3),dx,mu_grad_n,mu_tmp,mu,damped_mu, mu_min
 dx = 1.d-5
 elec_a = 0.d0
 elec_b = 0.d0
 mu_min = mu_erf
 do ipoint = 1, n_points_extra_final_grid
  weight = final_weight_at_r_vector_extra(ipoint)
  r(:) = final_grid_points_extra(:,ipoint)
  mu_tmp = mu_grad_n(r)
  if(damped_mu_of_r)then
   mu = damped_mu(mu_tmp,mu_min)
  else
   mu = max(mu_tmp,mu_of_r_min)
  endif
  mu_of_r_extra_grid_grad_n(ipoint,1) = mu
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight
 enddo 
END_PROVIDER 

 BEGIN_PROVIDER [double precision, average_mu_grad_n      ]
&BEGIN_PROVIDER [double precision, mu_of_r_grad_n , (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_grad_n , (3,n_points_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 double precision :: weight, rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs
 double precision :: elec_a,elec_b,grad_mu(3),mu_grad_n
 double precision :: mu, r(3),dx,damped_mu,mu_min,mu_tmp
 mu_min = mu_erf
 dx = 1.d-5
 average_mu_grad_n = 0.d0
 elec_a = 0.d0
 elec_b = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  r(:) = final_grid_points(:,ipoint)
  mu_tmp = mu_grad_n(r)
  if(damped_mu_of_r)then
   mu = damped_mu(mu_tmp,mu_min)
   mu = max(mu_tmp,mu_of_r_min)
  else
   mu = mu_tmp
  endif

  mu_of_r_grad_n(ipoint,1) = mu

  average_mu_grad_n    +=  mu_of_r_grad_n(ipoint,1) * weight * rho_hf
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight

  r(:) = final_grid_points(:,ipoint) 
  if(damped_mu_of_r)then
   call get_grad_damped_mu_grad_n(r,dx,mu_min,grad_mu)
  else
   call get_grad_mu_grad_n(r,dx,grad_mu)
   if(mu_tmp.lt.mu)then
    grad_mu = 0.d0
   endif
  endif
  grad_mu_of_r_grad_n(:,ipoint,1) = grad_mu(:)

 enddo 
 average_mu_grad_n    = average_mu_grad_n   / dble(elec_a+ elec_b)

END_PROVIDER 



double precision function mu_grad_n(r)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, allocatable :: aos_array(:),grad_aos_array(:,:)
 double precision, allocatable :: dm_a(:),dm_b(:), dm_a_grad(:,:), dm_b_grad(:,:)
 double precision :: grad_n, rho
 integer :: m
 allocate(dm_a(N_states),dm_b(N_states), dm_a_grad(3,N_states), dm_b_grad(3,N_states))
 allocate(aos_array(ao_num),grad_aos_array(3,ao_num))
 call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a,dm_b,  dm_a_grad, dm_b_grad, aos_array, grad_aos_array)
 grad_n = 0.D0
 do m = 1, 3
  grad_n += dm_a_grad(m,1)**2 + dm_b_grad(m,1)**2 + 2.d0 * dm_a_grad(m,1)*dm_b_grad(m,1)
 enddo
 rho = dm_a(1) + dm_b(1)
 grad_n = dsqrt(grad_n)
 if(dabs(rho).gt.1.d-20)then
  mu_grad_n = grad_n/(4.d0 * rho)
 else
  mu_grad_n = 0.d0
 endif
end

double precision function mu_grad_n_damped(r,mu_min)
 implicit none
 double precision, intent(in) :: r(3), mu_min
 double precision :: mu_grad_n, damped_mu, mu_tmp
 mu_tmp =  mu_grad_n(r)
 mu_grad_n_damped = damped_mu(mu_tmp,mu_min)
end

subroutine get_grad_mu_grad_n(r,dx,grad_mu)
 implicit none
 double precision, intent(in) :: r(3),dx
 double precision, intent(out):: grad_mu(3)
 double precision :: r1(3),mu_plus,mu_minus,mu_grad_n
 integer :: m
 do m = 1, 3 ! compute grad mu
  r1 = r
  r1(m) += dx 
  mu_plus = mu_grad_n(r1)
  r1 = r 
  r1(m) -= dx 
  mu_minus = mu_grad_n(r1)
  grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
 enddo
end

subroutine get_grad_damped_mu_grad_n(r,dx,mu_min,grad_mu)
 implicit none
 double precision, intent(in) :: r(3),dx,mu_min
 double precision, intent(out):: grad_mu(3)
 double precision :: r1(3),mu_plus,mu_minus,mu_grad_n_damped
 integer :: m
 do m = 1, 3 ! compute grad mu
  r1 = r
  r1(m) += dx 
  mu_plus = mu_grad_n_damped(r1,mu_min)
  r1 = r 
  r1(m) -= dx 
  mu_minus = mu_grad_n_damped(r1,mu_min)
  grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
 enddo
end
