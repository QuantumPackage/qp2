
 BEGIN_PROVIDER [double precision, mu_of_r_extra_grid_lda , (n_points_extra_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 double precision :: weight, rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs
 double precision :: elec_a,elec_b
 double precision :: r(3),dx,mu_lda,mu_tmp,mu,damped_mu, mu_min
 dx = 1.d-5
 elec_a = 0.d0
 elec_b = 0.d0
 mu_min = mu_erf
 do ipoint = 1, n_points_extra_final_grid
  weight = final_weight_at_r_vector_extra(ipoint)
  r(:) = final_grid_points_extra(:,ipoint)
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  mu_tmp = mu_lda(rho_a_hf,rho_b_hf)
  if(damped_mu_of_r)then
   mu = damped_mu(mu_tmp,mu_min)
  else
   mu = max(mu_tmp,mu_of_r_min)
  endif
  if(rescaled_on_top_mu)then
   mu = mu * dsqrt(0.25d0 * (rho_a_hf+rho_b_hf)**2/(rho_a_hf*rho_b_hf+1.d-12))
  else 
   mu=mu
  endif
  mu_of_r_extra_grid_lda(ipoint,1) = mu

  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight

 enddo 

END_PROVIDER 




 BEGIN_PROVIDER [double precision, average_mu_lda      ]
&BEGIN_PROVIDER [double precision, mu_of_r_lda , (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_lda , (3,n_points_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 double precision :: weight, rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs
 double precision :: elec_a,elec_b,grad_mu(3),mu_lda
 double precision :: mu, r(3),dx,damped_mu,mu_min,mu_tmp
 mu_min = mu_erf
 dx = 1.d-5
 average_mu_lda = 0.d0
 elec_a = 0.d0
 elec_b = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  rho_a_hf = one_e_dm_and_grad_alpha_in_r(4,ipoint,1)
  rho_b_hf = one_e_dm_and_grad_beta_in_r(4,ipoint,1)
  rho_hf = rho_a_hf + rho_b_hf
  mu_tmp = mu_lda(rho_a_hf,rho_b_hf)
  if(damped_mu_of_r)then
   mu = damped_mu(mu_tmp,mu_min)
   mu = max(mu_tmp,mu_of_r_min)
  else
   mu = mu_tmp 
  endif
  if(rescaled_on_top_mu)then
   mu = mu * dsqrt(0.25d0 * (rho_a_hf+rho_b_hf)**2/(rho_a_hf*rho_b_hf+1.d-12))
  else 
   mu=mu
  endif

  mu_of_r_lda(ipoint,1) = mu

  average_mu_lda    +=  mu_of_r_lda(ipoint,1) * weight * rho_hf
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight

  r(:) = final_grid_points(:,ipoint) 
  if(damped_mu_of_r)then
   call get_grad_damped_mu_lda(r,dx,mu_min,grad_mu)
  else
   call get_grad_mu_lda(r,dx,grad_mu)
   if(mu_tmp.lt.mu)then
    grad_mu = 0.d0
   endif
  endif
  grad_mu_of_r_lda(:,ipoint,1) = grad_mu(:)

 enddo 
 average_mu_lda    = average_mu_lda   / dble(elec_a+ elec_b)

END_PROVIDER 


double precision function mu_lda(rho_a,rho_b)
 implicit none 
 double precision, intent(in) :: rho_a,rho_b
 include 'constants.include.F'
 double precision :: g0,g0_UEG_mu_inf
 g0 = g0_UEG_mu_inf(rho_a,rho_b)
 mu_lda = - 1.d0 / (dlog(2.d0 * g0) * sqpi) 
end

double precision function damped_mu(mu,mu_min)
 implicit none
 double precision, intent(in) :: mu, mu_min
 damped_mu = mu_min * (1.d0 - derf(mu)) + derf(mu)*mu
end

double precision function mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
 implicit none
 double precision, intent(in) :: rho_b_hf,rho_a_hf,mu_min
 double precision :: mu_lda,mu,damped_mu
 mu = mu_lda(rho_a_hf,rho_b_hf)
 mu_lda_damped = damped_mu(mu,mu_min)
end

subroutine get_grad_mu_lda(r,dx,grad_mu)
  implicit none
  double precision, intent(in) :: r(3),dx
  double precision, intent(out):: grad_mu(3)
  double precision :: r1(3),rho_a_hf,rho_b_hf,mu_plus,mu_minus,mu_lda
  integer :: m
  do m = 1, 3 ! compute grad mu
   r1 = r
   r1(m) += dx 
   call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
   mu_plus = mu_lda(rho_a_hf,rho_b_hf)
   r1 = r 
   r1(m) -= dx 
   call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
   mu_minus = mu_lda(rho_a_hf,rho_b_hf)
   grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
  enddo
end

subroutine get_grad_damped_mu_lda(r,dx,mu_min,grad_mu)
  implicit none
  double precision, intent(in) :: r(3),dx,mu_min
  double precision, intent(out):: grad_mu(3)
  double precision :: r1(3),rho_a_hf,rho_b_hf,mu_plus,mu_minus,mu_lda_damped
  integer :: m
  do m = 1, 3 ! compute grad mu
   r1 = r
   r1(m) += dx 
   call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
   mu_plus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
   r1 = r 
   r1(m) -= dx 
   call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
   mu_minus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
   grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
  enddo
end

double precision function grad_mu_lda_comp(r,dx,mu_min,m)
  implicit none
  double precision, intent(in) :: r(3),dx,mu_min
  integer, intent(in)          :: m
  double precision :: r1(3),rho_a_hf,rho_b_hf,mu_plus,mu_minus,mu_lda_damped
   r1 = r
   r1(m) += dx 
   call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
   mu_plus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
   r1 = r 
   r1(m) -= dx 
   call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
   mu_minus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
   grad_mu_lda_comp = (mu_plus - mu_minus)/(2.d0 * dx)
end

subroutine get_lapl_mu_lda(r,dx,mu_min,lapl_mu)
 implicit none
  double precision, intent(in) :: r(3),dx,mu_min
  double precision, intent(out):: lapl_mu(3)
  double precision :: r1(3),rho_a_hf,rho_b_hf,mu_plus,mu_minus,mu_lda_damped,mu
  integer :: m
  do m = 1, 3
   r1 = r
   call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
   mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)

   r1 = r
   r1(m) += 2.d0 * dx 
   call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
   mu_plus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)

   r1 = r 
   r1(m) -= 2.d0 * dx 
   call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
   mu_minus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)

   lapl_mu(m) = (mu_plus + mu_minus - 2.d0 * mu)/(4.d0 * dx * dx)
  enddo
  lapl_mu = 0.d0
end

BEGIN_PROVIDER [ double precision, average_mu_rs_c_lda]
 implicit none
 average_mu_rs_c_lda = 0.5d0 * (average_mu_rs_c + average_mu_lda)
END_PROVIDER 



