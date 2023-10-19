BEGIN_PROVIDER [double precision, mu_erf_dft]
 implicit none
 BEGIN_DOC
! range separation parameter used in RS-DFT. 
! 
! It is set to mu_erf in order to be consistent with the module "hamiltonian"
 END_DOC
 mu_erf_dft = mu_erf

END_PROVIDER

BEGIN_PROVIDER [double precision, mu_of_r_dft, (n_points_final_grid)]
 implicit none
 integer :: i
 if(mu_dft_type == "Read")then
   call ezfio_get_mu_of_r_mu_of_r_disk(mu_of_r_dft)
 else
  do i = 1, n_points_final_grid
   if(mu_dft_type == "cst")then
    mu_of_r_dft(i) = mu_erf_dft
   else if(mu_dft_type == "hf")then
    mu_of_r_dft(i) = mu_of_r_hf(i)
   else if(mu_dft_type == "rsc")then
    mu_of_r_dft(i) = mu_rsc_of_r(i)
   else if(mu_dft_type == "grad_rho")then
    mu_of_r_dft(i) = mu_grad_rho(i)
   else 
    print*,'mu_dft_type is not of good type = ',mu_dft_type
    print*,'it must be of type Read, cst, hf, rsc'
    print*,'Stopping ...'
    stop
   endif
  enddo
 endif
END_PROVIDER 

BEGIN_PROVIDER [double precision, mu_rsc_of_r, (n_points_final_grid)]
 implicit none
 integer :: i
 double precision :: mu_rs_c,rho,r(3), dm_a, dm_b
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
  rho = dm_a + dm_b
  mu_rsc_of_r(i) = mu_rs_c(rho)
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, mu_grad_rho, (n_points_final_grid)]
 implicit none
 integer :: i
 double precision :: mu_grad_rho_func, r(3)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  mu_grad_rho(i) = mu_grad_rho_func(r)
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, mu_of_r_dft_average]
 implicit none
 integer :: i
 double precision :: mu_rs_c,rho,r(3), dm_a, dm_b
 mu_of_r_dft_average = 0.d0
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
  rho = dm_a + dm_b
  if(mu_of_r_dft(i).gt.1.d+3)cycle
  mu_of_r_dft_average += rho * mu_of_r_dft(i) * final_weight_at_r_vector(i)
 enddo
 mu_of_r_dft_average = mu_of_r_dft_average / dble(elec_alpha_num + elec_beta_num)
 print*,'mu_of_r_dft_average = ',mu_of_r_dft_average
END_PROVIDER 
