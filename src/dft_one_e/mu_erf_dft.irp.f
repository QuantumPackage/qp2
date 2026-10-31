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
 BEGIN_DOC
 ! mu(r) for RS-DFT at the n_points_final_grid points of the final grid.
 ! If mu_dft_type is "Read", read from disk (mu_of_r_disk). Otherwise,
 ! set from mu_erf_dft ("cst"), mu_of_r_hf ("hf"), mu_rsc_of_r ("rsc")
 ! or mu_grad_rho ("grad_rho"); stops otherwise.
 END_DOC
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
 BEGIN_DOC
 ! "rsc" range separation function mu(r) at the n_points_final_grid
 ! points of the final grid: for each point r, mu_rsc_of_r(i) =
 ! mu_rs_c(rho) with rho = dm_a + dm_b the total DFT density at r
 ! (dm_dft_alpha_beta_at_r). mu_rs_c is defined in
 ! dft_utils_func/mu_of_r_dft.irp.f.
 END_DOC
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
 BEGIN_DOC
 ! "grad_rho" range separation function mu(r) at the
 ! n_points_final_grid points of the final grid: mu_grad_rho(i) =
 ! mu_grad_rho_func(r), defined in dft_utils_func/mu_of_r_dft.irp.f
 ! (eta * |grad rho| / rho with eta = mu_erf; 1.d-10 where rho < 1.d-12).
 END_DOC
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
 BEGIN_DOC
 ! Sum over the final grid of rho(i) * mu_of_r_dft(i) *
 ! final_weight_at_r_vector(i) (the grid weights), skipping points
 ! where mu_of_r_dft(i) > 1.d+3, divided by elec_alpha_num +
 ! elec_beta_num. The value is printed; no other in-tree provider
 ! uses it.
 ! The 1.d+3 threshold is a cutoff to avoid unphysical values of
 ! mu_of_r_dft (author, 2026-10-31).
 END_DOC
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
