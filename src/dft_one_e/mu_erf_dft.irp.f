BEGIN_PROVIDER [double precision, mu_erf_dft]
 implicit none
 BEGIN_DOC
! range separation parameter used in RS-DFT. 
! 
! It is set to mu_erf in order to be consistent with the module "ao_two_e_erf_ints"
 END_DOC
 mu_erf_dft = mu_erf

END_PROVIDER

BEGIN_PROVIDER [double precision, mu_of_r_dft, (n_points_final_grid)]
 implicit none
 integer :: i
 do i = 1, n_points_final_grid
  if(mu_dft_type == "cst")then
   mu_of_r_dft(i) = mu_erf_dft
  else
   mu_of_r_dft(i) = mu_of_r_hf(i)
  endif
 enddo
END_PROVIDER 
