program pouet
 implicit none
 call test
end
subroutine test
 implicit none
! provide mos_times_cholesky_r1
! provide mos_times_cholesky_r2
 integer :: ipoint
 double precision :: accu,weight
 accu = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
!  accu += dabs(mu_of_r_hf(ipoint) - mu_of_r_hf_old(ipoint)) * weight
  accu += dabs(f_hf_sparse_cholesky(ipoint) - f_hf_cholesky(ipoint)) * weight
 enddo
 print*,'accu = ',accu
end
