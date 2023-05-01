program pt2_tc_cisd

  BEGIN_DOC
  !
  ! TODO : Reads psi_det in the EZFIO folder and prints out the left- and right-eigenvectors together 
  !        with the energy. Saves the left-right wave functions at the end. 
  !
  END_DOC

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  print*, ' nb of states = ', N_states
  print*, ' nb of det    = ', N_det

  call routine
end

subroutine routine
 implicit none
 integer :: i,h1,p1,h2,p2,s1,s2
 double precision :: h0i,hi0,e00,ei,delta_e
 double precision :: norm,e_corr,coef
 norm = 0.d0 
 e_corr = 0.d0
 call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,1), psi_det(1,1,1), N_int, e00) 
 do i = 2, N_det
  call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,i), psi_det(1,1,1), N_int, hi0) 
  call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,1), psi_det(1,1,i), N_int, h0i) 
  call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,i), psi_det(1,1,i), N_int, ei) 
  delta_e = e00 - ei
  coef = hi0/delta_e
  norm += coef*coef
  e_corr += dabs(coef* h0i)
 enddo
 print*,'e_corr = ',e_corr
 print*,'norm   = ',norm

end
