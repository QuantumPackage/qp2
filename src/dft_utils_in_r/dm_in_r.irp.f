 BEGIN_PROVIDER [double precision, one_e_dm_and_grad_alpha_in_r, (4,n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, one_e_dm_and_grad_beta_in_r,  (4,n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, one_e_grad_2_dm_alpha_at_r, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, one_e_grad_2_dm_beta_at_r, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, scal_prod_grad_one_e_dm_ab, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, one_e_stuff_for_pbe, (3,n_points_final_grid,N_states) ]
 BEGIN_DOC
 ! one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
 !
 ! one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
 !
 ! one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
 !
 ! one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
 !
 ! one_e_grad_2_dm_alpha_at_r(i,istate)      = (d\dx n_alpha(r_i,istate))^2 + (d\dy n_alpha(r_i,istate))^2 + (d\dz n_alpha(r_i,istate))^2
 !
 ! scal_prod_grad_one_e_dm_ab(i,istate)      = grad n_alpha(r_i) . grad n_beta(r_i)
 !
 ! where r_i is the ith point of the grid and istate is the state number
 !
 ! !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed
 END_DOC
 implicit none
 integer :: i,j,k,l,m,istate
 double precision :: contrib
 double precision :: r(3)
 double precision, allocatable :: aos_array(:),grad_aos_array(:,:)
 double precision, allocatable :: dm_a(:),dm_b(:), dm_a_grad(:,:), dm_b_grad(:,:)
 allocate(dm_a(N_states),dm_b(N_states), dm_a_grad(3,N_states), dm_b_grad(3,N_states))
 allocate(aos_array(ao_num),grad_aos_array(3,ao_num))
 do istate = 1, N_states
  do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)

   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a,dm_b,  dm_a_grad, dm_b_grad, aos_array, grad_aos_array)

   ! alpha/beta density 
   one_e_dm_and_grad_alpha_in_r(4,i,istate)  =  dm_a(istate)
   one_e_dm_and_grad_beta_in_r(4,i,istate)   =  dm_b(istate)

   ! alpha/beta density gradients 
   one_e_dm_and_grad_alpha_in_r(1,i,istate)  =  dm_a_grad(1,istate)
   one_e_dm_and_grad_alpha_in_r(2,i,istate)  =  dm_a_grad(2,istate)
   one_e_dm_and_grad_alpha_in_r(3,i,istate)  =  dm_a_grad(3,istate)

   one_e_dm_and_grad_beta_in_r(1,i,istate)  =  dm_b_grad(1,istate)
   one_e_dm_and_grad_beta_in_r(2,i,istate)  =  dm_b_grad(2,istate)
   one_e_dm_and_grad_beta_in_r(3,i,istate)  =  dm_b_grad(3,istate)

   ! alpha/beta squared of the gradients 
   one_e_grad_2_dm_alpha_at_r(i,istate) = dm_a_grad(1,istate) * dm_a_grad(1,istate)  & 
                                        + dm_a_grad(2,istate) * dm_a_grad(2,istate)  & 
                                        + dm_a_grad(3,istate) * dm_a_grad(3,istate)
   one_e_grad_2_dm_beta_at_r(i,istate)  = dm_b_grad(1,istate) * dm_b_grad(1,istate)  & 
                                        + dm_b_grad(2,istate) * dm_b_grad(2,istate)  & 
                                        + dm_b_grad(3,istate) * dm_b_grad(3,istate)

   ! scalar product between alpha and beta density gradient 
   scal_prod_grad_one_e_dm_ab(i,istate) = dm_a_grad(1,istate) * dm_b_grad(1,istate)  & 
                                        + dm_a_grad(2,istate) * dm_b_grad(2,istate)  & 
                                        + dm_a_grad(3,istate) * dm_b_grad(3,istate)

   ! some stuffs needed for GGA type potentials 
   one_e_stuff_for_pbe(1,i,istate) = 2.D0 * (dm_a_grad(1,istate) + dm_b_grad(1,istate) ) & 
                                                 * (dm_a(istate) + dm_b(istate))
   one_e_stuff_for_pbe(2,i,istate) = 2.D0 * (dm_a_grad(2,istate) + dm_b_grad(2,istate) ) & 
                                                 * (dm_a(istate) + dm_b(istate))
   one_e_stuff_for_pbe(3,i,istate) = 2.D0 * (dm_a_grad(3,istate) + dm_b_grad(3,istate) ) & 
                                                 * (dm_a(istate) + dm_b(istate))
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, elec_beta_num_grid_becke , (N_states) ]
&BEGIN_PROVIDER [double precision, elec_alpha_num_grid_becke , (N_states) ]
&BEGIN_PROVIDER [double precision, elec_num_grid_becke , (N_states) ]
 implicit none
 BEGIN_DOC
 ! number of electrons when the one-e alpha/beta densities are numerically integrated on the DFT grid
 !
 ! !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed
 END_DOC
 integer :: i,istate
 double precision :: r(3),weight
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)
 
   elec_alpha_num_grid_becke(istate) += one_e_dm_and_grad_alpha_in_r(4,i,istate) * weight
   elec_beta_num_grid_becke(istate)  += one_e_dm_and_grad_beta_in_r(4,i,istate) * weight
  enddo
  elec_num_grid_becke(istate) = elec_alpha_num_grid_becke(istate) + elec_beta_num_grid_becke(istate)
 enddo

END_PROVIDER


