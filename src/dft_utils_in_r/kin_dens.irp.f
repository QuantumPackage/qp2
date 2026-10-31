BEGIN_PROVIDER [double precision, kinetic_density_generalized,  (n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! kinetic_density_generalized(p) = generalized kinetic energy density
 ! at the pth grid point: 0.5 times the double sum over MOs (i,j) of
 ! the MO density matrix element one_e_dm_mo_for_dft(j,i,1) (alpha +
 ! beta, state 1) times the scalar product of the MO gradients
 END_DOC
 integer :: i,j,m,i_point
 kinetic_density_generalized = 0.d0
 do i_point = 1, n_points_final_grid
  do i = 1, mo_num
   do j = 1, mo_num 
    do m = 1, 3
     kinetic_density_generalized(i_point) += 0.5d0 * mos_grad_in_r_array_tranp(m,j,i_point) * mos_grad_in_r_array_tranp(m,i,i_point) * one_e_dm_mo_for_dft(j,i,1) 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 
