program print_tc_energy
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 call write_tc_energy
end

subroutine write_tc_energy()

  implicit none
  integer          :: i, j, k
  double precision :: hmono, htwoe, htot
  double precision :: E_TC, O_TC

  do k = 1, n_states

    E_TC = 0.d0
    do i = 1, N_det
      do j = 1, N_det
        call hmat_bi_ortho(psi_det(1,1,i), psi_det(1,1,j), N_int, hmono, htwoe, htot)
        E_TC = E_TC + psi_l_coef_bi_ortho(i,k) * psi_r_coef_bi_ortho(j,k) * htot
        !E_TC = E_TC + leigvec_tc_bi_orth(i,k) * reigvec_tc_bi_orth(j,k) * htot
      enddo
    enddo

    O_TC = 0.d0
    do i = 1, N_det
      O_TC = O_TC + psi_l_coef_bi_ortho(i,k) * psi_r_coef_bi_ortho(i,k)
      !O_TC = O_TC + leigvec_tc_bi_orth(i,k) * reigvec_tc_bi_orth(i,k)
    enddo

    print *, ' state :', k
    print *, " E_TC = ", E_TC
    print *, " O_TC = ", O_TC

  enddo

end

