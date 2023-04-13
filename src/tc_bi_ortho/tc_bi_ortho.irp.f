program tc_bi_ortho

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

  call routine_diag()
  call write_tc_energy()
  call save_tc_bi_ortho_wavefunction()
end

subroutine test
 implicit none
 integer :: i,j 
 double precision :: hmono,htwoe,hthree,htot
 use bitmasks
 print*,'reading the wave function '
 do i = 1, N_det
  call debug_det(psi_det(1,1,i),N_int)
  print*,i,psi_l_coef_bi_ortho(i,1)*psi_r_coef_bi_ortho(i,1)
  print*,i,psi_l_coef_bi_ortho(i,1),psi_r_coef_bi_ortho(i,1)
 enddo

end

subroutine routine_diag()

  implicit none
  integer          :: i, j, k
  double precision :: dE

  ! provide eigval_right_tc_bi_orth
  ! provide overlap_bi_ortho
  ! provide htilde_matrix_elmt_bi_ortho

  if(N_states .eq. 1) then

    print*,'eigval_right_tc_bi_orth = ',eigval_right_tc_bi_orth(1)
    print*,'e_tc_left_right         = ',e_tc_left_right
    print*,'e_tilde_bi_orth_00      = ',e_tilde_bi_orth_00
    print*,'e_pt2_tc_bi_orth        = ',e_pt2_tc_bi_orth
    print*,'e_pt2_tc_bi_orth_single = ',e_pt2_tc_bi_orth_single
    print*,'e_pt2_tc_bi_orth_double = ',e_pt2_tc_bi_orth_double
    print*,'***'
    print*,'e_corr_bi_orth          = ',e_corr_bi_orth
    print*,'e_corr_bi_orth_proj     = ',e_corr_bi_orth_proj
    print*,'e_corr_single_bi_orth   = ',e_corr_single_bi_orth
    print*,'e_corr_double_bi_orth   = ',e_corr_double_bi_orth
    print*,'Left/right eigenvectors'
    do i = 1,N_det
      write(*,'(I5,X,(100(F12.7,X)))')i,leigvec_tc_bi_orth(i,1),reigvec_tc_bi_orth(i,1),leigvec_tc_bi_orth(i,1)*reigvec_tc_bi_orth(i,1)
    enddo

  else

    print*,'eigval_right_tc_bi_orth : '
    do i = 1, N_states
      print*, i, eigval_right_tc_bi_orth(i)
    enddo

    print*,''
    print*,'******************************************************'
    print*,'TC Excitation energies (au)                     (eV)'
    do i = 2, N_states
      dE = eigval_right_tc_bi_orth(i) - eigval_right_tc_bi_orth(1)
      print*, i, dE, dE/0.0367502d0
    enddo
    print*,''

  endif

end



