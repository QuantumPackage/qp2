
subroutine write_tc_energy()

  implicit none
  integer                       :: i, j, k
  double precision              :: hmono, htwoe, hthree, htot
  double precision              :: E_TC, O_TC
  double precision              :: E_1e, E_2e, E_3e
  double precision, allocatable :: E_TC_tmp(:), E_1e_tmp(:), E_2e_tmp(:), E_3e_tmp(:)

  ! GS 
  ! ---

  allocate(E_TC_tmp(N_det), E_1e_tmp(N_det), E_2e_tmp(N_det), E_3e_tmp(N_det))

  !$OMP PARALLEL                                                                &
  !$OMP DEFAULT (NONE)                                                          &
  !$OMP PRIVATE(i, j, hmono, htwoe, hthree, htot)                               &
  !$OMP SHARED(N_det, psi_det, N_int, psi_l_coef_bi_ortho, psi_r_coef_bi_ortho, &
  !$OMP        E_TC_tmp, E_1e_tmp, E_2e_tmp, E_3e_tmp)
  !$OMP DO
  do i = 1, N_det
    E_TC_tmp(i) = 0.d0
    E_1e_tmp(i) = 0.d0
    E_2e_tmp(i) = 0.d0
    E_3e_tmp(i) = 0.d0
    do j = 1, N_det
      call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,i), psi_det(1,1,j), N_int, hmono, htwoe, hthree, htot)
      E_TC_tmp(i) = E_TC_tmp(i) + psi_l_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(j,1) * htot
      E_1e_tmp(i) = E_1e_tmp(i) + psi_l_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(j,1) * hmono
      E_2e_tmp(i) = E_2e_tmp(i) + psi_l_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(j,1) * htwoe
      E_3e_tmp(i) = E_3e_tmp(i) + psi_l_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(j,1) * hthree
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  E_1e = 0.d0
  E_2e = 0.d0
  E_3e = 0.d0
  E_TC = 0.d0
  O_TC = 0.d0
  do i = 1, N_det
    E_1e = E_1e + E_1e_tmp(i)
    E_2e = E_2e + E_2e_tmp(i)
    E_3e = E_3e + E_3e_tmp(i)
    E_TC = E_TC + E_TC_tmp(i)
    O_TC = O_TC + psi_l_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(i,1)
  enddo

  print *, ' state :', 1
  print *, " E_TC = ", E_TC / O_TC
  print *, " E_1e = ", E_1e / O_TC
  print *, " E_2e = ", E_2e / O_TC
  print *, " E_3e = ", E_3e / O_TC
  print *, " O_TC = ", O_TC

  call ezfio_set_tc_bi_ortho_tc_gs_energy(E_TC/O_TC)

  ! ---

  ! ES 
  ! ---

  do k = 2, n_states

    E_TC = 0.d0
    E_1e = 0.d0
    E_2e = 0.d0
    E_3e = 0.d0
    do i = 1, N_det
      do j = 1, N_det
        call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,i), psi_det(1,1,j), N_int, hmono, htwoe, hthree, htot)
        E_TC = E_TC + psi_l_coef_bi_ortho(i,k) * psi_r_coef_bi_ortho(j,k) * htot
        E_1e = E_1e + psi_l_coef_bi_ortho(i,k) * psi_r_coef_bi_ortho(j,k) * hmono
        E_2e = E_2e + psi_l_coef_bi_ortho(i,k) * psi_r_coef_bi_ortho(j,k) * htwoe
        E_3e = E_3e + psi_l_coef_bi_ortho(i,k) * psi_r_coef_bi_ortho(j,k) * hthree
      enddo
    enddo

    O_TC = 0.d0
    do i = 1, N_det
      O_TC = O_TC + psi_l_coef_bi_ortho(i,k) * psi_r_coef_bi_ortho(i,k)
    enddo

    print *, ' state :', k
    print *, " E_TC = ", E_TC / O_TC
    print *, " E_1e = ", E_1e / O_TC
    print *, " E_2e = ", E_2e / O_TC
    print *, " E_3e = ", E_3e / O_TC
    print *, " O_TC = ", O_TC

  enddo

  deallocate(E_TC_tmp, E_1e_tmp, E_2e_tmp, E_3e_tmp)

end

! ---

subroutine write_tc_var()

  implicit none
  integer          :: i, j, k
  double precision :: hmono, htwoe, hthree, htot_1j, htot_j1
  double precision :: SIGMA_TC

  do k = 1, n_states

    SIGMA_TC = 0.d0
    do j = 2, N_det
      call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,j), N_int, hmono, htwoe, hthree, htot_1j)
      call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,1), N_int, hmono, htwoe, hthree, htot_j1)
      SIGMA_TC = SIGMA_TC + htot_1j * htot_j1
    enddo

    print *, " state    : ", k
    print *, " SIGMA_TC = ", SIGMA_TC

  enddo

end

! ---

subroutine write_tc_gs_var_HF()

  implicit none
  integer          :: i, j, k
  double precision :: hmono, htwoe, hthree, htot
  double precision :: SIGMA_TC

  SIGMA_TC = 0.d0
  do j = 2, N_det
    call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,1), N_int, hmono, htwoe, hthree, htot)
    SIGMA_TC = SIGMA_TC + htot * htot
  enddo

  print *, " SIGMA_TC = ", SIGMA_TC

  call ezfio_set_tc_bi_ortho_tc_gs_var(SIGMA_TC)

end

! ---


