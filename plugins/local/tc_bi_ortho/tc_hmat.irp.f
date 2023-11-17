
! ---

BEGIN_PROVIDER [double precision, htilde_matrix_elmt_bi_ortho, (N_det,N_det)]

  BEGIN_DOC
  !
  ! htilde_matrix_elmt_bi_ortho(j,i) = <J| H^tilde |I> 
  !
  ! WARNING !!!!!!!!! IT IS NOT HERMITIAN !!!!!!!!!
  !
  END_DOC
 
  implicit none
  integer          :: i, j
  double precision :: htot

  call provide_all_three_ints_bi_ortho

  i = 1
  j = 1
  call htilde_mu_mat_opt_bi_ortho_tot(psi_det(1,1,j), psi_det(1,1,i), N_int, htot)

  !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) PRIVATE(i,j, htot) &
  !$OMP SHARED (N_det, psi_det, N_int, htilde_matrix_elmt_bi_ortho)
  do i = 1, N_det
    do j = 1, N_det
      ! < J |Htilde | I >
      call htilde_mu_mat_opt_bi_ortho_tot(psi_det(1,1,j), psi_det(1,1,i), N_int, htot)

      htilde_matrix_elmt_bi_ortho(j,i) = htot
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, htilde_matrix_elmt_bi_ortho_tranp, (N_det,N_det)]

  implicit none
  integer ::i,j

  do i = 1, N_det
    do j = 1, N_det
      htilde_matrix_elmt_bi_ortho_tranp(j,i) = htilde_matrix_elmt_bi_ortho(i,j)
    enddo
  enddo
END_PROVIDER 

! ---

