
 BEGIN_PROVIDER [double precision, htilde_matrix_elmt_bi_ortho, (N_det,N_det)]

  BEGIN_DOC
  ! htilde_matrix_elmt_bi_ortho(j,i) = <J| H^tilde |I> 
  !
  ! WARNING !!!!!!!!! IT IS NOT HERMITIAN !!!!!!!!!
  END_DOC
 
  implicit none
  integer          :: i, j
  double precision :: hmono,htwoe,hthree,htot

  PROVIDE N_int
 !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) PRIVATE(i,j,hmono, htwoe, hthree, htot) &
 !$OMP SHARED (N_det, psi_det, N_int,htilde_matrix_elmt_bi_ortho)
    do i = 1, N_det
      do j = 1, N_det
        ! < J | Htilde | I >
        call htilde_mu_mat_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)

        !print *, ' hmono  = ', hmono
        !print *, ' htwoe  = ', htwoe
        !print *, ' hthree = ', hthree
        htilde_matrix_elmt_bi_ortho(j,i) = htot
      enddo
    enddo
 !$OMP END PARALLEL DO
! print*,'htilde_matrix_elmt_bi_ortho = '
! do i = 1, min(100,N_det)
!  write(*,'(100(F16.10,X))')htilde_matrix_elmt_bi_ortho(1:min(100,N_det),i)
! enddo


END_PROVIDER 

 BEGIN_PROVIDER [double precision, htilde_matrix_elmt_bi_ortho_tranp, (N_det,N_det)]
 implicit none
 integer ::i,j
  do i = 1, N_det
    do j = 1, N_det
      htilde_matrix_elmt_bi_ortho_tranp(j,i) = htilde_matrix_elmt_bi_ortho(i,j)
    enddo
  enddo
END_PROVIDER 
