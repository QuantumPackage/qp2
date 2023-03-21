subroutine htc_bi_ortho_calc_tdav(v, u, N_st, sze)

  use bitmasks

  BEGIN_DOC
    ! Application of H_TC on a vector 
    !
    ! v(i,istate) = \sum_j u(j,istate) H_TC(i,j), with: 
    !   H_TC(i,j) = < Di | H_TC | Dj > 
    !
  END_DOC

  implicit none

  integer, intent(in)             :: N_st, sze
  double precision, intent(in)    :: u(sze,N_st)
  double precision, intent(inout) :: v(sze,N_st)

  integer                         :: i, j, istate
  double precision                :: htot

  PROVIDE N_int 
  PROVIDE psi_det


 ! TODO : transform it with the bi-linear representation in terms of alpha-beta. 

  i = 1
  j = 1
  call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,i), psi_det(1,1,j), N_int, htot)

  v = 0.d0
 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
 !$OMP SHARED(N_st, sze, N_int, psi_det, u, v)       &
 !$OMP PRIVATE(istate, i, j, htot)
  do istate = 1, N_st
    do i = 1, sze
      do j = 1, sze
        call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,i), psi_det(1,1,j), N_int, htot)
        v(i,istate) = v(i,istate) + htot * u(j,istate)
      enddo
    enddo 
  enddo
 !$OMP END PARALLEL DO

end 

subroutine htcdag_bi_ortho_calc_tdav(v, u, N_st, sze)

  use bitmasks

  BEGIN_DOC
    ! Application of (H_TC)^dagger on a vector 
    !
    ! v(i,istate) = \sum_j u(j,istate) H_TC(j,i), with: 
    !   H_TC(i,j) = < Di | H_TC | Dj > 
    !
  END_DOC

  implicit none

  integer, intent(in)             :: N_st, sze
  double precision, intent(in)    :: u(sze,N_st)
  double precision, intent(inout) :: v(sze,N_st)

  integer                         :: i, j, istate
  double precision                :: htot

  PROVIDE N_int
  PROVIDE psi_det

  i = 1
  j = 1
  call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,i), psi_det(1,1,j), N_int, htot)

  v = 0.d0

 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
 !$OMP SHARED(N_st, sze, N_int, psi_det, u, v)       &
 !$OMP PRIVATE(istate, i, j, htot)
  do istate = 1, N_st
    do i = 1, sze
      do j = 1, sze
        call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,j), psi_det(1,1,i), N_int, htot)
        v(i,istate) = v(i,istate) + htot * u(j,istate)
      enddo
    enddo
  enddo
 !$OMP END PARALLEL DO

end 

