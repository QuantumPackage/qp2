subroutine htc_bi_ortho_calc_tdav_slow(v, u, N_st, sze)

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
  call htilde_mu_mat_bi_ortho_tot_slow(psi_det(1,1,i), psi_det(1,1,j), N_int, htot)

  v = 0.d0
 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
 !$OMP SHARED(N_st, sze, N_int, psi_det, u, v)       &
 !$OMP PRIVATE(istate, i, j, htot)
  do istate = 1, N_st
    do i = 1, sze
      do j = 1, sze
        call htilde_mu_mat_bi_ortho_tot_slow(psi_det(1,1,i), psi_det(1,1,j), N_int, htot)
        v(i,istate) = v(i,istate) + htot * u(j,istate)
      enddo
    enddo 
  enddo
 !$OMP END PARALLEL DO

end 

subroutine htcdag_bi_ortho_calc_tdav_slow(v, u, N_st, sze)

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
  call htilde_mu_mat_bi_ortho_tot_slow(psi_det(1,1,i), psi_det(1,1,j), N_int, htot)

  v = 0.d0

 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
 !$OMP SHARED(N_st, sze, N_int, psi_det, u, v)       &
 !$OMP PRIVATE(istate, i, j, htot)
  do istate = 1, N_st
    do i = 1, sze
      do j = 1, sze
        call htilde_mu_mat_bi_ortho_tot_slow(psi_det(1,1,j), psi_det(1,1,i), N_int, htot)
        v(i,istate) = v(i,istate) + htot * u(j,istate)
      enddo
    enddo
  enddo
 !$OMP END PARALLEL DO

end 

subroutine i_H_tc_psi_phi(key,keys,coef_l,coef_r,Nint,Ndet,Ndet_max,Nstate,chi_H_i_array,i_H_phi_array)
  use bitmasks
  implicit none
  BEGIN_DOC
! Computes $\langle i|H|Phi \rangle  = \sum_J c^R_J \langle i | H | J \rangle$.
!
! AND      $\langle Chi|H|  i \rangle  = \sum_J c^L_J \langle J | H | i \rangle$.
!
! CONVENTION: i_H_phi_array(0) =   total matrix element,  
!
!             i_H_phi_array(1) =   one-electron matrix element,
!
!             i_H_phi_array(2) =   two-electron matrix element,
!
!             i_H_phi_array(3) = three-electron matrix element,
!
! Uses filter_connected_i_H_psi0 to get all the $|J \rangle$ to which $|i \rangle$
! is connected.
!
! The i_H_psi_minilist is much faster but requires to build the
! minilists.
  END_DOC
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef_l(Ndet_max,Nstate),coef_r(Ndet_max,Nstate)
  double precision, intent(out)  :: chi_H_i_array(0:3,Nstate),i_H_phi_array(0:3,Nstate)

  integer                        :: i, ii,j
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: hmono, htwoe, hthree, htot
  integer, allocatable           :: idx(:)

  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  allocate(idx(0:Ndet))

  chi_H_i_array = 0.d0
  i_H_phi_array = 0.d0

  call filter_connected_i_H_psi0(keys,key,Nint,Ndet,idx)
  if (Nstate == 1) then

    do ii=1,idx(0)
      i = idx(ii)
      ! computes <Chi|H_tc|i>
      !DIR$ FORCEINLINE
      call htilde_mu_mat_opt_bi_ortho(keys(1,1,i), key, Nint, hmono, htwoe, hthree, htot)
      chi_H_i_array(0,1) = chi_H_i_array(0,1) + coef_l(i,1)*htot
      chi_H_i_array(1,1) = chi_H_i_array(1,1) + coef_l(i,1)*hmono
      chi_H_i_array(2,1) = chi_H_i_array(2,1) + coef_l(i,1)*htwoe
      chi_H_i_array(3,1) = chi_H_i_array(3,1) + coef_l(i,1)*hthree
      ! computes <i|H_tc|Phi>
      !DIR$ FORCEINLINE
      call htilde_mu_mat_opt_bi_ortho(key,keys(1,1,i), Nint, hmono, htwoe, hthree, htot)
      i_H_phi_array(0,1) = i_H_phi_array(0,1) + coef_r(i,1)*htot
      i_H_phi_array(1,1) = i_H_phi_array(1,1) + coef_r(i,1)*hmono
      i_H_phi_array(2,1) = i_H_phi_array(2,1) + coef_r(i,1)*htwoe
      i_H_phi_array(3,1) = i_H_phi_array(3,1) + coef_r(i,1)*hthree
    enddo

  else

    do ii=1,idx(0)
      i = idx(ii)
      ! computes <Chi|H_tc|i>
      !DIR$ FORCEINLINE
      call htilde_mu_mat_opt_bi_ortho(keys(1,1,i), key, Nint, hmono, htwoe, hthree, htot)
      do j = 1, Nstate
        chi_H_i_array(0,j) = chi_H_i_array(0,j) + coef_l(i,j)*htot
        chi_H_i_array(1,j) = chi_H_i_array(1,j) + coef_l(i,j)*hmono
        chi_H_i_array(2,j) = chi_H_i_array(2,j) + coef_l(i,j)*htwoe
        chi_H_i_array(3,j) = chi_H_i_array(3,j) + coef_l(i,j)*hthree
      enddo
      ! computes <i|H_tc|Phi>
      !DIR$ FORCEINLINE
      call htilde_mu_mat_opt_bi_ortho(key,keys(1,1,i), Nint, hmono, htwoe, hthree, htot)
      do j = 1, Nstate
        i_H_phi_array(0,j) = i_H_phi_array(0,j) + coef_r(i,j)*htot
        i_H_phi_array(1,j) = i_H_phi_array(1,j) + coef_r(i,j)*hmono
        i_H_phi_array(2,j) = i_H_phi_array(2,j) + coef_r(i,j)*htwoe
        i_H_phi_array(3,j) = i_H_phi_array(3,j) + coef_r(i,j)*hthree
      enddo
    enddo

  endif

end

