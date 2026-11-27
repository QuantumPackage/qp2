program compute_deltamu_right

  implicit none

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  PROVIDE N_int 
  call delta_right()

end

! ---

subroutine delta_right()

  implicit none
  integer                       :: k
  integer                       :: N_det4
  integer, allocatable          :: order4(:)
  double precision, allocatable :: delta(:,:) 

  allocate( delta(N_det,N_states) )
  delta = 0.d0

  ! dset_order is a 32-bit template routine : down-cast the order
  N_det4 = int(N_det, 4)
  allocate(order4(N_det4))
  order4(1:N_det4) = psi_bilinear_matrix_order(1:N_det4)

  do k = 1, N_states
  !do k = 1, 1

    ! get < I_left | H_mu - H | psi_right > 
    !call get_h_bitc_right(psi_det, psi_r_coef_bi_ortho(:,k), N_det, N_int, delta(:,k))
    call get_delta_bitc_right(psi_det, psi_r_coef_bi_ortho(:,k), N_det, N_int, delta(:,k))

    ! order as QMCCHEM
    call dset_order(delta(:,k), order4, N_det4)

  enddo

!  call ezfio_set_dmc_dress_dmc_delta_h(delta)

  deallocate(delta)

  return
end

! ---

