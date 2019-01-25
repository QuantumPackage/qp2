use bitmasks


BEGIN_PROVIDER [ integer , N_det_delta_ij ]
  implicit none
  N_det_delta_ij = N_det
END_PROVIDER

BEGIN_PROVIDER [ double precision, delta_ij, (N_states, N_det, 2) ]
  implicit none
  if(.true.) then
    delta_ij(:,:N_det_delta_ij, :) = delta_ij_tmp(:,:,:)
  endif
  delta_ij(:,N_det_delta_ij+1:,:) = 0d0
END_PROVIDER

BEGIN_PROVIDER [ double precision, delta_ij_tmp, (N_states,N_det_delta_ij,2) ]
  use bitmasks
  implicit none

  integer                        :: i,j,k

  double precision, allocatable  :: dress(:), del(:,:), del_s2(:,:)
  double precision               :: E_CI_before(N_states)
  integer :: cnt = 0

  allocate(dress(N_states), del(N_states, N_det), del_s2(N_states, N_det))

  delta_ij_tmp = 0d0

  E_CI_before(:) = dress_E0_denominator(:) + nuclear_repulsion

  call write_double(6,dress_relative_error,"Convergence of the stochastic algorithm")

  call ZMQ_dress(E_CI_before, dress, del, del_s2, abs(dress_relative_error))
  delta_ij_tmp(:,1:N_det_delta_ij,1) = del(:,1:N_det_delta_ij)
  delta_ij_tmp(:,1:N_det_delta_ij,2) = del_s2(:,1:N_det_delta_ij)


  deallocate(dress, del, del_s2)
END_PROVIDER



