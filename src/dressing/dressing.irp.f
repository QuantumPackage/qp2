use bitmasks


BEGIN_PROVIDER [ integer , N_det_delta_ij ]
  implicit none
  BEGIN_DOC
  ! Number of determinant columns stored in delta_ij_tmp:
  ! equal to N_det, so delta_ij and delta_ij_tmp currently
  ! have the same number of columns.
  END_DOC
  N_det_delta_ij = N_det
END_PROVIDER

BEGIN_PROVIDER [ double precision, delta_ij, (N_states, N_det, 2) ]
  implicit none
  BEGIN_DOC
  ! Dressing contribution $\Delta_{ij}$, per state and
  ! determinant: the first N_det_delta_ij columns are copied
  ! from delta_ij_tmp, the remaining columns are zeroed.
  ! Component 1 feeds dressing_column_h, component 2
  ! dressing_column_s: the two components are a contribution to
  ! the Hamiltonian (component 1) and a contribution to $S^2$
  ! (component 2) (author, 2026-10-31).
  END_DOC
  if(.true.) then
    delta_ij(:,:N_det_delta_ij, :) = delta_ij_tmp(:,:,:)
  endif
  delta_ij(:,N_det_delta_ij+1:,:) = 0d0
END_PROVIDER

BEGIN_PROVIDER [ double precision, delta_ij_tmp, (N_states,N_det_delta_ij,2) ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Stochastic estimate of the dressing, computed by
  ! ZMQ_dress with reference energies
  ! dress_E0_denominator + nuclear_repulsion and relative error
  ! criterion dress_relative_error. Component 1 is del,
  ! component 2 is del_s2: a contribution to the Hamiltonian
  ! (component 1) and a contribution to $S^2$ (component 2)
  ! (author, 2026-10-31).
  END_DOC

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



