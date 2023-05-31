
! ---

BEGIN_PROVIDER [ double precision, CI_energy_nonsym_dressed, (N_states_diag) ]

  BEGIN_DOC
  ! N_states lowest eigenvalues of the CI matrix
  END_DOC

  implicit none
  integer       :: j
  character*(8) :: st

  call write_time(6)
  do j = 1, min(N_det, N_states_diag)
    CI_energy_nonsym_dressed(j) = CI_electronic_energy_nonsym_dressed(j) + nuclear_repulsion
  enddo

  do j = 1, min(N_det, N_states)
    write(st, '(I4)') j
    call write_double(6, CI_energy_nonsym_dressed(j), 'Energy of state '//trim(st))
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, CI_electronic_energy_nonsym_dressed, (N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_eigenvectors_nonsym_dressed, (N_det,N_states_diag) ]

  BEGIN_DOC
  ! Eigenvectors/values of the CI matrix
  END_DOC

  implicit none
  logical                       :: converged
  integer                       :: i, j, k
  integer                       :: i_other_state
  integer                       :: i_state
  logical,          allocatable :: good_state_array(:)
  integer,          allocatable :: index_good_state_array(:)
  double precision, allocatable :: eigenvectors(:,:), eigenvalues(:)

  PROVIDE threshold_nonsym_davidson nthreads_davidson

  ! Guess values for the "N_states" states of the CI_eigenvectors_nonsym_dressed
  do j = 1, min(N_states, N_det)
    do i = 1, N_det
      CI_eigenvectors_nonsym_dressed(i,j) = psi_coef(i,j)
    enddo
  enddo

  do j = min(N_states, N_det)+1, N_states_diag
    do i = 1, N_det
      CI_eigenvectors_nonsym_dressed(i,j) = 0.d0
    enddo
  enddo

  ! ---

  if(diag_algorithm == "Davidson") then

    ASSERT(n_states_diag .lt. n_states)

    do j = 1, min(N_states, N_det)
      do i = 1, N_det
        CI_eigenvectors_nonsym_dressed(i,j) = psi_coef(i,j)
      enddo
    enddo

    converged = .False.
    call davidson_diag_nonsym_h( psi_det, CI_eigenvectors_nonsym_dressed &
                               , size(CI_eigenvectors_nonsym_dressed, 1) &
                               , CI_electronic_energy_nonsym_dressed     &
                               , N_det, min(N_det, N_states), min(N_det, N_states_diag), N_int, 1, converged )

  else if(diag_algorithm == "Lapack") then

    allocate(eigenvectors(size(H_matrix_nonsym_dressed, 1),N_det))
    allocate(eigenvalues(N_det))

    call diag_nonsym_right( N_det, H_matrix_nonsym_dressed, size(H_matrix_nonsym_dressed, 1)       &
                          , eigenvectors, size(eigenvectors, 1), eigenvalues, size(eigenvalues, 1) )

    CI_electronic_energy_nonsym_dressed(:) = 0.d0

    ! Select the "N_states_diag" states of lowest energy
    do j = 1, min(N_det, N_states_diag)
      do i = 1, N_det
        CI_eigenvectors_nonsym_dressed(i,j) = eigenvectors(i,j)
      enddo
      CI_electronic_energy_nonsym_dressed(j) = eigenvalues(j)
    enddo

    deallocate(eigenvectors, eigenvalues)

    ! --- ---

  endif

  ! ---

END_PROVIDER

! ---

subroutine diagonalize_CI_nonsym_dressed()

  BEGIN_DOC
  !  Replace the coefficients of the CI states by the coefficients of the
  !  eigenstates of the CI matrix
  END_DOC

  implicit none
  integer :: i, j

  PROVIDE dressing_delta

  do j = 1, N_states
    do i = 1, N_det
      psi_coef(i,j) = CI_eigenvectors_nonsym_dressed(i,j)
    enddo
  enddo

  SOFT_TOUCH psi_coef

end subroutine diagonalize_CI_nonsym_dressed

! ---

BEGIN_PROVIDER [ double precision, H_matrix_nonsym_dressed, (N_det,N_det) ]

  BEGIN_DOC
  ! Dressed H with Delta_ij
  END_DOC

  implicit none
  integer          :: i, j, l, k
  double precision :: f

  H_matrix_nonsym_dressed(1:N_det,1:N_det) = h_matrix_all_dets(1:N_det,1:N_det)

  if(N_states == 1) then

!    !symmetric formula
!    l = dressed_column_idx(1)
!    f = 1.0d0/psi_coef(l,1)
!    do i=1,N_det
!      h_matrix_nonsym_dressed(i,l) +=  dressing_column_h(i,1) *f
!      h_matrix_nonsym_dressed(l,i) +=  dressing_column_h(i,1) *f
!    enddo

!    l = dressed_column_idx(1)
!    f = 1.0d0 / psi_coef(l,1)
!    do j = 1, N_det
!      H_matrix_nonsym_dressed(j,l) += f * dressing_delta(j,1) 
!    enddo

    k = 1
    l = 1
    f = overlap_states_inv(k,l)
    do j = 1, N_det
      do i = 1, N_det
        H_matrix_nonsym_dressed(i,j) = H_matrix_nonsym_dressed(i,j) + f * dressing_delta(i,k) * psi_coef(j,l)
      enddo
    enddo

  else

    do k = 1, N_states
      do l = 1, N_states
        f = overlap_states_inv(k,l)

        do j = 1, N_det
          do i = 1, N_det
            H_matrix_nonsym_dressed(i,j) = H_matrix_nonsym_dressed(i,j) + f * dressing_delta(i,k) * psi_coef(j,l)
          enddo
        enddo

      enddo
    enddo

  endif

END_PROVIDER

! ---

