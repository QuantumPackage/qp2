subroutine dress_zmq()
  implicit none
  double precision, allocatable  :: energy(:)
  allocate (energy(N_states))
  threshold_generators = 1d0

  read_wf = .True.
  SOFT_TOUCH read_wf threshold_generators

  if (.True.) then
    integer :: i,j
    do j=1,N_states
      do i=1,N_det
        psi_coef(i,j) = CI_eigenvectors(i,j)
      enddo
    enddo
    SOFT_TOUCH psi_coef
  endif
  call run_dressing(N_states,energy)
  deallocate(energy)
end

