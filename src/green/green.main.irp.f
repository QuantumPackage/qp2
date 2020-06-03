program green
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
  provide n_green_vec
  print*,'ref_bitmask_energy = ',ref_bitmask_energy
!  call psicoefprinttest
  call print_lanczos_eigvals
  call print_spec
end

subroutine psicoefprinttest
  implicit none
  integer :: i
  TOUCH psi_coef_complex
  print *, 'printing ndet', N_det
end
subroutine print_lanczos_eigvals
  implicit none
  integer :: i, iunit, j
  integer :: getunitandopen
  character(5) :: jstr

  do j=1,n_green_vec
    write(jstr,'(I0.3)') j
    iunit = getunitandopen('lanczos_eigval_alpha_beta.out.'//trim(jstr),'w')
    print *, 'printing lanczos eigenvalues, alpha, beta to "lanczos_eigval_alpha_beta.out.'//trim(jstr)//'"'
    do i=1,n_lanczos_iter
      write(iunit,'(I6,3(E25.15))') i, lanczos_eigvals(i,j), alpha_lanczos(i,j), beta_lanczos(i,j)
    enddo
    close(iunit)
  enddo
end
subroutine print_spec
  implicit none
  integer :: i, iunit, j
  integer :: getunitandopen
  character(5) :: jstr
  do j=1,n_green_vec
    write(jstr,'(I0.3)') j
    iunit = getunitandopen('omega_A.out.'//trim(jstr),'w')
    print *, 'printing frequency, spectral density to "omega_A.out.'//trim(jstr)//'"'
    do i=1,n_omega
      write(iunit,'(2(E25.15))') omega_list(i), spectral_lanczos(i,j)
    enddo
    close(iunit)
  enddo
end
