subroutine only_act_bitmask
 implicit none
  integer :: i,j,k
  do k = 1, N_generators_bitmask
   do j = 1, 6
    do i = 1, N_int
    generators_bitmask(i,1,j,k) = act_bitmask(i,1)
    generators_bitmask(i,2,j,k) = act_bitmask(i,2)
    enddo
   enddo
  enddo
  touch generators_bitmask
end

subroutine run_cisd
  implicit none
  integer :: i

  if(pseudo_sym)then
   call H_apply_cisd_sym
  else
   call H_apply_cisd
  endif
  print *,  'N_det = ', N_det
  print*,'******************************'
  print *,  'Energies  of the states:'
  do i = 1,N_states
    print *,  i, CI_energy(i)
  enddo
  if (N_states > 1) then
    print*,'******************************'
    print*,'Excitation energies '
    do i = 2, N_states
      print*, i ,CI_energy(i) - CI_energy(1)
    enddo
  endif
  psi_coef = ci_eigenvectors
  SOFT_TOUCH psi_coef
  call save_wavefunction
  call ezfio_set_cisd_energy(CI_energy)

end
