program save_bitcpsileft_for_qmcchem

  integer          :: iunit
  logical          :: exists
  double precision :: e_ref

  print *, ' '
  print *, ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '  call save_for_qmcchem before  '
  print *, ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, ' '

  call write_lr_spindeterminants()

  e_ref = 0.d0
  iunit = 13
  open(unit=iunit, file=trim(ezfio_filename)//'/simulation/e_ref', action='write')

    call ezfio_has_fci_energy_pt2(exists)
    if(.not.exists) then

      call ezfio_has_fci_energy(exists)
      if(.not.exists) then

        call ezfio_has_cisd_energy(exists)
        if(.not.exists) then

          call ezfio_has_tc_scf_bitc_energy(exists)
          if(exists) then
            call ezfio_get_tc_scf_bitc_energy(e_ref)
          endif

        else
          call ezfio_get_cisd_energy(e_ref)
        endif

      else
        call ezfio_get_fci_energy(e_ref)
      endif

      else
        call ezfio_get_fci_energy_pt2(e_ref)
    endif

    write(iunit,*) e_ref

  close(iunit)

end

! --

subroutine write_lr_spindeterminants()

  use bitmasks

  implicit none

  integer                       :: k, l
  double precision, allocatable :: buffer(:,:)

  PROVIDE psi_bitcleft_bilinear_matrix_values

  allocate(buffer(N_det,N_states))
  do l = 1, N_states
    do k = 1, N_det
      buffer(k,l) = psi_bitcleft_bilinear_matrix_values(k,l)
    enddo
  enddo
  call ezfio_set_spindeterminants_psi_left_coef_matrix_values(buffer)
  deallocate(buffer)

end subroutine write_lr_spindeterminants

! ---

