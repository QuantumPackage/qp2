program save_natorb
  implicit none
  BEGIN_DOC
! Save natural |MOs| into the |EZFIO|.
!
! This program reads the wave function stored in the |EZFIO| directory,
! extracts the corresponding natural orbitals and setd them as the new
! |MOs|.
!
! If this is a multi-state calculation, the density matrix that produces
! the natural orbitals is obtained from an average of the density
! matrices of each state with the corresponding
! :option:`determinants state_average_weight`
  END_DOC
  read_wf = .True.
  touch read_wf
  call save_natural_mos
  call save_ref_determinant
  call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('None')
  call ezfio_set_mo_two_e_ints_io_df_mo_integrals('None')
  call ezfio_set_mo_one_e_ints_io_mo_one_e_integrals('None')
  call ezfio_set_mo_one_e_ints_io_mo_integrals_kinetic('None')
  call ezfio_set_mo_one_e_ints_io_mo_integrals_e_n('None')
  call ezfio_set_mo_one_e_ints_io_mo_integrals_pseudo('None')
end

