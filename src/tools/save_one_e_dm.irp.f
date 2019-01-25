program save_one_e_dm
  implicit none
 BEGIN_DOC
! programs that computes the one body density on the mo basis for alpha and beta electrons

! from the wave function stored in the EZFIO folder, and then save it into the EZFIO folder aux_quantities.
!
! Then, the global variable data_one_e_dm_alpha_mo and data_one_e_dm_beta_mo will automatically read this density in a further calculation.
!
! This can be used to perform damping on the density in RS-DFT calculation (see the density_for_dft module).
 END_DOC
  read_wf = .True.
  touch read_wf
  call routine_save_one_e_dm

end

subroutine routine_save_one_e_dm
 implicit none
 BEGIN_DOC
 ! routine called by :c:func:`save_one_e_dm`
 END_DOC
 call ezfio_set_aux_quantities_data_one_e_dm_alpha_mo(one_e_dm_mo_alpha)
 call ezfio_set_aux_quantities_data_one_e_dm_beta_mo(one_e_dm_mo_beta)
end
