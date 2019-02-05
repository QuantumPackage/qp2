program save_one_e_dm
  implicit none
 BEGIN_DOC
! Program that computes the one body density on the |MO| basis
! for $\alpha$ and $\beta$ electrons from the wave function
! stored in the |EZFIO| directory, and then saves it into the
! :ref:`module_aux_quantities`.
!
! Then, the global variable :option:`aux_quantities data_one_e_dm_alpha_mo`
! and :option:`aux_quantities data_one_e_dm_beta_mo` will automatically
! read this density in the next calculation. This can be used to perform
! damping on the density in |RSDFT| calculations (see
! :ref:`module_density_for_dft`).
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
