program save_one_e_dm_cplx
  implicit none
 BEGIN_DOC
! Program that computes the one body density on the |MO| and |AO| basis
! for $\alpha$ and $\beta$ electrons from the wave function
! stored in the |EZFIO| directory, and then saves it into the
! :ref:`module_aux_quantities`.
!
! Then, the global variable :option:`aux_quantities data_one_e_dm_alpha_mo_complex`
! and :option:`aux_quantities data_one_e_dm_beta_mo_complex` (and the corresponding for |AO|)
 END_DOC
  read_wf = .True.
  touch read_wf
  call routine_save_one_e_dm_cplx

end

subroutine routine_save_one_e_dm_cplx
 implicit none
 BEGIN_DOC
 ! routine called by :c:func:`save_one_e_dm_cplx`
 END_DOC
 call ezfio_set_aux_quantities_data_one_e_dm_alpha_mo_complex(one_e_dm_mo_alpha_complex)
 call ezfio_set_aux_quantities_data_one_e_dm_beta_mo_complex(one_e_dm_mo_beta_complex)
 call ezfio_set_aux_quantities_data_one_e_dm_alpha_ao_complex(one_e_dm_ao_alpha_complex_nst)
 call ezfio_set_aux_quantities_data_one_e_dm_beta_ao_complex(one_e_dm_ao_beta_complex_nst)
end
