subroutine hcore_guess
  use sizet
  BEGIN_DOC
! Produce `H_core` MO orbital
  END_DOC
  implicit none
  character*(64)                 :: label
  label = 'Guess'
    call mo_as_eigvectors_of_mo_matrix(mo_one_e_integrals, &
                                     int(size(mo_one_e_integrals,1),size_t), &
                                     int(size(mo_one_e_integrals,2),size_t), &
                                     label, 1, .true.)
  call nullify_small_elements(int(ao_num, size_t),int(mo_num, size_t), &
      mo_coef, int(size(mo_coef,1), size_t), 1.d-12 )
  call save_mos
  TOUCH mo_coef mo_label
end
