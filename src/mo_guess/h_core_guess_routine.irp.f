subroutine hcore_guess
  BEGIN_DOC
! Produce `H_core` MO orbital
  END_DOC
  implicit none
  character*(64)                 :: label
  label = "Guess"
  call mo_as_eigvectors_of_mo_matrix(mo_one_e_integrals,          &
                                     size(mo_one_e_integrals,1),  &
                                     size(mo_one_e_integrals,2),label,1,.false.)
  call save_mos
  TOUCH mo_coef mo_label
end
