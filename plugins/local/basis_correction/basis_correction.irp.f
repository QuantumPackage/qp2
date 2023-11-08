program basis_correction
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf
  no_core_density = .True.
  touch no_core_density
  if(io_mo_two_e_integrals .ne. "Read")then
   provide ao_two_e_integrals_in_map
  endif
  provide mo_two_e_integrals_in_map
  call print_basis_correction
end

