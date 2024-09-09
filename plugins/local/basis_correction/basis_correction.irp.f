program basis_correction
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf
  no_core_density = .True.
  touch no_core_density
  call print_basis_correction
end

