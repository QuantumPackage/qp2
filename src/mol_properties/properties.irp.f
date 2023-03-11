program mol_properties

  implicit none

  BEGIN_DOC
  ! Run the propertie calculations
  END_DOC

  read_wf = .True.
  touch read_wf

  call print_mol_properties()

end
