program export_trexio
  implicit none
  BEGIN_DOC
  !     Exports the wave function in TREXIO format
  END_DOC
  logical :: update, full_path
  read_wf = .True.
  SOFT_TOUCH read_wf
  update = .False.
  full_path = .False.
  call export_trexio_subroutine(update, full_path)
end

