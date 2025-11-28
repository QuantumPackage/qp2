program export_trexio_
  implicit none
  logical :: update, full_path
  read_wf = .True.
  SOFT_TOUCH read_wf
  update = .False.
  full_path = .False.
  call export_trexio(update, full_path)
end

