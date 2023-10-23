program export_trexio_prog
  implicit none
  read_wf = .True.
  SOFT_TOUCH read_wf
  call export_trexio(.False.)
end

