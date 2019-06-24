program casscf_new
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  no_vvvv_integrals = .True.
  SOFT_TOUCH no_vvvv_integrals
  call run
end

subroutine run
  implicit none
  call run_cipsi
end
