program ccsd

  implicit none

  BEGIN_DOC
  ! CCSD program
  END_DOC

  read_wf = .True.
  touch read_wf

  if (.not. cc_ref_is_open_shell) then
    call run_ccsd_space_orb
  else 
    call run_ccsd_spin_orb
  endif

end
