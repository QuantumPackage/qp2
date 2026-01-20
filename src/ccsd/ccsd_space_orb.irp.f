program ccsd

  implicit none

  BEGIN_DOC
  ! Closed-shell CCSD
  END_DOC
  read_wf = .True.
  touch read_wf

  call run_ccsd_space_orb
  
end
