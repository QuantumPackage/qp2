! Prog

program ccsd

  implicit none

  BEGIN_DOC
  ! CCSD in spin orbitals
  END_DOC

  read_wf = .True.
  touch read_wf

  call run_ccsd_spin_orb
  
end
