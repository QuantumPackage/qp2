! Orbital optimization program

! This is an optimization program for molecular orbitals. It produces
! orbital rotations in order to lower the energy of a truncated wave
! function.  
! This program just optimize the orbitals for a fixed number of
! determinants. This optimization process must be repeated for different
! number of determinants.




! Main program : orb_opt_trust


program orb_opt
  BEGIN_DOC
  ! Orbital optimization for a fixed number of determinants: applies
  ! orbital rotations to the active orbitals (trust region algorithm
  ! in run_orb_opt_trust_v2) to lower the energy of the CI wave
  ! function, while the number of determinants is kept fixed. Sets
  ! read_wf so that the wave function is read from the EZFIO file.
  END_DOC
  read_wf = .true. ! must be True for the orbital optimization !!!
  TOUCH read_wf
  io_mo_two_e_integrals = 'None'
  TOUCH io_mo_two_e_integrals
  call run_orb_opt_trust_v2
end
