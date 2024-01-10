program cas_complete
  implicit none
  BEGIN_DOC
! Diagonalizes the Hamiltonian in the complete active space
  END_DOC

  call generate_cas_space
  call diagonalize_ci
  call save_wavefunction

end


