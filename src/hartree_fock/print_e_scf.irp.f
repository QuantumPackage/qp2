program print_e_scf
  call run
end

subroutine run

  use bitmasks
  implicit none

  call print_debug_scf_complex

  print*,'hf 1e,2e,total energy'
  print*,hf_one_electron_energy
  print*,hf_two_electron_energy
  print*,hf_energy
  
end


