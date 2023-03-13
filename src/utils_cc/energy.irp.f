subroutine det_energy(det,energy)

  implicit none

  integer(bit_kind), intent(in) :: det

  double precision, intent(out) :: energy

  call i_H_j(det,det,N_int,energy)

  energy = energy + nuclear_repulsion
  
end
