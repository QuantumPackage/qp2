subroutine det_energy(det,energy)

  implicit none

  integer(bit_kind), intent(in) :: det

  double precision, intent(out) :: energy
  double precision, external :: diag_H_mat_elem

  energy = diag_H_mat_elem(det,N_int) + nuclear_repulsion

end
