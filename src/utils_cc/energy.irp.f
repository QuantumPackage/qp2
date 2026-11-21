subroutine det_energy(det,energy)

  implicit none
  BEGIN_DOC
  ! Energy of the input determinant: the diagonal Hamiltonian matrix
  ! element $\langle i|H|i \rangle$ plus the nuclear repulsion.
  END_DOC

  integer(bit_kind), intent(in) :: det

  double precision, intent(out) :: energy
  double precision, external :: diag_H_mat_elem

  energy = diag_H_mat_elem(det,N_int) + nuclear_repulsion

end
