program print_energy
 implicit none
 BEGIN_DOC
 ! Prints the energy of the wave function stored in the |EZFIO| directory.
 END_DOC

 ! this has to be done in order to be sure that N_det, psi_det and
 ! psi_coef_sorted are the wave function stored in the |EZFIO| directory.
 read_wf = .True.
 touch read_wf
 PROVIDE N_states
 call run
end

subroutine run
 implicit none
 BEGIN_DOC
! E = \sum_{ij} c_i c_j H_{ij}  / \sum{i} c_i^2
 END_DOC
 integer :: i,j
 double precision :: i_H_psi_array(N_states)
 double precision :: E(N_states)
 double precision :: norm(N_states)

 E(1:N_states) = nuclear_repulsion
 norm(1:N_states) = 0.d0
 do i=1,N_det
  call i_H_psi(psi_det(1,1,i), psi_det, psi_coef, N_int, N_det, &
               size(psi_coef,1), N_states, i_H_psi_array)
  do j=1,N_states
    norm(j) += psi_coef(i,j)*psi_coef(i,j)
    E(j) += i_H_psi_array(j) * psi_coef(i,j)
  enddo
 enddo

 print *, 'Energy:'
 do i=1,N_states
   print *, E(i)/norm(i)
 enddo
end
