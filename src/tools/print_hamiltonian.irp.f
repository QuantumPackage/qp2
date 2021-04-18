program print_hamiltonian
 implicit none
 BEGIN_DOC
 ! Prints the Hamiltonian matrix defined in the space of determinants
 ! present in the |EZFIO| directory.
 END_DOC

 ! this has to be done in order to be sure that N_det, psi_det and
 ! psi_coef_sorted are the wave function stored in the |EZFIO| directory.
 read_wf = .True.
 touch read_wf
 call run
end

subroutine run
 implicit none
 integer :: i, j
 double precision :: hij

 do j=1,N_det
   do i=1,N_det
     call i_H_j(psi_det(1,1,i), psi_det(1,1,j), N_int, hij)
     if (dabs(hij) > 1.d-20) then
       print *, i, j, hij
     endif
   enddo
 enddo

end
