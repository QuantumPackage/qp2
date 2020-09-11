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
 if (is_complex) then
  call run_complex
 else
  call run
 endif
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

subroutine run_complex
 implicit none
 integer :: i, j
 complex*16 :: hij
 double precision :: s2

 print*,'i,j,Hij'
 do j=1,N_det
   do i=1,N_det
     call i_h_j_complex(psi_det(1,1,i), psi_det(1,1,j), N_int, hij)
     if (cdabs(hij) > 1.d-20) then
       print *, i, j, dble(hij), dimag(hij)
     endif
   enddo
 enddo
 print*,'i,j,S2ij'
 do j=1,N_det
   do i=1,N_det
     call get_s2(psi_det(1,1,i), psi_det(1,1,j), N_int, s2)
     if (dabs(s2) > 1.d-20) then
       print *, i, j, s2
     endif
   enddo
 enddo
!  use bitmasks
  integer                        :: degree

 print*,'i,j,degij'
 do j=1,N_det
   do i=1,N_det
     call get_excitation_degree(psi_det(1,1,i), psi_det(1,1,j), degree, N_int)
     if (degree.le.2) then
       print *, i, j, degree
     endif
   enddo
 enddo

end
