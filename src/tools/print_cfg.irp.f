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
 integer :: i,j

 do i=1,N_configuration
   print *, i, sum(popcnt(psi_configuration(:,1,i)))
 enddo

 print *, ''
 do i=0,elec_num
   print *, i, cfg_seniority_index(i)
 enddo
end
