program print_energy
 implicit none
 BEGIN_DOC
 ! Prints the energy of the wave function stored in the |EZFIO| directory.
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
 integer :: i
 double precision :: i_H_psi_array(N_states)
 double precision :: E(N_states)
 double precision :: norm(N_states)

 E(:) = nuclear_repulsion
 norm(:) = 0.d0
 do i=1,N_det
  call i_H_psi(psi_det(1,1,i), psi_det, psi_coef, N_int, N_det, &
               size(psi_coef,1), N_states, i_H_psi_array)
  norm(:) += psi_coef(i,:)**2
  E(:) += i_H_psi_array(:) * psi_coef(i,:)
 enddo

 print *, 'Energy:'
 do i=1,N_states
   print *, E(i)/norm(i)
 enddo
end

subroutine run_complex
  implicit none
  integer :: i
  complex*16 :: i_h_psi_array(n_states)
  double precision :: e(n_states)
  double precision :: norm(n_states)
 
  e(:) = nuclear_repulsion
  norm(:) = 0.d0
  do i=1,n_det
   call i_H_psi_complex(psi_det(1,1,i), psi_det, psi_coef_complex, N_int, N_det, &
                size(psi_coef_complex,1), N_states, i_H_psi_array)
   norm(:) += cdabs(psi_coef_complex(i,:))**2
   E(:) += dble(i_h_psi_array(:) * dconjg(psi_coef_complex(i,:)))
  enddo
 
  print *, 'Energy:'
  do i=1,N_states
    print *, E(i)/norm(i)
  enddo
end
