BEGIN_PROVIDER [double precision, shifting_constant, (N_states)]
 implicit none
 BEGIN_DOC
 ! shifting_constant = (E_{Hxc} - <\Psi | V_{Hxc} | \Psi>) / N_elec
 ! constant to add to the potential in order to obtain the variational energy as
 ! the eigenvalue of the effective long-range Hamiltonian
 ! (see original paper of Levy PRL 113, 113002 (2014), equation (17) )
 END_DOC
 integer :: istate
 do istate = 1, N_states
  shifting_constant(istate) = energy_x(istate) + energy_c(istate) + short_range_Hartree(istate) - Trace_v_Hxc(istate)
 enddo
 shifting_constant = shifting_constant / dble(elec_num)


END_PROVIDER
