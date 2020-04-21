 BEGIN_PROVIDER [double precision, psi_dft_energy_kinetic, (N_states) ]
&BEGIN_PROVIDER [double precision, psi_dft_energy_nuclear_elec, (N_states) ]
&BEGIN_PROVIDER [double precision, psi_dft_energy_h_core, (N_states) ]
 implicit none
 BEGIN_DOC
! kinetic, electron-nuclear and total h_core energy computed with the density matrix one_e_dm_mo_beta_for_dft+one_e_dm_mo_alpha_for_dft
 END_DOC
 integer :: i,j,istate
 double precision :: accu
 psi_dft_energy_kinetic = 0.d0
 psi_dft_energy_nuclear_elec = 0.d0
 do istate = 1, N_states
  do i = 1, mo_num
   do j = 1, mo_num
    psi_dft_energy_kinetic(istate)      += ( one_e_dm_mo_alpha_for_dft(j,i,istate)+one_e_dm_mo_beta_for_dft(j,i,istate)) * mo_kinetic_integrals(j,i)
    psi_dft_energy_nuclear_elec(istate) += ( one_e_dm_mo_alpha_for_dft(j,i,istate)+one_e_dm_mo_beta_for_dft(j,i,istate)) * mo_integrals_n_e(j,i)
   enddo
  enddo
 enddo
 do i = 1, N_states
  accu = 0.d0
  do j = 1, mo_num
   accu += one_e_dm_mo_alpha_for_dft(j,j,i) + one_e_dm_mo_beta_for_dft(j,j,i)
  enddo
  accu = (elec_alpha_num + elec_beta_num ) / accu
  psi_dft_energy_kinetic(i)      = psi_dft_energy_kinetic(i) * accu
  psi_dft_energy_nuclear_elec(i) = psi_dft_energy_nuclear_elec(i) * accu
  psi_dft_energy_h_core(i)       = psi_dft_energy_nuclear_elec(i)  +  psi_dft_energy_kinetic(i)
 enddo

END_PROVIDER
