program basis_correction
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf
  no_core_density = .True.
  touch no_core_density
  provide mo_two_e_integrals_in_map
  call print_basis_correction
!  call print_e_b
end

subroutine print_e_b
 implicit none
  print *, 'Hello world'
  print*,'ecmd_lda_mu_of_r             = ',ecmd_lda_mu_of_r
  print*,'ecmd_pbe_ueg_mu_of_r         = ',ecmd_pbe_ueg_mu_of_r
  print*,'ecmd_pbe_ueg_eff_xi_mu_of_r  = ',ecmd_pbe_ueg_eff_xi_mu_of_r
  print*,''
  print*,'psi_energy + E^B_LDA         = ',psi_energy + ecmd_lda_mu_of_r
  print*,'psi_energy + E^B_PBE_UEG     = ',psi_energy + ecmd_pbe_ueg_mu_of_r
  print*,'psi_energy + E^B_PBE_UEG_Xi  = ',psi_energy + ecmd_pbe_ueg_eff_xi_mu_of_r
  print*,''
  print*,'mu_average_prov              = ',mu_average_prov
end
