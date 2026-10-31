 BEGIN_PROVIDER [double precision, ao_potential_alpha_xc, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_potential_beta_xc, (ao_num, ao_num)]
 implicit none
  BEGIN_DOC
  ! AO-basis (ao_num x ao_num) exchange/correlation potential for
  ! alpha and beta electrons (both provided in this block); added to
  ! the KS Fock matrix in fock_matrix_ks.irp.f. If same_xc_func,
  ! copies potential_xc_alpha_ao / potential_xc_beta_ao at state
  ! index 1; otherwise sums potential_x_*_ao and potential_c_*_ao
  ! at state index 1.
  END_DOC
 integer :: i,j,k,l
 ao_potential_alpha_xc = 0.d0
 ao_potential_beta_xc = 0.d0
  if(same_xc_func)then
   do i = 1, ao_num
    do j = 1, ao_num
     ao_potential_alpha_xc(i,j) =  potential_xc_alpha_ao(i,j,1) 
     ao_potential_beta_xc(i,j)  =  potential_xc_beta_ao(i,j,1)  
    enddo
   enddo
  else
   do i = 1, ao_num
    do j = 1, ao_num
     ao_potential_alpha_xc(i,j) =  potential_c_alpha_ao(i,j,1) + potential_x_alpha_ao(i,j,1)
     ao_potential_beta_xc(i,j)  =  potential_c_beta_ao(i,j,1)  + potential_x_beta_ao(i,j,1)
    enddo
   enddo
  endif
END_PROVIDER

BEGIN_PROVIDER [double precision, e_exchange_dft]
 implicit none
  BEGIN_DOC
  ! DFT exchange energy: energy_x(1), i.e. the exchange energy at
  ! state index 1 for the functional named by exchange_functional,
  ! scaled by (1 - HF_exchange) (general provider energy_x in
  ! dft_one_e/e_xc_general.irp.f). Added to KS_energy in
  ! ks_enery.irp.f.
  END_DOC
  e_exchange_dft = energy_x(1)

END_PROVIDER

BEGIN_PROVIDER [double precision, e_correlation_dft]
 implicit none
  BEGIN_DOC
  ! DFT correlation energy: energy_c(1), i.e. the correlation
  ! energy at state index 1 for the functional named by
  ! correlation_functional (general provider energy_c in
  ! dft_one_e/e_xc_general.irp.f). Added to KS_energy in
  ! ks_enery.irp.f.
  END_DOC
  e_correlation_dft = energy_c(1)

END_PROVIDER
