BEGIN_PROVIDER [double precision, one_e_dm_mo_alpha_for_dft, (mo_num,mo_num, N_states)]
 implicit none
 BEGIN_DOC
! density matrix for alpha electrons in the MO basis used for all DFT calculations based on the density
 END_DOC
 double precision :: delta_alpha(mo_num,mo_num,N_states)
 if(density_for_dft .EQ. "damping_rs_dft")then
  delta_alpha = one_e_dm_mo_alpha - data_one_e_dm_alpha_mo
  one_e_dm_mo_alpha_for_dft = data_one_e_dm_alpha_mo + damping_for_rs_dft * delta_alpha
 else if (density_for_dft .EQ. "input_density")then
  one_e_dm_mo_alpha_for_dft = data_one_e_dm_alpha_mo
 else if (density_for_dft .EQ. "WFT")then
  provide mo_coef
  one_e_dm_mo_alpha_for_dft = one_e_dm_mo_alpha
 else if (density_for_dft .EQ. "KS")then
  provide mo_coef
  one_e_dm_mo_alpha_for_dft = one_body_dm_mo_alpha_one_det
 endif
 
 if(no_core_density .EQ. "no_core_dm")then
  integer :: i,j
  do i = 1, n_core_orb
   do j = 1, mo_num
    one_e_dm_mo_alpha_for_dft(j,i,:) = 0.d0
    one_e_dm_mo_alpha_for_dft(i,j,:) = 0.d0
   enddo
  enddo
 endif

END_PROVIDER

BEGIN_PROVIDER [double precision, one_e_dm_mo_beta_for_dft, (mo_num,mo_num, N_states)]
 implicit none
 BEGIN_DOC
! density matrix for beta  electrons in the MO basis used for all DFT calculations based on the density
 END_DOC
 double precision :: delta_beta(mo_num,mo_num,N_states)
 if(density_for_dft .EQ. "damping_rs_dft")then
  delta_beta = one_e_dm_mo_beta - data_one_e_dm_beta_mo
  one_e_dm_mo_beta_for_dft = data_one_e_dm_beta_mo + damping_for_rs_dft * delta_beta
 else if (density_for_dft .EQ. "input_density")then
  one_e_dm_mo_beta_for_dft = data_one_e_dm_beta_mo
 else if (density_for_dft .EQ. "WFT")then
  provide mo_coef
  one_e_dm_mo_beta_for_dft = one_e_dm_mo_beta
 else if (density_for_dft .EQ. "KS")then
  provide mo_coef
  one_e_dm_mo_beta_for_dft = one_body_dm_mo_beta_one_det
 endif

 if(no_core_density .EQ. "no_core_dm")then
  integer :: i,j
  do i = 1, n_core_orb
   do j = 1, mo_num
    one_e_dm_mo_beta_for_dft(j,i,:) = 0.d0
    one_e_dm_mo_beta_for_dft(i,j,:) = 0.d0
   enddo
  enddo
 endif
END_PROVIDER

BEGIN_PROVIDER [double precision, one_e_dm_mo_for_dft, (mo_num,mo_num, N_states)]
 implicit none
 one_e_dm_mo_for_dft = one_e_dm_mo_beta_for_dft + one_e_dm_mo_alpha_for_dft
END_PROVIDER

BEGIN_PROVIDER [double precision, one_e_dm_average_mo_for_dft, (mo_num,mo_num)]
 implicit none
 integer :: i
 one_e_dm_average_mo_for_dft = 0.d0
 do i = 1, N_states
  one_e_dm_average_mo_for_dft(:,:) +=  one_e_dm_mo_for_dft(:,:,i) * state_average_weight(i)
 enddo
END_PROVIDER

 BEGIN_PROVIDER [ double precision, one_e_dm_alpha_ao_for_dft, (ao_num,ao_num,N_states) ]
&BEGIN_PROVIDER [ double precision, one_e_dm_beta_ao_for_dft, (ao_num,ao_num,N_states) ]
 BEGIN_DOC
! one body density matrix on the AO basis based on one_e_dm_mo_alpha_for_dft
 END_DOC
 implicit none
 integer :: istate
 double precision :: mo_alpha,mo_beta

 one_e_dm_alpha_ao_for_dft = 0.d0
 one_e_dm_beta_ao_for_dft = 0.d0
 do istate = 1, N_states
  call mo_to_ao_no_overlap( one_e_dm_mo_alpha_for_dft(1,1,istate), &
                            size(one_e_dm_mo_alpha_for_dft,1),     &
                            one_e_dm_alpha_ao_for_dft(1,1,istate), &
                            size(one_e_dm_alpha_ao_for_dft,1) )
  call mo_to_ao_no_overlap( one_e_dm_mo_beta_for_dft(1,1,istate), &
                            size(one_e_dm_mo_beta_for_dft,1),     &
                            one_e_dm_beta_ao_for_dft(1,1,istate), &
                            size(one_e_dm_beta_ao_for_dft,1) )
 enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, one_body_dm_mo_alpha_one_det, (mo_num,mo_num, N_states)]
&BEGIN_PROVIDER [double precision, one_body_dm_mo_beta_one_det, (mo_num,mo_num, N_states)]
 implicit none
 BEGIN_DOC
! One body density matrix on the |MO| basis for a single determinant
 END_DOC
 integer :: i
 one_body_dm_mo_alpha_one_det = 0.d0
 one_body_dm_mo_beta_one_det = 0.d0
 do i =1, elec_alpha_num
  one_body_dm_mo_alpha_one_det(i,i, 1:N_states) = 1.d0
 enddo
 do i =1, elec_beta_num
  one_body_dm_mo_beta_one_det(i,i, 1:N_states) = 1.d0
 enddo
END_PROVIDER
