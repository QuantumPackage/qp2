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
 else if (density_for_dft .EQ. "input_density_ao")then
  call ao_to_mo(data_one_e_dm_alpha_mo,size(data_one_e_dm_alpha_mo,1),one_e_dm_mo_alpha_for_dft,size(one_e_dm_mo_alpha_for_dft,1))
 else if (density_for_dft .EQ. "WFT")then
  provide mo_coef
  one_e_dm_mo_alpha_for_dft = one_e_dm_mo_alpha
 else if (density_for_dft .EQ. "KS")then
  provide mo_coef
  one_e_dm_mo_alpha_for_dft = one_body_dm_mo_alpha_one_det
 else if (density_for_dft .EQ. "state_average_dens")then
  one_e_dm_mo_alpha_for_dft = 0.d0
  one_e_dm_mo_alpha_for_dft(:,:,1) = one_e_dm_mo_alpha_average(:,:)
 endif
 
 if(no_core_density)then
  integer :: ii,i,j
  do ii = 1, n_core_orb
   i = list_core(ii)
   do j = 1, mo_num
    one_e_dm_mo_alpha_for_dft(j,i,:) = 0.d0
    one_e_dm_mo_alpha_for_dft(i,j,:) = 0.d0
   enddo
  enddo
  if(normalize_dm)then
   double precision :: elec_alpha_frozen_num, elec_alpha_valence(N_states)
   elec_alpha_frozen_num = elec_alpha_num - n_core_orb 
   elec_alpha_valence = 0.d0
   integer :: istate 
   do istate = 1, N_states
    do i = 1, mo_num
     elec_alpha_valence(istate) += one_e_dm_mo_alpha_for_dft(i,i,istate)
    enddo
    elec_alpha_valence(istate) = elec_alpha_frozen_num/elec_alpha_valence(istate)
    if( dabs(elec_alpha_valence(istate)) .lt.1.d-12)then
     one_e_dm_mo_alpha_for_dft = 0.d0
    else
     one_e_dm_mo_alpha_for_dft(:,:,istate) = one_e_dm_mo_alpha_for_dft(:,:,istate) * elec_alpha_valence(istate)
    endif
   enddo 
   
  endif
 endif

END_PROVIDER

BEGIN_PROVIDER [double precision, one_e_dm_mo_beta_for_dft, (mo_num,mo_num, N_states)]
 implicit none
 BEGIN_DOC
! density matrix for beta  electrons in the MO basis used for all DFT calculations based on the density
 END_DOC
 double precision :: delta_beta(mo_num,mo_num,N_states)
 one_e_dm_mo_beta_for_dft = 0.d0
 if(density_for_dft .EQ. "damping_rs_dft")then
  delta_beta = one_e_dm_mo_beta - data_one_e_dm_beta_mo
  one_e_dm_mo_beta_for_dft = data_one_e_dm_beta_mo + damping_for_rs_dft * delta_beta
 else if (density_for_dft .EQ. "input_density")then
  one_e_dm_mo_beta_for_dft = data_one_e_dm_beta_mo
 else if (density_for_dft .EQ. "input_density_ao")then
  call ao_to_mo(data_one_e_dm_beta_mo,size(data_one_e_dm_beta_mo,1),one_e_dm_mo_beta_for_dft,size(one_e_dm_mo_beta_for_dft,1))
 else if (density_for_dft .EQ. "WFT")then
  provide mo_coef
  one_e_dm_mo_beta_for_dft = one_e_dm_mo_beta
 else if (density_for_dft .EQ. "KS")then
  provide mo_coef
  one_e_dm_mo_beta_for_dft = one_body_dm_mo_beta_one_det
 else if (density_for_dft .EQ. "state_average_dens")then
  one_e_dm_mo_beta_for_dft = 0.d0
  one_e_dm_mo_beta_for_dft(:,:,1) = one_e_dm_mo_beta_average(:,:)
 endif


 if(no_core_density)then
  integer :: ii,i,j
  do ii = 1, n_core_orb
   i = list_core(ii)
   do j = 1, mo_num
    one_e_dm_mo_beta_for_dft(j,i,:) = 0.d0
    one_e_dm_mo_beta_for_dft(i,j,:) = 0.d0
   enddo
  enddo
  double precision :: elec_beta_valence(N_states),elec_beta_frozen_num
  integer :: istate 
  if(normalize_dm)then
   elec_beta_frozen_num = elec_beta_num - n_core_orb 
   elec_beta_valence = 0.d0
   do istate = 1, N_states
    do i = 1, mo_num
     elec_beta_valence(istate) += one_e_dm_mo_beta_for_dft(i,i,istate)
    enddo
    if(dabs(elec_beta_valence(istate)).lt.1.d-12)then
     one_e_dm_mo_beta_for_dft = 0.d0
    else
     elec_beta_valence(istate) = elec_beta_frozen_num/elec_beta_valence(istate)
     one_e_dm_mo_beta_for_dft(:,:,istate) = one_e_dm_mo_beta_for_dft(:,:,istate) * elec_beta_valence(istate)
    endif
   enddo 
  endif
 endif
END_PROVIDER

BEGIN_PROVIDER [double precision, one_e_dm_mo_for_dft, (mo_num,mo_num, N_states)]
 implicit none
 one_e_dm_mo_for_dft = one_e_dm_mo_beta_for_dft + one_e_dm_mo_alpha_for_dft
END_PROVIDER

BEGIN_PROVIDER [double precision, one_e_dm_average_mo_for_dft, (mo_num,mo_num)]
 implicit none
 integer :: i
 one_e_dm_average_mo_for_dft = one_e_dm_average_alpha_mo_for_dft + one_e_dm_average_beta_mo_for_dft
END_PROVIDER


BEGIN_PROVIDER [double precision, one_e_dm_average_alpha_mo_for_dft, (mo_num,mo_num)]
 implicit none
 integer :: i
 one_e_dm_average_alpha_mo_for_dft = 0.d0
 do i = 1, N_states
  one_e_dm_average_alpha_mo_for_dft(:,:) +=  one_e_dm_mo_alpha_for_dft(:,:,i) * state_average_weight(i)
 enddo
END_PROVIDER


BEGIN_PROVIDER [double precision, one_e_dm_average_beta_mo_for_dft, (mo_num,mo_num)]
 implicit none
 integer :: i
 one_e_dm_average_beta_mo_for_dft = 0.d0
 do i = 1, N_states
  one_e_dm_average_beta_mo_for_dft(:,:) +=  one_e_dm_mo_beta_for_dft(:,:,i) * state_average_weight(i)
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

 if (density_for_dft .EQ. "input_density_ao")then
  one_e_dm_alpha_ao_for_dft = data_one_e_dm_alpha_ao 
  one_e_dm_beta_ao_for_dft = data_one_e_dm_beta_ao 
 else
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
 endif

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




BEGIN_PROVIDER [double precision, one_e_dm_mo_alpha_for_dft_no_core, (mo_num,mo_num, N_states)]
 implicit none
 BEGIN_DOC
! density matrix for alpha electrons in the MO basis without the core orbitals
 END_DOC
 one_e_dm_mo_alpha_for_dft_no_core = one_e_dm_mo_alpha_for_dft 

 integer :: ii,i,j
 do ii = 1, n_core_orb
  i = list_core(ii)
  do j = 1, mo_num
   one_e_dm_mo_alpha_for_dft_no_core(j,i,:) = 0.d0
   one_e_dm_mo_alpha_for_dft_no_core(i,j,:) = 0.d0
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, one_e_dm_mo_beta_for_dft_no_core, (mo_num,mo_num, N_states)]
 implicit none
 BEGIN_DOC
! density matrix for beta  electrons in the MO basis without the core orbitals 
 END_DOC
 one_e_dm_mo_beta_for_dft_no_core = one_e_dm_mo_beta_for_dft 
 integer :: ii,i,j
 do ii = 1, n_core_orb
  i = list_core(ii)
  do j = 1, mo_num
   one_e_dm_mo_beta_for_dft_no_core(j,i,:) = 0.d0
   one_e_dm_mo_beta_for_dft_no_core(i,j,:) = 0.d0
  enddo
 enddo
END_PROVIDER

 BEGIN_PROVIDER [ double precision, one_e_dm_alpha_ao_for_dft_no_core, (ao_num,ao_num,N_states) ]
&BEGIN_PROVIDER [ double precision, one_e_dm_beta_ao_for_dft_no_core, (ao_num,ao_num,N_states) ]
 BEGIN_DOC
! one body density matrix on the AO basis based on one_e_dm_mo_alpha_for_dft_no_core
 END_DOC
 implicit none
 integer :: istate
 double precision :: mo_alpha,mo_beta

 one_e_dm_alpha_ao_for_dft_no_core = 0.d0
 one_e_dm_beta_ao_for_dft_no_core = 0.d0
 do istate = 1, N_states
  call mo_to_ao_no_overlap( one_e_dm_mo_alpha_for_dft_no_core(1,1,istate), &
                            size(one_e_dm_mo_alpha_for_dft_no_core,1),     &
                            one_e_dm_alpha_ao_for_dft_no_core(1,1,istate), &
                            size(one_e_dm_alpha_ao_for_dft_no_core,1) )
  call mo_to_ao_no_overlap( one_e_dm_mo_beta_for_dft_no_core(1,1,istate), &
                            size(one_e_dm_mo_beta_for_dft_no_core,1),     &
                            one_e_dm_beta_ao_for_dft_no_core(1,1,istate), &
                            size(one_e_dm_beta_ao_for_dft_no_core,1) )
 enddo

END_PROVIDER

