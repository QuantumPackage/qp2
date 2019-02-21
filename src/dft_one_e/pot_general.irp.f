
 BEGIN_PROVIDER [double precision, potential_x_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! general providers for the alpha/beta exchange/correlation potentials on the AO basis
 END_DOC

  if(trim(exchange_functional)=="short_range_LDA")then
   potential_x_alpha_ao = potential_sr_x_alpha_ao_LDA
   potential_x_beta_ao = potential_sr_x_beta_ao_LDA
  else if(exchange_functional.EQ."short_range_PBE")then
   potential_x_alpha_ao = potential_sr_x_alpha_ao_PBE
   potential_x_beta_ao = potential_sr_x_beta_ao_PBE
  else if(trim(exchange_functional)=="LDA")then
   potential_x_alpha_ao = potential_x_alpha_ao_LDA
   potential_x_beta_ao = potential_x_beta_ao_LDA
  else if(exchange_functional.EQ."PBE")then
   potential_x_alpha_ao = potential_x_alpha_ao_PBE
   potential_x_beta_ao = potential_x_beta_ao_PBE
  else if(exchange_functional.EQ."my_functional")then
   potential_x_alpha_ao = potential_new_functional_x_alpha_ao
   potential_x_beta_ao  = potential_new_functional_x_beta_ao
  else if(exchange_functional.EQ."None")then
   potential_x_alpha_ao = 0.d0
   potential_x_beta_ao = 0.d0
  else
   print*, 'Exchange functional required does not exist ...'
   print*,'exchange_functional',exchange_functional
   stop
  endif

  if(trim(correlation_functional)=="short_range_LDA")then
   potential_c_alpha_ao = potential_sr_c_alpha_ao_LDA
   potential_c_beta_ao = potential_sr_c_beta_ao_LDA
  else if(trim(correlation_functional)=="LDA")then
   potential_c_alpha_ao = potential_c_alpha_ao_LDA
   potential_c_beta_ao = potential_c_beta_ao_LDA
  else if(correlation_functional.EQ."short_range_PBE")then
   potential_c_alpha_ao = potential_sr_c_alpha_ao_PBE
   potential_c_beta_ao = potential_sr_c_beta_ao_PBE
  else if(correlation_functional.EQ."PBE")then
   potential_c_alpha_ao = potential_c_alpha_ao_PBE
   potential_c_beta_ao = potential_c_beta_ao_PBE
  else if(correlation_functional.EQ."my_functional")then
   potential_c_alpha_ao = potential_new_functional_c_alpha_ao
   potential_c_beta_ao  = potential_new_functional_c_beta_ao
  else if(correlation_functional.EQ."None")then
   potential_c_alpha_ao = 0.d0
   potential_c_beta_ao = 0.d0
  else
   print*, 'Correlation functional required does not ecist ...'
   print*,'correlation_functional',correlation_functional
   stop
  endif


END_PROVIDER





 BEGIN_PROVIDER [double precision, potential_x_alpha_mo,(mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_mo,(mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_mo,(mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_mo,(mo_num,mo_num,N_states)]
 implicit none
 BEGIN_DOC
! general providers for the alpha/beta exchange/correlation potentials on the MO basis
 END_DOC
 integer :: istate
 do istate = 1, N_states
    call ao_to_mo(                                                   &
        potential_x_alpha_ao(1,1,istate),                                 &
        size(potential_x_alpha_ao,1),                                &
        potential_x_alpha_mo(1,1,istate),                                 &
        size(potential_x_alpha_mo,1)                                 &
        )

    call ao_to_mo(                                                   &
        potential_x_beta_ao(1,1,istate),                                  &
        size(potential_x_beta_ao,1),                                 &
        potential_x_beta_mo(1,1,istate),                                  &
        size(potential_x_beta_mo,1)                                  &
        )


    call ao_to_mo(                                                   &
        potential_c_alpha_ao(1,1,istate),                                 &
        size(potential_c_alpha_ao,1),                                &
        potential_c_alpha_mo(1,1,istate),                                 &
        size(potential_c_alpha_mo,1)                                 &
        )

    call ao_to_mo(                                                   &
        potential_c_beta_ao(1,1,istate),                                  &
        size(potential_c_beta_ao,1),                                 &
        potential_c_beta_mo(1,1,istate),                                  &
        size(potential_c_beta_mo,1)                                  &
        )

 enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, Trace_v_xc, (N_states)]
&BEGIN_PROVIDER [double precision, Trace_v_H, (N_states)]
&BEGIN_PROVIDER [double precision, Trace_v_Hxc, (N_states)]
 implicit none
 integer :: i,j,istate
 double precision :: dm
 BEGIN_DOC
! Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
! Trace_v_Hxc = \sum_{i,j} v^{H}_{ij} (rho_{ij}_\alpha + rho_{ij}_\beta)
! Trace_v_Hxc = \sum_{i,j} rho_{ij} v^{Hxc}_{ij}
 END_DOC
 do istate = 1, N_states
  Trace_v_xc(istate) = 0.d0
  Trace_v_H(istate) = 0.d0
  do i = 1, mo_num
   do j = 1, mo_num
     Trace_v_xc(istate) += (potential_x_alpha_mo(j,i,istate) + potential_c_alpha_mo(j,i,istate)) * one_e_dm_mo_alpha_for_dft(j,i,istate)
     Trace_v_xc(istate) += (potential_x_beta_mo(j,i,istate)  + potential_c_beta_mo(j,i,istate) ) * one_e_dm_mo_beta_for_dft(j,i,istate)
     dm = one_e_dm_mo_alpha_for_dft(j,i,istate) + one_e_dm_mo_beta_for_dft(j,i,istate)
     Trace_v_H(istate) += dm * short_range_Hartree_operator(j,i,istate)
   enddo
  enddo
  Trace_v_Hxc(istate) = Trace_v_xc(istate) + Trace_v_H(istate)
 enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, Trace_v_xc_new, (N_states)]
 implicit none
 integer :: i,j,istate
 double precision :: dm
 BEGIN_DOC
! Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
 END_DOC
 do istate = 1, N_states
  Trace_v_xc_new(istate) = 0.d0
  do i = 1, mo_num
   do j = 1, mo_num
     Trace_v_xc_new(istate) += (potential_xc_alpha_mo(j,i,istate) ) * one_e_dm_mo_alpha_for_dft(j,i,istate)
     Trace_v_xc_new(istate) += (potential_xc_beta_mo(j,i,istate)  ) * one_e_dm_mo_beta_for_dft(j,i,istate)
   enddo
  enddo
 enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, potential_xc_alpha_mo,(mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_xc_beta_mo,(mo_num,mo_num,N_states)]
 implicit none
 integer :: istate
 
 do istate = 1, N_states
    call ao_to_mo(                                                   &
        potential_xc_alpha_ao(1,1,istate),                                 &
        size(potential_xc_alpha_ao,1),                                &
        potential_xc_alpha_mo(1,1,istate),                                 &
        size(potential_xc_alpha_mo,1)                                 &
        )

    call ao_to_mo(                                                   &
        potential_xc_beta_ao(1,1,istate),                                  &
        size(potential_xc_beta_ao,1),                                 &
        potential_xc_beta_mo(1,1,istate),                                  &
        size(potential_xc_beta_mo,1)                                  &
        )
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, potential_xc_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_xc_beta_ao,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! general providers for the alpha/beta exchange/correlation potentials on the AO basis
 END_DOC

  if(trim(exchange_functional)=="short_range_LDA")then
   potential_xc_alpha_ao = potential_sr_xc_alpha_ao_LDA
   potential_xc_beta_ao  = potential_sr_xc_beta_ao_LDA
  else if(trim(exchange_functional)=="LDA")then
   potential_xc_alpha_ao = potential_xc_alpha_ao_LDA
   potential_xc_beta_ao  = potential_xc_beta_ao_LDA
  else if(exchange_functional.EQ."None")then
   potential_xc_alpha_ao = 0.d0
   potential_xc_beta_ao = 0.d0
  else if(trim(exchange_functional)=="short_range_PBE")then
   potential_xc_alpha_ao = potential_sr_xc_alpha_ao_PBE
   potential_xc_beta_ao  = potential_sr_xc_beta_ao_PBE
  else if(trim(exchange_functional)=="PBE")then
   potential_xc_alpha_ao = potential_xc_alpha_ao_PBE
   potential_xc_beta_ao  = potential_xc_beta_ao_PBE
  else if(exchange_functional.EQ."None")then
   potential_xc_alpha_ao = 0.d0
   potential_xc_beta_ao = 0.d0
  else
   print*, 'Exchange functional required does not exist ...'
   print*,'exchange_functional',exchange_functional
   stop
  endif

END_PROVIDER

