 BEGIN_PROVIDER [double precision, effective_one_e_potential, (mo_num, mo_num,N_states)]
&BEGIN_PROVIDER [double precision, effective_one_e_potential_without_kin, (mo_num, mo_num,N_states)]
 implicit none
 integer :: i,j,istate
 effective_one_e_potential = 0.d0
 BEGIN_DOC
! Effective_one_e_potential(i,j) = $\rangle i_{MO}| v_{H}^{sr} |j_{MO}\rangle  + \rangle i_{MO}| h_{core} |j_{MO}\rangle  + \rangle i_{MO}|v_{xc} |j_{MO}\rangle$
!
! on the |MO| basis
! 
! Taking the expectation value does not provide any energy, but
!
! effective_one_e_potential(i,j) is the potential coupling DFT and WFT parts 
!
! and it is used in any RS-DFT based calculations  
 END_DOC
 do istate = 1, N_states
  do j = 1, mo_num
   do i = 1, mo_num

    effective_one_e_potential(i,j,istate) = short_range_Hartree_operator(i,j,istate) + mo_integrals_n_e(i,j) + mo_kinetic_integrals(i,j)   &
                                   + 0.5d0 * (potential_x_alpha_mo(i,j,istate) + potential_c_alpha_mo(i,j,istate)                          &
                                   +          potential_x_beta_mo(i,j,istate)  + potential_c_beta_mo(i,j,istate)   )

    effective_one_e_potential_without_kin(i,j,istate) = short_range_Hartree_operator(i,j,istate) + mo_integrals_n_e(i,j)                   &
                                   + 0.5d0 * (potential_x_alpha_mo(i,j,istate) + potential_c_alpha_mo(i,j,istate)                          &
                                   +          potential_x_beta_mo(i,j,istate)  + potential_c_beta_mo(i,j,istate)   )
   enddo
  enddo
 enddo
END_PROVIDER


 BEGIN_PROVIDER [double precision, effective_one_e_potential_sa, (mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, effective_one_e_potential_without_kin_sa, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! State-averaged potential in MO basis
 END_DOC

 integer :: istate

 effective_one_e_potential_sa(:,:) = 0.d0
 effective_one_e_potential_without_kin_sa(:,:) = 0.d0
 do istate = 1, N_states
   effective_one_e_potential_sa(:,:) += effective_one_e_potential(:,:,istate) * state_average_weight(istate)
   effective_one_e_potential_without_kin_sa(:,:) += effective_one_e_potential_without_kin(:,:,istate) * state_average_weight(istate)
 enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, ao_effective_one_e_potential, (ao_num, ao_num,N_states)]
&BEGIN_PROVIDER [double precision, ao_effective_one_e_potential_without_kin, (ao_num, ao_num,N_states)]
 implicit none
 integer :: i,j,istate
 ao_effective_one_e_potential = 0.d0
 ao_effective_one_e_potential_without_kin = 0.d0
 BEGIN_DOC
! Effective_one_e_potential(i,j) = $\rangle i_{AO}| v_{H}^{sr} |j_{AO}\rangle  + \rangle i_{AO}| h_{core} |j_{AO}\rangle  + \rangle i_{AO}|v_{xc} |j_{AO}\rangle$
!
! on the |MO| basis
! 
! Taking the expectation value does not provide any energy, but
!
! ao_effective_one_e_potential(i,j) is the potential coupling DFT and WFT parts 
!
! and it is used in any RS-DFT based calculations  
 END_DOC
 do istate = 1, N_states
    call mo_to_ao(effective_one_e_potential(1,1,istate),mo_num,ao_effective_one_e_potential(1,1,istate),ao_num)
    call mo_to_ao(effective_one_e_potential_without_kin(1,1,istate),mo_num,ao_effective_one_e_potential_without_kin(1,1,istate),ao_num)
 enddo
END_PROVIDER


 BEGIN_PROVIDER [double precision, ao_effective_one_e_potential_sa, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_effective_one_e_potential_without_kin_sa, (ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! State-averaged potential in AO basis
 END_DOC

 integer :: istate

 ao_effective_one_e_potential_sa(:,:) = 0.d0
 ao_effective_one_e_potential_without_kin_sa(:,:) = 0.d0
 do istate = 1, N_states
   ao_effective_one_e_potential_sa(:,:) += ao_effective_one_e_potential(:,:,istate) * state_average_weight(istate)
   ao_effective_one_e_potential_without_kin_sa(:,:) += ao_effective_one_e_potential_without_kin(:,:,istate) * state_average_weight(istate)
 enddo

END_PROVIDER 

