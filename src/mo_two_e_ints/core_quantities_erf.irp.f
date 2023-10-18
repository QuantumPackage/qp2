BEGIN_PROVIDER [double precision, core_energy_erf]
 implicit none
 BEGIN_DOC
! energy from the core : contains all core-core contributionswith the erf interaction
 END_DOC
 integer :: i,j,k,l
 core_energy_erf = 0.d0
 do i = 1, n_core_orb
  j = list_core(i)
  core_energy_erf +=  mo_two_e_int_erf_jj(j,j)
  do k = i+1, n_core_orb
   l = list_core(k)
   core_energy_erf += 2.d0 * (2.d0 * mo_two_e_int_erf_jj(j,l) - mo_two_e_int_erf_jj_exchange(j,l))
  enddo
 enddo
 core_energy_erf += nuclear_repulsion

END_PROVIDER

BEGIN_PROVIDER [double precision, core_fock_operator_erf, (mo_num,mo_num)]
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: get_mo_two_e_integral_erf
 BEGIN_DOC
! this is the contribution to the Fock operator from the core electrons with the erf interaction
 END_DOC
 core_fock_operator_erf = 0.d0
 do i = 1, n_act_orb
  j = list_act(i)
  do k = 1, n_act_orb
   l = list_act(k)
   do m = 1, n_core_orb
    n = list_core(m)
    core_fock_operator_erf(j,l) += 2.d0 * get_mo_two_e_integral_erf(j,n,l,n,mo_integrals_erf_map) - get_mo_two_e_integral_erf(j,n,n,l,mo_integrals_erf_map)
   enddo
  enddo
 enddo
END_PROVIDER
