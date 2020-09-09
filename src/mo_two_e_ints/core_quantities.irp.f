BEGIN_PROVIDER [double precision, core_energy]
 implicit none
 BEGIN_DOC
! energy from the core : contains all core-core contributions
 END_DOC
 integer :: i,j,k,l
 core_energy = 0.d0
 do i = 1, n_core_orb
  j = list_core(i)
  core_energy += 2.d0 * mo_one_e_integrals_diag(j) + mo_two_e_integrals_jj(j,j)
  do k = i+1, n_core_orb
   l = list_core(k)
   core_energy += 2.d0 * (2.d0 * mo_two_e_integrals_jj(j,l) - mo_two_e_integrals_jj_exchange(j,l))
  enddo
 enddo
 core_energy += nuclear_repulsion

END_PROVIDER

BEGIN_PROVIDER [double precision, core_fock_operator, (mo_num,mo_num)]
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: get_two_e_integral
 BEGIN_DOC
! this is the contribution to the Fock operator from the core electrons
 END_DOC
 core_fock_operator = 0.d0
 do i = 1, n_act_orb
  j = list_act(i)
  do k = 1, n_act_orb
   l = list_act(k)
   do m = 1, n_core_orb
    n = list_core(m)
    core_fock_operator(j,l) += 2.d0 * get_two_e_integral(j,n,l,n,mo_integrals_map) - get_two_e_integral(j,n,n,l,mo_integrals_map)
   enddo
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [complex*16, core_fock_operator_complex, (mo_num,mo_num)]
 implicit none
 integer :: i,j,k,l,m,n
 complex*16 :: get_two_e_integral_complex
 BEGIN_DOC
! this is the contribution to the Fock operator from the core electrons
 END_DOC
 core_fock_operator_complex = (0.d0,0.d0)
 do i = 1, n_act_orb
  j = list_act(i)
  do k = 1, n_act_orb
   l = list_act(k)
   do m = 1, n_core_orb
    n = list_core(m)
    core_fock_operator_complex(j,l) += 2.d0 * &
                    get_two_e_integral_complex(j,n,l,n,mo_integrals_map,mo_integrals_map_2) - &
                    get_two_e_integral_complex(j,n,n,l,mo_integrals_map,mo_integrals_map_2)
   enddo
  enddo
 enddo
END_PROVIDER
