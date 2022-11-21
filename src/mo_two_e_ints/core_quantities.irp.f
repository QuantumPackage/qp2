BEGIN_PROVIDER [double precision, core_energy]
 implicit none
 BEGIN_DOC
! energy from the core : contains all core-core contributions
 END_DOC
 integer :: i,j,k,l
 core_energy = 0.d0
 do i = 1, n_core_orb
  j = list_core(i)
  core_energy += 2.d0 * mo_one_e_integrals(j,j) + mo_two_e_integrals_jj(j,j)
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

BEGIN_PROVIDER [ double precision, h_core_ri, (mo_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Core Hamiltonian with 3-index exchange integrals:
 !
 ! $\tilde{h}{pq} = h_{pq} - \frac{1}{2}\sum_{k} g(pk,kq)$
 END_DOC

 integer :: i,j, k

 do j=1,mo_num
   do i=1,mo_num
     h_core_ri(i,j) = mo_one_e_integrals(i,j)
   enddo
   do k=1,mo_num
     do i=1,mo_num
       h_core_ri(i,j) = h_core_ri(i,j) - 0.5d0 * big_array_exchange_integrals(k,i,j)
     enddo
   enddo
 enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, h_act_ri, (mo_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Active Hamiltonian with 3-index exchange integrals:
 !
 ! $\tilde{h}{pq} = h_{pq} - \frac{1}{2}\sum_{k} g(pk,kq)$
 END_DOC

 integer :: i,j, k
 integer :: p,q, r
 ! core-core contribution
 h_act_ri = core_fock_operator
 !print *,' Bef----hact(1,14)=',h_act_ri(4,14)
 ! act-act contribution
 do p=1,n_act_orb
   j=list_act(p)
   do q=1,n_act_orb
     i=list_act(q)
     h_act_ri(i,j) = mo_one_e_integrals(i,j)
   enddo
   do r=1,n_act_orb
     k=list_act(r)
     do q=1,n_act_orb
       i=list_act(q)
       h_act_ri(i,j) = h_act_ri(i,j) - 0.5 * big_array_exchange_integrals(k,i,j)
     enddo
   enddo
 enddo
 ! core-act contribution
 !do p=1,n_act_orb
 !  j=list_core(p)
 !  do k=1,n_core_orb
 !    do q=1,n_act_orb
 !      i=list_act(q)
 !      h_act_ri(i,j) = h_act_ri(i,j) - 0.5 * big_array_exchange_integrals(k,i,j)
 !    enddo
 !  enddo
 !enddo
 !print *,' Aft----hact(1,14)=',h_act_ri(4,14), mo_one_e_integrals(4,14)
END_PROVIDER

