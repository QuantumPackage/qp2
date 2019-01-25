 BEGIN_PROVIDER [ double precision, psi_energy_h_core, (N_states) ]
  implicit none
  integer :: i
  integer :: j,k
  double precision :: tmp(mo_num,mo_num),mono_ints(mo_num,mo_num)
  BEGIN_DOC
! psi_energy_h_core = $\langle \Psi | h_{core} |\Psi \rangle$
!
! computed using the :c:data:`one_e_dm_mo_alpha` +
! :c:data:`one_e_dm_mo_beta` and :c:data:`mo_one_e_integrals`
  END_DOC
  psi_energy_h_core = 0.d0
  do i = 1, N_states
   do j = 1, mo_num
    do k = 1, mo_num
     psi_energy_h_core(i) += mo_one_e_integrals(k,j) * (one_e_dm_mo_alpha(k,j,i) + one_e_dm_mo_beta(k,j,i))
    enddo
   enddo
  enddo
 double precision :: accu
 accu = 0.d0
 do i = 1, N_states
  do j = 1, mo_num
   accu += one_e_dm_mo_alpha(j,j,i) + one_e_dm_mo_beta(j,j,i)
  enddo
  accu = (elec_alpha_num + elec_beta_num ) / accu
  psi_energy_h_core(i) = psi_energy_h_core(i) * accu
 enddo
END_PROVIDER
