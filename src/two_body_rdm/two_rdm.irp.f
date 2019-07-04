 BEGIN_PROVIDER [double precision, two_rdm_alpha_beta_mo, (mo_num,mo_num,mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, two_rdm_alpha_alpha_mo, (mo_num,mo_num,mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, two_rdm_beta_beta_mo, (mo_num,mo_num,mo_num,mo_num,N_states)]
  implicit none
  BEGIN_DOC
  !  two_rdm_alpha_beta(i,j,k,l) = <Psi| a^{dagger}_{j,alpha} a^{dagger}_{l,beta} a_{k,beta} a_{i,alpha} | Psi>
  !                     1 1 2 2  = chemist notations
  !  note that no 1/2 factor is introduced in order to take into acccount for the spin symmetry
  !
  END_DOC
  integer                        :: dim1,dim2,dim3,dim4
  double precision               :: cpu_0,cpu_1
  dim1 = mo_num
  dim2 = mo_num
  dim3 = mo_num
  dim4 = mo_num
  two_rdm_alpha_beta_mo = 0.d0
  two_rdm_alpha_alpha_mo= 0.d0
  two_rdm_beta_beta_mo  = 0.d0
  print*,'providing two_rdm_alpha_beta ...'
  call wall_time(cpu_0)
  call all_two_rdm_dm_nstates(two_rdm_alpha_alpha_mo,two_rdm_beta_beta_mo,two_rdm_alpha_beta_mo,dim1,dim2,dim3,dim4,psi_coef,size(psi_coef,2),size(psi_coef,1))
  call wall_time(cpu_1)
  print*,'two_rdm_alpha_beta provided in',dabs(cpu_1-cpu_0)
  
END_PROVIDER


 BEGIN_PROVIDER [double precision, two_rdm_alpha_beta_mo_physicist, (mo_num,mo_num,mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, two_rdm_alpha_alpha_mo_physicist, (mo_num,mo_num,mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, two_rdm_beta_beta_mo_physicist, (mo_num,mo_num,mo_num,mo_num,N_states)]
  implicit none
  BEGIN_DOC
  !  two_rdm_alpha_beta_mo_physicist,(i,j,k,l) = <Psi| a^{dagger}_{k,alpha} a^{dagger}_{l,beta} a_{j,beta} a_{i,alpha} | Psi>
  !                                   1 2 1 2  = physicist notations
  !  note that no 1/2 factor is introduced in order to take into acccount for the spin symmetry
  !
  END_DOC
  integer                        :: i,j,k,l,istate
  double precision               :: cpu_0,cpu_1
  two_rdm_alpha_beta_mo_physicist = 0.d0
  print*,'providing two_rdm_alpha_beta_mo_physicist ...'
  call wall_time(cpu_0)
  do istate = 1, N_states
    do i = 1, mo_num
      do j = 1, mo_num
        do k = 1, mo_num
          do l = 1, mo_num
            !                               1 2 1 2                                 1 1 2 2
            two_rdm_alpha_beta_mo_physicist(l,k,i,j,istate) = two_rdm_alpha_beta_mo(i,l,j,k,istate)
            two_rdm_alpha_alpha_mo_physicist(l,k,i,j,istate) = two_rdm_alpha_alpha_mo(i,l,j,k,istate)
            two_rdm_beta_beta_mo_physicist(l,k,i,j,istate) = two_rdm_beta_beta_mo(i,l,j,k,istate)
          enddo
        enddo
      enddo
    enddo
  enddo
  call wall_time(cpu_1)
  print*,'two_rdm_alpha_beta_mo_physicist provided in',dabs(cpu_1-cpu_0)
  
END_PROVIDER

