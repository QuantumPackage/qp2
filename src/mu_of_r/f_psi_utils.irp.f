BEGIN_PROVIDER [double precision, core_inact_act_V_kl_contracted, (n_core_inact_act_orb,n_core_inact_act_orb,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! core_inact_act_V_kl_contracted(k,l,ipoint) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
!
! This is needed to build the function f_{\Psi^B}(X_1,X_2) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) 
! 
 END_DOC
 integer :: ipoint,k,l
 do k = 1, n_core_inact_act_orb
  do l = 1, n_core_inact_act_orb
   do ipoint = 1, n_points_final_grid
    core_inact_act_V_kl_contracted(k,l,ipoint) = full_occ_v_kl_cntrctd(ipoint,k,l)
   enddo
  enddo
 enddo
 free full_occ_v_kl_cntrctd 

END_PROVIDER 

BEGIN_PROVIDER [double precision, full_occ_2_rdm_cntrctd, (n_core_inact_act_orb,n_core_inact_act_orb,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! full_occ_2_rdm_cntrctd(k,l,ipoint,istate) = \sum_{ij} \Gamma_{ij}^{kl}  phi_i(r_ipoint) phi_j(r_ipoint) 
!
! where \Gamma_{ij}^{kl}(istate)  = <Psi_{istate}| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi_{istate}>
!
! This is needed to build the function f_{\Psi^B}(X_1,X_2) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) 
! 
 END_DOC
 integer :: ipoint,k,l,istate
 do istate = 1, N_states
  do k = 1, n_core_inact_act_orb
   do l = 1, n_core_inact_act_orb
    do ipoint = 1, n_points_final_grid
     full_occ_2_rdm_cntrctd(k,l,ipoint,istate) = full_occ_2_rdm_cntrctd_trans(ipoint,k,l,istate)
    enddo
   enddo
  enddo
 enddo
 free full_occ_2_rdm_cntrctd_trans 
END_PROVIDER 




BEGIN_PROVIDER [double precision, full_occ_2_rdm_cntrctd_trans, (n_points_final_grid,n_core_inact_act_orb,n_core_inact_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! full_occ_2_rdm_cntrctd_trans(ipoint,k,l,istate) = \sum_{ij} \Gamma_{ij}^{kl}  phi_i(r_ipoint) phi_j(r_ipoint) 
!
! where \Gamma_{ij}^{kl}(istate)  = <Psi_{istate}| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi_{istate}>
!
! This is needed to build the function f_{\Psi^B}(X_1,X_2) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) 
! 
 END_DOC
 integer :: i,j,k,l,istate
 integer :: ipoint
 double precision, allocatable :: mos_array_r(:),r(:)
 provide full_occ_2_rdm_ab_mo
 double precision :: wall0,wall1
 print*,'Providing  full_occ_2_rdm_cntrctd_trans ..... '
 call wall_time(wall0)
 full_occ_2_rdm_cntrctd_trans = 0.d0
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,k,l,i,j,istate) & 
 !$OMP SHARED  (n_core_inact_act_orb, n_points_final_grid, full_occ_2_rdm_cntrctd_trans, final_grid_points,full_occ_2_rdm_ab_mo,core_inact_act_mos_in_r_array,N_states )
 !$OMP DO              
 do istate = 1, N_states
  do l = 1, n_core_inact_act_orb ! 2 
   do k = 1, n_core_inact_act_orb ! 1 
    do ipoint = 1, n_points_final_grid
     do j = 1, n_core_inact_act_orb
      do i = 1, n_core_inact_act_orb
                                             !                                 1 2 1 2 
       full_occ_2_rdm_cntrctd_trans(ipoint,k,l,istate) += full_occ_2_rdm_ab_mo(i,j,k,l,istate) * core_inact_act_mos_in_r_array(j,ipoint) * core_inact_act_mos_in_r_array(i,ipoint)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide full_occ_2_rdm_cntrctd_trans = ',wall1 - wall0

END_PROVIDER 



BEGIN_PROVIDER [double precision, full_occ_v_kl_cntrctd, (n_points_final_grid,n_core_inact_act_orb,n_core_inact_act_orb)]
 implicit none
 BEGIN_DOC
! full_occ_v_kl_cntrctd(ipoint,k,l) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
!
! This is needed to build the function f_{\Psi^B}(X_1,X_2) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) 
! 
 END_DOC
 integer :: i,j,k,l,kk,ll,ii,jj
 integer :: ipoint
 double precision, allocatable :: integrals_array(:,:), mos_array_r(:),r(:), integrals_basis(:,:)
 ! just not to mess with parallelization
 allocate(integrals_array(mo_num,mo_num))
  k = 1
  l = 1
  call get_mo_two_e_integrals_ij(k,l,mo_num,integrals_array,mo_integrals_map) 
  provide basis_mos_in_r_array
 deallocate(integrals_array)
 double precision :: wall0,wall1
 call wall_time(wall0)

 full_occ_v_kl_cntrctd = 0.d0
 print*,'Providing  full_occ_v_kl_cntrctd ..... '
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,kk,ll,k,l,i,j,ii,jj,integrals_array,integrals_basis) & 
 !$OMP SHARED (mo_num, n_points_final_grid, n_basis_orb, list_basis, full_occ_v_kl_cntrctd, mo_integrals_map,final_grid_points,basis_mos_in_r_array, n_core_inact_act_orb, list_core_inact_act)
 allocate(integrals_array(mo_num,mo_num), integrals_basis(n_basis_orb,n_basis_orb))
 !$OMP DO              
  do l = 1, n_core_inact_act_orb! 2 
   ll = list_core_inact_act(l)
   do k = 1, n_core_inact_act_orb ! 1 
    kk = list_core_inact_act(k)
    call get_mo_two_e_integrals_ij(kk,ll,mo_num,integrals_array,mo_integrals_map) 
    do j = 1, n_basis_orb 
     jj = list_basis(j)
     do i = 1, n_basis_orb 
      ii = list_basis(i)
      integrals_basis(i,j) = integrals_array(ii,jj)
     enddo
    enddo
    do ipoint = 1, n_points_final_grid
     do j = 1, n_basis_orb ! condition on mo_num in order to ensure the correct CBS limit 
      do i = 1, n_basis_orb ! 
                                                           !1 2                     1 2 
       full_occ_v_kl_cntrctd(ipoint,k,l) += integrals_basis(i,j) * basis_mos_in_r_array(j,ipoint) * basis_mos_in_r_array(i,ipoint)
      enddo
     enddo
    enddo
   enddo
  enddo
 !$OMP END DO
 deallocate(integrals_array,integrals_basis)
 !$OMP END PARALLEL

 call wall_time(wall1)
 print*,'Time to provide full_occ_v_kl_cntrctd = ',wall1 - wall0
END_PROVIDER 

