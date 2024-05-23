BEGIN_PROVIDER [integer, list_couple_hf_orb_r1, (2,n_couple_orb_r1)]
 implicit none
 integer :: ii,i,mm,m,itmp
 itmp = 0
  do ii = 1, n_occ_val_orb_for_hf(1)
   i = list_valence_orb_for_hf(ii,1)
   do mm = 1, n_basis_orb ! electron 1 
    m = list_basis(mm)
    itmp += 1
    list_couple_hf_orb_r1(1,itmp) = i
    list_couple_hf_orb_r1(2,itmp) = m
   enddo
  enddo
END_PROVIDER 


BEGIN_PROVIDER [integer, list_couple_hf_orb_r2, (2,n_couple_orb_r2)]
 implicit none
 integer :: ii,i,mm,m,itmp
 itmp = 0
  do ii = 1, n_occ_val_orb_for_hf(2)
   i = list_valence_orb_for_hf(ii,2)
   do mm = 1, n_basis_orb ! electron 1 
    m = list_basis(mm)
    itmp += 1
    list_couple_hf_orb_r2(1,itmp) = i
    list_couple_hf_orb_r2(2,itmp) = m
   enddo
  enddo
END_PROVIDER 


BEGIN_PROVIDER [integer, n_couple_orb_r1] 
 implicit none
 BEGIN_DOC
 ! number of couples of alpha occupied times any basis orbital
 END_DOC
 n_couple_orb_r1 = n_occ_val_orb_for_hf(1) * n_basis_orb
END_PROVIDER 

BEGIN_PROVIDER [integer, n_couple_orb_r2] 
 implicit none
 BEGIN_DOC
 ! number of couples of beta occupied times any basis orbital
 END_DOC
 n_couple_orb_r2 = n_occ_val_orb_for_hf(2) * n_basis_orb
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mos_times_cholesky_r1, (cholesky_mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! V1_AR = \sum_{I}V_AI Phi_IR where "R" specifies the index of the grid point and A the number of cholesky point
 ! 
 ! here Phi_IR is phi_i(R)xphi_b(R) for r1 and V_AI = (ib|A) chollesky vector
 END_DOC
 double precision, allocatable :: mos_ib_r1(:,:),mo_chol_r1(:,:)
 double precision, allocatable :: test(:,:)
 double precision :: mo_i_r1,mo_b_r1
 integer :: ii,i,mm,m,itmp,ipoint,ll
 allocate(mos_ib_r1(n_couple_orb_r1,n_points_final_grid))
 allocate(mo_chol_r1(cholesky_mo_num,n_couple_orb_r1))

 do ipoint = 1, n_points_final_grid
  itmp = 0
  do ii = 1, n_occ_val_orb_for_hf(1)
   i = list_valence_orb_for_hf(ii,1)
   mo_i_r1 = mos_in_r_array_omp(i,ipoint)
   do mm = 1, n_basis_orb ! electron 1 
    m = list_basis(mm)
    mo_b_r1 = mos_in_r_array_omp(m,ipoint)
    itmp += 1
    mos_ib_r1(itmp,ipoint) = mo_i_r1 * mo_b_r1
   enddo
  enddo
 enddo

 itmp = 0
 do ii = 1, n_occ_val_orb_for_hf(1)
  i = list_valence_orb_for_hf(ii,1)
  do mm = 1, n_basis_orb ! electron 1 
   m = list_basis(mm)
   itmp += 1
   do ll = 1, cholesky_mo_num
    mo_chol_r1(ll,itmp) = cholesky_mo_transp(ll,m,i)
   enddo
   enddo
  enddo

 call get_AB_prod(mo_chol_r1,cholesky_mo_num,n_couple_orb_r1,mos_ib_r1,n_points_final_grid,mos_times_cholesky_r1)
   

END_PROVIDER 

BEGIN_PROVIDER [ double precision, mos_times_cholesky_r2, (cholesky_mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! V1_AR = \sum_{I}V_AI Phi_IR where "R" specifies the index of the grid point and A the number of cholesky point
 ! 
 ! here Phi_IR is phi_i(R)xphi_b(R) for r2 and V_AI = (ib|A) chollesky vector
 END_DOC
 double precision, allocatable :: mos_ib_r2(:,:),mo_chol_r2(:,:)
 double precision, allocatable :: test(:,:)
 double precision :: mo_i_r2,mo_b_r2
 integer :: ii,i,mm,m,itmp,ipoint,ll
 allocate(mos_ib_r2(n_couple_orb_r2,n_points_final_grid))
 allocate(mo_chol_r2(cholesky_mo_num,n_couple_orb_r2))

 do ipoint = 1, n_points_final_grid
  itmp = 0
  do ii = 1, n_occ_val_orb_for_hf(2)
   i = list_valence_orb_for_hf(ii,2)
   mo_i_r2 = mos_in_r_array_omp(i,ipoint)
   do mm = 1, n_basis_orb ! electron 1 
    m = list_basis(mm)
    mo_b_r2 = mos_in_r_array_omp(m,ipoint)
    itmp += 1
    mos_ib_r2(itmp,ipoint) = mo_i_r2 * mo_b_r2
   enddo
  enddo
 enddo

 itmp = 0
 do ii = 1, n_occ_val_orb_for_hf(2)
  i = list_valence_orb_for_hf(ii,2)
  do mm = 1, n_basis_orb ! electron 1 
   m = list_basis(mm)
   itmp += 1
   do ll = 1, cholesky_mo_num
    mo_chol_r2(ll,itmp) = cholesky_mo_transp(ll,m,i)
   enddo
   enddo
  enddo

 call get_AB_prod(mo_chol_r2,cholesky_mo_num,n_couple_orb_r2,mos_ib_r2,n_points_final_grid,mos_times_cholesky_r2)

END_PROVIDER 


BEGIN_PROVIDER [ double precision, f_hf_cholesky, (n_points_final_grid)]
 implicit none
 integer :: ipoint,m,k
 !!f(R) =  \sum_{I} \sum_{J} Phi_I(R) Phi_J(R) V_IJ
 !!     =  \sum_{I}\sum_{J}\sum_A Phi_I(R) Phi_J(R) V_AI V_AJ
 !!     =  \sum_A \sum_{I}Phi_I(R)V_AI \sum_{J}V_AJ Phi_J(R)
 !!     =  \sum_A V_AR G_AR 
 !! V_AR = \sum_{I}Phi_IR V_AI = \sum_{I}Phi^t_RI V_AI
 double precision :: u_dot_v,wall0,wall1
 if(elec_alpha_num == elec_beta_num)then
  print*,'providing f_hf_cholesky ...'
  call wall_time(wall0)
  provide mos_times_cholesky_r1
  !$OMP PARALLEL DO &
  !$OMP DEFAULT (NONE)  &
  !$OMP PRIVATE (ipoint,m) & 
  !$OMP ShARED (mos_times_cholesky_r1,cholesky_mo_num,f_hf_cholesky,n_points_final_grid) 
   do ipoint = 1, n_points_final_grid
    f_hf_cholesky(ipoint) = 0.d0
    do m = 1, cholesky_mo_num
     f_hf_cholesky(ipoint) =  f_hf_cholesky(ipoint) + &
       mos_times_cholesky_r1(m,ipoint) * mos_times_cholesky_r1(m,ipoint)
    enddo
    f_hf_cholesky(ipoint) *= 2.D0
   enddo
  !$OMP END PARALLEL DO
  
  call wall_time(wall1)
  print*,'Time to provide f_hf_cholesky = ',wall1-wall0
  free mos_times_cholesky_r1
 else
  print*,'providing f_hf_cholesky ...'
  call wall_time(wall0)
  provide mos_times_cholesky_r2 mos_times_cholesky_r1
  !$OMP PARALLEL DO &
  !$OMP DEFAULT (NONE)  &
  !$OMP PRIVATE (ipoint,m) & 
  !$OMP ShARED (mos_times_cholesky_r2,mos_times_cholesky_r1,cholesky_mo_num,f_hf_cholesky,n_points_final_grid) 
  do ipoint = 1, n_points_final_grid
   f_hf_cholesky(ipoint) = 0.D0
    do m = 1, cholesky_mo_num
     f_hf_cholesky(ipoint) =  f_hf_cholesky(ipoint) + &
            mos_times_cholesky_r2(m,ipoint)*mos_times_cholesky_r1(m,ipoint)
    enddo
    f_hf_cholesky(ipoint) *= 2.D0
  enddo
  !$OMP END PARALLEL DO
  call wall_time(wall1)
  print*,'Time to provide f_hf_cholesky = ',wall1-wall0
  free mos_times_cholesky_r2 mos_times_cholesky_r1
 endif
END_PROVIDER 

BEGIN_PROVIDER [ double precision, f_hf_cholesky_sparse, (n_points_final_grid)]
 implicit none
 integer :: ipoint,m,mm,i,ii,p
 !!f(R) =  \sum_{I} \sum_{J} Phi_I(R) Phi_J(R) V_IJ
 !!     =  \sum_{I}\sum_{J}\sum_A Phi_I(R) Phi_J(R) V_AI V_AJ
 !!     =  \sum_A \sum_{I}Phi_I(R)V_AI \sum_{J}V_AJ Phi_J(R)
 !!     =  \sum_A V_AR G_AR 
 !! V_AR = \sum_{I}Phi_IR V_AI = \sum_{I}Phi^t_RI V_AI
 double precision :: u_dot_v,wall0,wall1,accu_1, accu_2,mo_i_r1,mo_b_r1
 double precision :: thresh_1,thresh_2
 double precision, allocatable :: accu_vec(:)
 thresh_2 = ao_cholesky_threshold * 100.d0
 thresh_1 = dsqrt(thresh_2)
 provide cholesky_mo_transp
 if(elec_alpha_num == elec_beta_num)then
  call wall_time(wall0)
  !$OMP PARALLEL DEFAULT(NONE)                                      &
  !$OMP PRIVATE (accu_vec,ipoint,p,ii,i,mm,m,mo_i_r1,mo_b_r1) & 
  !$OMP ShARED (n_occ_val_orb_for_hf,list_valence_orb_for_hf,list_basis,mos_in_r_array_omp,thresh_1,thresh_2) & 
  !$OMP ShARED (cholesky_mo_num,f_hf_cholesky_sparse,n_points_final_grid,cholesky_mo_transp,n_basis_orb) 
  allocate(accu_vec(cholesky_mo_num))
  !$OMP DO 
   do ipoint = 1, n_points_final_grid
    f_hf_cholesky_sparse(ipoint) = 0.d0
     accu_vec = 0.d0
     do ii = 1, n_occ_val_orb_for_hf(1)
      i = list_valence_orb_for_hf(ii,1)
      mo_i_r1 = mos_in_r_array_omp(i,ipoint)
      if(dabs(mo_i_r1).lt.thresh_1)cycle
      do mm = 1, n_basis_orb ! electron 1 
       m = list_basis(mm)
       mo_b_r1 = mos_in_r_array_omp(m,ipoint)
       if(dabs(mo_i_r1*mo_b_r1).lt.thresh_2)cycle
       do p = 1, cholesky_mo_num
        accu_vec(p) += mo_i_r1 * mo_b_r1 * cholesky_mo_transp(p,m,i)
       enddo
      enddo
     enddo
     do p = 1, cholesky_mo_num
      f_hf_cholesky_sparse(ipoint) += accu_vec(p) * accu_vec(p)
     enddo
    f_hf_cholesky_sparse(ipoint) *= 2.D0
   enddo
  !$OMP END DO
  deallocate(accu_vec)
  !$OMP END PARALLEL
  
  call wall_time(wall1)
  print*,'Time to provide f_hf_cholesky_sparse = ',wall1-wall0
 else
  call wall_time(wall0)
  !$OMP PARALLEL DO &
  !$OMP DEFAULT (NONE)  &
  !$OMP PRIVATE (accu_2,accu_1,ipoint,p,ii,i,mm,m,mo_i_r1,mo_b_r1) & 
  !$OMP ShARED (n_occ_val_orb_for_hf,list_valence_orb_for_hf,list_basis,mos_in_r_array_omp) & 
  !$OMP ShARED (cholesky_mo_num,f_hf_cholesky_sparse,n_points_final_grid,cholesky_mo,n_basis_orb) 
   do ipoint = 1, n_points_final_grid
    f_hf_cholesky_sparse(ipoint) = 0.d0
    do p = 1, cholesky_mo_num
     accu_2 = 0.d0
     do ii = 1, n_occ_val_orb_for_hf(2)
      i = list_valence_orb_for_hf(ii,2)
      mo_i_r1 = mos_in_r_array_omp(i,ipoint)
      do mm = 1, n_basis_orb ! electron 1 
       m = list_basis(mm)
       mo_b_r1 = mos_in_r_array_omp(m,ipoint)
       accu_2 += mo_i_r1 * mo_b_r1 * cholesky_mo(m,i,p)
      enddo
     enddo
     accu_1 = accu_2
     do ii = n_occ_val_orb_for_hf(2)+1,n_occ_val_orb_for_hf(1)
      i = list_valence_orb_for_hf(ii,1)
      mo_i_r1 = mos_in_r_array_omp(i,ipoint)
      do mm = 1, n_basis_orb ! electron 1 
       m = list_basis(mm)
       mo_b_r1 = mos_in_r_array_omp(m,ipoint)
       accu_1 += mo_i_r1 * mo_b_r1 * cholesky_mo(m,i,p)
      enddo
     enddo
     f_hf_cholesky_sparse(ipoint) += accu_1 * accu_2
    enddo
    f_hf_cholesky_sparse(ipoint) *= 2.D0
   enddo
  !$OMP END PARALLEL DO
  call wall_time(wall1)
  print*,'Time to provide f_hf_cholesky_sparse = ',wall1-wall0
 endif
END_PROVIDER 

BEGIN_PROVIDER [ double precision, on_top_hf_grid, (n_points_final_grid)]
 implicit none
 integer :: ipoint,i,ii
 double precision :: dm_a, dm_b,wall0,wall1
 print*,'providing on_top_hf_grid ...'
 provide mos_in_r_array_omp
 call wall_time(wall0)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,dm_a,dm_b,ii,i) & 
 !$OMP ShARED (n_points_final_grid,n_occ_val_orb_for_hf,mos_in_r_array_omp,list_valence_orb_for_hf,on_top_hf_grid) 
 do ipoint = 1, n_points_final_grid
  dm_a = 0.d0
  do ii = 1, n_occ_val_orb_for_hf(1)
   i = list_valence_orb_for_hf(ii,1)
   dm_a += mos_in_r_array_omp(i,ipoint)*mos_in_r_array_omp(i,ipoint)
  enddo
  dm_b = 0.d0
  do ii = 1, n_occ_val_orb_for_hf(2)
   i = list_valence_orb_for_hf(ii,2)
   dm_b += mos_in_r_array_omp(i,ipoint)*mos_in_r_array_omp(i,ipoint)
  enddo
   on_top_hf_grid(ipoint) = 2.D0 * dm_a*dm_b
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide on_top_hf_grid = ',wall1-wall0
END_PROVIDER 

