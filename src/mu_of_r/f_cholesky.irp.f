BEGIN_PROVIDER [integer, list_couple_orb_r1, (2,n_couple_orb_r1)]
 implicit none
 integer :: ii,i,mm,m,itmp
 itmp = 0
  do ii = 1, n_occ_val_orb_for_hf(1)
   i = list_valence_orb_for_hf(ii,1)
   do mm = 1, n_basis_orb ! electron 1 
    m = list_basis(mm)
    itmp += 1
    list_couple_orb_r1(1,itmp) = i
    list_couple_orb_r1(2,itmp) = m
   enddo
  enddo
END_PROVIDER 


BEGIN_PROVIDER [integer, list_couple_orb_r2, (2,n_couple_orb_r2)]
 implicit none
 integer :: ii,i,mm,m,itmp
 itmp = 0
  do ii = 1, n_occ_val_orb_for_hf(2)
   i = list_valence_orb_for_hf(ii,2)
   do mm = 1, n_basis_orb ! electron 1 
    m = list_basis(mm)
    itmp += 1
    list_couple_orb_r2(1,itmp) = i
    list_couple_orb_r2(2,itmp) = m
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
 allocate(test(cholesky_mo_num,n_points_final_grid))
 test = 0.d0
 do ipoint = 1, n_points_final_grid
  do itmp = 1, n_couple_orb_r1
    i = list_couple_orb_r1(1,itmp)
    m = list_couple_orb_r1(2,itmp) 
    mo_i_r1 = mos_in_r_array_omp(i,ipoint)
    mo_b_r1 = mos_in_r_array_omp(m,ipoint)
    do mm = 1, cholesky_mo_num
     test(mm,ipoint) += mo_i_r1 * mo_b_r1 * mo_chol_r1(mm,itmp)
    enddo
   enddo
  enddo
 double precision :: accu
 accu = 0.d0
 do ipoint = 1, n_points_final_grid
  do mm = 1, cholesky_mo_num
   accu += dabs(mos_times_cholesky_r1(mm,ipoint) - test(mm,ipoint) )
   if(dabs(mos_times_cholesky_r1(mm,ipoint) - test(mm,ipoint)).gt.1.d-10)then
    print*,'problem ! ',dabs(mos_times_cholesky_r1(mm,ipoint) - test(mm,ipoint)) & 
                       ,     mos_times_cholesky_r1(mm,ipoint) , test(mm,ipoint) 
   endif
  enddo
 enddo
 print*,'accu = ',accu
   

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
 allocate(test(cholesky_mo_num,n_points_final_grid))
 test = 0.d0
 do ipoint = 1, n_points_final_grid
  do itmp = 1, n_couple_orb_r2
    i = list_couple_orb_r2(1,itmp)
    m = list_couple_orb_r2(2,itmp) 
    mo_i_r2 = mos_in_r_array_omp(i,ipoint)
    mo_b_r2 = mos_in_r_array_omp(m,ipoint)
    do mm = 1, cholesky_mo_num
     test(mm,ipoint) += mo_i_r2 * mo_b_r2 * mo_chol_r2(mm,itmp)
    enddo
   enddo
  enddo
 double precision :: accu
 accu = 0.d0
 do ipoint = 1, n_points_final_grid
  do mm = 1, cholesky_mo_num
   accu += dabs(mos_times_cholesky_r2(mm,ipoint) - test(mm,ipoint) )
   if(dabs(mos_times_cholesky_r2(mm,ipoint) - test(mm,ipoint)).gt.1.d-10)then
    print*,'problem ! ',dabs(mos_times_cholesky_r2(mm,ipoint) - test(mm,ipoint)) & 
                       ,     mos_times_cholesky_r2(mm,ipoint) , test(mm,ipoint) 
   endif
  enddo
 enddo
 print*,'accu = ',accu

END_PROVIDER 


BEGIN_PROVIDER [ double precision, f_hf_cholesky, (n_points_final_grid)]
 implicit none
 integer :: ipoint
 !!f(R) =  \sum_{I} \sum_{J} Phi_I(R) Phi_J(R) V_IJ
 !!     =  \sum_{I}\sum_{J}\sum_A Phi_I(R) Phi_J(R) V_AI V_AJ
 !!     =  \sum_A \sum_{I}Phi_I(R)V_AI \sum_{J}V_AJ Phi_J(R)
 !!     =  \sum_A V_AR G_AR 
 !! V_AR = \sum_{I}Phi_IR V_AI = \sum_{I}Phi^t_RI V_AI
 double precision :: u_dot_v
 do ipoint = 1, n_points_final_grid
  f_hf_cholesky(ipoint) = 2.D0 * u_dot_v(mos_times_cholesky_r2(1,ipoint),mos_times_cholesky_r1(1,ipoint),cholesky_mo_num)
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, on_top_hf_grid, (n_points_final_grid)]
 implicit none
 integer :: ipoint,i,ii
 double precision :: dm_a, dm_b
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
END_PROVIDER 

