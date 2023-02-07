
! ---

BEGIN_PROVIDER [double precision, ao_non_hermit_term_chemist, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !                            1 1 2 2      1 2                                1 2 
  !
  ! ao_non_hermit_term_chemist(k,i,l,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the AO basis
  !
  END_DOC

  implicit none
  integer                       :: i, j, k, l, ipoint, m
  double precision              :: weight1, r(3)
  double precision              :: wall1, wall0
  double precision, allocatable :: b_mat(:,:,:,:), ac_mat(:,:,:,:)

  provide v_ij_erf_rk_cst_mu x_v_ij_erf_rk_cst_mu

  call wall_time(wall0)
  allocate(b_mat(n_points_final_grid,ao_num,ao_num,3), ac_mat(ao_num,ao_num,ao_num,ao_num))

 !$OMP PARALLEL                         &
 !$OMP DEFAULT (NONE)                   &
 !$OMP PRIVATE (i,k,m,ipoint,r,weight1) & 
 !$OMP SHARED (aos_in_r_array_transp,aos_grad_in_r_array_transp_bis,b_mat)& 
 !$OMP SHARED (ao_num,n_points_final_grid,final_grid_points,final_weight_at_r_vector)
 !$OMP DO SCHEDULE (static)
  do m = 1, 3
    do i = 1, ao_num
      do k = 1, ao_num
        do ipoint = 1, n_points_final_grid
          r(1) = final_grid_points(1,ipoint)
          r(2) = final_grid_points(2,ipoint)
          r(3) = final_grid_points(3,ipoint)
          weight1 = final_weight_at_r_vector(ipoint)
          b_mat(ipoint,k,i,m) = 0.5d0 * aos_in_r_array_transp(ipoint,k) * r(m) * weight1 * aos_grad_in_r_array_transp_bis(ipoint,i,m) 
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  ! (A)                b_mat(ipoint,k,i,m) X v_ij_erf_rk_cst_mu(j,l,r1)
  ! 1/2 \int dr1 x1 phi_k(1) d/dx1 phi_i(1) \int dr2 (1 - erf(mu_r12))/r12  phi_j(2) phi_l(2)
  ac_mat = 0.d0
  do m = 1, 3
    !           A   B^T  dim(A,1)       dim(B,2)       dim(A,2)        alpha * A                LDA 

    call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0             &
              , v_ij_erf_rk_cst_mu(1,1,1), ao_num*ao_num, b_mat(1,1,1,m), n_points_final_grid &
              , 1.d0, ac_mat, ao_num*ao_num) 

  enddo

 !$OMP PARALLEL                       &
 !$OMP DEFAULT (NONE)                 &
 !$OMP PRIVATE (i,k,m,ipoint,weight1) & 
 !$OMP SHARED (aos_in_r_array_transp,aos_grad_in_r_array_transp_bis,b_mat,ao_num,n_points_final_grid,final_weight_at_r_vector)
 !$OMP DO SCHEDULE (static)
  do m = 1, 3
    do i = 1, ao_num
      do k = 1, ao_num
        do ipoint = 1, n_points_final_grid
          weight1 = final_weight_at_r_vector(ipoint)
          b_mat(ipoint,k,i,m) = 0.5d0 * aos_in_r_array_transp(ipoint,k) * weight1 * aos_grad_in_r_array_transp_bis(ipoint,i,m) 
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

 ! (B)                b_mat(ipoint,k,i,m) X x_v_ij_erf_rk_cst_mu(j,l,r1,m)
 ! 1/2 \int dr1 phi_k(1) d/dx1 phi_i(1) \int dr2 x2(1 - erf(mu_r12))/r12  phi_j(2) phi_l(2)
  do m = 1, 3
   !           A   B^T  dim(A,1)       dim(B,2)       dim(A,2)        alpha * A                LDA 

    call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, -1.d0                &
              , x_v_ij_erf_rk_cst_mu(1,1,1,m), ao_num*ao_num, b_mat(1,1,1,m), n_points_final_grid &
              , 1.d0, ac_mat, ao_num*ao_num) 
  enddo

 !$OMP PARALLEL          &
 !$OMP DEFAULT (NONE)    &
 !$OMP PRIVATE (i,k,j,l) & 
 !$OMP SHARED (ac_mat,ao_non_hermit_term_chemist,ao_num)
 !$OMP DO SCHEDULE (static)
  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          !                          (ki|lj)           (ki|lj)           (lj|ki)
          ao_non_hermit_term_chemist(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i)    
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time dgemm ', wall1 - wall0

END_PROVIDER 

! ---

! TODO :: optimization :: transform into DGEM

BEGIN_PROVIDER [double precision, mo_non_hermit_term_chemist, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !                            1 1 2 2      1 2                                1 2 
  !
  ! mo_non_hermit_term_chemist(k,i,l,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the MO basis
  END_DOC

  implicit none
  integer                       :: i, j, k, l, m, n, p, q
  double precision, allocatable :: mo_tmp_1(:,:,:,:), mo_tmp_2(:,:,:,:)
 
  allocate(mo_tmp_1(mo_num,ao_num,ao_num,ao_num))
  mo_tmp_1 = 0.d0

  do m = 1, ao_num
    do p = 1, ao_num
      do n = 1, ao_num
        do q = 1, ao_num
          do k = 1, mo_num
            !       (k n|p m)    = sum_q c_qk * (q n|p m)
            mo_tmp_1(k,n,p,m) += mo_coef_transp(k,q) * ao_non_hermit_term_chemist(q,n,p,m)
          enddo
        enddo
      enddo
    enddo
  enddo
  free ao_non_hermit_term_chemist 

  allocate(mo_tmp_2(mo_num,mo_num,ao_num,ao_num))
  mo_tmp_2 = 0.d0

  do m = 1, ao_num
    do p = 1, ao_num
      do n = 1, ao_num
        do i = 1, mo_num
          do k = 1, mo_num
            !       (k i|p m) = sum_n c_ni * (k n|p m)
            mo_tmp_2(k,i,p,m) += mo_coef_transp(i,n) * mo_tmp_1(k,n,p,m)
          enddo
        enddo
      enddo
    enddo
  enddo
  deallocate(mo_tmp_1)

  allocate(mo_tmp_1(mo_num,mo_num,mo_num,ao_num))
  mo_tmp_1 = 0.d0

  do m = 1, ao_num
    do p = 1, ao_num
      do l = 1, mo_num
        do i = 1, mo_num
          do k = 1, mo_num
            mo_tmp_1(k,i,l,m) += mo_coef_transp(l,p) * mo_tmp_2(k,i,p,m)
          enddo
        enddo
      enddo
    enddo
  enddo
  deallocate(mo_tmp_2)

  mo_non_hermit_term_chemist = 0.d0
  do m = 1, ao_num
    do j = 1, mo_num
      do l = 1, mo_num
        do i = 1, mo_num
          do k = 1, mo_num
            mo_non_hermit_term_chemist(k,i,l,j) += mo_coef_transp(j,m) * mo_tmp_1(k,i,l,m)
          enddo
        enddo
      enddo
    enddo
  enddo
  deallocate(mo_tmp_1)

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, mo_non_hermit_term, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !                    1 2 1 2      1 2                                1 2 
  !
  ! mo_non_hermit_term(k,l,i,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the MO basis
  END_DOC

  implicit none
  integer :: i, j, k, l

  do j = 1, mo_num
    do i = 1, mo_num
      do l = 1, mo_num
        do k = 1, mo_num
          mo_non_hermit_term(k,l,i,j) = mo_non_hermit_term_chemist(k,i,l,j)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

