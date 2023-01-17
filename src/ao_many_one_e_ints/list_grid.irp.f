 BEGIN_PROVIDER [ integer, n_pts_grid_ao_prod, (ao_num, ao_num)]
&BEGIN_PROVIDER [ integer, max_n_pts_grid_ao_prod]
 implicit none
 integer :: i,j,ipoint
 double precision :: overlap, r(3),thr, overlap_abs_gauss_r12_ao,overlap_gauss_r12_ao
 double precision :: sigma,dist,center_ij(3),fact_gauss, alpha, center(3)
 n_pts_grid_ao_prod = 0
 thr = 1.d-11
 print*,' expo_good_j_mu_1gauss = ',expo_good_j_mu_1gauss
 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, r, overlap, thr,fact_gauss, alpha, center,dist,sigma,center_ij) &
 !$OMP SHARED  (n_points_final_grid, ao_num, ao_overlap_abs_grid,n_pts_grid_ao_prod,expo_good_j_mu_1gauss,&
 !$OMP          final_grid_points,ao_prod_center,ao_prod_sigma,ao_nucl)
 !$OMP DO
 do i = 1, ao_num
! do i = 3,3
  do j = 1, ao_num
! do i = 22,22
!  do j = 9,9
   center_ij(1:3) = ao_prod_center(1:3,j,i)
   sigma = ao_prod_sigma(j,i)
   sigma *= sigma
   sigma = 0.5d0 /sigma
!   if(dabs(ao_overlap_abs_grid(j,i)).lt.1.d-10)cycle
   do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    dist  = (center_ij(1) - r(1))*(center_ij(1) - r(1))
    dist += (center_ij(2) - r(2))*(center_ij(2) - r(2))
    dist += (center_ij(3) - r(3))*(center_ij(3) - r(3))
    dist = dsqrt(dist)
    call gaussian_product(sigma, center_ij, expo_good_j_mu_1gauss, r, fact_gauss, alpha, center)
!    print*,''
!    print*,j,i,ao_overlap_abs_grid(j,i),ao_overlap_abs(j,i)
!    print*,r
!    print*,dist,sigma
!    print*,fact_gauss
    if( fact_gauss*ao_overlap_abs_grid(j,i).lt.1.d-11)cycle
    if(ao_nucl(i) == ao_nucl(j))then
     overlap = overlap_abs_gauss_r12_ao(r, expo_good_j_mu_1gauss, i, j)
    else
     overlap = overlap_gauss_r12_ao(r, expo_good_j_mu_1gauss, i, j)
    endif
!    print*,overlap
    if(dabs(overlap).lt.thr)cycle
    n_pts_grid_ao_prod(j,i) += 1
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 integer :: list(ao_num)
 do i = 1, ao_num
  list(i) = maxval(n_pts_grid_ao_prod(:,i))
 enddo
 max_n_pts_grid_ao_prod = maxval(list) 
END_PROVIDER 
