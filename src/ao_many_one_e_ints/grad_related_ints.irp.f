BEGIN_PROVIDER [ double precision, v_ij_erf_rk_cst_mu, ( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) (erf(mu(R) |r - R| - 1)/|r - R|
 END_DOC
 integer :: i,j,ipoint
 double precision :: mu,r(3),NAI_pol_mult_erf_ao
 double precision :: int_mu, int_coulomb
 provide mu_erf final_grid_points 
 double precision :: wall0, wall1
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,int_coulomb) & 
 !$OMP SHARED (ao_num,n_points_final_grid,v_ij_erf_rk_cst_mu,final_grid_points,mu_erf)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_erf
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = NAI_pol_mult_erf_ao(i,j,mu,r)
     int_coulomb = NAI_pol_mult_erf_ao(i,j,1.d+9,r)
     v_ij_erf_rk_cst_mu(j,i,ipoint)= (int_mu - int_coulomb )
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
     v_ij_erf_rk_cst_mu(j,i,ipoint)= v_ij_erf_rk_cst_mu(i,j,ipoint)
    enddo
   enddo
  enddo

 call wall_time(wall1)
 print*,'wall time for v_ij_erf_rk_cst_mu  ',wall1 - wall0
END_PROVIDER 

BEGIN_PROVIDER [ double precision, v_ij_erf_rk_cst_mu_transp, (n_points_final_grid, ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) (erf(mu(R) |r - R| - 1)/|r - R|
 END_DOC
 integer :: i,j,ipoint
 double precision :: mu,r(3),NAI_pol_mult_erf_ao
 double precision :: int_mu, int_coulomb
 provide mu_erf final_grid_points 
 double precision :: wall0, wall1
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,int_coulomb) & 
 !$OMP SHARED (ao_num,n_points_final_grid,v_ij_erf_rk_cst_mu_transp,final_grid_points,mu_erf)
 !$OMP DO SCHEDULE (dynamic)
 do i = 1, ao_num
  do j = i, ao_num
   do ipoint = 1, n_points_final_grid
     mu = mu_erf 
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = NAI_pol_mult_erf_ao(i,j,mu,r)
     int_coulomb = NAI_pol_mult_erf_ao(i,j,1.d+9,r)
     v_ij_erf_rk_cst_mu_transp(ipoint,j,i)= (int_mu - int_coulomb )
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do i = 1, ao_num
  do j = 1, i-1
   do ipoint = 1, n_points_final_grid
     v_ij_erf_rk_cst_mu_transp(ipoint,j,i)= v_ij_erf_rk_cst_mu_transp(ipoint,i,j)
    enddo
   enddo
  enddo

 call wall_time(wall1)
 print*,'wall time for v_ij_erf_rk_cst_mu_transp  ',wall1 - wall0
END_PROVIDER 


BEGIN_PROVIDER [ double precision, x_v_ij_erf_rk_cst_mu_tmp, (3,ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr x * phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1)/|r - R|
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),ints(3),ints_coulomb(3)
 double precision :: wall0, wall1
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,ints,m,ints_coulomb) & 
 !$OMP SHARED (ao_num,n_points_final_grid,x_v_ij_erf_rk_cst_mu_tmp,final_grid_points,mu_erf)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = i, ao_num
    mu = mu_erf
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    call NAI_pol_x_mult_erf_ao(i,j,mu,r,ints)
    call NAI_pol_x_mult_erf_ao(i,j,1.d+9,r,ints_coulomb)
    do m = 1, 3
     x_v_ij_erf_rk_cst_mu_tmp(m,j,i,ipoint) = ( ints(m) - ints_coulomb(m))
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
    do m = 1, 3
     x_v_ij_erf_rk_cst_mu_tmp(m,j,i,ipoint)= x_v_ij_erf_rk_cst_mu_tmp(m,i,j,ipoint)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall1)
 print*,'wall time for x_v_ij_erf_rk_cst_mu_tmp',wall1 - wall0


END_PROVIDER 

BEGIN_PROVIDER [ double precision, x_v_ij_erf_rk_cst_mu, (ao_num, ao_num,n_points_final_grid,3)]
 implicit none
 BEGIN_DOC
! int dr x * phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1)/|r - R|
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),ints,ints_coulomb
 double precision :: wall0, wall1
 call wall_time(wall0)
 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, ao_num
    do m = 1, 3
     x_v_ij_erf_rk_cst_mu(j,i,ipoint,m)= x_v_ij_erf_rk_cst_mu_tmp(m,j,i,ipoint)
    enddo
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for x_v_ij_erf_rk_cst_mu',wall1 - wall0

END_PROVIDER 


BEGIN_PROVIDER [ double precision, x_v_ij_erf_rk_cst_mu_transp, (ao_num, ao_num,3,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr x * phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1)/|r - R|
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),ints,ints_coulomb
 double precision :: wall0, wall1
 call wall_time(wall0)
 do ipoint = 1, n_points_final_grid
  do m = 1, 3
   do i = 1, ao_num
    do j = 1, ao_num
     x_v_ij_erf_rk_cst_mu_transp(j,i,m,ipoint)= x_v_ij_erf_rk_cst_mu_tmp(m,j,i,ipoint)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall1)
 print*,'wall time for x_v_ij_erf_rk_cst_mu_transp',wall1 - wall0


END_PROVIDER 

BEGIN_PROVIDER [ double precision, x_v_ij_erf_rk_cst_mu_transp_bis, (n_points_final_grid,ao_num, ao_num,3)]
 implicit none
 BEGIN_DOC
! int dr x * phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1)/|r - R|
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),ints,ints_coulomb
 double precision :: wall0, wall1
 call wall_time(wall0)
 do m = 1, 3
  do i = 1, ao_num
   do j = 1, ao_num
    do ipoint = 1, n_points_final_grid
     x_v_ij_erf_rk_cst_mu_transp_bis(ipoint,j,i,m)= x_v_ij_erf_rk_cst_mu_tmp(m,j,i,ipoint)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall1)
 print*,'wall time for x_v_ij_erf_rk_cst_mu_transp',wall1 - wall0


END_PROVIDER 


BEGIN_PROVIDER [ double precision, d_dx_v_ij_erf_rk_cst_mu_tmp, (3,n_points_final_grid,ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! d_dx_v_ij_erf_rk_cst_mu_tmp(m,R,j,i) = int dr phi_j(r)) (erf(mu(R) |r - R|) - 1)/|r - R| d/dx (phi_i(r) 
!
! with m == 1 -> d/dx , m == 2 -> d/dy , m == 3 -> d/dz
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),ints(3),ints_coulomb(3)
 double precision :: wall0, wall1
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,ints,m,ints_coulomb) & 
 !$OMP SHARED (ao_num,n_points_final_grid,d_dx_v_ij_erf_rk_cst_mu_tmp,final_grid_points,mu_erf)
 !$OMP DO SCHEDULE (dynamic)
 do i = 1, ao_num
  do j = 1, ao_num
   do ipoint = 1, n_points_final_grid
    mu = mu_erf 
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    call phi_j_erf_mu_r_dxyz_phi(j,i,mu, r, ints)
    call phi_j_erf_mu_r_dxyz_phi(j,i,1.d+9, r, ints_coulomb)
    do m = 1, 3
     d_dx_v_ij_erf_rk_cst_mu_tmp(m,ipoint,j,i) = ( ints(m) - ints_coulomb(m))
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 call wall_time(wall1)
 print*,'wall time for d_dx_v_ij_erf_rk_cst_mu_tmp',wall1 - wall0


END_PROVIDER 

BEGIN_PROVIDER [ double precision, d_dx_v_ij_erf_rk_cst_mu, (n_points_final_grid,ao_num, ao_num,3)]
 implicit none
 BEGIN_DOC
! d_dx_v_ij_erf_rk_cst_mu_tmp(j,i,R,m) = int dr phi_j(r)) (erf(mu(R) |r - R|) - 1)/|r - R| d/dx (phi_i(r) 
!
! with m == 1 -> d/dx , m == 2 -> d/dy , m == 3 -> d/dz
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),ints,ints_coulomb
 double precision :: wall0, wall1
 call wall_time(wall0)
 do i = 1, ao_num
  do j = 1, ao_num
   do m = 1, 3
    do ipoint = 1, n_points_final_grid
     d_dx_v_ij_erf_rk_cst_mu(ipoint,j,i,m)= d_dx_v_ij_erf_rk_cst_mu_tmp(m,ipoint,j,i)
    enddo
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for d_dx_v_ij_erf_rk_cst_mu',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, x_d_dx_v_ij_erf_rk_cst_mu_tmp, (3,n_points_final_grid,ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! x_d_dx_v_ij_erf_rk_cst_mu_tmp(m,j,i,R) = int dr x phi_j(r)) (erf(mu(R) |r - R|) - 1)/|r - R| d/dx (phi_i(r) 
!
! with m == 1 -> d/dx , m == 2 -> d/dy , m == 3 -> d/dz
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),ints(3),ints_coulomb(3)
 double precision :: wall0, wall1
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,ints,m,ints_coulomb) & 
 !$OMP SHARED (ao_num,n_points_final_grid,x_d_dx_v_ij_erf_rk_cst_mu_tmp,final_grid_points,mu_erf)
 !$OMP DO SCHEDULE (dynamic)
 do i = 1, ao_num
  do j = 1, ao_num
   do ipoint = 1, n_points_final_grid
    mu = mu_erf
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    call phi_j_erf_mu_r_xyz_dxyz_phi(j,i,mu, r, ints)
    call phi_j_erf_mu_r_xyz_dxyz_phi(j,i,1.d+9, r, ints_coulomb)
    do m = 1, 3
     x_d_dx_v_ij_erf_rk_cst_mu_tmp(m,ipoint,j,i) = ( ints(m) - ints_coulomb(m))
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 call wall_time(wall1)
 print*,'wall time for x_d_dx_v_ij_erf_rk_cst_mu_tmp',wall1 - wall0


END_PROVIDER 

BEGIN_PROVIDER [ double precision, x_d_dx_v_ij_erf_rk_cst_mu, (n_points_final_grid,ao_num, ao_num,3)]
 implicit none
 BEGIN_DOC
! x_d_dx_v_ij_erf_rk_cst_mu_tmp(j,i,R,m) = int dr x phi_j(r)) (erf(mu(R) |r - R|) - 1)/|r - R| d/dx (phi_i(r) 
!
! with m == 1 -> d/dx , m == 2 -> d/dy , m == 3 -> d/dz
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),ints,ints_coulomb
 double precision :: wall0, wall1
 call wall_time(wall0)
 do i = 1, ao_num
  do j = 1, ao_num
   do ipoint = 1, n_points_final_grid
    do m = 1, 3
     x_d_dx_v_ij_erf_rk_cst_mu(ipoint,j,i,m)= x_d_dx_v_ij_erf_rk_cst_mu_tmp(m,ipoint,j,i)
    enddo
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for x_d_dx_v_ij_erf_rk_cst_mu',wall1 - wall0

END_PROVIDER 

