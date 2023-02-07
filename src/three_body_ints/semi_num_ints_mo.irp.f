
BEGIN_PROVIDER [ double precision, mo_v_ij_erf_rk_cst_mu_naive, ( mo_num, mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1 )/(2|r - R|) on the MO basis
 END_DOC
 integer :: i,j,k,l,ipoint
 do ipoint = 1, n_points_final_grid
  mo_v_ij_erf_rk_cst_mu_naive(:,:,ipoint) = 0.d0
  do i = 1, mo_num
   do j = 1, mo_num
    do k = 1, ao_num
     do l = 1, ao_num
      mo_v_ij_erf_rk_cst_mu_naive(j,i,ipoint) += mo_coef(l,j) * 0.5d0 * v_ij_erf_rk_cst_mu(l,k,ipoint) * mo_coef(k,i)
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_v_ij_erf_rk_cst_mu, ( mo_num, mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1)/(2|r - R|) on the MO basis
 END_DOC
 integer :: ipoint
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (ipoint) & 
 !$OMP SHARED (n_points_final_grid,v_ij_erf_rk_cst_mu,mo_v_ij_erf_rk_cst_mu)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   call ao_to_mo(v_ij_erf_rk_cst_mu(1,1,ipoint),size(v_ij_erf_rk_cst_mu,1),mo_v_ij_erf_rk_cst_mu(1,1,ipoint),size(mo_v_ij_erf_rk_cst_mu,1))
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 mo_v_ij_erf_rk_cst_mu = mo_v_ij_erf_rk_cst_mu * 0.5d0
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_v_ij_erf_rk_cst_mu_transp, ( n_points_final_grid,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1)/(2|r - R|) on the MO basis
 END_DOC
 integer :: ipoint,i,j
 do i = 1, mo_num
  do j = 1, mo_num
   do ipoint = 1, n_points_final_grid
    mo_v_ij_erf_rk_cst_mu_transp(ipoint,j,i) = mo_v_ij_erf_rk_cst_mu(j,i,ipoint)
   enddo
  enddo
 enddo
 FREE mo_v_ij_erf_rk_cst_mu
END_PROVIDER 


BEGIN_PROVIDER [ double precision, mo_x_v_ij_erf_rk_cst_mu_naive, ( mo_num, mo_num,3,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr  x * phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1 )/|r - R| on the MO basis
 END_DOC
 integer :: i,j,k,l,ipoint,m
 do ipoint = 1, n_points_final_grid
  mo_x_v_ij_erf_rk_cst_mu_naive(:,:,:,ipoint) = 0.d0
  do i = 1, mo_num
   do j = 1, mo_num
    do m = 1, 3
     do k = 1, ao_num
      do l = 1, ao_num
       mo_x_v_ij_erf_rk_cst_mu_naive(j,i,m,ipoint) += mo_coef(l,j) * 0.5d0 * x_v_ij_erf_rk_cst_mu_transp(l,k,m,ipoint) * mo_coef(k,i)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_x_v_ij_erf_rk_cst_mu, ( mo_num, mo_num,3,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr x * phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1)/2|r - R| on the MO basis
 END_DOC
 integer :: ipoint,m
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (ipoint,m) & 
 !$OMP SHARED (n_points_final_grid,x_v_ij_erf_rk_cst_mu_transp,mo_x_v_ij_erf_rk_cst_mu)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
  do m = 1, 3
   call ao_to_mo(x_v_ij_erf_rk_cst_mu_transp(1,1,m,ipoint),size(x_v_ij_erf_rk_cst_mu_transp,1),mo_x_v_ij_erf_rk_cst_mu(1,1,m,ipoint),size(mo_x_v_ij_erf_rk_cst_mu,1))
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 mo_x_v_ij_erf_rk_cst_mu = 0.5d0 * mo_x_v_ij_erf_rk_cst_mu

END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_x_v_ij_erf_rk_cst_mu_transp, (n_points_final_grid,3, mo_num, mo_num)]
 implicit none
 integer :: i,j,m,ipoint
 do i = 1, mo_num
  do j = 1, mo_num
   do m = 1, 3
    do ipoint = 1, n_points_final_grid
     mo_x_v_ij_erf_rk_cst_mu_transp(ipoint,m,j,i) = mo_x_v_ij_erf_rk_cst_mu(j,i,m,ipoint)
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [ double precision, x_W_ij_erf_rk, ( n_points_final_grid,3,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! W_mn^X(R) = \int dr phi_m(r) phi_n(r) (1 - erf(mu |r-R|)) (x-X)
 END_DOC
 include 'constants.include.F'
 integer :: ipoint,m,i,j
 double precision :: xyz,cst
 double precision :: wall0, wall1

 cst = 0.5d0 * inv_sq_pi
 print*,'providing x_W_ij_erf_rk ...'
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (ipoint,m,i,j,xyz) & 
 !$OMP SHARED (x_W_ij_erf_rk,n_points_final_grid,mo_x_v_ij_erf_rk_cst_mu_transp,mo_v_ij_erf_rk_cst_mu_transp,mo_num,final_grid_points) 
 !$OMP DO SCHEDULE (dynamic)
 do i = 1, mo_num
  do j = 1, mo_num
   do m = 1, 3
    do ipoint = 1, n_points_final_grid
     xyz = final_grid_points(m,ipoint)
     x_W_ij_erf_rk(ipoint,m,j,i)  =  mo_x_v_ij_erf_rk_cst_mu_transp(ipoint,m,j,i) - xyz * mo_v_ij_erf_rk_cst_mu_transp(ipoint,j,i)
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 FREE mo_v_ij_erf_rk_cst_mu_transp 
 FREE mo_x_v_ij_erf_rk_cst_mu_transp
 call wall_time(wall1)
 print*,'time to provide x_W_ij_erf_rk = ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, sqrt_weight_at_r, (n_points_final_grid)]
 implicit none
 integer :: ipoint
 do ipoint = 1, n_points_final_grid
  sqrt_weight_at_r(ipoint) = dsqrt(final_weight_at_r_vector(ipoint))
 enddo
END_PROVIDER 

!BEGIN_PROVIDER [ double precision, mos_in_r_array_transp_sq_weight, (n_points_final_grid,mo_num)]


!BEGIN_PROVIDER [ double precision, gauss_ij_rk_transp, (ao_num, ao_num, n_points_final_grid) ]
! implicit none
! integer :: i,j,ipoint
! do ipoint = 1, n_points_final_grid
!  do j = 1, ao_num
!   do i = 1, ao_num
!    gauss_ij_rk_transp(i,j,ipoint) = gauss_ij_rk(ipoint,i,j)
!   enddo
!  enddo
! enddo
!END_PROVIDER 
!
!
!BEGIN_PROVIDER [ double precision, mo_gauss_ij_rk, ( mo_num, mo_num,n_points_final_grid)]
! implicit none
! integer :: ipoint
! !$OMP PARALLEL                  &
! !$OMP DEFAULT (NONE)            &
! !$OMP PRIVATE (ipoint) & 
! !$OMP SHARED (n_points_final_grid,gauss_ij_rk_transp,mo_gauss_ij_rk)
! !$OMP DO SCHEDULE (dynamic)
! do ipoint = 1, n_points_final_grid
!  call ao_to_mo(gauss_ij_rk_transp(1,1,ipoint),size(gauss_ij_rk_transp,1),mo_gauss_ij_rk(1,1,ipoint),size(mo_gauss_ij_rk,1))
! enddo
! !$OMP END DO
! !$OMP END PARALLEL
!
!END_PROVIDER 
!
!BEGIN_PROVIDER [ double precision, mo_gauss_ij_rk_transp, (n_points_final_grid, mo_num, mo_num)]
! implicit none
! integer :: i,j,ipoint
! do ipoint = 1, n_points_final_grid
!  do j = 1, mo_num
!   do i = 1, mo_num
!    mo_gauss_ij_rk_transp(ipoint,i,j) = mo_gauss_ij_rk(i,j,ipoint)
!   enddo
!  enddo
! enddo
!
!END_PROVIDER 
!
