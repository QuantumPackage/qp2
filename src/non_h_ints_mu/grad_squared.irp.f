
! ---

!       TODO : strong optmization : write the loops in a different way
!            : for each couple of AO, the gaussian product are done once for all 

BEGIN_PROVIDER [ double precision, gradu_squared_u_ij_mu, (ao_num, ao_num, n_points_final_grid) ]

  BEGIN_DOC
  !
  ! if J(r1,r2) = u12:
  !
  ! gradu_squared_u_ij_mu = -0.50 x \int r2 [ (grad_1 u12)^2 + (grad_2 u12^2)] \phi_i(2) \phi_j(2)
  !                       = -0.25 x \int r2       (1 - erf(mu*r12))^2          \phi_i(2) \phi_j(2)
  ! and
  !   (1 - erf(mu*r12))^2 = \sum_i coef_gauss_1_erf_x_2(i) * exp(-expo_gauss_1_erf_x_2(i) * r12^2)
  !
  ! if J(r1,r2) = u12 x v1 x v2
  !
  ! gradu_squared_u_ij_mu = -0.50 x \int r2 \phi_i(2) \phi_j(2) [ v1^2  v2^2 ((grad_1 u12)^2 + (grad_2 u12^2)]) + u12^2  v2^2 (grad_1 v1)^2 + 2 u12 v1 v2^2 (grad_1 u12) . (grad_1 v1) ]
  !                       = -0.25 x        v1^2      \int r2 \phi_i(2) \phi_j(2) [1 - erf(mu r12)]^2 v2^2
  !                       + -0.50 x    (grad_1 v1)^2 \int r2 \phi_i(2) \phi_j(2)      u12^2          v2^2
  !                       + -1.00 x v1 (grad_1 v1)   \int r2 \phi_i(2) \phi_j(2)   (grad_1 u12)      v2^2 
  !                       =                 v1^2        x   int2_grad1u2_grad2u2_j1b2
  !                       + -0.5 x      (grad_1 v1)^2   x   int2_u2_j1b2
  !                       + -1.0 X V1 x (grad_1 v1)   \cdot [ int2_u_grad1u_j1b2 x r - int2_u_grad1u_x_j1b ]
  !
  !
  END_DOC
 
  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: x, y, z, r(3), delta, coef
  double precision           :: tmp_v, tmp_x, tmp_y, tmp_z
  double precision           :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing gradu_squared_u_ij_mu ...'
  call wall_time(time0)

  PROVIDE j1b_type

  if(j1b_type .eq. 3) then

    do ipoint = 1, n_points_final_grid

      x     = final_grid_points(1,ipoint)
      y     = final_grid_points(2,ipoint)
      z     = final_grid_points(3,ipoint)
      tmp_v = v_1b       (ipoint)
      tmp_x = v_1b_grad(1,ipoint)
      tmp_y = v_1b_grad(2,ipoint)
      tmp_z = v_1b_grad(3,ipoint)

      tmp1 = tmp_v * tmp_v
      tmp2 = -0.5d0 * (tmp_x * tmp_x + tmp_y * tmp_y + tmp_z * tmp_z)
      tmp3 = tmp_v * tmp_x
      tmp4 = tmp_v * tmp_y
      tmp5 = tmp_v * tmp_z

      tmp6 = -x * tmp3
      tmp7 = -y * tmp4
      tmp8 = -z * tmp5

      do j = 1, ao_num
        do i = 1, ao_num

          tmp9 = int2_u_grad1u_j1b2(i,j,ipoint)

          gradu_squared_u_ij_mu(i,j,ipoint) = tmp1 * int2_grad1u2_grad2u2_j1b2(i,j,ipoint)            &
                                            + tmp2 * int2_u2_j1b2             (i,j,ipoint)            &
                                            + tmp6 * tmp9 + tmp3 * int2_u_grad1u_x_j1b2(1,i,j,ipoint) &
                                            + tmp7 * tmp9 + tmp4 * int2_u_grad1u_x_j1b2(2,i,j,ipoint) &
                                            + tmp8 * tmp9 + tmp5 * int2_u_grad1u_x_j1b2(3,i,j,ipoint)
        enddo
      enddo
    enddo

  else

    gradu_squared_u_ij_mu = 0.d0
    do ipoint = 1, n_points_final_grid
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      do j = 1, ao_num
        do i = 1, ao_num
          do igauss = 1, n_max_fit_slat
            delta = expo_gauss_1_erf_x_2(igauss)
            coef  = coef_gauss_1_erf_x_2(igauss)
            gradu_squared_u_ij_mu(i,j,ipoint) += -0.25d0 * coef * overlap_gauss_r12_ao(r, delta, i, j)
          enddo
        enddo
      enddo
    enddo

  endif

  call wall_time(time1)
  print*, ' Wall time for gradu_squared_u_ij_mu = ', time1 - time0

END_PROVIDER 

! ---

!BEGIN_PROVIDER [double precision, tc_grad_square_ao, (ao_num, ao_num, ao_num, ao_num)]
!
!  BEGIN_DOC
!  !
!  ! tc_grad_square_ao(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_1 u(r1,r2)|^2 | ij>
!  !
!  END_DOC
!
!  implicit none
!  integer                       :: ipoint, i, j, k, l
!  double precision              :: weight1, ao_ik_r, ao_i_r
!  double precision, allocatable :: ac_mat(:,:,:,:)
!
!  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num))
!  ac_mat = 0.d0
!
!  do ipoint = 1, n_points_final_grid
!    weight1 = final_weight_at_r_vector(ipoint)
!
!    do i = 1, ao_num
!      ao_i_r = weight1 * aos_in_r_array_transp(ipoint,i)
!
!      do k = 1, ao_num
!        ao_ik_r = ao_i_r * aos_in_r_array_transp(ipoint,k)
!
!        do j = 1, ao_num
!          do l = 1, ao_num
!            ac_mat(k,i,l,j) += ao_ik_r * gradu_squared_u_ij_mu(l,j,ipoint)
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo
!
!  do j = 1, ao_num
!    do l = 1, ao_num
!      do i = 1, ao_num
!        do k = 1, ao_num
!          tc_grad_square_ao(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i)
!          !write(11,*) tc_grad_square_ao(k,i,l,j)
!        enddo
!      enddo
!    enddo
!  enddo
!
!  deallocate(ac_mat)
!
!END_PROVIDER 

! ---

 BEGIN_PROVIDER [ double precision, grad_1_squared_u_ij_mu_new, (n_points_final_grid, ao_num, ao_num)]
 implicit none
 integer :: ipoint,i,j,m,igauss
 BEGIN_DOC
 ! grad_1_squared_u_ij_mu(j,i,ipoint) = -1/2 \int dr2 phi_j(r2) phi_i(r2) |\grad_r1 u(r1,r2,\mu)|^2 
 ! |\grad_r1 u(r1,r2,\mu)|^2 = 1/4 * (1 - erf(mu*r12))^2 
 ! ! (1 - erf(mu*r12))^2 = \sum_i coef_gauss_1_erf_x_2(i) * exp(-expo_gauss_1_erf_x_2(i) * r12^2)
 END_DOC
  include 'constants.include.F'
 double precision :: r(3),delta,coef
 double precision :: overlap_gauss_r12_ao,time0,time1
 integer :: num_a,num_b,power_A(3), power_B(3),l,k
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta,analytical_j
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new, coef_i, coef_j, k_ab,center_new(3),p_new,c_tmp,coef_last        ! constant factor
 double precision  :: coefxy, coefx, coefy, coefz,coefxyz
 integer           :: d(3),lx,ly,lz,iorder_tmp(3),dim1
 double precision  :: overlap,overlap_x,overlap_y,overlap_z,thr
 dim1=100
 thr = 0.d0
 print*,'providing grad_1_squared_u_ij_mu_new ...'
 grad_1_squared_u_ij_mu_new = 0.d0
 call wall_time(time0)
 !TODO : strong optmization : write the loops in a different way
 !     : for each couple of AO, the gaussian product are done once for all 
 d = 0
 do i = 1, ao_num
  do j = 1, ao_num
   ! \int dr2 phi_j(r2) phi_i(r2) (1 - erf(mu*r12))^2 
   ! = \sum_i coef_gauss_1_erf_x_2(i) \int dr2 phi_j(r2) phi_i(r2) exp(-expo_gauss_1_erf_x_2(i) * (r_1 - r_2)^2)
   if(ao_overlap_abs(j,i).lt.1.d-12)then
    cycle
   endif
   num_A = ao_nucl(i)
   power_A(1:3)= ao_power(i,1:3)
   A_center(1:3) = nucl_coord(num_A,1:3)
   num_B = ao_nucl(j)
   power_B(1:3)= ao_power(j,1:3)
   B_center(1:3) = nucl_coord(num_B,1:3)
   do l=1,ao_prim_num(i)
    coef_i = ao_coef_normalized_ordered_transp(l,i)  
    alpha = ao_expo_ordered_transp(l,i)     
    do k=1,ao_prim_num(j)
     beta = ao_expo_ordered_transp(k,j)     
     coef_j = ao_coef_normalized_ordered_transp(k,j) 

     ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order                                
     ! from gaussian_A * gaussian_B
     call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , & 
                                      beta,alpha,power_B,power_A,B_center,A_center,n_pt_max_integrals)
     c_tmp = coef_i*coef_j*fact_a_new
     if(dabs(c_tmp).lt.thr)cycle
     do ipoint = 1, n_points_final_grid
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      do igauss = 1, n_max_fit_slat
       delta = expo_gauss_1_erf_x_2(igauss)
       coef  = coef_gauss_1_erf_x_2(igauss)
       coef_last = c_tmp * coef
       if(dabs(coef_last).lt.thr)cycle
       do lx = 0, iorder_a_new(1)
        coefx = A_new(lx,1)
        coefx *= coef_last 
!        if(dabs(coefx).lt.thr)cycle
        iorder_tmp(1) = lx
        do ly = 0, iorder_a_new(2)
         coefy = A_new(ly,2)
         coefxy= coefx*coefy
!         if(dabs(coefxy).lt.thr)cycle
         iorder_tmp(2) = ly
         do lz = 0, iorder_a_new(3)
          coefz = A_new(lz,3)
          coefxyz = coefz * coefxy
!          if(dabs(coefxyz).lt.thr)cycle
          iorder_tmp(3) = lz
!          call gaussian_product(alpha_new,A_center_new,delta,r,k_ab,p_new,center_new)
!          if(dabs(coef_last*k_ab).lt.thr)cycle
          call overlap_gaussian_xyz(A_center_new,r,alpha_new,delta,iorder_tmp,d,overlap_x,overlap_y,overlap_z,overlap,dim1)
          grad_1_squared_u_ij_mu_new(ipoint,j,i) += -0.25 * coefxyz * overlap
         enddo ! igauss
        enddo ! ipoint 
       enddo ! lz
      enddo ! ly
     enddo ! lx 
    enddo ! k
   enddo ! l
  enddo ! j
 enddo ! i
 call wall_time(time1)
 print*,'Wall time for grad_1_squared_u_ij_mu_new = ',time1 - time0
 END_PROVIDER 


BEGIN_PROVIDER [double precision, tc_grad_square_ao, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! tc_grad_square_ao(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_1 u(r1,r2)|^2 | ij>
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l
  double precision              :: weight1, ao_ik_r, ao_i_r
  double precision, allocatable :: ac_mat(:,:,:,:), bc_mat(:,:,:,:)

  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num))
  ac_mat = 0.d0
  allocate(bc_mat(ao_num,ao_num,ao_num,ao_num))
  bc_mat = 0.d0

  do ipoint = 1, n_points_final_grid
    weight1 = final_weight_at_r_vector(ipoint)

    do i = 1, ao_num
      ao_i_r = weight1 * aos_in_r_array_transp(ipoint,i)

      do k = 1, ao_num
        ao_ik_r = ao_i_r * aos_in_r_array_transp(ipoint,k)

        do j = 1, ao_num
          do l = 1, ao_num
            ac_mat(k,i,l,j) += ao_ik_r * ( u12sq_j1bsq(l,j,ipoint) + u12_grad1_u12_j1b_grad1_j1b(l,j,ipoint) )
            bc_mat(k,i,l,j) += ao_ik_r * grad12_j12(l,j,ipoint)
          enddo
        enddo
      enddo
    enddo
  enddo

  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          tc_grad_square_ao(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i) + bc_mat(k,i,l,j)
        enddo
      enddo
    enddo
  enddo

  deallocate(ac_mat)
  deallocate(bc_mat)

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, grad12_j12, (ao_num, ao_num, n_points_final_grid) ]
 
  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: r(3), delta, coef
  double precision           :: tmp1
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing grad12_j12 ...'
  call wall_time(time0)

  PROVIDE j1b_type

  if(j1b_type .eq. 3) then

    do ipoint = 1, n_points_final_grid
      tmp1 = v_1b(ipoint)
      tmp1 = tmp1 * tmp1
      do j = 1, ao_num
        do i = 1, ao_num
          grad12_j12(i,j,ipoint) = tmp1 * int2_grad1u2_grad2u2_j1b2(i,j,ipoint)
        enddo
      enddo
    enddo

  else

    grad12_j12 = 0.d0
    do ipoint = 1, n_points_final_grid
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      do j = 1, ao_num
        do i = 1, ao_num
          do igauss = 1, n_max_fit_slat
            delta = expo_gauss_1_erf_x_2(igauss)
            coef  = coef_gauss_1_erf_x_2(igauss)
            grad12_j12(i,j,ipoint) += -0.25d0 * coef * overlap_gauss_r12_ao(r, delta, i, j)
          enddo
        enddo
      enddo
    enddo

  endif

  call wall_time(time1)
  print*, ' Wall time for grad12_j12 = ', time1 - time0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, u12sq_j1bsq, (ao_num, ao_num, n_points_final_grid) ]

  implicit none
  integer                    :: ipoint, i, j
  double precision           :: tmp_x, tmp_y, tmp_z
  double precision           :: tmp1
  double precision           :: time0, time1

  print*, ' providing u12sq_j1bsq ...'
  call wall_time(time0)

  do ipoint = 1, n_points_final_grid
    tmp_x = v_1b_grad(1,ipoint)
    tmp_y = v_1b_grad(2,ipoint)
    tmp_z = v_1b_grad(3,ipoint)
    tmp1  = -0.5d0 * (tmp_x * tmp_x + tmp_y * tmp_y + tmp_z * tmp_z)
    do j = 1, ao_num
      do i = 1, ao_num
        u12sq_j1bsq(i,j,ipoint) = tmp1 * int2_u2_j1b2(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(time1)
  print*, ' Wall time for u12sq_j1bsq = ', time1 - time0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, u12_grad1_u12_j1b_grad1_j1b, (ao_num, ao_num, n_points_final_grid) ]
 
  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: x, y, z
  double precision           :: tmp_v, tmp_x, tmp_y, tmp_z
  double precision           :: tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing u12_grad1_u12_j1b_grad1_j1b ...'
  call wall_time(time0)

  do ipoint = 1, n_points_final_grid

    x     = final_grid_points(1,ipoint)
    y     = final_grid_points(2,ipoint)
    z     = final_grid_points(3,ipoint)
    tmp_v = v_1b       (ipoint)
    tmp_x = v_1b_grad(1,ipoint)
    tmp_y = v_1b_grad(2,ipoint)
    tmp_z = v_1b_grad(3,ipoint)

    tmp3 = tmp_v * tmp_x
    tmp4 = tmp_v * tmp_y
    tmp5 = tmp_v * tmp_z

    tmp6 = -x * tmp3
    tmp7 = -y * tmp4
    tmp8 = -z * tmp5

    do j = 1, ao_num
      do i = 1, ao_num

        tmp9 = int2_u_grad1u_j1b2(i,j,ipoint)

        u12_grad1_u12_j1b_grad1_j1b(i,j,ipoint) = tmp6 * tmp9 + tmp3 * int2_u_grad1u_x_j1b2(1,i,j,ipoint) &
                                                + tmp7 * tmp9 + tmp4 * int2_u_grad1u_x_j1b2(2,i,j,ipoint) &
                                                + tmp8 * tmp9 + tmp5 * int2_u_grad1u_x_j1b2(3,i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(time1)
  print*, ' Wall time for u12_grad1_u12_j1b_grad1_j1b = ', time1 - time0

END_PROVIDER 

! ---

