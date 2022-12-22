! ---

subroutine overlap_gauss_xyz_r12_ao(D_center,delta,i,j,gauss_ints)

 implicit none
 BEGIN_DOC
! gauss_ints(m) = \int dr AO_i(r) AO_j(r) x/y/z e^{-delta |r-D_center|^2}
!
! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in)  :: D_center(3), delta
 double precision, intent(out) :: gauss_ints(3)

 integer :: num_a,num_b,power_A(3), power_B(3),l,k,m
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta,gauss_ints_tmp(3)
 gauss_ints = 0.d0
 if(ao_overlap_abs(j,i).lt.1.d-12)then
  return
 endif
 num_A = ao_nucl(i)
 power_A(1:3)= ao_power(i,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j)
 power_B(1:3)= ao_power(j,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 do l=1,ao_prim_num(i)
  alpha = ao_expo_ordered_transp(l,i)
  do k=1,ao_prim_num(j)
   beta = ao_expo_ordered_transp(k,j)
   call overlap_gauss_xyz_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,gauss_ints_tmp)
   do m = 1, 3
    gauss_ints(m) += gauss_ints_tmp(m) *  ao_coef_normalized_ordered_transp(l,i)             &
                                       *  ao_coef_normalized_ordered_transp(k,j)
   enddo
  enddo
 enddo

end



double precision function overlap_gauss_xyz_r12_ao_specific(D_center,delta,i,j,mx)
 implicit none
 BEGIN_DOC
! \int dr AO_i(r) AO_j(r) x/y/z e^{-delta |r-D_center|^2}
!
! with mx == 1 ==> x, mx == 2 ==> y, mx == 3 ==> z
 END_DOC
 integer, intent(in) :: i,j,mx
 double precision, intent(in)  :: D_center(3), delta

 integer :: num_a,num_b,power_A(3), power_B(3),l,k
 double precision :: gauss_int
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta
 double precision :: overlap_gauss_xyz_r12_specific
 overlap_gauss_xyz_r12_ao_specific = 0.d0
 if(ao_overlap_abs(j,i).lt.1.d-12)then
  return
 endif
 num_A = ao_nucl(i)
 power_A(1:3)= ao_power(i,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j)
 power_B(1:3)= ao_power(j,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 do l=1,ao_prim_num(i)
  alpha = ao_expo_ordered_transp(l,i)
  do k=1,ao_prim_num(j)
   beta = ao_expo_ordered_transp(k,j)
   gauss_int = overlap_gauss_xyz_r12_specific(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,mx)
   overlap_gauss_xyz_r12_ao_specific = gauss_int *  ao_coef_normalized_ordered_transp(l,i)             &
                                                 *  ao_coef_normalized_ordered_transp(k,j)
  enddo
 enddo
end


subroutine overlap_gauss_r12_all_ao(D_center,delta,aos_ints)
 implicit none
 double precision, intent(in) :: D_center(3), delta
 double precision, intent(out):: aos_ints(ao_num,ao_num)

 integer :: num_a,num_b,power_A(3), power_B(3),l,k,i,j
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta,analytical_j
 aos_ints = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   if(ao_overlap_abs(j,i).lt.1.d-12)cycle
   num_A = ao_nucl(i)
   power_A(1:3)= ao_power(i,1:3)
   A_center(1:3) = nucl_coord(num_A,1:3)
   num_B = ao_nucl(j)
   power_B(1:3)= ao_power(j,1:3)
   B_center(1:3) = nucl_coord(num_B,1:3)
   do l=1,ao_prim_num(i)
    alpha = ao_expo_ordered_transp(l,i)
    do k=1,ao_prim_num(j)
     beta = ao_expo_ordered_transp(k,j)
     analytical_j = overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
     aos_ints(j,i) += analytical_j *  ao_coef_normalized_ordered_transp(l,i)             &
                                   *  ao_coef_normalized_ordered_transp(k,j)
    enddo
   enddo
  enddo
 enddo
end

! ---

! TODO :: PUT CYCLES IN LOOPS
double precision function overlap_gauss_r12_ao(D_center, delta, i, j)

  BEGIN_DOC
  ! \int dr AO_i(r) AO_j(r) e^{-delta |r-D_center|^2}
  END_DOC

  implicit none
  integer,          intent(in) :: i, j
  double precision, intent(in) :: D_center(3), delta

  integer                      :: power_A(3), power_B(3), l, k
  double precision             :: A_center(3), B_center(3), alpha, beta, coef, coef1, analytical_j

  double precision, external   :: overlap_gauss_r12

  overlap_gauss_r12_ao = 0.d0

  if(ao_overlap_abs(j,i).lt.1.d-12) then
    return
  endif

  power_A(1:3) = ao_power(i,1:3)
  power_B(1:3) = ao_power(j,1:3)

  A_center(1:3) = nucl_coord(ao_nucl(i),1:3)
  B_center(1:3) = nucl_coord(ao_nucl(j),1:3)

  do l = 1, ao_prim_num(i)
    alpha = ao_expo_ordered_transp           (l,i)
    coef1 = ao_coef_normalized_ordered_transp(l,i)

    do k = 1, ao_prim_num(j)
      beta = ao_expo_ordered_transp(k,j)
      coef = coef1 * ao_coef_normalized_ordered_transp(k,j)

      if(dabs(coef) .lt. 1d-12) cycle

      analytical_j = overlap_gauss_r12(D_center, delta, A_center, B_center, power_A, power_B, alpha, beta)

      overlap_gauss_r12_ao += coef * analytical_j
    enddo
  enddo

end function overlap_gauss_r12_ao

! --

double precision function overlap_abs_gauss_r12_ao(D_center, delta, i, j)

  BEGIN_DOC
  ! \int dr AO_i(r) AO_j(r) e^{-delta |r-D_center|^2}
  END_DOC

  implicit none
  integer,          intent(in) :: i, j
  double precision, intent(in) :: D_center(3), delta

  integer                      :: power_A(3), power_B(3), l, k
  double precision             :: A_center(3), B_center(3), alpha, beta, coef, coef1, analytical_j

  double precision, external   :: overlap_abs_gauss_r12

  overlap_abs_gauss_r12_ao = 0.d0

  if(ao_overlap_abs(j,i).lt.1.d-12) then
    return
  endif

  power_A(1:3) = ao_power(i,1:3)
  power_B(1:3) = ao_power(j,1:3)

  A_center(1:3) = nucl_coord(ao_nucl(i),1:3)
  B_center(1:3) = nucl_coord(ao_nucl(j),1:3)

  do l = 1, ao_prim_num(i)
    alpha = ao_expo_ordered_transp           (l,i)
    coef1 = ao_coef_normalized_ordered_transp(l,i)

    do k = 1, ao_prim_num(j)
      beta = ao_expo_ordered_transp(k,j)
      coef = coef1 * ao_coef_normalized_ordered_transp(k,j)

      if(dabs(coef) .lt. 1d-12) cycle

      analytical_j = overlap_abs_gauss_r12(D_center, delta, A_center, B_center, power_A, power_B, alpha, beta)

      overlap_abs_gauss_r12_ao += dabs(coef * analytical_j)
    enddo
  enddo

end function overlap_gauss_r12_ao

! --

subroutine overlap_gauss_r12_ao_v(D_center, LD_D, delta, i, j, resv, LD_resv, n_points)

  BEGIN_DOC
  !
  ! \int dr AO_i(r) AO_j(r) e^{-delta |r-D_center|^2}
  !
  ! n_points: nb of integrals <= min(LD_D, LD_resv)
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: i, j, LD_D, LD_resv, n_points
  double precision, intent(in)  :: D_center(LD_D,3), delta
  double precision, intent(out) :: resv(LD_resv)

  integer                       :: ipoint
  integer                       :: power_A(3), power_B(3), l, k
  double precision              :: A_center(3), B_center(3), alpha, beta, coef, coef1 
  double precision, allocatable :: analytical_j(:)

  resv(:) = 0.d0
  if(ao_overlap_abs(j,i) .lt. 1.d-12) then
    return
  endif

  power_A(1:3) = ao_power(i,1:3)
  power_B(1:3) = ao_power(j,1:3)

  A_center(1:3) = nucl_coord(ao_nucl(i),1:3)
  B_center(1:3) = nucl_coord(ao_nucl(j),1:3)

  allocate(analytical_j(n_points))

  do l = 1, ao_prim_num(i)
    alpha = ao_expo_ordered_transp           (l,i)
    coef1 = ao_coef_normalized_ordered_transp(l,i)

    do k = 1, ao_prim_num(j)
      beta = ao_expo_ordered_transp(k,j)
      coef = coef1 * ao_coef_normalized_ordered_transp(k,j)

      if(dabs(coef) .lt. 1d-12) cycle

      call overlap_gauss_r12_v(D_center, LD_D, delta, A_center, B_center, power_A, power_B, alpha, beta, analytical_j, n_points, n_points)

      do ipoint = 1, n_points
        resv(ipoint) = resv(ipoint) + coef * analytical_j(ipoint)
      enddo

    enddo
  enddo

  deallocate(analytical_j)

end subroutine overlap_gauss_r12_ao_v

! ---

double precision function overlap_gauss_r12_ao_with1s(B_center, beta, D_center, delta, i, j)

  BEGIN_DOC
  !
  ! \int dr AO_i(r) AO_j(r) e^{-beta |r-B_center^2|} e^{-delta |r-D_center|^2}
  !
  END_DOC

  implicit none
  integer,          intent(in) :: i, j
  double precision, intent(in) :: B_center(3), beta, D_center(3), delta

  integer                      :: power_A1(3), power_A2(3), l, k
  double precision             :: A1_center(3), A2_center(3), alpha1, alpha2, coef1, coef12, analytical_j
  double precision             :: G_center(3), gama, fact_g, gama_inv

  double precision, external   :: overlap_gauss_r12, overlap_gauss_r12_ao

  if(beta .lt. 1d-10) then
    overlap_gauss_r12_ao_with1s = overlap_gauss_r12_ao(D_center, delta, i, j)
    return
  endif

  overlap_gauss_r12_ao_with1s = 0.d0

  if(ao_overlap_abs(j,i) .lt. 1.d-12) then
    return
  endif

  ! e^{-beta |r-B_center^2|} e^{-delta |r-D_center|^2} = fact_g e^{-gama |r - G|^2}

  gama        = beta + delta
  gama_inv    = 1.d0 / gama
  G_center(1) = (beta * B_center(1) + delta * D_center(1)) * gama_inv
  G_center(2) = (beta * B_center(2) + delta * D_center(2)) * gama_inv
  G_center(3) = (beta * B_center(3) + delta * D_center(3)) * gama_inv
  fact_g      = beta * delta * gama_inv * ( (B_center(1) - D_center(1)) * (B_center(1) - D_center(1)) &
                                          + (B_center(2) - D_center(2)) * (B_center(2) - D_center(2)) &
                                          + (B_center(3) - D_center(3)) * (B_center(3) - D_center(3)) )
  if(fact_g .gt. 10d0) return
  fact_g = dexp(-fact_g)

  ! ---

  power_A1(1:3) = ao_power(i,1:3)
  power_A2(1:3) = ao_power(j,1:3)

  A1_center(1:3) = nucl_coord(ao_nucl(i),1:3)
  A2_center(1:3) = nucl_coord(ao_nucl(j),1:3)

  do l = 1, ao_prim_num(i)
    alpha1 = ao_expo_ordered_transp                    (l,i)
    coef1  = fact_g * ao_coef_normalized_ordered_transp(l,i)
    if(dabs(coef1) .lt. 1d-12) cycle

    do k = 1, ao_prim_num(j)
      alpha2 = ao_expo_ordered_transp                   (k,j)
      coef12 = coef1 * ao_coef_normalized_ordered_transp(k,j)
      if(dabs(coef12) .lt. 1d-12) cycle

      analytical_j = overlap_gauss_r12(G_center, gama, A1_center, A2_center, power_A1, power_A2, alpha1, alpha2)

      overlap_gauss_r12_ao_with1s += coef12 * analytical_j
    enddo
  enddo

end function overlap_gauss_r12_ao_with1s

! ---

subroutine overlap_gauss_r12_ao_with1s_v(B_center, beta, D_center, LD_D, delta, i, j, resv, LD_resv, n_points)

  BEGIN_DOC
  !
  ! \int dr AO_i(r) AO_j(r) e^{-beta |r-B_center^2|} e^{-delta |r-D_center|^2}
  ! using an array of D_centers.
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: i, j, n_points, LD_D, LD_resv
  double precision, intent(in)  :: B_center(3), beta, D_center(LD_D,3), delta
  double precision, intent(out) :: resv(LD_resv)

  integer                       :: ipoint
  integer                       :: power_A1(3), power_A2(3), l, k
  double precision              :: A1_center(3), A2_center(3), alpha1, alpha2, coef1
  double precision              :: coef12, coef12f
  double precision              :: gama, gama_inv
  double precision              :: bg, dg, bdg
  double precision, allocatable :: fact_g(:), G_center(:,:), analytical_j(:)

  if(ao_overlap_abs(j,i) .lt. 1.d-12) then
    return
  endif

  ASSERT(beta .gt. 0.d0)

  if(beta .lt. 1d-10) then
    call overlap_gauss_r12_ao_v(D_center, LD_D, delta, i, j, resv, LD_resv, n_points)
    return
  endif

  resv(:) = 0.d0

  ! e^{-beta |r-B_center^2|} e^{-delta |r-D_center|^2} = fact_g e^{-gama |r - G|^2}

  gama     = beta + delta
  gama_inv = 1.d0 / gama

  power_A1(1:3) = ao_power(i,1:3)
  power_A2(1:3) = ao_power(j,1:3)

  A1_center(1:3) = nucl_coord(ao_nucl(i),1:3)
  A2_center(1:3) = nucl_coord(ao_nucl(j),1:3)

  allocate(fact_g(n_points), G_center(n_points,3), analytical_j(n_points))

  bg  = beta  * gama_inv
  dg  = delta * gama_inv
  bdg = bg * delta 

  do ipoint = 1, n_points

    G_center(ipoint,1) = bg * B_center(1) + dg * D_center(ipoint,1)
    G_center(ipoint,2) = bg * B_center(2) + dg * D_center(ipoint,2)
    G_center(ipoint,3) = bg * B_center(3) + dg * D_center(ipoint,3)
    fact_g(ipoint) = bdg * ( (B_center(1) - D_center(ipoint,1)) * (B_center(1) - D_center(ipoint,1)) &
                           + (B_center(2) - D_center(ipoint,2)) * (B_center(2) - D_center(ipoint,2)) &
                           + (B_center(3) - D_center(ipoint,3)) * (B_center(3) - D_center(ipoint,3)) )

    if(fact_g(ipoint) < 10d0) then
      fact_g(ipoint) = dexp(-fact_g(ipoint))
    else
      fact_g(ipoint) = 0.d0
    endif

  enddo

  do l = 1, ao_prim_num(i)
    alpha1 = ao_expo_ordered_transp           (l,i)
    coef1  = ao_coef_normalized_ordered_transp(l,i)

    do k = 1, ao_prim_num(j)
      alpha2 = ao_expo_ordered_transp                   (k,j)
      coef12 = coef1 * ao_coef_normalized_ordered_transp(k,j)
      if(dabs(coef12) .lt. 1d-12) cycle

      call overlap_gauss_r12_v(G_center, n_points, gama, A1_center, A2_center, power_A1, power_A2, alpha1, alpha2, analytical_j, n_points, n_points)

      do ipoint = 1, n_points
        coef12f = coef12 * fact_g(ipoint)
        resv(ipoint) += coef12f * analytical_j(ipoint)
      enddo
    enddo
  enddo

  deallocate(fact_g, G_center, analytical_j)

end subroutine overlap_gauss_r12_ao_with1s_v

! ---

