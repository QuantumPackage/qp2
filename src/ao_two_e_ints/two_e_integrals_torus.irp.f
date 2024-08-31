
! ---

double precision function ao_two_e_integral_torus(i, j, k, l)

  BEGIN_DOC
  !
  ! TODO
  !
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  !
  END_DOC

  implicit none

  include 'utils/constants.include.F'

  integer, intent(in)        :: i, j, k, l
                             
  integer                    :: p, q, r, s
  integer                    :: num_i, num_j, num_k, num_l, dim1
  integer                    :: I_power(3), J_power(3), K_power(3), L_power(3)
  integer                    :: iorder_p(3), iorder_q(3)
  double precision           :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision           :: integral
  double precision           :: P_new(0:max_dim,3), P_center(3), fact_p, pp
  double precision           :: Q_new(0:max_dim,3), Q_center(3), fact_q, qq
  double precision           :: coef1, coef2, coef3, coef4
  double precision           :: p_inv, q_inv

  double precision, external :: ERI
  double precision, external :: general_primitive_integral_torus


  PROVIDE torus_length
  PROVIDE n_pt_max_integrals

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)
  ao_two_e_integral_torus = 0.d0

  if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k) then
    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
      I_center(p) = nucl_coord(num_i,p)
      J_center(p) = nucl_coord(num_j,p)
      K_center(p) = nucl_coord(num_k,p)
      L_center(p) = nucl_coord(num_l,p)
    enddo

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)

      do q = 1, ao_prim_num(j)
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)

        call give_explicit_poly_and_gaussian_torus(P_new, P_center, pp, fact_p, iorder_p, &
                                                   ao_expo_ordered_transp(p,i),           &
                                                   ao_expo_ordered_transp(q,j),           &
                                                   I_power, J_power, I_center, J_center, torus_length, dim1)


        p_inv = 1.d0 / pp
        do r = 1, ao_prim_num(k)
          coef3 = coef2 * ao_coef_normalized_ordered_transp(r,k)

          do s = 1, ao_prim_num(l)
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)

            call give_explicit_poly_and_gaussian_torus(Q_new, Q_center, qq, fact_q, iorder_q, &
                                                       ao_expo_ordered_transp(r,k),           &
                                                       ao_expo_ordered_transp(s,l),           &
                                                       K_power, L_power, K_center, L_center, torus_length, dim1)

            q_inv = 1.d0 / qq

            integral = general_primitive_integral_torus(torus_length, dim1, &
                P_new, P_center, fact_p, pp, p_inv, iorder_p,               &
                Q_new, Q_center, fact_q, qq, q_inv, iorder_q)

            ao_two_e_integral_torus += coef4 * integral
          enddo ! s
        enddo ! r
      enddo ! q
    enddo ! p

  else

    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
    enddo

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)

      do q = 1, ao_prim_num(j)
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)

        do r = 1, ao_prim_num(k)
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)

          do s = 1, ao_prim_num(l)
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)

            integral = ERI(ao_expo_ordered_transp(p,i), ao_expo_ordered_transp(q,j), &
                           ao_expo_ordered_transp(r,k), ao_expo_ordered_transp(s,l), &
                           I_power(1), J_power(1), K_power(1), L_power(1),           &
                           I_power(2), J_power(2), K_power(2), L_power(2),           &
                           I_power(3), J_power(3), K_power(3), L_power(3))

            ao_two_e_integral_torus = ao_two_e_integral_torus + coef4 * integral
          enddo ! s
        enddo ! r
      enddo ! q
    enddo ! p

  endif

end

! ---

double precision function general_primitive_integral_torus(torus_L, dim,   &
                              P_new, P_center, fact_p, p, p_inv, iorder_p, &
                              Q_new, Q_center, fact_q, q, q_inv, iorder_q)

  implicit none

  include 'utils/constants.include.F'

  BEGIN_DOC
  !
  ! TODO
  !
  ! Computes the integral <pq|rs> where p,q,r,s are Gaussian primitives
  !
  END_DOC

  integer,          intent(in) :: dim
  integer,          intent(in) :: iorder_p(3), iorder_q(3)
  double precision, intent(in) :: torus_L(3)
  double precision, intent(in) :: P_new(0:max_dim,3), P_center(3), fact_p, p, p_inv
  double precision, intent(in) :: Q_new(0:max_dim,3), Q_center(3), fact_q, q, q_inv

  integer                      :: ix,iy,iz,jx,jy,jz,i
  integer                      :: n_Ix,n_Iy,n_Iz,nx,ny,nz
  integer                      :: n_pt_tmp,n_pt_out, iorder
  integer                      :: ib, ic
  double precision             :: rho, dist
  double precision             :: dx(0:max_dim), Ix_pol(0:max_dim)
  double precision             :: dy(0:max_dim), Iy_pol(0:max_dim)
  double precision             :: dz(0:max_dim), Iz_pol(0:max_dim)
  double precision             :: a, b, c, d, e, f, accu, pq, const
  double precision             :: pq_inv, p10_1, p10_2, p01_1, p01_2,pq_inv_2
  double precision             :: d1(0:max_dim),d_poly(0:max_dim),d1_screened(0:max_dim)
  double precision             :: dist_tmp_x, dist_tmp_y, dist_tmp_z


  double precision, external   :: rint_sum

  general_primitive_integral_torus = 0.d0

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: dx,Ix_pol,dy,Iy_pol,dz,Iz_pol
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: d1, d_poly

  ! Gaussian Product
  ! ----------------

  pq = p_inv*0.5d0*q_inv
  pq_inv = 0.5d0/(p+q)
  p10_1 = q*pq  ! 1/(2p)
  p01_1 = p*pq  ! 1/(2q)
  pq_inv_2 = pq_inv+pq_inv
  p10_2 = pq_inv_2 * p10_1*q !0.5d0*q/(pq + p*p)
  p01_2 = pq_inv_2 * p01_1*p !0.5d0*p/(q*q + pq)


  accu = 0.d0
  iorder = iorder_p(1)+iorder_q(1)+iorder_p(1)+iorder_q(1)
  do ix=0,iorder
    Ix_pol(ix) = 0.d0
  enddo
  n_Ix = 0
  do ix = 0, iorder_p(1)
    if (abs(P_new(ix,1)) < thresh) cycle
    a = P_new(ix,1)
    do jx = 0, iorder_q(1)
      d = a*Q_new(jx,1)
      if (abs(d) < thresh) cycle
      !DIR$ FORCEINLINE
      call give_polynom_mult_center_x_torus(torus_L(1), P_center(1), Q_center(1), ix, jx, p, q, iorder, &
                                            pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dx, nx)
      !DIR$ FORCEINLINE
      call add_poly_multiply(dx, nx, d, Ix_pol, n_Ix)
    enddo
  enddo
  if (n_Ix == -1) then
    return
  endif
  iorder = iorder_p(2)+iorder_q(2)+iorder_p(2)+iorder_q(2)
  do ix=0, iorder
    Iy_pol(ix) = 0.d0
  enddo
  n_Iy = 0
  do iy = 0, iorder_p(2)
    if (abs(P_new(iy,2)) > thresh) then
      b = P_new(iy,2)
      do jy = 0, iorder_q(2)
        e = b*Q_new(jy,2)
        if (abs(e) < thresh) cycle
        !DIR$ FORCEINLINE
        call give_polynom_mult_center_x_torus(torus_L(2), P_center(2), Q_center(2), iy, jy, p, q, iorder, &
                                              pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dy, ny)
        !DIR$ FORCEINLINE
        call add_poly_multiply(dy,ny,e,Iy_pol,n_Iy)
      enddo
    endif
  enddo
  if (n_Iy == -1) then
    return
  endif

  iorder = iorder_p(3)+iorder_q(3)+iorder_p(3)+iorder_q(3)
  do ix=0,iorder
    Iz_pol(ix) = 0.d0
  enddo
  n_Iz = 0
  do iz = 0, iorder_p(3)
    if (abs(P_new(iz,3)) > thresh) then
      c = P_new(iz,3)
      do jz = 0, iorder_q(3)
        f = c*Q_new(jz,3)
        if (abs(f) < thresh) cycle
        !DIR$ FORCEINLINE
        call give_polynom_mult_center_x_torus(torus_L(3), P_center(3), Q_center(3), iz, jz, p, q, iorder, &
                                              pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dz, nz)
        !DIR$ FORCEINLINE
        call add_poly_multiply(dz,nz,f,Iz_pol,n_Iz)
      enddo
    endif
  enddo
  if (n_Iz == -1) then
    return
  endif

  rho = p * q * pq_inv_2

  ! old
  !dist = (P_center(1) - Q_center(1))*(P_center(1) - Q_center(1)) +  &
  !       (P_center(2) - Q_center(2))*(P_center(2) - Q_center(2)) +  &
  !       (P_center(3) - Q_center(3))*(P_center(3) - Q_center(3))
  ! new
  call ssd_euc_torus(P_center(1), Q_center(1), torus_L(1), dist_tmp_x)
  call ssd_euc_torus(P_center(2), Q_center(2), torus_L(2), dist_tmp_y)
  call ssd_euc_torus(P_center(3), Q_center(3), torus_L(3), dist_tmp_z)
  dist = dist_tmp_x * dist_tmp_x + dist_tmp_y * dist_tmp_y + dist_tmp_z * dist_tmp_z

  const = dist * rho

  n_pt_tmp = n_Ix+n_Iy
  do i=0,n_pt_tmp
    d_poly(i)=0.d0
  enddo

  if (ior(n_Ix,n_Iy) >= 0) then
    do ib=0,n_Ix
      do ic = 0,n_Iy
        d_poly(ib+ic) = d_poly(ib+ic) + Iy_pol(ic) * Ix_pol(ib)
      enddo
    enddo

    do n_pt_tmp = n_Ix+n_Iy, 0, -1
      if (d_poly(n_pt_tmp) /= 0.d0) exit
    enddo
  endif

  if (n_pt_tmp == -1) then
    return
  endif
  n_pt_out = n_pt_tmp+n_Iz
  do i=0,n_pt_out
    d1(i)=0.d0
  enddo

  if (ior(n_pt_tmp,n_Iz) >= 0) then
    ! Bottleneck here
    do ib=0,n_pt_tmp
      do ic = 0,n_Iz
        d1(ib+ic) = d1(ib+ic) + Iz_pol(ic) * d_poly(ib)
      enddo
    enddo

    do n_pt_out = n_pt_tmp+n_Iz, 0, -1
      if (d1(n_pt_out) /= 0.d0) exit
    enddo
  endif

  accu = accu + rint_sum(n_pt_out, const, d1)

  general_primitive_integral_torus = fact_p * fact_q * accu *pi_5_2*p_inv*q_inv/dsqrt(p+q)

end

! ---

subroutine give_polynom_mult_center_x_torus(Lx, P_center, Q_center, a_x, d_x, p, q, n_pt_in, pq_inv, pq_inv_2, &
                                            p10_1, p01_1, p10_2, p01_2, d, n_pt_out)

  implicit none

  include 'utils/constants.include.F'

  BEGIN_DOC
  ! subroutine that returns the explicit polynom in term of the "t"
  ! variable of the following polynomw :
  !
  ! $I_{x_1}(a_x,d_x,p,q) \, I_{x_1}(a_y,d_y,p,q) \ I_{x_1}(a_z,d_z,p,q)$
  END_DOC

  integer,          intent(in)  :: n_pt_in
  integer,          intent(in)  :: a_x, d_x
  double precision, intent(in)  :: Lx
  double precision, intent(in)  :: P_center, Q_center
  double precision, intent(in)  :: p, q, pq_inv, p10_1, p01_1, p10_2, p01_2, pq_inv_2

  integer,          intent(out) :: n_pt_out
  double precision, intent(out) :: d(0:max_dim)

  integer                       :: n_pt1, dim, i
  double precision              :: B10(0:2), B01(0:2), B00(0:2),C00(0:2),D00(0:2)

  double precision              :: accu, tmp



  accu = 0.d0


  B10(0) = p10_1
  B10(1) = 0.d0
  B10(2) = - p10_2

  B01(0) = p01_1
  B01(1) = 0.d0
  B01(2) = - p01_2

  B00(0) = 0.d0
  B00(1) = 0.d0
  B00(2) = pq_inv

  do i = 0, n_pt_in
    d(i) = 0.d0
  enddo

  n_pt1 = n_pt_in

  ! ---

  C00(0) = 0.d0
  C00(1) = 0.d0
  

  ! old
  !C00(2) = -q*(P_center-Q_center) * pq_inv_2

  ! torus
  !call ssd_torus(P_center, Q_center, Lx, tmp)
  call ssd_euc_torus(P_center, Q_center, Lx, tmp)
  C00(2) = -q * tmp * pq_inv_2

  ! ---

  D00(0) = 0.d0
  D00(1) = 0.d0


  ! old
  !D00(2) = -p*(Q_center-P_center) * pq_inv_2

  ! torus
  !call ssd_torus(Q_center, P_center, Lx, tmp)
  call ssd_euc_torus(Q_center, P_center, Lx, tmp)
  D00(2) = -p * tmp * pq_inv_2

  ! ---

  !DIR$ FORCEINLINE
  call I_x1_pol_mult(a_x, d_x, B10, B01, B00, C00, D00, d, n_pt1, n_pt_in)

  n_pt_out = n_pt1
  if(n_pt1 < 0) then
    n_pt_out = -1
    do i = 0, n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif

end

