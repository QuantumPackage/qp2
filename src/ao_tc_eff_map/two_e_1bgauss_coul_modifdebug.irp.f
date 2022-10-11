double precision function j1b_gauss_coul_modifdebug(i, j, k, l)

  include 'utils/constants.include.F'

  implicit none

  integer, intent(in) :: i, j, k, l

  integer             :: p, q, r, s, ii
  integer             :: num_i, num_j, num_k, num_l, num_ii 
  integer             :: I_power(3), J_power(3), K_power(3), L_power(3)
  integer             :: iorder_p(3), iorder_q(3)
  integer             :: shift_P(3), shift_Q(3)
  integer             :: dim1

  double precision    :: coef1, coef2, coef3, coef4
  double precision    :: expo1, expo2, expo3, expo4
  double precision    :: p_inv, q_inv
  double precision    :: P_new_tmp(0:max_dim,3), P_center_tmp(3), fact_p_tmp, pp_tmp
  double precision    :: Q_new_tmp(0:max_dim,3), Q_center_tmp(3), fact_q_tmp, qq_tmp
  double precision    :: P_new(0:max_dim,3), P_center(3), fact_p, pp
  double precision    :: Q_new(0:max_dim,3), Q_center(3), fact_q, qq
  double precision    :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision    :: expoii, factii, Centerii(3)
  double precision    :: ff, gg, cx, cy, cz

  double precision    :: general_primitive_integral_coul
  double precision    :: general_primitive_integral_coul_shifted
  double precision    :: ao_two_e_integral
  
  PROVIDE j1b_gauss_pen

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)

  do p = 1, 3
    I_power(p)  = ao_power(i,p)
    J_power(p)  = ao_power(j,p)
    K_power(p)  = ao_power(k,p)
    L_power(p)  = ao_power(l,p)
    I_center(p) = nucl_coord(num_i,p)
    J_center(p) = nucl_coord(num_j,p)
    K_center(p) = nucl_coord(num_k,p)
    L_center(p) = nucl_coord(num_l,p)
  enddo

  j1b_gauss_coul_modifdebug = 0.d0

!  do ii = 1, nucl_num
!    expoii = j1b_gauss_pen(ii)
!    j1b_gauss_coul_modifdebug += expoii * ao_two_e_integral(i, j, k, l)
!  enddo


  ! -------------------------------------------------------------------------------------------------------------------
  !
  !                               [ 1 / r12 ] \sum_A a_A exp(-aA r1A^2)
  !
  ! -------------------------------------------------------------------------------------------------------------------

  shift_P = (/ 0, 0, 0 /)
  shift_Q = (/ 0, 0, 0 /)
 
  do p = 1, ao_prim_num(i)
    coef1 = ao_coef_normalized_ordered_transp(p, i)
    expo1 = ao_expo_ordered_transp(p, i)

    do q = 1, ao_prim_num(j)
      coef2 = coef1 * ao_coef_normalized_ordered_transp(q, j)
      expo2 = ao_expo_ordered_transp(q, j)

      call give_explicit_poly_and_gaussian( P_new_tmp, P_center_tmp, pp_tmp, fact_p_tmp, iorder_p, expo1, expo2 &
                                          , I_power, J_power, I_center, J_center, dim1                          )

      do r = 1, ao_prim_num(k)
        coef3 = coef2 * ao_coef_normalized_ordered_transp(r, k)
        expo3 = ao_expo_ordered_transp(r, k)

        do s = 1, ao_prim_num(l)
          coef4 = coef3 * ao_coef_normalized_ordered_transp(s, l)
          expo4 = ao_expo_ordered_transp(s, l)
 
          call give_explicit_poly_and_gaussian( Q_new, Q_center, qq, fact_q, iorder_q, expo3, expo4 &
                                              , K_power, L_power, K_center, L_center, dim1          )
          q_inv = 1.d0 / qq

          cx = 0.d0
          do ii = 1, nucl_num
            expoii        = j1b_gauss_pen(ii)
            Centerii(1:3) = nucl_coord(ii, 1:3)

            call gaussian_product(pp_tmp, P_center_tmp, expoii, Centerii, factii, pp, P_center)
            fact_p = fact_p_tmp * factii
            p_inv  = 1.d0 / pp
            P_new(:,:) = 0.d0
            call pol_modif_center( P_center_tmp, P_center, iorder_p, P_new_tmp, P_new)

            ! ----------------------------------------------------------------------------------------------------
            ! x term:

            cx = cx + expoii * general_primitive_integral_coul_shifted( dim1 &
                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P &
                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q )

            ! ----------------------------------------------------------------------------------------------------

          enddo

          j1b_gauss_coul_modifdebug = j1b_gauss_coul_modifdebug + coef4 * cx 
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------


  ! -------------------------------------------------------------------------------------------------------------------
  !
  !                               [ 1 / r12 ] \sum_A a_A exp(-aA r2A^2)
  !
  ! -------------------------------------------------------------------------------------------------------------------

  shift_P = (/ 0, 0, 0 /)
  shift_Q = (/ 0, 0, 0 /)
 
  do p = 1, ao_prim_num(i)
    coef1 = ao_coef_normalized_ordered_transp(p, i)
    expo1 = ao_expo_ordered_transp(p, i)

    do q = 1, ao_prim_num(j)
      coef2 = coef1 * ao_coef_normalized_ordered_transp(q, j)
      expo2 = ao_expo_ordered_transp(q, j)

      call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p, expo1, expo2 &
                                          , I_power, J_power, I_center, J_center, dim1 )
      p_inv  = 1.d0 / pp

      do r = 1, ao_prim_num(k)
        coef3 = coef2 * ao_coef_normalized_ordered_transp(r, k)
        expo3 = ao_expo_ordered_transp(r, k)

        do s = 1, ao_prim_num(l)
          coef4 = coef3 * ao_coef_normalized_ordered_transp(s, l)
          expo4 = ao_expo_ordered_transp(s, l)
 
          call give_explicit_poly_and_gaussian( Q_new_tmp, Q_center_tmp, qq_tmp, fact_q_tmp, iorder_q, expo3, expo4 &
                                              , K_power, L_power, K_center, L_center, dim1 )

          cx = 0.d0
          do ii = 1, nucl_num
            expoii        = j1b_gauss_pen(ii)
            Centerii(1:3) = nucl_coord(ii, 1:3)

            call gaussian_product(qq_tmp, Q_center_tmp, expoii, Centerii, factii, qq, Q_center)
            fact_q = fact_q_tmp * factii
            Q_inv  = 1.d0 / qq
            call pol_modif_center( Q_center_tmp, Q_center, iorder_q, Q_new_tmp, Q_new)

            ! ----------------------------------------------------------------------------------------------------
            ! x term:

            cx = cx + expoii * general_primitive_integral_coul_shifted( dim1 &
                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P &
                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q )

            ! ----------------------------------------------------------------------------------------------------

          enddo

          j1b_gauss_coul_modifdebug = j1b_gauss_coul_modifdebug + coef4 * cx 
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  return
end function j1b_gauss_coul_modifdebug




double precision function general_primitive_integral_coul(dim,            &
      P_new,P_center,fact_p,p,p_inv,iorder_p,                        &
      Q_new,Q_center,fact_q,q,q_inv,iorder_q)
  implicit none
  BEGIN_DOC
  ! Computes the integral <pq|rs> where p,q,r,s are Gaussian primitives
  END_DOC
  integer,intent(in)             :: dim
  include 'utils/constants.include.F'
  double precision, intent(in)   :: P_new(0:max_dim,3),P_center(3),fact_p,p,p_inv
  double precision, intent(in)   :: Q_new(0:max_dim,3),Q_center(3),fact_q,q,q_inv
  integer, intent(in)            :: iorder_p(3)
  integer, intent(in)            :: iorder_q(3)

  double precision               :: r_cut,gama_r_cut,rho,dist
  double precision               :: dx(0:max_dim),Ix_pol(0:max_dim),dy(0:max_dim),Iy_pol(0:max_dim),dz(0:max_dim),Iz_pol(0:max_dim)
  integer                        :: n_Ix,n_Iy,n_Iz,nx,ny,nz
  double precision               :: bla
  integer                        :: ix,iy,iz,jx,jy,jz,i
  double precision               :: a,b,c,d,e,f,accu,pq,const
  double precision               :: pq_inv, p10_1, p10_2, p01_1, p01_2,pq_inv_2
  integer                        :: n_pt_tmp,n_pt_out, iorder
  double precision               :: d1(0:max_dim),d_poly(0:max_dim),rint,d1_screened(0:max_dim)

  general_primitive_integral_coul = 0.d0

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
      call give_polynom_mult_center_x(P_center(1),Q_center(1),ix,jx,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dx,nx)
      !DIR$ FORCEINLINE
      call add_poly_multiply(dx,nx,d,Ix_pol,n_Ix)
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
        call   give_polynom_mult_center_x(P_center(2),Q_center(2),iy,jy,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dy,ny)
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
        call   give_polynom_mult_center_x(P_center(3),Q_center(3),iz,jz,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dz,nz)
        !DIR$ FORCEINLINE
        call add_poly_multiply(dz,nz,f,Iz_pol,n_Iz)
      enddo
    endif
  enddo
  if (n_Iz == -1) then
    return
  endif

  rho = p*q *pq_inv_2
  dist =  (P_center(1) - Q_center(1))*(P_center(1) - Q_center(1)) +  &
      (P_center(2) - Q_center(2))*(P_center(2) - Q_center(2)) +      &
      (P_center(3) - Q_center(3))*(P_center(3) - Q_center(3))
  const = dist*rho

  n_pt_tmp = n_Ix+n_Iy
  do i=0,n_pt_tmp
    d_poly(i)=0.d0
  enddo

  !DIR$ FORCEINLINE
  call multiply_poly(Ix_pol,n_Ix,Iy_pol,n_Iy,d_poly,n_pt_tmp)
  if (n_pt_tmp == -1) then
    return
  endif
  n_pt_out = n_pt_tmp+n_Iz
  do i=0,n_pt_out
    d1(i)=0.d0
  enddo

  !DIR$ FORCEINLINE
  call multiply_poly(d_poly ,n_pt_tmp ,Iz_pol,n_Iz,d1,n_pt_out)
  double precision               :: rint_sum
  accu = accu + rint_sum(n_pt_out,const,d1)

  general_primitive_integral_coul = fact_p * fact_q * accu *pi_5_2*p_inv*q_inv/dsqrt(p+q)
end function general_primitive_integral_coul
