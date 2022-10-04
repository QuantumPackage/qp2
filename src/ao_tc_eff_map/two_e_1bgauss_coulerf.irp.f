double precision function j1b_gauss_coulerf(i, j, k, l)

  BEGIN_DOC
  ! 
  !  integral in the AO basis:
  !     i(r1) j(r1) f(r12) k(r2) l(r2)
  !
  !  with:
  !     f(r12) = - [ (0.5 - 0.5 erf(mu r12)) / r12 ] (r1-r2) \cdot \sum_A (-2 a_A) [ r1A exp(-aA r1A^2) - r2A exp(-aA r2A^2) ]
  !            = [ (1 - erf(mu r12) / r12 ] \sum_A a_A [ (r1-RA)^2 exp(-aA r1A^2)
  !                                                    + (r2-RA)^2 exp(-aA r2A^2) 
  !                                                    - (r1-RA) \cdot (r2-RA) exp(-aA r1A^2)
  !                                                    - (r1-RA) \cdot (r2-RA) exp(-aA r2A^2) ]
  !
  END_DOC

  include 'utils/constants.include.F'

  implicit none

  integer, intent(in) :: i, j, k, l

  integer             :: p, q, r, s
  integer             :: num_i, num_j, num_k, num_l, num_ii 
  integer             :: I_power(3), J_power(3), K_power(3), L_power(3)
  integer             :: iorder_p(3), iorder_q(3)
  integer             :: shift_P(3), shift_Q(3)
  integer             :: dim1

  double precision    :: coef1, coef2, coef3, coef4
  double precision    :: expo1, expo2, expo3, expo4
  double precision    :: P1_new(0:max_dim,3), P1_center(3), fact_p1, pp1, p1_inv
  double precision    :: Q1_new(0:max_dim,3), Q1_center(3), fact_q1, qq1, q1_inv
  double precision    :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision    :: ff, gg, cx, cy, cz

  double precision    :: j1b_gauss_coulerf_schwartz
  
  PROVIDE j1b_gauss_pen

  dim1 = n_pt_max_integrals

  if( ao_prim_num(i) * ao_prim_num(j) * ao_prim_num(k) * ao_prim_num(l) > 1024 ) then
    j1b_gauss_coulerf = j1b_gauss_coulerf_schwartz(i, j, k, l)
    return
  endif

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

  j1b_gauss_coulerf = 0.d0

  do p = 1, ao_prim_num(i)
    coef1 = ao_coef_normalized_ordered_transp(p, i)
    expo1 = ao_expo_ordered_transp(p, i)

    do q = 1, ao_prim_num(j)
      coef2 = coef1 * ao_coef_normalized_ordered_transp(q, j)
      expo2 = ao_expo_ordered_transp(q, j)

      call give_explicit_poly_and_gaussian( P1_new, P1_center, pp1, fact_p1, iorder_p, expo1, expo2 &
                                          , I_power, J_power, I_center, J_center, dim1 )
      p1_inv = 1.d0 / pp1

      do r = 1, ao_prim_num(k)
        coef3 = coef2 * ao_coef_normalized_ordered_transp(r, k)
        expo3 = ao_expo_ordered_transp(r, k)

        do s = 1, ao_prim_num(l)
          coef4 = coef3 * ao_coef_normalized_ordered_transp(s, l)
          expo4 = ao_expo_ordered_transp(s, l)
 
          call give_explicit_poly_and_gaussian( Q1_new, Q1_center, qq1, fact_q1, iorder_q, expo3, expo4 &
                                              , K_power, L_power, K_center, L_center, dim1 )
          q1_inv = 1.d0 / qq1

          call get_cxcycz( dim1, cx, cy, cz                                  &
                         , P1_center, P1_new, pp1, fact_p1, p1_inv, iorder_p &
                         , Q1_center, Q1_new, qq1, fact_q1, q1_inv, iorder_q )

          j1b_gauss_coulerf = j1b_gauss_coulerf + coef4 * ( cx + cy + cz )
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

  return
end function j1b_gauss_coulerf

