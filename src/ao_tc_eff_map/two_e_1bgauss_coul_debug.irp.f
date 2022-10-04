double precision function j1b_gauss_coul_debug(i, j, k, l)

  BEGIN_DOC
  ! 
  !  integral in the AO basis:
  !     i(r1) j(r1) f(r12) k(r2) l(r2)
  !
  !  with:
  !     f(r12) = - [ 0.5 / r12 ] (r1-r2) \cdot \sum_A (-2 a_A) [ r1A exp(-aA r1A^2) - r2A exp(-aA r2A^2) ]
  !            = [ 1 / r12 ] \sum_A a_A [ (r1-RA)^2 exp(-aA r1A^2)
  !                                     + (r2-RA)^2 exp(-aA r2A^2) 
  !                                     - (r1-RA) \cdot (r2-RA) exp(-aA r1A^2) 
  !                                     - (r1-RA) \cdot (r2-RA) exp(-aA r2A^2) ]
  !
  END_DOC

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

  double precision    :: general_primitive_integral_coul_shifted
  
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

  j1b_gauss_coul_debug = 0.d0


  ! -------------------------------------------------------------------------------------------------------------------
  !
  !                               [ 1 / r12 ] \sum_A a_A [ (r1-RA)^2 exp(-aA r1A^2)
  !
  ! -------------------------------------------------------------------------------------------------------------------

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

            ! pol centerd on P_center_tmp ==> centerd on P_center
            call pol_modif_center( P_center_tmp, P_center, iorder_p, P_new_tmp, P_new)

            ! ----------------------------------------------------------------------------------------------------
            ! x term:

            ff = P_center(1) - Centerii(1) 

            shift_P = (/ 2, 0, 0 /)
            cx = cx + expoii * general_primitive_integral_coul_shifted( dim1 &
                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P &
                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q )

            shift_P = (/ 1, 0, 0 /)
            cx = cx + expoii * 2.d0 * ff * general_primitive_integral_coul_shifted( dim1 &
                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P             &
                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q             )

            shift_P = (/ 0, 0, 0 /)
            cx = cx + expoii * ff * ff * general_primitive_integral_coul_shifted( dim1 &
                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P           &
                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q           )

            ! ----------------------------------------------------------------------------------------------------

          enddo

          j1b_gauss_coul_debug = j1b_gauss_coul_debug + coef4 * cx 
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------


!  ! -------------------------------------------------------------------------------------------------------------------
!  !
!  !                               [ 1 / r12 ] \sum_A a_A [ (r2-RA)^2 exp(-aA r2A^2)
!  !
!  ! -------------------------------------------------------------------------------------------------------------------
!
!  shift_P = (/ 0, 0, 0 /)
! 
!  do p = 1, ao_prim_num(i)
!    coef1 = ao_coef_normalized_ordered_transp(p, i)
!    expo1 = ao_expo_ordered_transp(p, i)
!
!    do q = 1, ao_prim_num(j)
!      coef2 = coef1 * ao_coef_normalized_ordered_transp(q, j)
!      expo2 = ao_expo_ordered_transp(q, j)
!
!      call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p, expo1, expo2 &
!                                          , I_power, J_power, I_center, J_center, dim1          )
!      p_inv = 1.d0 / pp
!
!      do r = 1, ao_prim_num(k)
!        coef3 = coef2 * ao_coef_normalized_ordered_transp(r, k)
!        expo3 = ao_expo_ordered_transp(r, k)
!
!        do s = 1, ao_prim_num(l)
!          coef4 = coef3 * ao_coef_normalized_ordered_transp(s, l)
!          expo4 = ao_expo_ordered_transp(s, l)
! 
!          call give_explicit_poly_and_gaussian( Q_new_tmp, Q_center_tmp, qq_tmp, fact_q_tmp, iorder_q, expo3, expo4 &
!                                              , K_power, L_power, K_center, L_center, dim1                          )
!
!          cx = 0.d0
!          do ii = 1, nucl_num
!            expoii        = j1b_gauss_pen(ii)
!            Centerii(1:3) = nucl_coord(ii, 1:3)
!
!            call gaussian_product(qq_tmp, Q_center_tmp, expoii, Centerii, factii, qq, Q_center)
!
!            fact_q = fact_q_tmp * factii
!            q_inv  = 1.d0 / qq
!
!            ! pol centerd on Q_center_tmp ==> centerd on Q_center
!            call pol_modif_center( Q_center_tmp, Q_center, iorder_q, Q_new_tmp, Q_new)
!
!            ! ----------------------------------------------------------------------------------------------------
!            ! x term:
!
!            ff = Q_center(1) - Centerii(1) 
!
!            shift_Q = (/ 2, 0, 0 /)
!            cx = cx + expoii * general_primitive_integral_coul_shifted( dim1 &
!                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P &
!                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q )
!
!            shift_Q = (/ 1, 0, 0 /)
!            cx = cx + expoii * 2.d0 * ff * general_primitive_integral_coul_shifted( dim1 &
!                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P             &
!                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q             )
!
!            shift_Q = (/ 0, 0, 0 /)
!            cx = cx + expoii * ff * ff * general_primitive_integral_coul_shifted( dim1 &
!                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P           &
!                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q           )
!
!            ! ----------------------------------------------------------------------------------------------------
!
!          enddo
!
!          j1b_gauss_coul_debug = j1b_gauss_coul_debug + coef4 * cx
!        enddo ! s
!      enddo  ! r
!    enddo   ! q
!  enddo    ! p
!
!  ! -------------------------------------------------------------------------------------------------------------------
!  ! -------------------------------------------------------------------------------------------------------------------


  ! -------------------------------------------------------------------------------------------------------------------
  !
  !                          - [ 1 / r12 ] \sum_A a_A [ (r1-RA) \cdot (r2-RA) exp(-aA r1A^2) ]
  !
  ! -------------------------------------------------------------------------------------------------------------------
 
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

            ! pol centerd on P_center_tmp ==> centerd on P_center
            call pol_modif_center( P_center_tmp, P_center, iorder_p, P_new_tmp, P_new)

            ! ----------------------------------------------------------------------------------------------------
            ! x term:

            ff = P_center(1) - Centerii(1) 
            gg = Q_center(1) - Centerii(1) 

            shift_P = (/ 1, 0, 0 /)
            shift_Q = (/ 1, 0, 0 /)
            cx = cx + expoii * general_primitive_integral_coul_shifted( dim1 &
                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P &
                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q )

            shift_P = (/ 1, 0, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cx = cx + expoii * gg * general_primitive_integral_coul_shifted( dim1 &
                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P      &
                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q      )

            shift_P = (/ 0, 0, 0 /)
            shift_Q = (/ 1, 0, 0 /)
            cx = cx + expoii * ff * general_primitive_integral_coul_shifted( dim1 &
                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P      &
                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q      )

            shift_P = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cx = cx + expoii * ff * gg * general_primitive_integral_coul_shifted( dim1 &
                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P           &
                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q           )

            ! ----------------------------------------------------------------------------------------------------

          enddo

          j1b_gauss_coul_debug = j1b_gauss_coul_debug - coef4 * cx
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------



!  ! -------------------------------------------------------------------------------------------------------------------
!  !
!  !                          - [ 1 / r12 ] \sum_A a_A [ (r1-RA) \cdot (r2-RA) exp(-aA r2A^2) ]
!  !
!  ! -------------------------------------------------------------------------------------------------------------------
! 
!  do p = 1, ao_prim_num(i)
!    coef1 = ao_coef_normalized_ordered_transp(p, i)
!    expo1 = ao_expo_ordered_transp(p, i)
!
!    do q = 1, ao_prim_num(j)
!      coef2 = coef1 * ao_coef_normalized_ordered_transp(q, j)
!      expo2 = ao_expo_ordered_transp(q, j)
!
!      call give_explicit_poly_and_gaussian( P_new, P_center, pp, fact_p, iorder_p, expo1, expo2 &
!                                          , I_power, J_power, I_center, J_center, dim1          )
!      p_inv = 1.d0 / pp
!
!      do r = 1, ao_prim_num(k)
!        coef3 = coef2 * ao_coef_normalized_ordered_transp(r, k)
!        expo3 = ao_expo_ordered_transp(r, k)
!
!        do s = 1, ao_prim_num(l)
!          coef4 = coef3 * ao_coef_normalized_ordered_transp(s, l)
!          expo4 = ao_expo_ordered_transp(s, l)
! 
!          call give_explicit_poly_and_gaussian( Q_new_tmp, Q_center_tmp, qq_tmp, fact_q_tmp, iorder_q, expo3, expo4 &
!                                              , K_power, L_power, K_center, L_center, dim1                          )
!
!          cx = 0.d0
!          do ii = 1, nucl_num
!            expoii        = j1b_gauss_pen(ii)
!            Centerii(1:3) = nucl_coord(ii, 1:3)
!
!            call gaussian_product(qq_tmp, Q_center_tmp, expoii, Centerii, factii, qq, Q_center)
!
!            fact_q = fact_q_tmp * factii
!            q_inv  = 1.d0 / qq
!
!            ! pol centerd on Q_center_tmp ==> centerd on Q_center
!            call pol_modif_center( Q_center_tmp, Q_center, iorder_q, Q_new_tmp, Q_new)
!
!            ! ----------------------------------------------------------------------------------------------------
!            ! x term:
!
!            ff = P_center(1) - Centerii(1) 
!            gg = Q_center(1) - Centerii(1) 
!
!            shift_P = (/ 1, 0, 0 /)
!            shift_Q = (/ 1, 0, 0 /)
!            cx = cx + expoii * general_primitive_integral_coul_shifted( dim1 &
!                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P &
!                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q )
!
!            shift_P = (/ 1, 0, 0 /)
!            shift_Q = (/ 0, 0, 0 /)
!            cx = cx + expoii * gg * general_primitive_integral_coul_shifted( dim1 &
!                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P      &
!                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q      )
!
!            shift_P = (/ 0, 0, 0 /)
!            shift_Q = (/ 1, 0, 0 /)
!            cx = cx + expoii * ff * general_primitive_integral_coul_shifted( dim1 &
!                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P      &
!                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q      )
!
!            shift_P = (/ 0, 0, 0 /)
!            shift_Q = (/ 0, 0, 0 /)
!            cx = cx + expoii * ff * gg * general_primitive_integral_coul_shifted( dim1 &
!                     , P_new, P_center, fact_p, pp, p_inv, iorder_p, shift_P           &
!                     , Q_new, Q_center, fact_q, qq, q_inv, iorder_q, shift_Q           )
!
!            ! ----------------------------------------------------------------------------------------------------
!
!          enddo
!
!          j1b_gauss_coul_debug = j1b_gauss_coul_debug - coef4 * cx
!
!        enddo ! s
!      enddo  ! r
!    enddo   ! q
!  enddo    ! p
!
!  ! -------------------------------------------------------------------------------------------------------------------
!  ! -------------------------------------------------------------------------------------------------------------------

  return
end function j1b_gauss_coul_debug

