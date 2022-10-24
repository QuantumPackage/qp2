
! ---

double precision function ao_two_e_integral_cosgtos(i, j, k, l)

  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in) :: i, j, k, l

  integer             :: p, q, r, s
  integer             :: num_i, num_j, num_k, num_l, dim1, I_power(3), J_power(3), K_power(3), L_power(3)
  integer             :: iorder_p1(3), iorder_p2(3), iorder_p3(3), iorder_p4(3), iorder_q1(3), iorder_q2(3)
  double precision    :: coef1, coef2, coef3, coef4
  complex*16          :: I_center(3), J_center(3), K_center(3), L_center(3)
  complex*16          :: expo1, expo2, expo3, expo4
  complex*16          :: P1_new(0:max_dim,3), P1_center(3), fact_p1, pp1, p1_inv
  complex*16          :: P2_new(0:max_dim,3), P2_center(3), fact_p2, pp2, p2_inv
  complex*16          :: P3_new(0:max_dim,3), P3_center(3), fact_p3, pp3, p3_inv
  complex*16          :: P4_new(0:max_dim,3), P4_center(3), fact_p4, pp4, p4_inv
  complex*16          :: Q1_new(0:max_dim,3), Q1_center(3), fact_q1, qq1, q1_inv
  complex*16          :: Q2_new(0:max_dim,3), Q2_center(3), fact_q2, qq2, q2_inv
  complex*16          :: integral1, integral2, integral3, integral4
  complex*16          :: integral5, integral6, integral7, integral8
  complex*16          :: integral_tot

  double precision    :: ao_two_e_integral_cosgtos_schwartz_accel
  complex*16          :: ERI_cosgtos
  complex*16          :: general_primitive_integral_cosgtos

  if(ao_prim_num(i) * ao_prim_num(j) * ao_prim_num(k) * ao_prim_num(l) > 1024) then

    !print *, ' with shwartz acc '
    ao_two_e_integral_cosgtos = ao_two_e_integral_cosgtos_schwartz_accel(i, j, k, l)

  else
    !print *, ' without shwartz acc '

    dim1 = n_pt_max_integrals

    num_i = ao_nucl(i)
    num_j = ao_nucl(j)
    num_k = ao_nucl(k)
    num_l = ao_nucl(l)

    ao_two_e_integral_cosgtos = 0.d0

    if(num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k) then
      !print *, ' not the same center'

      do p = 1, 3
        I_power(p)  = ao_power(i,p)
        J_power(p)  = ao_power(j,p)
        K_power(p)  = ao_power(k,p)
        L_power(p)  = ao_power(l,p)
        I_center(p) = nucl_coord(num_i,p) * (1.d0, 0.d0) 
        J_center(p) = nucl_coord(num_j,p) * (1.d0, 0.d0)
        K_center(p) = nucl_coord(num_k,p) * (1.d0, 0.d0)
        L_center(p) = nucl_coord(num_l,p) * (1.d0, 0.d0)
      enddo

      do p = 1, ao_prim_num(i)
        coef1 = ao_coef_norm_ord_transp_cosgtos(p,i)
        expo1 = ao_expo_ord_transp_cosgtos(p,i) 

        do q = 1, ao_prim_num(j)
          coef2 = coef1 * ao_coef_norm_ord_transp_cosgtos(q,j)
          expo2 = ao_expo_ord_transp_cosgtos(q,j) 

          call give_explicit_cpoly_and_cgaussian( P1_new, P1_center, pp1, fact_p1, iorder_p1               &
                                                , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
          p1_inv = (1.d0,0.d0) / pp1

          call give_explicit_cpoly_and_cgaussian( P2_new, P2_center, pp2, fact_p2, iorder_p2                      &
                                                , conjg(expo1), expo2, I_power, J_power, I_center, J_center, dim1 )
          p2_inv = (1.d0,0.d0) / pp2

          call give_explicit_cpoly_and_cgaussian( P3_new, P3_center, pp3, fact_p3, iorder_p3                      &
                                                , expo1, conjg(expo2), I_power, J_power, I_center, J_center, dim1 )
          p3_inv = (1.d0,0.d0) / pp3

          call give_explicit_cpoly_and_cgaussian( P4_new, P4_center, pp4, fact_p4, iorder_p4                             &
                                                , conjg(expo1), conjg(expo2), I_power, J_power, I_center, J_center, dim1 )
          p4_inv = (1.d0,0.d0) / pp4

          !integer :: ii
          !do ii = 1, 3
          !  print *, 'fact_p1', fact_p1
          !  print *, 'fact_p2', fact_p2
          !  print *, 'fact_p3', fact_p3
          !  print *, 'fact_p4', fact_p4
          !  !print *, pp1, p1_inv
          !  !print *, pp2, p2_inv
          !  !print *, pp3, p3_inv
          !  !print *, pp4, p4_inv
          !enddo 
          !  if( abs(aimag(P1_center(ii))) .gt. 0.d0 ) then
          !    print *, ' P_1 is complex !!'
          !    print *, P1_center
          !    print *, expo1, expo2
          !    print *, conjg(expo1), conjg(expo2)
          !    stop
          !  endif
          !  if( abs(aimag(P2_center(ii))) .gt. 0.d0 ) then
          !    print *, ' P_2 is complex !!'
          !    print *, P2_center
          !    print *, ' old expos:'
          !    print *, expo1, expo2
          !    print *, conjg(expo1), conjg(expo2)
          !    print *, ' new expo:'
          !    print *, pp2, p2_inv
          !    print *, ' factor:'
          !    print *, fact_p2
          !    print *, ' old centers:'
          !    print *, I_center, J_center
          !    print *, ' powers:'
          !    print *, I_power, J_power
          !    stop
          !  endif
          !  if( abs(aimag(P3_center(ii))) .gt. 0.d0 ) then
          !    print *, ' P_3 is complex !!'
          !    print *, P3_center
          !    print *, expo1, expo2
          !    print *, conjg(expo1), conjg(expo2)
          !    stop
          !  endif
          !  if( abs(aimag(P4_center(ii))) .gt. 0.d0 ) then
          !    print *, ' P_4 is complex !!'
          !    print *, P4_center
          !    print *, expo1, expo2
          !    print *, conjg(expo1), conjg(expo2)
          !    stop
          !  endif
          !enddo

          do r = 1, ao_prim_num(k)
            coef3 = coef2 * ao_coef_norm_ord_transp_cosgtos(r,k)
            expo3 = ao_expo_ord_transp_cosgtos(r,k) 

            do s = 1, ao_prim_num(l)
              coef4 = coef3 * ao_coef_norm_ord_transp_cosgtos(s,l)
              expo4 = ao_expo_ord_transp_cosgtos(s,l) 

              call give_explicit_cpoly_and_cgaussian( Q1_new, Q1_center, qq1, fact_q1, iorder_q1               &
                                                    , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
              q1_inv = (1.d0,0.d0) / qq1

              call give_explicit_cpoly_and_cgaussian( Q2_new, Q2_center, qq2, fact_q2, iorder_q2               &
                                                    , conjg(expo3), expo4, K_power, L_power, K_center, L_center, dim1 )
              q2_inv = (1.d0,0.d0) / qq2

              !do ii = 1, 3
              !  !print *, qq1, q1_inv
              !  !print *, qq2, q2_inv
              !  print *, 'fact_q1', fact_q1
              !  print *, 'fact_q2', fact_q2
              !enddo 
              !  if( abs(aimag(Q1_center(ii))) .gt. 0.d0 ) then
              !    print *, ' Q_1 is complex !!'
              !    print *, Q1_center
              !    print *, expo3, expo4
              !    print *, conjg(expo3), conjg(expo4)
              !    stop
              !  endif
              !  if( abs(aimag(Q2_center(ii))) .gt. 0.d0 ) then
              !    print *, ' Q_2 is complex !!'
              !    print *, Q2_center
              !    print *, expo3, expo4
              !    print *, conjg(expo3), conjg(expo4)
              !    stop
              !  endif
              !enddo


              integral1 = general_primitive_integral_cosgtos( dim1, P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 &
                                                                  , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q1 )

              integral2 = general_primitive_integral_cosgtos( dim1, P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 &
                                                                  , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q2 )

              integral3 = general_primitive_integral_cosgtos( dim1, P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 &
                                                                  , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q1 )

              integral4 = general_primitive_integral_cosgtos( dim1, P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 &
                                                                  , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q2 )

              integral5 = general_primitive_integral_cosgtos( dim1, P3_new, P3_center, fact_p3, pp3, p3_inv, iorder_p3 &
                                                                  , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q1 )

              integral6 = general_primitive_integral_cosgtos( dim1, P3_new, P3_center, fact_p3, pp3, p3_inv, iorder_p3 &
                                                                  , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q2 )

              integral7 = general_primitive_integral_cosgtos( dim1, P4_new, P4_center, fact_p4, pp4, p4_inv, iorder_p4 &
                                                                  , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q1 )

              integral8 = general_primitive_integral_cosgtos( dim1, P4_new, P4_center, fact_p4, pp4, p4_inv, iorder_p4 &
                                                                  , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q2 )

              integral_tot = integral1 + integral2 + integral3 + integral4 + integral5 + integral6 + integral7 + integral8

              !integral_tot = integral1
              !print*, integral_tot

              ao_two_e_integral_cosgtos = ao_two_e_integral_cosgtos + coef4 * 2.d0 * real(integral_tot)
            enddo ! s
          enddo  ! r
        enddo   ! q
      enddo    ! p

    else
      !print *, ' the same center'

      do p = 1, 3
        I_power(p) = ao_power(i,p)
        J_power(p) = ao_power(j,p)
        K_power(p) = ao_power(k,p)
        L_power(p) = ao_power(l,p)
      enddo

      do p = 1, ao_prim_num(i)
        coef1 = ao_coef_norm_ord_transp_cosgtos(p,i)
        expo1 = ao_expo_ord_transp_cosgtos(p,i) 

        do q = 1, ao_prim_num(j)
          coef2 = coef1 * ao_coef_norm_ord_transp_cosgtos(q,j)
          expo2 = ao_expo_ord_transp_cosgtos(q,j) 

          do r = 1, ao_prim_num(k)
            coef3 = coef2 * ao_coef_norm_ord_transp_cosgtos(r,k)
            expo3 = ao_expo_ord_transp_cosgtos(r,k) 

            do s = 1, ao_prim_num(l)
              coef4 = coef3 * ao_coef_norm_ord_transp_cosgtos(s,l)
              expo4 = ao_expo_ord_transp_cosgtos(s,l) 

              integral1 = ERI_cosgtos( expo1, expo2, expo3, expo4                     &
                                     , I_power(1), J_power(1), K_power(1), L_power(1) &
                                     , I_power(2), J_power(2), K_power(2), L_power(2) &
                                     , I_power(3), J_power(3), K_power(3), L_power(3) )

              integral2 = ERI_cosgtos( expo1, expo2, conjg(expo3), expo4              &
                                     , I_power(1), J_power(1), K_power(1), L_power(1) &
                                     , I_power(2), J_power(2), K_power(2), L_power(2) &
                                     , I_power(3), J_power(3), K_power(3), L_power(3) )

              integral3 = ERI_cosgtos( conjg(expo1), expo2, expo3, expo4              &
                                     , I_power(1), J_power(1), K_power(1), L_power(1) &
                                     , I_power(2), J_power(2), K_power(2), L_power(2) &
                                     , I_power(3), J_power(3), K_power(3), L_power(3) )

              integral4 = ERI_cosgtos( conjg(expo1), expo2, conjg(expo3), expo4       &
                                     , I_power(1), J_power(1), K_power(1), L_power(1) &
                                     , I_power(2), J_power(2), K_power(2), L_power(2) &
                                     , I_power(3), J_power(3), K_power(3), L_power(3) )

              integral5 = ERI_cosgtos( expo1, conjg(expo2), expo3, expo4              &
                                     , I_power(1), J_power(1), K_power(1), L_power(1) &
                                     , I_power(2), J_power(2), K_power(2), L_power(2) &
                                     , I_power(3), J_power(3), K_power(3), L_power(3) )

              integral6 = ERI_cosgtos( expo1, conjg(expo2), conjg(expo3), expo4       &
                                     , I_power(1), J_power(1), K_power(1), L_power(1) &
                                     , I_power(2), J_power(2), K_power(2), L_power(2) &
                                     , I_power(3), J_power(3), K_power(3), L_power(3) )

              integral7 = ERI_cosgtos( conjg(expo1), conjg(expo2), expo3, expo4       &
                                     , I_power(1), J_power(1), K_power(1), L_power(1) &
                                     , I_power(2), J_power(2), K_power(2), L_power(2) &
                                     , I_power(3), J_power(3), K_power(3), L_power(3) )

              integral8 = ERI_cosgtos( conjg(expo1), conjg(expo2), conjg(expo3), expo4 &
                                     , I_power(1), J_power(1), K_power(1), L_power(1)  &
                                     , I_power(2), J_power(2), K_power(2), L_power(2)  &
                                     , I_power(3), J_power(3), K_power(3), L_power(3)  )

              integral_tot = integral1 + integral2 + integral3 + integral4 + integral5 + integral6 + integral7 + integral8


              ao_two_e_integral_cosgtos = ao_two_e_integral_cosgtos + coef4 * 2.d0 * real(integral_tot)
            enddo ! s
          enddo  ! r
        enddo   ! q
      enddo    ! p

    endif
  endif

end function ao_two_e_integral_cosgtos

! ---

double precision function ao_two_e_integral_cosgtos_schwartz_accel(i, j, k, l)

  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in)           :: i, j, k, l

  integer                       :: p, q, r, s
  integer                       :: num_i, num_j, num_k, num_l, dim1, I_power(3), J_power(3), K_power(3), L_power(3)
  integer                       :: iorder_p1(3), iorder_p2(3), iorder_p3(3), iorder_p4(3), iorder_q1(3), iorder_q2(3)
  double precision              :: coef1, coef2, coef3, coef4
  complex*16                    :: I_center(3), J_center(3), K_center(3), L_center(3)
  complex*16                    :: expo1, expo2, expo3, expo4
  complex*16                    :: P1_new(0:max_dim,3), P1_center(3), fact_p1, pp1, p1_inv
  complex*16                    :: P2_new(0:max_dim,3), P2_center(3), fact_p2, pp2, p2_inv
  complex*16                    :: P3_new(0:max_dim,3), P3_center(3), fact_p3, pp3, p3_inv
  complex*16                    :: P4_new(0:max_dim,3), P4_center(3), fact_p4, pp4, p4_inv
  complex*16                    :: Q1_new(0:max_dim,3), Q1_center(3), fact_q1, qq1, q1_inv
  complex*16                    :: Q2_new(0:max_dim,3), Q2_center(3), fact_q2, qq2, q2_inv
  complex*16                    :: integral1, integral2, integral3, integral4
  complex*16                    :: integral5, integral6, integral7, integral8
  complex*16                    :: integral_tot

  double precision, allocatable :: schwartz_kl(:,:)
  double precision              :: thr
  double precision              :: schwartz_ij

  complex*16                    :: ERI_cosgtos
  complex*16                    :: general_primitive_integral_cosgtos

  ao_two_e_integral_cosgtos_schwartz_accel = 0.d0

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)


  thr = ao_integrals_threshold*ao_integrals_threshold

  allocate( schwartz_kl(0:ao_prim_num(l),0:ao_prim_num(k)) )

  if(num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k) then

    do p = 1, 3
      I_power(p)  = ao_power(i,p)
      J_power(p)  = ao_power(j,p)
      K_power(p)  = ao_power(k,p)
      L_power(p)  = ao_power(l,p)
      I_center(p) = nucl_coord(num_i,p) * (1.d0, 0.d0)
      J_center(p) = nucl_coord(num_j,p) * (1.d0, 0.d0)
      K_center(p) = nucl_coord(num_k,p) * (1.d0, 0.d0)
      L_center(p) = nucl_coord(num_l,p) * (1.d0, 0.d0)
    enddo

    schwartz_kl(0,0) = 0.d0
    do r = 1, ao_prim_num(k)
      coef1 = ao_coef_norm_ord_transp_cosgtos(r,k) * ao_coef_norm_ord_transp_cosgtos(r,k)
      expo1 = ao_expo_ord_transp_cosgtos(r,k) 

      schwartz_kl(0,r) = 0.d0
      do s = 1, ao_prim_num(l)
        coef2 = coef1 * ao_coef_norm_ord_transp_cosgtos(s,l) * ao_coef_norm_ord_transp_cosgtos(s,l)
        expo2 = ao_expo_ord_transp_cosgtos(s,l) 

        call give_explicit_cpoly_and_cgaussian( P1_new, P1_center, pp1, fact_p1, iorder_p1               &
                                              , expo1, expo2, K_power, L_power, K_center, L_center, dim1 )
        p1_inv = (1.d0,0.d0) / pp1

        call give_explicit_cpoly_and_cgaussian( P2_new, P2_center, pp2, fact_p2, iorder_p2                      &
                                              , conjg(expo1), expo2, K_power, L_power, K_center, L_center, dim1 )
        p2_inv = (1.d0,0.d0) / pp2

        call give_explicit_cpoly_and_cgaussian( P3_new, P3_center, pp3, fact_p3, iorder_p3                      &
                                              , expo1, conjg(expo2), K_power, L_power, K_center, L_center, dim1 )
        p3_inv = (1.d0,0.d0) / pp3

        call give_explicit_cpoly_and_cgaussian( P4_new, P4_center, pp4, fact_p4, iorder_p4                             &
                                              , conjg(expo1), conjg(expo2), K_power, L_power, K_center, L_center, dim1 )
        p4_inv = (1.d0,0.d0) / pp4

        integral1 = general_primitive_integral_cosgtos( dim1, P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 &
                                                            , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 )

        integral2 = general_primitive_integral_cosgtos( dim1, P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 &
                                                            , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 )

        integral3 = general_primitive_integral_cosgtos( dim1, P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 &
                                                            , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 )

        integral4 = general_primitive_integral_cosgtos( dim1, P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 &
                                                            , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 )

        integral5 = general_primitive_integral_cosgtos( dim1, P3_new, P3_center, fact_p3, pp3, p3_inv, iorder_p3 &
                                                            , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 )

        integral6 = general_primitive_integral_cosgtos( dim1, P3_new, P3_center, fact_p3, pp3, p3_inv, iorder_p3 &
                                                            , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 )

        integral7 = general_primitive_integral_cosgtos( dim1, P4_new, P4_center, fact_p4, pp4, p4_inv, iorder_p4 &
                                                            , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 )

        integral8 = general_primitive_integral_cosgtos( dim1, P4_new, P4_center, fact_p4, pp4, p4_inv, iorder_p4 &
                                                            , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 )

        integral_tot = integral1 + integral2 + integral3 + integral4 + integral5 + integral6 + integral7 + integral8


        schwartz_kl(s,r) = coef2 * 2.d0 * real(integral_tot)

        schwartz_kl(0,r) = max(schwartz_kl(0,r), schwartz_kl(s,r))
      enddo

      schwartz_kl(0,0) = max(schwartz_kl(0,r), schwartz_kl(0,0))
    enddo


    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_norm_ord_transp_cosgtos(p,i)
      expo1 = ao_expo_ord_transp_cosgtos(p,i) 

      do q = 1, ao_prim_num(j)
        coef2 = coef1 * ao_coef_norm_ord_transp_cosgtos(q,j)
        expo2 = ao_expo_ord_transp_cosgtos(q,j) 

        call give_explicit_cpoly_and_cgaussian( P1_new, P1_center, pp1, fact_p1, iorder_p1               &
                                              , expo1, expo2, I_power, J_power, I_center, J_center, dim1 )
        p1_inv = (1.d0,0.d0) / pp1

        call give_explicit_cpoly_and_cgaussian( P2_new, P2_center, pp2, fact_p2, iorder_p2                      &
                                              , conjg(expo1), expo2, I_power, J_power, I_center, J_center, dim1 )
        p2_inv = (1.d0,0.d0) / pp2

        call give_explicit_cpoly_and_cgaussian( P3_new, P3_center, pp3, fact_p3, iorder_p3                      &
                                              , expo1, conjg(expo2), I_power, J_power, I_center, J_center, dim1 )
        p3_inv = (1.d0,0.d0) / pp3

        call give_explicit_cpoly_and_cgaussian( P4_new, P4_center, pp4, fact_p4, iorder_p4                             &
                                              , conjg(expo1), conjg(expo2), I_power, J_power, I_center, J_center, dim1 )
        p4_inv = (1.d0,0.d0) / pp4

        integral1 = general_primitive_integral_cosgtos( dim1, P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 &
                                                            , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 )

        integral2 = general_primitive_integral_cosgtos( dim1, P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 &
                                                            , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 )

        integral3 = general_primitive_integral_cosgtos( dim1, P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 &
                                                            , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 )

        integral4 = general_primitive_integral_cosgtos( dim1, P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 &
                                                            , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 )

        integral5 = general_primitive_integral_cosgtos( dim1, P3_new, P3_center, fact_p3, pp3, p3_inv, iorder_p3 &
                                                            , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 )

        integral6 = general_primitive_integral_cosgtos( dim1, P3_new, P3_center, fact_p3, pp3, p3_inv, iorder_p3 &
                                                            , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 )

        integral7 = general_primitive_integral_cosgtos( dim1, P4_new, P4_center, fact_p4, pp4, p4_inv, iorder_p4 &
                                                            , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 )

        integral8 = general_primitive_integral_cosgtos( dim1, P4_new, P4_center, fact_p4, pp4, p4_inv, iorder_p4 &
                                                            , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 )

        integral_tot = integral1 + integral2 + integral3 + integral4 + integral5 + integral6 + integral7 + integral8

        schwartz_ij = coef2 * coef2 * 2.d0 * real(integral_tot)

        if(schwartz_kl(0,0)*schwartz_ij < thr) cycle

        do r = 1, ao_prim_num(k)
          if(schwartz_kl(0,r)*schwartz_ij < thr) cycle

          coef3 = coef2 * ao_coef_norm_ord_transp_cosgtos(r,k)
          expo3 = ao_expo_ord_transp_cosgtos(r,k)

          do s = 1, ao_prim_num(l)
            if(schwartz_kl(s,r)*schwartz_ij < thr) cycle

            coef4 = coef3 * ao_coef_norm_ord_transp_cosgtos(s,l)
            expo4 = ao_expo_ord_transp_cosgtos(s,l)

            call give_explicit_cpoly_and_cgaussian( Q1_new, Q1_center, qq1, fact_q1, iorder_q1               &
                                                  , expo3, expo4, K_power, L_power, K_center, L_center, dim1 )
            q1_inv = (1.d0,0.d0) / qq1

            call give_explicit_cpoly_and_cgaussian( Q2_new, Q2_center, qq2, fact_q2, iorder_q2               &
                                                  , conjg(expo3), expo4, K_power, L_power, K_center, L_center, dim1 )
            q2_inv = (1.d0,0.d0) / qq2

            integral1 = general_primitive_integral_cosgtos( dim1, P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 &
                                                                , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q1 )

            integral2 = general_primitive_integral_cosgtos( dim1, P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p1 &
                                                                , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q2 )

            integral3 = general_primitive_integral_cosgtos( dim1, P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 &
                                                                , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q1 )

            integral4 = general_primitive_integral_cosgtos( dim1, P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p2 &
                                                                , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q2 )


            integral5 = general_primitive_integral_cosgtos( dim1, P3_new, P3_center, fact_p3, pp3, p3_inv, iorder_p3 &
                                                                , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q1 )

            integral6 = general_primitive_integral_cosgtos( dim1, P3_new, P3_center, fact_p3, pp3, p3_inv, iorder_p3 &
                                                                , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q2 )

            integral7 = general_primitive_integral_cosgtos( dim1, P4_new, P4_center, fact_p4, pp4, p4_inv, iorder_p4 &
                                                                , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q1 )

            integral8 = general_primitive_integral_cosgtos( dim1, P4_new, P4_center, fact_p4, pp4, p4_inv, iorder_p4 &
                                                                , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q2 )

            integral_tot = integral1 + integral2 + integral3 + integral4 + integral5 + integral6 + integral7 + integral8

            ao_two_e_integral_cosgtos_schwartz_accel = ao_two_e_integral_cosgtos_schwartz_accel &
                                                     + coef4 * 2.d0 * real(integral_tot)
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  else

    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
    enddo

    schwartz_kl(0,0) = 0.d0
    do r = 1, ao_prim_num(k)
      coef1 = ao_coef_norm_ord_transp_cosgtos(r,k) * ao_coef_norm_ord_transp_cosgtos(r,k)
      expo1 = ao_expo_ord_transp_cosgtos(r,k)

      schwartz_kl(0,r) = 0.d0
      do s = 1, ao_prim_num(l)
        coef2 = coef1 * ao_coef_norm_ord_transp_cosgtos(s,l) * ao_coef_norm_ord_transp_cosgtos(s,l)
        expo2 = ao_expo_ord_transp_cosgtos(s,l)

        integral1 = ERI_cosgtos( expo1, expo2, expo1, expo2                     &
                               , K_power(1), L_power(1), K_power(1), L_power(1) &
                               , K_power(2), L_power(2), K_power(2), L_power(2) &
                               , K_power(3), L_power(3), K_power(3), L_power(3) )
        integral2 = ERI_cosgtos( expo1, expo2, conjg(expo1), expo2              &
                               , K_power(1), L_power(1), K_power(1), L_power(1) &
                               , K_power(2), L_power(2), K_power(2), L_power(2) &
                               , K_power(3), L_power(3), K_power(3), L_power(3) )

        integral3 = ERI_cosgtos( conjg(expo1), expo2, expo1, expo2              &
                               , K_power(1), L_power(1), K_power(1), L_power(1) &
                               , K_power(2), L_power(2), K_power(2), L_power(2) &
                               , K_power(3), L_power(3), K_power(3), L_power(3) )
        integral4 = ERI_cosgtos( conjg(expo1), expo2, conjg(expo1), expo2       &
                               , K_power(1), L_power(1), K_power(1), L_power(1) &
                               , K_power(2), L_power(2), K_power(2), L_power(2) &
                               , K_power(3), L_power(3), K_power(3), L_power(3) )

        integral5 = ERI_cosgtos( expo1, conjg(expo2), expo1, expo2              &
                               , K_power(1), L_power(1), K_power(1), L_power(1) &
                               , K_power(2), L_power(2), K_power(2), L_power(2) &
                               , K_power(3), L_power(3), K_power(3), L_power(3) )
        integral6 = ERI_cosgtos( expo1, conjg(expo2), conjg(expo1), expo2       &
                               , K_power(1), L_power(1), K_power(1), L_power(1) &
                               , K_power(2), L_power(2), K_power(2), L_power(2) &
                               , K_power(3), L_power(3), K_power(3), L_power(3) )

        integral7 = ERI_cosgtos( conjg(expo1), conjg(expo2), expo1, expo2       &
                               , K_power(1), L_power(1), K_power(1), L_power(1) &
                               , K_power(2), L_power(2), K_power(2), L_power(2) &
                               , K_power(3), L_power(3), K_power(3), L_power(3) )
        integral8 = ERI_cosgtos( conjg(expo1), conjg(expo2), conjg(expo1), expo2 &
                               , K_power(1), L_power(1), K_power(1), L_power(1)  &
                               , K_power(2), L_power(2), K_power(2), L_power(2)  &
                               , K_power(3), L_power(3), K_power(3), L_power(3)  )

        integral_tot = integral1 + integral2 + integral3 + integral4 + integral5 + integral6 + integral7 + integral8


        schwartz_kl(s,r) = coef2 * 2.d0 * real(integral_tot)

        schwartz_kl(0,r) = max(schwartz_kl(0,r), schwartz_kl(s,r))
      enddo
      schwartz_kl(0,0) = max(schwartz_kl(0,r), schwartz_kl(0,0))
    enddo

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_norm_ord_transp_cosgtos(p,i)
      expo1 = ao_expo_ord_transp_cosgtos(p,i)

      do q = 1, ao_prim_num(j)
        coef2 = coef1 * ao_coef_norm_ord_transp_cosgtos(q,j)
        expo2 = ao_expo_ord_transp_cosgtos(q,j)

        integral1 = ERI_cosgtos( expo1, expo2, expo1, expo2                     &
                               , I_power(1), J_power(1), I_power(1), J_power(1) &
                               , I_power(2), J_power(2), I_power(2), J_power(2) &
                               , I_power(3), J_power(3), I_power(3), J_power(3) )

        integral2 = ERI_cosgtos( expo1, expo2, conjg(expo1), expo2              &
                               , I_power(1), J_power(1), I_power(1), J_power(1) &
                               , I_power(2), J_power(2), I_power(2), J_power(2) &
                               , I_power(3), J_power(3), I_power(3), J_power(3) )

        integral3 = ERI_cosgtos( conjg(expo1), expo2, expo1, expo2              &
                               , I_power(1), J_power(1), I_power(1), J_power(1) &
                               , I_power(2), J_power(2), I_power(2), J_power(2) &
                               , I_power(3), J_power(3), I_power(3), J_power(3) )

        integral4 = ERI_cosgtos( conjg(expo1), expo2, conjg(expo1), expo2       &
                               , I_power(1), J_power(1), I_power(1), J_power(1) &
                               , I_power(2), J_power(2), I_power(2), J_power(2) &
                               , I_power(3), J_power(3), I_power(3), J_power(3) )

        integral5 = ERI_cosgtos( expo1, conjg(expo2), expo1, expo2              &
                               , I_power(1), J_power(1), I_power(1), J_power(1) &
                               , I_power(2), J_power(2), I_power(2), J_power(2) &
                               , I_power(3), J_power(3), I_power(3), J_power(3) )

        integral6 = ERI_cosgtos( expo1, conjg(expo2), conjg(expo1), expo2       &
                               , I_power(1), J_power(1), I_power(1), J_power(1) &
                               , I_power(2), J_power(2), I_power(2), J_power(2) &
                               , I_power(3), J_power(3), I_power(3), J_power(3) )

        integral7 = ERI_cosgtos( conjg(expo1), conjg(expo2), expo1, expo2       &
                               , I_power(1), J_power(1), I_power(1), J_power(1) &
                               , I_power(2), J_power(2), I_power(2), J_power(2) &
                               , I_power(3), J_power(3), I_power(3), J_power(3) )

        integral8 = ERI_cosgtos( conjg(expo1), conjg(expo2), conjg(expo1), expo2 &
                               , I_power(1), J_power(1), I_power(1), J_power(1)  &
                               , I_power(2), J_power(2), I_power(2), J_power(2)  &
                               , I_power(3), J_power(3), I_power(3), J_power(3)  )

        integral_tot = integral1 + integral2 + integral3 + integral4 + integral5 + integral6 + integral7 + integral8

        schwartz_ij = coef2 * coef2 * 2.d0 * real(integral_tot)

        if(schwartz_kl(0,0)*schwartz_ij < thr) cycle
        do r = 1, ao_prim_num(k)
          if(schwartz_kl(0,r)*schwartz_ij < thr) cycle

          coef3 = coef2 * ao_coef_norm_ord_transp_cosgtos(r,k)
          expo3 = ao_expo_ord_transp_cosgtos(r,k)

          do s = 1, ao_prim_num(l)
            if(schwartz_kl(s,r)*schwartz_ij < thr) cycle

            coef4 = coef3 * ao_coef_norm_ord_transp_cosgtos(s,l)
            expo4 = ao_expo_ord_transp_cosgtos(s,l)

            integral1 = ERI_cosgtos( expo1, expo2, expo3, expo4                     &
                                   , I_power(1), J_power(1), K_power(1), L_power(1) &
                                   , I_power(2), J_power(2), K_power(2), L_power(2) &
                                   , I_power(3), J_power(3), K_power(3), L_power(3) )

            integral2 = ERI_cosgtos( expo1, expo2, conjg(expo3), expo4              &
                                   , I_power(1), J_power(1), K_power(1), L_power(1) &
                                   , I_power(2), J_power(2), K_power(2), L_power(2) &
                                   , I_power(3), J_power(3), K_power(3), L_power(3) )

            integral3 = ERI_cosgtos( conjg(expo1), expo2, expo3, expo4              &
                                   , I_power(1), J_power(1), K_power(1), L_power(1) &
                                   , I_power(2), J_power(2), K_power(2), L_power(2) &
                                   , I_power(3), J_power(3), K_power(3), L_power(3) )

            integral4 = ERI_cosgtos( conjg(expo1), expo2, conjg(expo3), expo4       &
                                   , I_power(1), J_power(1), K_power(1), L_power(1) &
                                   , I_power(2), J_power(2), K_power(2), L_power(2) &
                                   , I_power(3), J_power(3), K_power(3), L_power(3) )

            integral5 = ERI_cosgtos( expo1, conjg(expo2), expo3, expo4              &
                                   , I_power(1), J_power(1), K_power(1), L_power(1) &
                                   , I_power(2), J_power(2), K_power(2), L_power(2) &
                                   , I_power(3), J_power(3), K_power(3), L_power(3) )

            integral6 = ERI_cosgtos( expo1, conjg(expo2), conjg(expo3), expo4       &
                                   , I_power(1), J_power(1), K_power(1), L_power(1) &
                                   , I_power(2), J_power(2), K_power(2), L_power(2) &
                                   , I_power(3), J_power(3), K_power(3), L_power(3) )

            integral7 = ERI_cosgtos( conjg(expo1), conjg(expo2), expo3, expo4       &
                                   , I_power(1), J_power(1), K_power(1), L_power(1) &
                                   , I_power(2), J_power(2), K_power(2), L_power(2) &
                                   , I_power(3), J_power(3), K_power(3), L_power(3) )

            integral8 = ERI_cosgtos( conjg(expo1), conjg(expo2), conjg(expo3), expo4 &
                                   , I_power(1), J_power(1), K_power(1), L_power(1)  &
                                   , I_power(2), J_power(2), K_power(2), L_power(2)  &
                                   , I_power(3), J_power(3), K_power(3), L_power(3)  )

            integral_tot = integral1 + integral2 + integral3 + integral4 + integral5 + integral6 + integral7 + integral8

            ao_two_e_integral_cosgtos_schwartz_accel = ao_two_e_integral_cosgtos_schwartz_accel &
                                      + coef4 * 2.d0 * real(integral_tot)
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  endif

  deallocate(schwartz_kl)

end function ao_two_e_integral_cosgtos_schwartz_accel

! ---

BEGIN_PROVIDER [ double precision, ao_two_e_integral_cosgtos_schwartz, (ao_num,ao_num)  ]

  BEGIN_DOC
  !  Needed to compute Schwartz inequalities
  END_DOC

  implicit none
  integer          :: i, k
  double precision :: ao_two_e_integral_cosgtos

  ao_two_e_integral_cosgtos_schwartz(1,1) = ao_two_e_integral_cosgtos(1, 1, 1, 1)

 !$OMP PARALLEL DO PRIVATE(i,k)                                       &
 !$OMP             DEFAULT(NONE)                                      &
 !$OMP             SHARED(ao_num, ao_two_e_integral_cosgtos_schwartz) &
 !$OMP             SCHEDULE(dynamic)
  do i = 1, ao_num
    do k = 1, i
      ao_two_e_integral_cosgtos_schwartz(i,k) = dsqrt(ao_two_e_integral_cosgtos(i, i, k, k))
      ao_two_e_integral_cosgtos_schwartz(k,i) = ao_two_e_integral_cosgtos_schwartz(i,k)
    enddo
  enddo
 !$OMP END PARALLEL DO

END_PROVIDER

! ---

complex*16 function general_primitive_integral_cosgtos( dim, P_new, P_center, fact_p, p, p_inv, iorder_p &
                                                           , Q_new, Q_center, fact_q, q, q_inv, iorder_q )

  BEGIN_DOC
  !
  ! Computes the integral <pq|rs> where p,q,r,s are cos-cGTOS primitives
  !
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in) :: dim
  integer,    intent(in) :: iorder_p(3), iorder_q(3)
  complex*16, intent(in) :: P_new(0:max_dim,3), P_center(3), fact_p, p, p_inv
  complex*16, intent(in) :: Q_new(0:max_dim,3), Q_center(3), fact_q, q, q_inv

  integer                :: i, j, nx, ny, nz, n_Ix, n_Iy, n_Iz, iorder, n_pt_tmp, n_pt_out
  double precision       :: tmp_mod
  double precision       :: ppq_re, ppq_im, ppq_mod, sq_ppq_re, sq_ppq_im
  complex*16             :: pq, pq_inv, pq_inv_2, p01_1, p01_2, p10_1, p10_2, ppq, sq_ppq
  complex*16             :: rho, dist, const
  complex*16             :: accu, tmp_p, tmp_q
  complex*16             :: dx(0:max_dim), Ix_pol(0:max_dim), dy(0:max_dim), Iy_pol(0:max_dim), dz(0:max_dim), Iz_pol(0:max_dim)
  complex*16             :: d1(0:max_dim), d_poly(0:max_dim)

  complex*16             :: crint_sum


  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: dx, Ix_pol, dy, Iy_pol, dz, Iz_pol
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: d1, d_poly

  general_primitive_integral_cosgtos = (0.d0, 0.d0)

  pq       = (0.5d0, 0.d0) * p_inv * q_inv
  pq_inv   = (0.5d0, 0.d0) / (p + q)
  pq_inv_2 = pq_inv + pq_inv
  p10_1    = q * pq               ! 1/(2p)
  p01_1    = p * pq               ! 1/(2q)
  p10_2    = pq_inv_2 * p10_1 * q ! 0.5d0*q/(pq + p*p)
  p01_2    = pq_inv_2 * p01_1 * p ! 0.5d0*p/(q*q + pq)

  ! get \sqrt(p + q)
  !ppq     = p + q
  !ppq_re  = REAL (ppq)
  !ppq_im  = AIMAG(ppq)
  !ppq_mod = dsqrt(ppq_re*ppq_re + ppq_im*ppq_im)
  !sq_ppq_re = sq_op5 * dsqrt(ppq_re + ppq_mod)
  !sq_ppq_im = 0.5d0 * ppq_im / sq_ppq_re 
  !sq_ppq    = sq_ppq_re + (0.d0, 1.d0) * sq_ppq_im
  sq_ppq = zsqrt(p + q)

  ! ---

  iorder = iorder_p(1) + iorder_q(1) + iorder_p(1) + iorder_q(1)

  do i = 0, iorder
    Ix_pol(i) = (0.d0, 0.d0)
  enddo

  n_Ix = 0
  do i = 0, iorder_p(1)

    tmp_p   = P_new(i,1)
    tmp_mod = dsqrt(REAL(tmp_p)*REAL(tmp_p) + AIMAG(tmp_p)*AIMAG(tmp_p))
    if(tmp_mod < thresh) cycle

    do j = 0, iorder_q(1)

      tmp_q   = tmp_p * Q_new(j,1)
      tmp_mod = dsqrt(REAL(tmp_q)*REAL(tmp_q) + AIMAG(tmp_q)*AIMAG(tmp_q))
      if(tmp_mod < thresh) cycle

      !DIR$ FORCEINLINE
      call give_cpolynom_mult_center_x(P_center(1), Q_center(1), i, j, p, q, iorder, pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dx, nx)
      !DIR$ FORCEINLINE
      call add_cpoly_multiply(dx, nx, tmp_q, Ix_pol, n_Ix)
    enddo
  enddo
  if(n_Ix == -1) then
    return
  endif

  ! ---

  iorder = iorder_p(2) + iorder_q(2) + iorder_p(2) + iorder_q(2)

  do i = 0, iorder
    Iy_pol(i) = (0.d0, 0.d0)
  enddo

  n_Iy = 0
  do i = 0, iorder_p(2)

    tmp_p   = P_new(i,2)
    tmp_mod = dsqrt(REAL(tmp_p)*REAL(tmp_p) + AIMAG(tmp_p)*AIMAG(tmp_p))
    if(tmp_mod < thresh) cycle

    do j = 0, iorder_q(2)

      tmp_q   = tmp_p * Q_new(j,2)
      tmp_mod = dsqrt(REAL(tmp_q)*REAL(tmp_q) + AIMAG(tmp_q)*AIMAG(tmp_q))
      if(tmp_mod < thresh) cycle

      !DIR$ FORCEINLINE
      call give_cpolynom_mult_center_x(P_center(2), Q_center(2), i, j, p, q, iorder, pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dy, ny)
      !DIR$ FORCEINLINE
      call add_cpoly_multiply(dy, ny, tmp_q, Iy_pol, n_Iy)
    enddo
  enddo

  if(n_Iy == -1) then
    return
  endif

  ! ---

  iorder = iorder_p(3) + iorder_q(3) + iorder_p(3) + iorder_q(3)

  do i = 0, iorder
    Iz_pol(i) = (0.d0, 0.d0)
  enddo

  n_Iz = 0
  do i = 0, iorder_p(3)

    tmp_p   = P_new(i,3)
    tmp_mod = dsqrt(REAL(tmp_p)*REAL(tmp_p) + AIMAG(tmp_p)*AIMAG(tmp_p))
    if(tmp_mod < thresh) cycle

    do j = 0, iorder_q(3)

      tmp_q   = tmp_p * Q_new(j,3)
      tmp_mod = dsqrt(REAL(tmp_q)*REAL(tmp_q) + AIMAG(tmp_q)*AIMAG(tmp_q))
      if(tmp_mod < thresh) cycle

      !DIR$ FORCEINLINE
      call give_cpolynom_mult_center_x(P_center(3), Q_center(3), i, j, p, q, iorder, pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dz, nz)
      !DIR$ FORCEINLINE
      call add_cpoly_multiply(dz, nz, tmp_q, Iz_pol, n_Iz)
    enddo
  enddo

  if(n_Iz == -1) then
    return
  endif

  ! ---

  rho = p * q * pq_inv_2
  dist = (P_center(1) - Q_center(1)) * (P_center(1) - Q_center(1)) &
       + (P_center(2) - Q_center(2)) * (P_center(2) - Q_center(2)) &
       + (P_center(3) - Q_center(3)) * (P_center(3) - Q_center(3))
  const = dist * rho

  n_pt_tmp = n_Ix + n_Iy
  do i = 0, n_pt_tmp
    d_poly(i) = (0.d0, 0.d0)
  enddo

  !DIR$ FORCEINLINE
  call multiply_cpoly(Ix_pol, n_Ix, Iy_pol, n_Iy, d_poly, n_pt_tmp)
  if(n_pt_tmp == -1) then
    return
  endif
  n_pt_out = n_pt_tmp + n_Iz
  do i = 0, n_pt_out
    d1(i) = (0.d0, 0.d0)
  enddo

  !DIR$ FORCEINLINE
  call multiply_cpoly(d_poly, n_pt_tmp, Iz_pol, n_Iz, d1, n_pt_out)

  accu = crint_sum(n_pt_out, const, d1)
!  print *, n_pt_out, real(d1(0:n_pt_out))
!  print *, real(accu)

  general_primitive_integral_cosgtos = fact_p * fact_q * accu * pi_5_2 * p_inv * q_inv / sq_ppq

end function general_primitive_integral_cosgtos

! ---

complex*16 function ERI_cosgtos(alpha, beta, delta, gama, a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y, a_z, b_z, c_z, d_z)

  BEGIN_DOC
  !  ATOMIC PRIMTIVE two-electron integral between the 4 primitives ::
  !         primitive_1 = x1**(a_x) y1**(a_y) z1**(a_z) exp(-alpha * r1**2)
  !         primitive_2 = x1**(b_x) y1**(b_y) z1**(b_z) exp(- beta * r1**2)
  !         primitive_3 = x2**(c_x) y2**(c_y) z2**(c_z) exp(-delta * r2**2)
  !         primitive_4 = x2**(d_x) y2**(d_y) z2**(d_z) exp(- gama * r2**2)
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in) :: a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y, a_z, b_z, c_z, d_z
  complex*16, intent(in) :: delta, gama, alpha, beta

  integer                :: a_x_2, b_x_2, c_x_2, d_x_2, a_y_2, b_y_2, c_y_2, d_y_2, a_z_2, b_z_2, c_z_2, d_z_2
  integer                :: i, j, k, l, n_pt
  integer                :: nx, ny, nz
  double precision       :: ppq_re, ppq_im, ppq_mod, sq_ppq_re, sq_ppq_im
  complex*16             :: p, q, ppq, sq_ppq, coeff, I_f

  ERI_cosgtos = (0.d0, 0.d0)

  ASSERT (REAL(alpha) >= 0.d0)
  ASSERT (REAL(beta ) >= 0.d0)
  ASSERT (REAL(delta) >= 0.d0)
  ASSERT (REAL(gama ) >= 0.d0)

  nx = a_x + b_x + c_x + d_x
  if(iand(nx,1) == 1) then
    ERI_cosgtos = (0.d0, 0.d0)
    return
  endif

  ny = a_y + b_y + c_y + d_y
  if(iand(ny,1) == 1) then
    ERI_cosgtos = (0.d0, 0.d0)
    return
  endif

  nz = a_z + b_z + c_z + d_z
  if(iand(nz,1) == 1) then
    ERI_cosgtos = (0.d0, 0.d0)
    return
  endif

  n_pt = shiftl(nx+ny+nz, 1)

  p = alpha + beta
  q = delta + gama

  ! get \sqrt(p + q)
  !ppq     = p + q
  !ppq_re  = REAL (ppq)
  !ppq_im  = AIMAG(ppq)
  !ppq_mod = dsqrt(ppq_re*ppq_re + ppq_im*ppq_im)
  !sq_ppq_re = sq_op5 * dsqrt(ppq_re + ppq_mod)
  !sq_ppq_im = 0.5d0 * ppq_im / sq_ppq_re 
  !sq_ppq    = sq_ppq_re + (0.d0, 1.d0) * sq_ppq_im
  sq_ppq = zsqrt(p + q)

  coeff = pi_5_2 / (p * q * sq_ppq)
  if(n_pt == 0) then
    ERI_cosgtos = coeff
    return
  endif

  call integrale_new_cosgtos(I_f, a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y, a_z, b_z, c_z, d_z, p, q, n_pt)

  ERI_cosgtos = I_f * coeff

end function ERI_cosgtos

! ---

subroutine integrale_new_cosgtos(I_f, a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y, a_z, b_z, c_z, d_z, p, q, n_pt)

  BEGIN_DOC
  ! Calculates the integral of the polynomial :
  !
  ! $I_{x_1}(a_x+b_x, c_x+d_x, p, q) \, I_{x_1}(a_y+b_y, c_y+d_y, p, q) \, I_{x_1}(a_z+b_z, c_z+d_z, p, q)$
  ! in $( 0 ; 1)$
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)  :: n_pt
  integer,    intent(in)  :: a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y, a_z, b_z, c_z, d_z
  complex*16, intent(out) :: I_f

  integer                 :: i, j, ix, iy, iz, jx, jy, jz, sx, sy, sz
  complex*16              :: p, q
  complex*16              :: pq_inv, p10_1, p10_2, p01_1, p01_2, pq_inv_2
  complex*16              :: B00(n_pt_max_integrals), B10(n_pt_max_integrals), B01(n_pt_max_integrals)
  complex*16              :: t1(n_pt_max_integrals), t2(n_pt_max_integrals)


  ASSERT (n_pt > 1)

  j = shiftr(n_pt, 1)

  pq_inv   = (0.5d0, 0.d0) / (p + q)
  p10_1    = (0.5d0, 0.d0) / p
  p01_1    = (0.5d0, 0.d0) / q
  p10_2    = (0.5d0, 0.d0) * q /(p * q + p * p)
  p01_2    = (0.5d0, 0.d0) * p /(q * q + q * p)
  pq_inv_2 = pq_inv + pq_inv

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: t1, t2, B10, B01, B00
  ix = a_x + b_x
  jx = c_x + d_x
  iy = a_y + b_y
  jy = c_y + d_y
  iz = a_z + b_z
  jz = c_z + d_z
  sx = ix  +  jx
  sy = iy  +  jy
  sz = iz  +  jz

  do i = 1, n_pt
    B10(i) = p10_1 - gauleg_t2(i, j) * p10_2
    B01(i) = p01_1 - gauleg_t2(i, j) * p01_2
    B00(i) = gauleg_t2(i, j) * pq_inv
  enddo

  if(sx > 0) then
    call I_x1_new_cosgtos(ix, jx, B10, B01, B00, t1, n_pt)
  else
    do i = 1, n_pt
      t1(i) = (1.d0, 0.d0)
    enddo
  endif

  if(sy > 0) then
    call I_x1_new_cosgtos(iy, jy, B10, B01, B00, t2, n_pt)
    do i = 1, n_pt
      t1(i) = t1(i) * t2(i)
    enddo
  endif

  if(sz > 0) then
    call I_x1_new_cosgtos(iz, jz, B10, B01, B00, t2, n_pt)
    do i = 1, n_pt
      t1(i) = t1(i) * t2(i)
    enddo
  endif

  I_f = (0.d0, 0.d0)
  do i = 1, n_pt
    I_f += gauleg_w(i, j) * t1(i)
  enddo

end subroutine integrale_new_cosgtos

! ---

recursive subroutine I_x1_new_cosgtos(a, c, B_10, B_01, B_00, res, n_pt)

  BEGIN_DOC
  !  recursive function involved in the two-electron integral
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)  :: a, c, n_pt
  complex*16, intent(in)  :: B_10(n_pt_max_integrals), B_01(n_pt_max_integrals), B_00(n_pt_max_integrals)
  complex*16, intent(out) :: res(n_pt_max_integrals)

  integer                 :: i
  complex*16              :: res2(n_pt_max_integrals)

  if(c < 0) then

    do i = 1, n_pt
      res(i) = (0.d0, 0.d0)
    enddo

  else if (a == 0) then

    call I_x2_new_cosgtos(c, B_10, B_01, B_00, res, n_pt)

  else if (a == 1) then

    call I_x2_new_cosgtos(c-1, B_10, B_01, B_00, res, n_pt)
    do i = 1, n_pt
      res(i) = dble(c) * B_00(i) * res(i)
    enddo

  else

    call I_x1_new_cosgtos(a-2, c  , B_10, B_01, B_00, res , n_pt)
    call I_x1_new_cosgtos(a-1, c-1, B_10, B_01, B_00, res2, n_pt)
    do i = 1, n_pt
      res(i) = dble(a-1) * B_10(i) * res(i) + dble(c) * B_00(i) * res2(i)
    enddo

  endif

end subroutine I_x1_new_cosgtos

! ---

recursive subroutine I_x2_new_cosgtos(c, B_10, B_01, B_00, res, n_pt)

  BEGIN_DOC
  !  recursive function involved in the two-electron integral
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)  :: c, n_pt
  complex*16, intent(in)  :: B_10(n_pt_max_integrals), B_01(n_pt_max_integrals), B_00(n_pt_max_integrals)
  complex*16, intent(out) :: res(n_pt_max_integrals)

  integer                 :: i

  if(c == 1) then

    do i = 1, n_pt
      res(i) = (0.d0, 0.d0)
    enddo

  elseif(c == 0) then

    do i = 1, n_pt
      res(i) = (1.d0, 0.d0)
    enddo

  else

    call I_x1_new_cosgtos(0, c-2, B_10, B_01, B_00, res, n_pt)
    do i = 1, n_pt
      res(i) = dble(c-1) * B_01(i) * res(i)
    enddo

  endif

end subroutine I_x2_new_cosgtos

! ---

subroutine give_cpolynom_mult_center_x( P_center, Q_center, a_x, d_x, p, q, n_pt_in &
                                      , pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, d, n_pt_out)

  BEGIN_DOC
  ! subroutine that returns the explicit polynom in term of the "t"
  ! variable of the following polynoms :
  !
  ! $I_{x_1}(a_x,d_x,p,q) \, I_{x_1}(a_y,d_y,p,q) \ I_{x_1}(a_z,d_z,p,q)$
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)  :: n_pt_in, a_x, d_x
  complex*16, intent(in)  :: P_center, Q_center, p, q, pq_inv, p10_1, p01_1, p10_2, p01_2, pq_inv_2
  integer,    intent(out) :: n_pt_out
  complex*16, intent(out) :: d(0:max_dim)

  integer                 :: n_pt1, i
  complex*16              :: B10(0:2), B01(0:2), B00(0:2), C00(0:2), D00(0:2)

  ASSERT (n_pt_in >= 0)

  B10(0) = p10_1
  B10(1) = (0.d0, 0.d0)
  B10(2) = -p10_2

  B01(0) = p01_1
  B01(1) = (0.d0, 0.d0)
  B01(2) = -p01_2

  B00(0) = (0.d0, 0.d0)
  B00(1) = (0.d0, 0.d0)
  B00(2) = pq_inv

  C00(0) = (0.d0, 0.d0)
  C00(1) = (0.d0, 0.d0)
  C00(2) =  -q * (P_center - Q_center) * pq_inv_2

  D00(0) = (0.d0, 0.d0)
  D00(1) = (0.d0, 0.d0)
  D00(2) =  -p * (Q_center - P_center) * pq_inv_2

  do i = 0, n_pt_in
    d(i) = (0.d0, 0.d0)
  enddo

  n_pt1 = n_pt_in

  !DIR$ FORCEINLINE
  call I_x1_pol_mult_cosgtos(a_x, d_x, B10, B01, B00, C00, D00, d, n_pt1, n_pt_in)
  n_pt_out = n_pt1

!  print *, ' '
!  print *, a_x, d_x
!  print *, real(B10), real(B01), real(B00), real(C00), real(D00)
!  print *, n_pt1, real(d(0:n_pt1))
!  print *, ' '

  if(n_pt1 < 0) then
    n_pt_out = -1
    do i = 0, n_pt_in
      d(i) = (0.d0, 0.d0)
    enddo
    return
  endif

end subroutine give_cpolynom_mult_center_x

! ---

subroutine I_x1_pol_mult_cosgtos(a, c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)

  BEGIN_DOC
  ! Recursive function involved in the two-electron integral
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)    :: n_pt_in, a, c
  complex*16, intent(in)    :: B_10(0:2), B_01(0:2), B_00(0:2), C_00(0:2), D_00(0:2)
  integer,    intent(inout) :: nd
  complex*16, intent(inout) :: d(0:max_dim)

  if( (c >= 0) .and. (nd >= 0) ) then

    if(a == 1) then
      call I_x1_pol_mult_a1_cosgtos(c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)
    else if(a == 2) then
      call I_x1_pol_mult_a2_cosgtos(c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)
    else if(a > 2) then
      call I_x1_pol_mult_recurs_cosgtos(a, c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)
    else  ! a == 0

      if(c == 0)then
        nd = 0
        d(0) = (1.d0, 0.d0)
        return
      endif

      call I_x2_pol_mult_cosgtos(c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)
    endif

  else

    nd = -1

  endif

end subroutine I_x1_pol_mult_cosgtos

! ---

recursive subroutine I_x1_pol_mult_recurs_cosgtos(a, c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)

  BEGIN_DOC
  ! Recursive function involved in the two-electron integral
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)    :: n_pt_in, a, c
  complex*16, intent(in)    :: B_10(0:2), B_01(0:2), B_00(0:2), C_00(0:2), D_00(0:2)
  integer,    intent(inout) :: nd
  complex*16, intent(inout) :: d(0:max_dim)

  integer                   :: nx, ix, iy, ny
  complex*16                :: X(0:max_dim)
  complex*16                :: Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X,Y

  ASSERT (a > 2)

  !DIR$ LOOP COUNT(8)
  do ix = 0, n_pt_in
    X(ix) = (0.d0, 0.d0)
  enddo

  nx = 0
  if(a == 3) then
    call I_x1_pol_mult_a1_cosgtos(c, B_10, B_01, B_00, C_00, D_00, X, nx, n_pt_in)
  elseif(a == 4) then
    call I_x1_pol_mult_a2_cosgtos(c, B_10, B_01, B_00, C_00, D_00, X, nx, n_pt_in)
  else
    ASSERT (a >= 5)
    call I_x1_pol_mult_recurs_cosgtos(a-2, c, B_10, B_01, B_00, C_00, D_00, X, nx, n_pt_in)
  endif

  !DIR$ LOOP COUNT(8)
  do ix = 0, nx
    X(ix) *= dble(a-1)
  enddo

  !DIR$ FORCEINLINE
  call multiply_cpoly(X, nx, B_10, 2, d, nd)
  nx = nd

  !DIR$ LOOP COUNT(8)
  do ix = 0, n_pt_in
    X(ix) = (0.d0, 0.d0)
  enddo

  if(c > 0) then

    if(a == 3) then
      call I_x1_pol_mult_a2_cosgtos(c-1, B_10, B_01, B_00, C_00, D_00, X, nx, n_pt_in)
    else
      ASSERT(a >= 4)
      call I_x1_pol_mult_recurs_cosgtos(a-1, c-1, B_10, B_01, B_00, C_00, D_00, X, nx, n_pt_in)
    endif

    if(c > 1) then
      !DIR$ LOOP COUNT(8)
      do ix = 0, nx
        X(ix) *= dble(c)
      enddo
    endif
    !DIR$ FORCEINLINE
    call multiply_cpoly(X, nx, B_00, 2, d, nd)

  endif

  ny = 0

  !DIR$ LOOP COUNT(8)
  do ix = 0, n_pt_in
    Y(ix) = (0.d0, 0.d0)
  enddo

  ASSERT (a > 2)

  if(a == 3) then
    call I_x1_pol_mult_a2_cosgtos(c, B_10, B_01, B_00, C_00, D_00, Y, ny, n_pt_in)
  else
    ASSERT(a >= 4)
    call I_x1_pol_mult_recurs_cosgtos(a-1, c, B_10, B_01, B_00, C_00, D_00, Y, ny, n_pt_in)
  endif

  !DIR$ FORCEINLINE
  call multiply_cpoly(Y, ny, C_00, 2, d, nd)

end subroutine I_x1_pol_mult_recurs_cosgtos

! ---

recursive subroutine I_x1_pol_mult_a1_cosgtos(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)

  BEGIN_DOC
  ! Recursive function involved in the two-electron integral
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)    :: n_pt_in, c
  complex*16, intent(in)    :: B_10(0:2), B_01(0:2), B_00(0:2), C_00(0:2), D_00(0:2)
  integer,    intent(inout) :: nd
  complex*16, intent(inout) :: d(0:max_dim)

  integer                   :: nx, ix, iy, ny
  complex*16                :: X(0:max_dim)
  complex*16                :: Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X,Y

  if( (c < 0) .or. (nd < 0) ) then
    nd = -1
    return
  endif

  nx = nd
  !DIR$ LOOP COUNT(8)
  do ix = 0, n_pt_in
    X(ix) = (0.d0, 0.d0)
  enddo
  call I_x2_pol_mult_cosgtos(c-1, B_10, B_01, B_00, C_00, D_00, X, nx, n_pt_in)

  if(c > 1) then
    !DIR$ LOOP COUNT(8)
    do ix = 0, nx
      X(ix) *= dble(c)
    enddo
  endif

  !DIR$ FORCEINLINE
  call multiply_cpoly(X, nx, B_00, 2, d, nd)

  ny = 0

  !DIR$ LOOP COUNT(8)
  do ix = 0, n_pt_in
    Y(ix) = (0.d0, 0.d0)
  enddo
  call I_x2_pol_mult_cosgtos(c, B_10, B_01, B_00, C_00, D_00, Y, ny, n_pt_in)

  !DIR$ FORCEINLINE
  call multiply_cpoly(Y, ny, C_00, 2, d, nd)

end subroutine I_x1_pol_mult_a1_cosgtos

! ---

recursive subroutine I_x1_pol_mult_a2_cosgtos(c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)

  BEGIN_DOC
  ! Recursive function involved in the two-electron integral
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)    :: n_pt_in, c
  complex*16, intent(in)    :: B_10(0:2), B_01(0:2), B_00(0:2), C_00(0:2), D_00(0:2)
  integer,    intent(inout) :: nd
  complex*16, intent(inout) :: d(0:max_dim)

  integer                   :: nx, ix, iy, ny
  complex*16                :: X(0:max_dim)
  complex*16                :: Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X,Y

  !DIR$ LOOP COUNT(8)
  do ix = 0, n_pt_in
    X(ix) = (0.d0, 0.d0)
  enddo

  nx = 0
  call I_x2_pol_mult_cosgtos(c, B_10, B_01, B_00, C_00, D_00, X, nx, n_pt_in)

  !DIR$ FORCEINLINE
  call multiply_cpoly(X, nx, B_10, 2, d, nd)

  nx = nd
  !DIR$ LOOP COUNT(8)
  do ix = 0, n_pt_in
    X(ix) = (0.d0, 0.d0)
  enddo

  !DIR$ FORCEINLINE
  call I_x1_pol_mult_a1_cosgtos(c-1, B_10, B_01, B_00, C_00, D_00, X, nx, n_pt_in)

  if (c>1) then
    !DIR$ LOOP COUNT(8)
    do ix = 0, nx
      X(ix) *= dble(c)
    enddo
  endif

  !DIR$ FORCEINLINE
  call multiply_cpoly(X, nx, B_00, 2, d, nd)

  ny = 0
  !DIR$ LOOP COUNT(8)
  do ix = 0, n_pt_in
    Y(ix) = 0.d0
  enddo
  !DIR$ FORCEINLINE
  call I_x1_pol_mult_a1_cosgtos(c, B_10, B_01, B_00, C_00, D_00, Y, ny, n_pt_in)

  !DIR$ FORCEINLINE
  call multiply_cpoly(Y, ny, C_00, 2, d, nd)

end subroutine I_x1_pol_mult_a2_cosgtos

! ---

recursive subroutine I_x2_pol_mult_cosgtos(c, B_10, B_01, B_00, C_00, D_00, d, nd, dim)

  BEGIN_DOC
  ! Recursive function involved in the two-electron integral
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)    :: dim, c
  complex*16, intent(in)    :: B_10(0:2), B_01(0:2), B_00(0:2), C_00(0:2), D_00(0:2)
  integer,    intent(inout) :: nd
  complex*16, intent(inout) :: d(0:max_dim)

  integer                   :: i
  integer                   :: nx, ix, ny
  complex*16                :: X(0:max_dim), Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X, Y

  select case (c)

    case (0)
      nd = 0
      d(0) = (1.d0, 0.d0)
      return

    case (:-1)
      nd = -1
      return

    case (1)
      nd = 2
      d(0) = D_00(0)
      d(1) = D_00(1)
      d(2) = D_00(2)
      return

    case (2)
      nd = 2
      d(0) = B_01(0)
      d(1) = B_01(1)
      d(2) = B_01(2)

      ny = 2
      Y(0) = D_00(0)
      Y(1) = D_00(1)
      Y(2) = D_00(2)

      !DIR$ FORCEINLINE
      call multiply_cpoly(Y, ny, D_00, 2, d, nd)
      return

    case default

      !DIR$ LOOP COUNT(6)
      do ix = 0, c+c
        X(ix) = (0.d0, 0.d0)
      enddo
      nx = 0
      call I_x2_pol_mult_cosgtos(c-2, B_10, B_01, B_00, C_00, D_00, X, nx, dim)

      !DIR$ LOOP COUNT(6)
      do ix = 0, nx
        X(ix) *= dble(c-1)
      enddo

      !DIR$ FORCEINLINE
      call multiply_cpoly(X, nx, B_01, 2, d, nd)

      ny = 0
      !DIR$ LOOP COUNT(6)
      do ix = 0, c+c
        Y(ix) = 0.d0
      enddo
      call I_x2_pol_mult_cosgtos(c-1, B_10, B_01, B_00, C_00, D_00, Y, ny, dim)

      !DIR$ FORCEINLINE
      call multiply_cpoly(Y, ny, D_00, 2, d, nd)

  end select

end subroutine I_x2_pol_mult_cosgtos

! ---


