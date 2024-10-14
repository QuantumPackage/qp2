
! ---

subroutine deb_ao_2eint_cgtos(i, j, k, l)

  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in)        :: i, j, k, l
                             
  integer                    :: p, q, r, s
  integer                    :: num_i, num_j, num_k, num_l, dim1, I_power(3), J_power(3), K_power(3), L_power(3)
  integer                    :: iorder_p1(3), iorder_p2(3), iorder_q1(3), iorder_q2(3)
  complex*16                 :: I_center(3), J_center(3), K_center(3), L_center(3)
  complex*16                 :: expo1, expo2, expo3, expo4
  complex*16                 :: P1_center(3), pp1
  complex*16                 :: P2_center(3), pp2
  complex*16                 :: Q1_center(3), qq1
  complex*16                 :: Q2_center(3), qq2



  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)

  if(num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k) then

    !print*, ao_prim_num(i), ao_prim_num(j), ao_prim_num(k), ao_prim_num(l)

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
      expo1 = ao_expo_cgtos_ord_transp(p,i) 
      !print*, "expo1 = ", expo1
      !print*, "center1 = ", I_center

      do q = 1, ao_prim_num(j)
        expo2 = ao_expo_cgtos_ord_transp(q,j) 
        !print*, "expo2 = ", expo2
        !print*, "center2 = ", J_center

        pp1 = expo1 + expo2
        P1_center(1:3) = (expo1 * I_center(1:3) + expo2 * J_center(1:3)) / pp1
        iorder_p1(1:3) = I_power(1:3) + J_power(1:3)

        pp2 = conjg(expo1) + expo2
        P2_center(1:3) = (conjg(expo1) * I_center(1:3) + expo2 * J_center(1:3)) / pp2
        iorder_p2(1:3) = I_power(1:3) + J_power(1:3)

        do r = 1, ao_prim_num(k)
          expo3 = ao_expo_cgtos_ord_transp(r,k) 
          !print*, "expo3 = ", expo3
          !print*, "center3 = ", K_center

          do s = 1, ao_prim_num(l)
            expo4 = ao_expo_cgtos_ord_transp(s,l) 
            !print*, "expo4 = ", expo4
            !print*, "center4 = ", L_center

            qq1 = expo3 + expo4
            Q1_center(1:3) = (expo3 * K_center(1:3) + expo4 * L_center(1:3)) / qq1
            iorder_q1(1:3) = K_power(1:3) + L_power(1:3)

            qq2 = conjg(expo3) + expo4
            Q2_center(1:3) = (conjg(expo3) * K_center(1:3) + expo4 * L_center(1:3)) / qq2
            iorder_q2(1:3) = K_power(1:3) + L_power(1:3)

            call deb_cboys(P1_center, pp1, iorder_p1, Q1_center, qq1, iorder_q1)
            call deb_cboys(P1_center, pp1, iorder_p1, Q2_center, qq2, iorder_q2)
            call deb_cboys(P2_center, pp2, iorder_p2, Q1_center, qq1, iorder_q1)
            call deb_cboys(P2_center, pp2, iorder_p2, Q2_center, qq2, iorder_q2)
            call deb_cboys(conjg(P2_center), conjg(pp2), iorder_p2, Q1_center, qq1, iorder_q1)
            call deb_cboys(conjg(P2_center), conjg(pp2), iorder_p2, Q2_center, qq2, iorder_q2)
            call deb_cboys(conjg(P1_center), conjg(pp1), iorder_p1, Q1_center, qq1, iorder_q1)
            call deb_cboys(conjg(P1_center), conjg(pp1), iorder_p1, Q2_center, qq2, iorder_q2)
          enddo ! s
        enddo ! r
      enddo ! q
    enddo ! p

  endif ! same centers

  return
end

! ---

subroutine deb_cboys(P_center, p, iorder_p, Q_center, q, iorder_q)


  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in) :: iorder_p(3), iorder_q(3)
  complex*16, intent(in) :: P_center(3), p
  complex*16, intent(in) :: Q_center(3), q

  integer                :: iorder, n
  complex*16             :: dist, rho
  complex*16             :: int1, int2

  complex*16, external   :: crint_2


  dist = (P_center(1) - Q_center(1)) * (P_center(1) - Q_center(1)) &
       + (P_center(2) - Q_center(2)) * (P_center(2) - Q_center(2)) &
       + (P_center(3) - Q_center(3)) * (P_center(3) - Q_center(3))
  rho = dist * p * q / (p + q)

  if(abs(rho) .lt. 1d-15) return

  iorder = 2*iorder_p(1)+2*iorder_q(1) + 2*iorder_p(2)+2*iorder_q(2) + 2*iorder_p(3)+2*iorder_q(3)
  n = shiftr(iorder, 1)
  
  !write(33,*) n, real(rho), aimag(rho)
  !print*, n, real(rho), aimag(rho)

  int1 = crint_2(n, rho)
  call crint_quad_12(n, rho, 1000000, int2)

  if(abs(int1 - int2) .gt. 1d-5) then
    print*, ' important error found: '
    print*, p!, P_center
    print*, q!, Q_center
    print*, dist
    print*, " n, tho = ", n, real(rho), aimag(rho)
    print*, real(int1), real(int2), dabs(real(int1-int2))
    print*, aimag(int1), aimag(int2), dabs(aimag(int1-int2))
    stop
  endif

end

! ---

