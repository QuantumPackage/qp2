
! ---

complex*16 function crint_1(n, rho)

  implicit none
  include 'constants.include.F'

  integer,    intent(in) :: n
  complex*16, intent(in) :: rho

  integer                :: i, mmax
  double precision       :: rho_mod
  double precision       :: tmp
  complex*16             :: rho_inv, rho_exp

  complex*16             :: crint_smallz

  rho_mod = zabs(rho)

  if(rho_mod < 3.5d0) then

    if(rho_mod .lt. 0.35d0) then

      select case(n)
      case(0)
        crint_1 = (((((((((1.3122532963802805073d-08 * rho &
                - 1.450385222315046877d-07) * rho &
                + 1.458916900093370682d-06) * rho &
                - 0.132275132275132275d-04) * rho &
                + 0.106837606837606838d-03) * rho &
                - 0.757575757575757576d-03) * rho &
                + 0.462962962962962963d-02) * rho &
                - 0.238095238095238095d-01) * rho &
                + 0.10000000000000000000d0) * rho &
                - 0.33333333333333333333d0) * rho &
                + 1.0d0
      case(1)
        crint_1 = (((((((((1.198144314086343d-08 * rho &
                - 1.312253296380281d-07) * rho &
                + 1.305346700083542d-06) * rho &
                - 1.167133520074696d-05) * rho &
                + 9.259259259259259d-05) * rho &
                - 6.410256410256410d-04) * rho &
                + 3.787878787878788d-03) * rho &
                - 1.851851851851852d-02) * rho &
                + 7.142857142857142d-02) * rho &
                - 2.000000000000000d-01) * rho &
                + 3.333333333333333d-01
      case(2)
        crint_1 = (((((((((1.102292768959436d-08 * rho &
                - 1.198144314086343d-07) * rho &
                + 1.181027966742252d-06) * rho &
                - 1.044277360066834d-05) * rho &
                + 8.169934640522875d-05) * rho &
                - 5.555555555555556d-04) * rho &
                + 3.205128205128205d-03) * rho &
                - 1.515151515151515d-02) * rho &
                + 5.555555555555555d-02) * rho &
                - 1.428571428571428d-01) * rho &
                + 2.000000000000000d-01 
      case(3)
        crint_1 = (((((((((1.020641452740218d-08 * rho &
                - 1.102292768959436d-07) * rho &
                + 1.078329882677709d-06) * rho &
                - 9.448223733938020d-06) * rho &
                + 7.309941520467836d-05) * rho &
                - 4.901960784313725d-04) * rho &
                + 2.777777777777778d-03) * rho &
                - 1.282051282051282d-02) * rho &
                + 4.545454545454546d-02) * rho &
                - 1.111111111111111d-01) * rho &
                + 1.428571428571428d-01 
      case default
        tmp = dble(n + n + 1)
        crint_1 = (((((((((2.755731922398589d-07 * rho / (tmp + 20.d0) &
                - 2.755731922398589d-06 / (tmp + 18.d0)) * rho &
                + 2.480158730158730d-05 / (tmp + 16.d0)) * rho &
                - 1.984126984126984d-04 / (tmp + 14.d0)) * rho &
                + 1.388888888888889d-03 / (tmp + 12.d0)) * rho &
                - 8.333333333333333d-03 / (tmp + 10.d0)) * rho &
                + 4.166666666666666d-02 / (tmp +  8.d0)) * rho &
                - 1.666666666666667d-01 / (tmp +  6.d0)) * rho &
                + 5.000000000000000d-01 / (tmp +  4.d0)) * rho &
                - 1.000000000000000d+00 / (tmp +  2.d0)) * rho &
                + 1.0d0 / tmp
      end select 

    else

      crint_1 = crint_smallz(n, rho)

    endif

  else

    rho_exp = 0.5d0 * zexp(-rho)
    rho_inv = (1.d0, 0.d0) / rho

    call zboysfun00_1(rho, crint_1)
    do i = 1, n
      crint_1 = ((dble(i) - 0.5d0) * crint_1 - rho_exp) * rho_inv
    enddo

  endif

  return
end

! ---

subroutine crint_1_vec(n_max, rho, vals)

  implicit none
  include 'constants.include.F'

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: vals(0:n_max)

  integer                :: n
  double precision       :: rho_mod
  double precision       :: tmp
  complex*16             :: rho_inv, rho_exp

  complex*16             :: crint_smallz

  rho_mod = zabs(rho)

  if(rho_mod < 3.5d0) then

    if(rho_mod .lt. 0.35d0) then

      vals(0) = (((((((((1.3122532963802805073d-08 * rho &
              - 1.450385222315046877d-07) * rho &
              + 1.458916900093370682d-06) * rho &
              - 0.132275132275132275d-04) * rho &
              + 0.106837606837606838d-03) * rho &
              - 0.757575757575757576d-03) * rho &
              + 0.462962962962962963d-02) * rho &
              - 0.238095238095238095d-01) * rho &
              + 0.10000000000000000000d0) * rho &
              - 0.33333333333333333333d0) * rho &
              + 1.0d0

      if(n_max > 0) then

        vals(1) = (((((((((1.198144314086343d-08 * rho &
                - 1.312253296380281d-07) * rho &
                + 1.305346700083542d-06) * rho &
                - 1.167133520074696d-05) * rho &
                + 9.259259259259259d-05) * rho &
                - 6.410256410256410d-04) * rho &
                + 3.787878787878788d-03) * rho &
                - 1.851851851851852d-02) * rho &
                + 7.142857142857142d-02) * rho &
                - 2.000000000000000d-01) * rho &
                + 3.333333333333333d-01

        if(n_max > 1) then

          vals(2) = (((((((((1.102292768959436d-08 * rho &
                  - 1.198144314086343d-07) * rho &
                  + 1.181027966742252d-06) * rho &
                  - 1.044277360066834d-05) * rho &
                  + 8.169934640522875d-05) * rho &
                  - 5.555555555555556d-04) * rho &
                  + 3.205128205128205d-03) * rho &
                  - 1.515151515151515d-02) * rho &
                  + 5.555555555555555d-02) * rho &
                  - 1.428571428571428d-01) * rho &
                  + 2.000000000000000d-01 

          if(n_max > 2) then

            vals(3) = (((((((((1.020641452740218d-08 * rho &
                    - 1.102292768959436d-07) * rho &
                    + 1.078329882677709d-06) * rho &
                    - 9.448223733938020d-06) * rho &
                    + 7.309941520467836d-05) * rho &
                    - 4.901960784313725d-04) * rho &
                    + 2.777777777777778d-03) * rho &
                    - 1.282051282051282d-02) * rho &
                    + 4.545454545454546d-02) * rho &
                    - 1.111111111111111d-01) * rho &
                    + 1.428571428571428d-01 

            do n = 4, n_max
              tmp = dble(n + n + 1)
              vals(n) = (((((((((2.755731922398589d-07 * rho / (tmp + 20.d0) &
                      - 2.755731922398589d-06 / (tmp + 18.d0)) * rho &
                      + 2.480158730158730d-05 / (tmp + 16.d0)) * rho &
                      - 1.984126984126984d-04 / (tmp + 14.d0)) * rho &
                      + 1.388888888888889d-03 / (tmp + 12.d0)) * rho &
                      - 8.333333333333333d-03 / (tmp + 10.d0)) * rho &
                      + 4.166666666666666d-02 / (tmp +  8.d0)) * rho &
                      - 1.666666666666667d-01 / (tmp +  6.d0)) * rho &
                      + 5.000000000000000d-01 / (tmp +  4.d0)) * rho &
                      - 1.000000000000000d+00 / (tmp +  2.d0)) * rho &
                      + 1.0d0 / tmp
            enddo

          endif ! n_max > 2
        endif ! n_max > 1
      endif ! n_max > 0

    else

      call crint_smallz_vec(n_max, rho, vals)

    endif

  else

    rho_exp = 0.5d0 * zexp(-rho)
    rho_inv = (1.d0, 0.d0) / rho

    call zboysfun00_1(rho, vals(0))
    do n = 1, n_max
      vals(n) = ((dble(n) - 0.5d0) * vals(n-1) - rho_exp) * rho_inv
    enddo

  endif

  return
end

! ---

complex*16 function crint_smallz(n, rho)

  BEGIN_DOC
  ! Standard version of rint
  END_DOC

  implicit none
  integer,    intent(in)      :: n
  complex*16, intent(in)      :: rho

  integer,          parameter :: kmax = 40
  double precision, parameter :: eps = 1.d-10

  integer                     :: k
  double precision            :: delta_mod
  complex*16                  :: rho_k, ct, delta_k

  ct           = 0.5d0 * zexp(-rho) * gamma(dble(n) + 0.5d0)
  crint_smallz = ct / gamma(dble(n) + 1.5d0)

  rho_k = (1.d0, 0.d0)
  do k = 1, kmax

    rho_k        = rho_k * rho
    delta_k      = ct * rho_k / gamma(dble(n+k) + 1.5d0)
    crint_smallz = crint_smallz + delta_k

    delta_mod = dsqrt(real(delta_k)*real(delta_k) + aimag(delta_k)*aimag(delta_k))
    if(delta_mod .lt. eps) return
  enddo

  if(delta_mod > eps) then
    write(*,*) ' pb in crint_smallz !'
    write(*,*) ' n, rho = ', n, rho
    write(*,*) ' value = ', crint_smallz
    write(*,*) ' delta_mod = ', delta_mod
    !stop 1
  endif

end

! ---

complex*16 function crint_2(n, rho)

  implicit none

  integer,    intent(in) :: n
  complex*16, intent(in) :: rho

  double precision       :: tmp, abs_rho
  complex*16             :: vals(0:n)

  complex*16, external   :: crint_smallz

  abs_rho = zabs(rho)

  if(abs_rho < 3.5d0) then

    if(abs_rho .lt. 0.35d0) then

      select case(n)
      case(0)
        crint_2 = (((((((((1.3122532963802805073d-08 * rho &
                - 1.450385222315046877d-07) * rho &
                + 1.458916900093370682d-06) * rho &
                - 0.132275132275132275d-04) * rho &
                + 0.106837606837606838d-03) * rho &
                - 0.757575757575757576d-03) * rho &
                + 0.462962962962962963d-02) * rho &
                - 0.238095238095238095d-01) * rho &
                + 0.10000000000000000000d0) * rho &
                - 0.33333333333333333333d0) * rho &
                + 1.0d0
      case(1)
        crint_2 = (((((((((1.198144314086343d-08 * rho &
                - 1.312253296380281d-07) * rho &
                + 1.305346700083542d-06) * rho &
                - 1.167133520074696d-05) * rho &
                + 9.259259259259259d-05) * rho &
                - 6.410256410256410d-04) * rho &
                + 3.787878787878788d-03) * rho &
                - 1.851851851851852d-02) * rho &
                + 7.142857142857142d-02) * rho &
                - 2.000000000000000d-01) * rho &
                + 3.333333333333333d-01
      case(2)
        crint_2 = (((((((((1.102292768959436d-08 * rho &
                - 1.198144314086343d-07) * rho &
                + 1.181027966742252d-06) * rho &
                - 1.044277360066834d-05) * rho &
                + 8.169934640522875d-05) * rho &
                - 5.555555555555556d-04) * rho &
                + 3.205128205128205d-03) * rho &
                - 1.515151515151515d-02) * rho &
                + 5.555555555555555d-02) * rho &
                - 1.428571428571428d-01) * rho &
                + 2.000000000000000d-01 
      case(3)
        crint_2 = (((((((((1.020641452740218d-08 * rho &
                - 1.102292768959436d-07) * rho &
                + 1.078329882677709d-06) * rho &
                - 9.448223733938020d-06) * rho &
                + 7.309941520467836d-05) * rho &
                - 4.901960784313725d-04) * rho &
                + 2.777777777777778d-03) * rho &
                - 1.282051282051282d-02) * rho &
                + 4.545454545454546d-02) * rho &
                - 1.111111111111111d-01) * rho &
                + 1.428571428571428d-01 
      case default
        tmp = dble(n + n + 1)
        crint_2 = (((((((((2.755731922398589d-07 * rho / (tmp + 20.d0) &
                - 2.755731922398589d-06 / (tmp + 18.d0)) * rho &
                + 2.480158730158730d-05 / (tmp + 16.d0)) * rho &
                - 1.984126984126984d-04 / (tmp + 14.d0)) * rho &
                + 1.388888888888889d-03 / (tmp + 12.d0)) * rho &
                - 8.333333333333333d-03 / (tmp + 10.d0)) * rho &
                + 4.166666666666666d-02 / (tmp +  8.d0)) * rho &
                - 1.666666666666667d-01 / (tmp +  6.d0)) * rho &
                + 5.000000000000000d-01 / (tmp +  4.d0)) * rho &
                - 1.000000000000000d+00 / (tmp +  2.d0)) * rho &
                + 1.0d0 / tmp
      end select 

    else

      crint_2 = crint_smallz(n, rho)

    endif

  else

    if(real(rho) .ge. 0.d0) then

      call zboysfun(n, rho, vals)
      crint_2 = vals(n)

    else

      call zboysfunnrp(n, rho, vals)
      crint_2 = vals(n) * zexp(-rho)

    endif
  endif

  return
end

! ---

subroutine zboysfun(n_max, x, vals)

  BEGIN_DOC
  !
  ! Computes values of the Boys function for n = 0, 1, ..., n_max
  ! for a complex valued argument
  !
  ! Input: x --- argument, complex*16, Re(x) >= 0
  ! Output: vals  --- values of the Boys function, n = 0, 1, ..., n_max
  !
  ! Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
  ! https://doi.org/10.1063/5.0062444
  !
  END_DOC

  implicit none

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: x
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n
  complex*16              :: yy, x_inv

  call zboysfun00_2(x, vals(0))

  yy = 0.5d0 * zexp(-x)
  x_inv = (1.d0, 0.d0) / x
  do n = 1, n_max
    vals(n) = ((dble(n) - 0.5d0) * vals(n-1) - yy) * x_inv
  enddo

  return
end

! ---

subroutine zboysfunnrp(n_max, x, vals)

  BEGIN_DOC
  !
  ! Computes values of e^z F(n,z) for n = 0, 1, ..., n_max
  ! (where F(n,z) are the Boys functions)
  ! for a complex valued argument WITH NEGATIVE REAL PART
  !
  ! Input: x  --- argument, complex *16 Re(x)<=0
  ! Output: vals  --- values of e^z F(n,z), n = 0, 1, ..., n_max
  !
  ! Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
  ! https://doi.org/10.1063/5.0062444
  !
  END_DOC

  implicit none

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: x
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n
  complex*16              :: x_inv

  call zboysfun00nrp(x, vals(0))

  x_inv = (1.d0, 0.d0) / x
  do n = 1, n_max
    vals(n) = ((dble(n) - 0.5d0) * vals(n-1) - 0.5d0) * x_inv
  enddo

  return
end

! ---

complex*16 function crint_sum(n_pt_out, rho, d1)

  implicit none

  integer,    intent(in)  :: n_pt_out
  complex*16, intent(in)  :: rho, d1(0:n_pt_out)
                          
  integer                 :: i
  integer                 :: n_max

  complex*16, allocatable :: vals(:)

  n_max = shiftr(n_pt_out, 1)
  allocate(vals(0:n_max))

  call crint_1_vec(n_max, rho, vals)
  !call crint_2_vec(n_max, rho, vals)
  ! FOR DEBUG
  !call crint_quad_12_vec(n_max, rho, vals)

  crint_sum = d1(0) * vals(0)
  do i = 2, n_pt_out, 2
    crint_sum += d1(i) * vals(shiftr(i, 1))
  enddo

  deallocate(vals)

  return
end

! ---

subroutine crint_2_vec(n_max, rho, vals)

  implicit none

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n
  double precision        :: tmp, abs_rho
  complex*16              :: erho


  abs_rho = abs(rho)

  if(abs_rho < 3.5d0) then

    if(abs_rho .lt. 0.35d0) then

      vals(0) = (((((((((1.3122532963802805073d-08 * rho &
              - 1.450385222315046877d-07) * rho &
              + 1.458916900093370682d-06) * rho &
              - 0.132275132275132275d-04) * rho &
              + 0.106837606837606838d-03) * rho &
              - 0.757575757575757576d-03) * rho &
              + 0.462962962962962963d-02) * rho &
              - 0.238095238095238095d-01) * rho &
              + 0.10000000000000000000d0) * rho &
              - 0.33333333333333333333d0) * rho &
              + 1.0d0

      if(n_max > 0) then

        vals(1) = (((((((((1.198144314086343d-08 * rho &
                - 1.312253296380281d-07) * rho &
                + 1.305346700083542d-06) * rho &
                - 1.167133520074696d-05) * rho &
                + 9.259259259259259d-05) * rho &
                - 6.410256410256410d-04) * rho &
                + 3.787878787878788d-03) * rho &
                - 1.851851851851852d-02) * rho &
                + 7.142857142857142d-02) * rho &
                - 2.000000000000000d-01) * rho &
                + 3.333333333333333d-01

        if(n_max > 1) then

          vals(2) = (((((((((1.102292768959436d-08 * rho &
                  - 1.198144314086343d-07) * rho &
                  + 1.181027966742252d-06) * rho &
                  - 1.044277360066834d-05) * rho &
                  + 8.169934640522875d-05) * rho &
                  - 5.555555555555556d-04) * rho &
                  + 3.205128205128205d-03) * rho &
                  - 1.515151515151515d-02) * rho &
                  + 5.555555555555555d-02) * rho &
                  - 1.428571428571428d-01) * rho &
                  + 2.000000000000000d-01 

          if(n_max > 2) then

            vals(3) = (((((((((1.020641452740218d-08 * rho &
                    - 1.102292768959436d-07) * rho &
                    + 1.078329882677709d-06) * rho &
                    - 9.448223733938020d-06) * rho &
                    + 7.309941520467836d-05) * rho &
                    - 4.901960784313725d-04) * rho &
                    + 2.777777777777778d-03) * rho &
                    - 1.282051282051282d-02) * rho &
                    + 4.545454545454546d-02) * rho &
                    - 1.111111111111111d-01) * rho &
                    + 1.428571428571428d-01 

            do n = 4, n_max
              tmp = dble(n + n + 1)
              vals(n) = (((((((((2.755731922398589d-07 * rho / (tmp + 20.d0) &
                      - 2.755731922398589d-06 / (tmp + 18.d0)) * rho &
                      + 2.480158730158730d-05 / (tmp + 16.d0)) * rho &
                      - 1.984126984126984d-04 / (tmp + 14.d0)) * rho &
                      + 1.388888888888889d-03 / (tmp + 12.d0)) * rho &
                      - 8.333333333333333d-03 / (tmp + 10.d0)) * rho &
                      + 4.166666666666666d-02 / (tmp +  8.d0)) * rho &
                      - 1.666666666666667d-01 / (tmp +  6.d0)) * rho &
                      + 5.000000000000000d-01 / (tmp +  4.d0)) * rho &
                      - 1.000000000000000d+00 / (tmp +  2.d0)) * rho &
                      + 1.0d0 / tmp
            enddo

          endif ! n_max > 2
        endif ! n_max > 1
      endif ! n_max > 0

    else

      call crint_smallz_vec(n_max, rho, vals)

    endif

  else

    if(real(rho) .ge. 0.d0) then

      call zboysfun(n_max, rho, vals)

    else

      call zboysfunnrp(n_max, rho, vals)
      erho = zexp(-rho)
      do n = 0, n_max
        vals(n) = vals(n) * erho
      enddo

    endif

  endif

  return
end

! ---

subroutine crint_smallz_vec(n_max, rho, vals)

  BEGIN_DOC
  ! Standard version of rint
  END_DOC

  implicit none
  integer,    intent(in)      :: n_max
  complex*16, intent(in)      :: rho
  complex*16, intent(out)     :: vals(0:n_max)

  integer,          parameter :: kmax = 40
  double precision, parameter :: eps = 1.d-10

  integer                     :: k, n
  complex*16                  :: ct, delta_k
  complex*16                  :: rhoe
  complex*16, allocatable     :: rho_k(:)


  allocate(rho_k(0:kmax))

  rho_k(0) = (1.d0, 0.d0)
  do k = 1, kmax
    rho_k(k) = rho_k(k-1) * rho
  enddo

  rhoe = 0.5d0 * zexp(-rho)

  do n = 0, n_max

    ct = rhoe * gamma(dble(n) + 0.5d0)
    vals(n) = ct / gamma(dble(n) + 1.5d0)
  
    do k = 1, kmax
      delta_k = ct * rho_k(k) / gamma(dble(n+k) + 1.5d0)
      vals(n) += delta_k
      if(abs(delta_k) .lt. eps) then
        exit
      endif
    enddo
  
    if(abs(delta_k) > eps) then
      write(*,*) ' pb in crint_smallz_vec !'
      write(*,*) ' n, rho = ', n, rho
      write(*,*) ' value = ', vals(n)
      write(*,*) ' |delta_k| = ', abs(delta_k)
    endif
  enddo

  deallocate(rho_k)

  return
end

! ---

subroutine crint_quad_1(n, rho, n_quad, crint_quad)

  implicit none

  integer,    intent(in)  :: n, n_quad
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: crint_quad

  integer                 :: i_quad
  double precision        :: tmp_inv, tmp0, tmp1, tmp2
  double precision        :: coef(0:3) = (/14.d0, 32.d0, 12.d0, 32.d0 /)

  tmp_inv = 1.d0 / dble(n_quad)

  crint_quad = 7.d0 * zexp(-rho)

  tmp0 = 0.d0
  select case (n)

    case (0)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1)
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (1)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1) * tmp1
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (2)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1) * tmp1 * tmp1
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (3)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1) * tmp1 * tmp1 * tmp1
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (4)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        tmp2 = tmp1 * tmp1
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1) * tmp2 * tmp2
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case default
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1) * tmp1**n
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv
  end select

end

! ---

subroutine crint_quad_2(n, rho, n_quad, crint_quad)

  implicit none

  integer,    intent(in)  :: n, n_quad
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: crint_quad

  integer                 :: i_quad
  double precision        :: tmp_inv, tmp0, tmp1, tmp2
  double precision        :: coef(0:3) = (/14.d0, 32.d0, 12.d0, 32.d0 /)
  complex*16              :: rhoc, rhoe

  tmp_inv = 1.d0 / dble(n_quad)

  crint_quad = 7.d0 * zexp(-rho)

  tmp0 = 0.d0
  rhoc = zexp(-rho*tmp_inv)
  rhoe = (1.d0, 0.d0)
  select case (n)

    case (0)
      !do i_quad = 1, n_quad - 1
      !  tmp0 = tmp0 + tmp_inv
      !  rhoe = rhoe * rhoc
      !  tmp1 = (rhoe - 1.d0) / dsqrt(tmp0)
      !  crint_quad = crint_quad + coef(iand(i_quad, 3)) * tmp1
      !enddo
      !crint_quad = 1.d0 + crint_quad * 0.022222222222222223d0 * tmp_inv
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe / dsqrt(tmp0)
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (1)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (2)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (3)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0 * tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (4)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        tmp2 = tmp1 * tmp1 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp2
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case default
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0**n / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

  end select

end

! ---

subroutine crint_quad_12(n, rho, n_quad, crint_quad)

  implicit none

  integer,    intent(in)  :: n, n_quad
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: crint_quad

  integer                 :: i_quad
  double precision        :: tmp_inv, tmp0, tmp1, tmp2
  double precision        :: coef(0:3) = (/14.d0, 32.d0, 12.d0, 32.d0 /)
  complex*16              :: rhoc, rhoe

  tmp_inv = 1.d0 / dble(n_quad)

  crint_quad = 7.d0 * zexp(-rho)

  tmp0 = 0.d0
  rhoc = zexp(-rho*tmp_inv)
  rhoe = (1.d0, 0.d0)
  select case (n)

    case (0)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1)
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (1)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (2)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (3)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0 * tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (4)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        tmp2 = tmp1 * tmp1 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp2
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case default
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0**n / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

  end select

end

! ---

subroutine crint_quad_12_vec(n_max, rho, vals)

  implicit none

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n

  do n = 0, n_max
    call crint_quad_12(n, rho, 10000000, vals(n))
  enddo

  return
end

! ---

