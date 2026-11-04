
! ---

complex*16 function crint_1(n, rho)

  BEGIN_DOC
  ! Complex Boys function $F_n(\rho) = \int_0^1 t^{2n} e^{-\rho t^2} dt$.
  !
  ! Series (crint_smallz) if $|\rho| < 3.5$ or $n > |\rho|$ (no cancellation
  ! there), else $F_0$ from erf and upward recursion (stable for
  ! $n \le |\rho|$).
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,    intent(in) :: n
  complex*16, intent(in) :: rho

  integer                :: i
  double precision       :: rho_mod
  complex*16             :: rho_inv, rho_exp

  complex*16             :: crint_smallz

  rho_mod = abs(rho)

  if((rho_mod < 3.5d0) .or. (dble(n) > rho_mod)) then

    crint_1 = crint_smallz(n, rho)

  else

    rho_exp = 0.5d0 * exp(-rho)
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

  BEGIN_DOC
  ! Stores crint_1(n, rho) for n = 0, 1, ..., n_max in vals(0:n_max),
  ! where crint_1(n, rho) = $\int_0^1 t^{2n} \exp(-\rho t^2) dt$ for
  ! a complex argument rho.
  !
  ! $|\rho| < 3.5$: series for $F_{n_{max}}$ and downward recursion
  ! (crint_smallz_vec). Otherwise: $F_0$ from erf and upward recursion for
  ! $n \le |\rho|$, and for $|\rho| < n \le n_{max}$ the series for
  ! $F_{n_{max}}$ and downward recursion. Each recursion is used only where it
  ! is stable, and the series only where its terms decrease from the start.
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n, n_up
  double precision        :: rho_mod
  complex*16              :: rho_inv, rho_exp

  complex*16              :: crint_smallz

  rho_mod = abs(rho)

  if(rho_mod < 3.5d0) then

    call crint_smallz_vec(n_max, rho, vals)

  else

    rho_exp = 0.5d0 * exp(-rho)
    rho_inv = (1.d0, 0.d0) / rho

    n_up = min(n_max, int(rho_mod))

    call zboysfun00_1(rho, vals(0))
    do n = 1, n_up
      vals(n) = ((dble(n) - 0.5d0) * vals(n-1) - rho_exp) * rho_inv
    enddo

    if(n_up < n_max) then
      vals(n_max) = crint_smallz(n_max, rho)
      do n = n_max, n_up + 2, -1
        vals(n-1) = (2.d0 * (rho * vals(n) + rho_exp)) / dble(n + n - 1)
      enddo
    endif

  endif

  return
end

! ---

complex*16 function crint_smallz(n, rho)

  BEGIN_DOC
  ! Complex Boys function by the series
  ! $F_n(\rho) = e^{-\rho} \sum_{k \ge 0}
  !   \frac{(2\rho)^k}{(2n+1)(2n+3)\cdots(2n+2k+1)}$.
  ! Accurate when $|\rho| < 3.5$ or $n > |\rho|$.
  END_DOC

  implicit none
  integer,    intent(in)      :: n
  complex*16, intent(in)      :: rho

  integer,          parameter :: kmax = 200
  double precision, parameter :: eps2 = 1.d-34 ! (relative tolerance 1e-17)^2

  integer                     :: k
  double precision            :: d2, s2
  complex*16                  :: rho2, delta_k

  rho2         = rho + rho
  delta_k      = (1.d0, 0.d0) / dble(n + n + 1)
  crint_smallz = delta_k

  do k = 1, kmax
    delta_k      = delta_k * (rho2 * (1.d0 / dble(n + n + k + k + 1)))
    crint_smallz = crint_smallz + delta_k
    d2 = real(delta_k)**2 + aimag(delta_k)**2
    s2 = real(crint_smallz)**2 + aimag(crint_smallz)**2
    if(d2 < eps2 * s2) exit
  enddo

  if(d2 >= eps2 * s2) then
    write(*,*) ' pb in crint_smallz !'
    write(*,*) ' n, rho = ', n, rho
    write(*,*) ' value = ', crint_smallz
    write(*,*) ' |delta_k| = ', abs(delta_k)
  endif

  crint_smallz = crint_smallz * exp(-rho)

end

! ---

complex*16 function crint_2(n, rho)

  implicit none

  integer,    intent(in) :: n
  complex*16, intent(in) :: rho

  double precision       :: tmp, abs_rho
  complex*16             :: vals(0:n)

  complex*16, external   :: crint_smallz

  abs_rho = abs(rho)

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
      crint_2 = vals(n) * exp(-rho)

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

  yy = 0.5d0 * exp(-x)
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
  ! Input: x  --- argument, complex \*16 Re(x)<=0
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

complex*16 function crint_sum_shift(n_pt_out, rho, d1, log_fact)

  BEGIN_DOC
  ! $e^{\rm log\_fact} \sum_{i\ \rm even} d_1(i) F_{i/2}(\rho)$, without
  ! forming $e^{\rm log\_fact}$ and $F_n(\rho)$ separately: for
  ! $\Re\rho < 0$ the scaled functions $e^{\rho} F_n(\rho)$ are multiplied
  ! by $e^{{\rm log\_fact} - \rho}$. With complex exponents $\Re\rho$ can be
  ! below $-709$, where $F_n(\rho)$ overflows and $e^{\rm log\_fact}$
  ! underflows, although their product is moderate.
  END_DOC

  implicit none

  integer,    intent(in)  :: n_pt_out
  complex*16, intent(in)  :: rho, d1(0:n_pt_out), log_fact

  integer                 :: i, n_max
  complex*16              :: accu

  complex*16, allocatable :: vals(:)

  n_max = shiftr(n_pt_out, 1)
  allocate(vals(0:n_max))

  if(real(rho) >= 0.d0) then
    call crint_1_vec(n_max, rho, vals)
  else
    call crint_1_vec_scaled(n_max, rho, vals)
  endif

  accu = d1(0) * vals(0)
  do i = 2, n_pt_out, 2
    accu += d1(i) * vals(shiftr(i, 1))
  enddo

  if(real(rho) >= 0.d0) then
    crint_sum_shift = exp(log_fact) * accu
  else
    crint_sum_shift = exp(log_fact - rho) * accu
  endif

  deallocate(vals)

  return
end

! ---

subroutine crint_1_vec_scaled(n_max, rho, vals)

  BEGIN_DOC
  ! Scaled complex Boys functions $e^{\rho} F_n(\rho)$, $n = 0, \dots,
  ! n_{max}$, for $\Re\rho < 0$ (then $|e^{\rho} F_n(\rho)| \le 1$).
  !
  ! $\Re\rho > -700$: crint_1_vec times $e^{\rho}$. Otherwise the asymptotic
  ! series $e^{\rho} F_n(\rho) = -\frac{1}{2\rho} \sum_{k \ge 0}
  ! \frac{(a-1)(a-2)\cdots(a-k)}{\rho^k}$, $a = n + 1/2$, which neglects
  ! $e^{\rho}\Gamma(a)/(2\rho^a)$, below $10^{-300}$ there.
  END_DOC

  implicit none

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n, k
  double precision        :: a
  double precision, parameter :: eps2 = 1.d-34
  complex*16              :: rho_inv, term, accu

  if(real(rho) > -700.d0) then

    call crint_1_vec(n_max, rho, vals)
    vals(0:n_max) = vals(0:n_max) * exp(rho)

  else

    rho_inv = (1.d0, 0.d0) / rho
    do n = 0, n_max
      a = dble(n) + 0.5d0
      term = (1.d0, 0.d0)
      accu = (1.d0, 0.d0)
      do k = 1, 100
        term = term * ((a - dble(k)) * rho_inv)
        accu = accu + term
        if(abs(term)**2 < eps2 * abs(accu)**2) exit
      enddo
      vals(n) = -0.5d0 * rho_inv * accu
    enddo

  endif

  return
end

! ---

subroutine crint_2_vec(n_max, rho, vals)

  BEGIN_DOC
  ! Stores crint_2(n, rho) for n = 0, 1, ..., n_max in vals(0:n_max),
  ! where crint_2(n, rho) = $\int_0^1 t^{2n} \exp(-\rho t^2) dt$ for
  ! a complex argument rho. For |rho| < 0.35 a Taylor polynomial is
  ! used, for 0.35 <= |rho| < 3.5 crint_smallz_vec; for |rho| >= 3.5,
  ! zboysfun if Re(rho) >= 0,
  ! and zboysfunnrp (which returns $\exp(\rho)$ times the integral)
  ! followed by the $\exp(-\rho)$ rescaling if Re(rho) < 0.
  END_DOC
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
      erho = exp(-rho)
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
  ! Complex Boys functions $F_n(\rho)$, $n = 0, \dots, n_{max}$, for
  ! $|\rho| < 3.5$: series for $F_{n_{max}}$ (crint_smallz) and downward
  ! recursion
  ! $F_{n-1} = (2 \rho F_n + e^{-\rho}) / (2n-1)$.
  END_DOC

  implicit none
  integer,    intent(in)      :: n_max
  complex*16, intent(in)      :: rho
  complex*16, intent(out)     :: vals(0:n_max)

  integer                     :: n
  complex*16                  :: rho2, rhoe

  complex*16                  :: crint_smallz

  rho2 = rho + rho
  rhoe = exp(-rho)

  vals(n_max) = crint_smallz(n_max, rho)
  do n = n_max, 1, -1
    vals(n-1) = (rho2 * vals(n) + rhoe) / dble(n + n - 1)
  enddo

  return
end

! ---

subroutine crint_quad_1(n, rho, n_quad, crint_quad)

  BEGIN_DOC
  ! Approximates crint_1(n, rho) with a fixed n_quad-point quadrature
  ! of $\int_0^1 t^{2n} \exp(-\rho t^2) dt$: nodes t = i/n_quad,
  ! weights cycling 32, 12, 32, 14 with period 4 (indexed by i mod 4),
  ! weight 7 at the endpoint t = 1, factor 2/(45 n_quad).
  END_DOC
  implicit none

  integer,    intent(in)  :: n, n_quad
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: crint_quad

  integer                 :: i_quad
  double precision        :: tmp_inv, tmp0, tmp1, tmp2
  double precision        :: coef(0:3) = (/14.d0, 32.d0, 12.d0, 32.d0 /)

  tmp_inv = 1.d0 / dble(n_quad)

  crint_quad = 7.d0 * exp(-rho)

  tmp0 = 0.d0
  select case (n)

    case (0)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * exp(-rho*tmp1)
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (1)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * exp(-rho*tmp1) * tmp1
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (2)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * exp(-rho*tmp1) * tmp1 * tmp1
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (3)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * exp(-rho*tmp1) * tmp1 * tmp1 * tmp1
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (4)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        tmp2 = tmp1 * tmp1
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * exp(-rho*tmp1) * tmp2 * tmp2
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case default
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * exp(-rho*tmp1) * tmp1**n
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv
  end select

end

! ---

subroutine crint_quad_2(n, rho, n_quad, crint_quad)

  BEGIN_DOC
  ! Approximates crint_2(n, rho) with a fixed n_quad-point quadrature
  ! of $\frac{1}{2}\int_0^1 t^{n-1/2} \exp(-\rho t) dt$, the form
  ! obtained from the crint_2 integrand with the change of variable
  ! t = s^2: nodes t = i/n_quad, weights cycling 32, 12, 32, 14 with
  ! period 4, weight 7 at the endpoint t = 1, factor 1/(45 n_quad).
  END_DOC
  implicit none

  integer,    intent(in)  :: n, n_quad
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: crint_quad

  integer                 :: i_quad
  double precision        :: tmp_inv, tmp0, tmp1, tmp2
  double precision        :: coef(0:3) = (/14.d0, 32.d0, 12.d0, 32.d0 /)
  complex*16              :: rhoc, rhoe

  tmp_inv = 1.d0 / dble(n_quad)

  crint_quad = 7.d0 * exp(-rho)

  tmp0 = 0.d0
  rhoc = exp(-rho*tmp_inv)
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

  BEGIN_DOC
  ! Approximates crint_1(n, rho) with a fixed n_quad-point quadrature:
  ! the rule of crint_quad_1 (nodes t = i/n_quad for
  ! $\int_0^1 t^{2n} \exp(-\rho t^2) dt$) is used for n = 0 and the
  ! rule of crint_quad_2 for n >= 1.
  END_DOC
  implicit none

  integer,    intent(in)  :: n, n_quad
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: crint_quad

  integer                 :: i_quad
  double precision        :: tmp_inv, tmp0, tmp1, tmp2
  double precision        :: coef(0:3) = (/14.d0, 32.d0, 12.d0, 32.d0 /)
  complex*16              :: rhoc, rhoe

  tmp_inv = 1.d0 / dble(n_quad)

  crint_quad = 7.d0 * exp(-rho)

  tmp0 = 0.d0
  rhoc = exp(-rho*tmp_inv)
  rhoe = (1.d0, 0.d0)
  select case (n)

    case (0)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * exp(-rho*tmp1)
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

  BEGIN_DOC
  ! Stores crint_quad_12(n, rho) for n = 0, 1, ..., n_max in
  ! vals(0:n_max), with n_quad fixed to 10000000.
  END_DOC
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

