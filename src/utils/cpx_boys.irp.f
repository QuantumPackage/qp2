
! ---

complex*16 function crint_1(n, rho)

  implicit none
  include 'constants.include.F'

  integer,    intent(in) :: n
  complex*16, intent(in) :: rho

  integer                :: i, mmax
  double precision       :: rho_mod, rho_re, rho_im 
  double precision       :: sq_rho_re, sq_rho_im
  double precision       :: n_tmp
  complex*16             :: sq_rho, rho_inv, rho_exp

  complex*16             :: crint_smallz, cpx_erf_1
  complex*16             :: cpx_erf_2

  rho_re  = real (rho)
  rho_im  = aimag(rho)
  rho_mod = dsqrt(rho_re*rho_re + rho_im*rho_im)

  if(rho_mod < 10.d0) then
    ! small z
    if(rho_mod .lt. 1.d-15) then
      crint_1 = 1.d0 / dble(n + n + 1)
    else
      crint_1 = crint_smallz(n, rho)
    endif

  else
    ! large z

    if(rho_mod .gt. 40.d0) then

      n_tmp = dble(n) + 0.5d0
      crint_1 = 0.5d0 * gamma(n_tmp) / (rho**n_tmp)

    else

      ! get \sqrt(rho)
      sq_rho_re = sq_op5 * dsqrt(rho_re + rho_mod)
      sq_rho_im = 0.5d0 * rho_im / sq_rho_re 
      sq_rho    = sq_rho_re + (0.d0, 1.d0) * sq_rho_im

      rho_exp = 0.5d0 * zexp(-rho)
      rho_inv = (1.d0, 0.d0) / rho

      !print*, sq_rho_re, sq_rho_im
      !print*, cpx_erf_1(sq_rho_re, sq_rho_im)
      !print*, cpx_erf_2(sq_rho_re, sq_rho_im)

      crint_1 = 0.5d0 * sqpi * cpx_erf_1(sq_rho_re, sq_rho_im) / sq_rho
      mmax = n
      if(mmax .gt. 0) then
        do i = 0, mmax-1
          crint_1 = ((dble(i) + 0.5d0) * crint_1 - rho_exp) * rho_inv
        enddo
      endif

    endif

  endif

end

! ---

complex*16 function crint_quad(n, rho)

  implicit none

  integer,    intent(in)  :: n
  complex*16, intent(in)  :: rho

  integer                 :: i_quad, n_quad
  double precision        :: tmp_inv, tmp

  n_quad = 1000000000
  tmp_inv = 1.d0 / dble(n_quad)

  !crint_quad = 0.5d0 * zexp(-rho)
  !do i_quad = 1, n_quad - 1
  !  tmp = tmp_inv * dble(i_quad)
  !  tmp = tmp * tmp
  !  crint_quad += zexp(-rho*tmp) * tmp**n
  !enddo
  !crint_quad = crint_quad * tmp_inv

  !crint_quad = 0.5d0 * zexp(-rho)
  !do i_quad = 1, n_quad - 1
  !  tmp = tmp_inv * dble(i_quad)
  !  crint_quad += zexp(-rho*tmp) * tmp**n / dsqrt(tmp)
  !enddo
  !crint_quad = crint_quad * 0.5d0 * tmp_inv

  ! Composite Boole's Rule
  crint_quad = 7.d0 * zexp(-rho)
  do i_quad = 1, n_quad - 1
    tmp = tmp_inv * dble(i_quad)
    tmp = tmp * tmp
    if(modulo(i_quad, 4) .eq. 0) then
      crint_quad += 14.d0 * zexp(-rho*tmp) * tmp**n
    else if(modulo(i_quad, 2) .eq. 0) then
      crint_quad += 12.d0 * zexp(-rho*tmp) * tmp**n
    else
      crint_quad += 32.d0 * zexp(-rho*tmp) * tmp**n
    endif
  enddo
  crint_quad = crint_quad * 2.d0 * tmp_inv / 45.d0

  ! Composite Simpson's 3/8 rule
  !crint_quad = zexp(-rho)
  !do i_quad = 1, n_quad - 1
  !  tmp = tmp_inv * dble(i_quad)
  !  tmp = tmp * tmp
  !  if(modulo(i_quad, 3) .eq. 0) then
  !    crint_quad += 2.d0 * zexp(-rho*tmp) * tmp**n
  !  else
  !    crint_quad += 3.d0 * zexp(-rho*tmp) * tmp**n
  !  endif
  !enddo
  !crint_quad = crint_quad * 3.d0 * tmp_inv / 8.d0

end

! ---

complex*16 function crint_sum_1(n_pt_out, rho, d1)

  implicit none
  include 'constants.include.F'

  integer,    intent(in)  :: n_pt_out
  complex*16, intent(in)  :: rho, d1(0:n_pt_out)

  integer                 :: n, i, mmax
  double precision        :: rho_mod, rho_re, rho_im 
  double precision        :: sq_rho_re, sq_rho_im
  complex*16              :: sq_rho, F0
  complex*16              :: rho_tmp, rho_inv, rho_exp
  complex*16, allocatable :: Fm(:)

  complex*16              :: crint_smallz, cpx_erf_1


  rho_re  = real (rho)
  rho_im  = aimag(rho)
  rho_mod = dsqrt(rho_re*rho_re + rho_im*rho_im)

!  ! debug
!  double precision        :: d1_real(0:n_pt_out)
!  double precision        :: rint_sum
!  do i = 0, n_pt_out
!    d1_real(i) = real(d1(i))
!  enddo
!  crint_sum_1 = rint_sum(n_pt_out, rho_re, d1_real)
!  return

  if(rho_mod < 10.d0) then
    ! small z

    if(rho_mod .lt. 1.d-15) then

      crint_sum_1 = d1(0)
      do i = 2, n_pt_out, 2
        n = shiftr(i, 1)
        crint_sum_1 = crint_sum_1 + d1(i) / dble(n+n+1)
      enddo

    else

      crint_sum_1 = d1(0) * crint_smallz(0, rho)
      do i = 2, n_pt_out, 2
        n = shiftr(i, 1)
        crint_sum_1 = crint_sum_1 + d1(i) * crint_smallz(n, rho)
      enddo

    endif

  else
    ! large z

    if(rho_mod .gt. 40.d0) then

      rho_inv = (1.d0, 0.d0) / rho
      rho_tmp = 0.5d0 * sqpi * zsqrt(rho_inv)

      crint_sum_1 = rho_tmp * d1(0)
      do i = 2, n_pt_out, 2
        n = shiftr(i, 1)
        rho_tmp   = rho_tmp * (dble(n) + 0.5d0) * rho_inv
        crint_sum_1 = crint_sum_1 + rho_tmp * d1(i)
      enddo

    else

      ! get \sqrt(rho)
      sq_rho_re = sq_op5 * dsqrt(rho_re + rho_mod)
      sq_rho_im = 0.5d0 * rho_im / sq_rho_re 
      sq_rho    = sq_rho_re + (0.d0, 1.d0) * sq_rho_im

      F0        = 0.5d0 * sqpi * cpx_erf_1(sq_rho_re, sq_rho_im) / sq_rho
      crint_sum_1 = F0 * d1(0)

      rho_exp = 0.5d0 * zexp(-rho)
      rho_inv = (1.d0, 0.d0) / rho

      mmax = shiftr(n_pt_out, 1)
      if(mmax .gt. 0) then

        allocate(Fm(mmax))
        Fm(1:mmax) = (0.d0, 0.d0)

        do n = 0, mmax-1
          F0      = ((dble(n) + 0.5d0) * F0 - rho_exp) * rho_inv
          Fm(n+1) = F0
        enddo
    
        do i = 2, n_pt_out, 2
          n = shiftr(i, 1)
          crint_sum_1 = crint_sum_1 + Fm(n) * d1(i)
        enddo

        deallocate(Fm)
      endif

    endif ! rho_mod
  endif ! rho_mod

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
  double precision, parameter :: eps = 1.d-13

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
    write(*,*) ' n, rho    = ', n, rho
    write(*,*) ' delta_mod = ', delta_mod
    !stop 1
  endif

end

! ---

complex*16 function crint_2(n, rho)

  implicit none

  integer,    intent(in) :: n
  complex*16, intent(in) :: rho

  double precision :: tmp
  complex*16 :: rho2
  complex*16 :: vals(0:n)
  complex*16, external :: crint_smallz

  if(abs(rho) < 10.d0) then

    if(abs(rho) .lt. 1d-6) then
      tmp = 2.d0 * dble(n)
      rho2 = rho * rho
      crint_2 = 1.d0 / (tmp + 1.d0)         &
              - rho / (tmp + 3.d0)          &
              + 0.5d0 * rho2 / (tmp + 5.d0) &
              - 0.16666666666666666d0 * rho * rho2 / (tmp + 7.d0)
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
  END_DOC

  implicit none

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: x
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n
  complex*16              :: yy, x_inv

  call zboysfun00(x, vals(0))

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

complex*16 function crint_sum_2(n_pt_out, rho, d1)

  implicit none

  integer,    intent(in)  :: n_pt_out
  complex*16, intent(in)  :: rho, d1(0:n_pt_out)
                          
  integer                 :: n, i
  integer                 :: n_max

  complex*16, allocatable :: vals(:)

  !complex*16, external    :: crint_2
  !crint_sum_2 = (0.d0, 0.d0)
  !do i = 0, n_pt_out, 2
  !  n = shiftr(i, 1)
  !  crint_sum_2 = crint_sum_2 + d1(i) * crint_2(n, rho)
  !enddo

  n_max = shiftr(n_pt_out, 1)

  allocate(vals(0:n_max))
  call crint_2_vec(n_max, rho, vals)

  crint_sum_2 = d1(0) * vals(0)
  do i = 2, n_pt_out, 2
    n = shiftr(i, 1)
    crint_sum_2 += d1(i) * vals(n)
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
  complex*16              :: rho2, rho3, erho


  abs_rho = abs(rho)

  if(abs_rho < 10.d0) then

    if(abs_rho .lt. 1d-6) then

      ! use finite expansion for very small rho

      ! rho^2 / 2
      rho2 = 0.5d0 * rho * rho
      ! rho^3 / 6
      rho3 = 0.3333333333333333d0 * rho * rho2

      vals(0) = 1.d0 - 0.3333333333333333d0 * rho + 0.2d0 * rho2 - 0.14285714285714285d0 * rho3
      do n = 1, n_max
        tmp = 2.d0 * dble(n)
        vals(n) = 1.d0 / (tmp + 1.d0) - rho  / (tmp + 3.d0) &
                + rho2 / (tmp + 5.d0) - rho3 / (tmp + 7.d0)
      enddo

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
  double precision, parameter :: eps = 1.d-13

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
  
    !if(abs(delta_k) > eps) then
    !  write(*,*) ' pb in crint_smallz_vec !'
    !  write(*,*) ' n, rho    = ', n, rho
    !  write(*,*) ' |delta_k| = ', abs(delta_k)
    !endif
  enddo

  deallocate(rho_k)

  return
end

! ---

