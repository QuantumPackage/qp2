
! ---

complex*16 function cpx_erf(x, y)

  BEGIN_DOC
  !
  ! compute erf(z) for z = x + i y
  !
  ! REF: Abramowitz and Stegun
  !
  END_DOC

  implicit none
  
  double precision, intent(in) :: x, y

  double precision             :: yabs
  complex*16                   :: erf_tmp1, erf_tmp2, erf_tmp3, erf_tot

  double precision             :: erf_F 
  complex*16                   :: erf_E, erf_G, erf_H

  yabs = dabs(y)

  if(yabs .lt. 1.d-15) then

    cpx_erf = (1.d0, 0.d0) * derf(x)
    return

  else

    erf_tmp1 = (1.d0, 0.d0) * derf(x)
    erf_tmp2 = erf_E(x, yabs) + erf_F(x, yabs)
    erf_tmp3 = zexp(-(0.d0, 2.d0) * x * yabs) * ( erf_G(x, yabs) + erf_H(x, yabs) )
    erf_tot  = erf_tmp1 + erf_tmp2 - erf_tmp3

  endif

  if(y .gt. 0.d0) then
    cpx_erf = erf_tot
  else
    cpx_erf = CONJG(erf_tot)
  endif

end function cpx_erf

! ---

complex*16 function erf_E(x, yabs)
 
  implicit none
  include 'constants.include.F'

  double precision, intent(in) :: x, yabs

  if( (dabs(x).gt.6.d0) .or. (x==0.d0) ) then
    erf_E = (0.d0, 0.d0)
    return
  endif

  if(dabs(x) .lt. 1.d-7) then

    erf_E = -inv_pi * (0.d0, 1.d0) * yabs

  else

    erf_E = 0.5d0 * inv_pi * dexp(-x*x) &
          * ((1.d0, 0.d0) - zexp(-(2.d0, 0.d0) * x * yabs)) / x

  endif

end function erf_E

! ---

double precision function erf_F(x, yabs)
 
  implicit none
  include 'constants.include.F'

  double precision, intent(in) :: x, yabs

  integer, parameter           :: Nmax = 13

  integer                      :: i
  double precision             :: tmp1, tmp2, x2, ct 


  if(dabs(x) .gt. 5.8d0) then

    erf_F = 0.d0

  else

    x2 = x * x
    ct = x * inv_pi

    erf_F = 0.d0
    do i = 1, Nmax

      tmp1  = 0.25d0 * dble(i) * dble(i) + x2
      tmp2  = dexp(-tmp1) / tmp1
      erf_F = erf_F + tmp2

      if(dabs(tmp2) .lt. 1d-15) exit
    enddo
    erf_F = ct * erf_F

  endif 

end function erf_F

! ---

complex*16 function erf_G(x, yabs)

  implicit none
  include 'constants.include.F'

  double precision, intent(in) :: x, yabs

  integer, parameter           :: Nmax = 13

  integer                      :: i, tmpi, imin, imax
  double precision             :: tmp0, tmp1, x2, idble
  complex*16                   :: tmp2

  if(x .eq. 0.d0) then
    erf_G = (0.d0, 0.d0)
    return
  endif

  tmpi = int(2.d0 * yabs)
  imin = max(1, tmpi-Nmax)
  imax = tmpi + Nmax

  x2 = x * x

  erf_G = 0.d0
  do i = imin, imax

    idble = dble(i)
    tmp0  = 0.5d0 * idble
    tmp1  = tmp0 * tmp0 + x2
    tmp2  = dexp( idble * yabs - tmp1 - dlog(tmp1) - dlog_2pi) * (x - (0.d0, 1.d0)*tmp0)

    erf_G = erf_G + tmp2

  enddo

end function erf_G

! ---

complex*16 function erf_H(x, yabs)
 
  implicit none
  include 'constants.include.F'

  double precision, intent(in) :: x, yabs

  integer, parameter           :: Nmax = 13

  integer                      :: i
  double precision             :: tmp0, tmp1, tmp_mod, x2, ct, idble
  complex*16                   :: tmp2

  if(x .eq. 0.d0) then
    erf_H = (0.d0, 0.d0)
    return
  endif


  if( (dabs(x) .lt. 10d0) .and. (yabs .lt. 6.1d0) ) then

    x2 = x * x
    ct = 0.5d0 * inv_pi

    erf_H = 0.d0
    do i = 1, Nmax

      idble = dble(i)
      tmp0  = 0.5d0 * idble
      tmp1  = tmp0 * tmp0 + x2
      tmp2  = dexp(-tmp1-idble*yabs) * (x + (0.d0, 1.d0)*tmp0) / tmp1
      erf_H = erf_H + tmp2

      tmp_mod = dsqrt(REAL(tmp2)*REAL(tmp2) + AIMAG(tmp2)*AIMAG(tmp2))
      if(tmp_mod .lt. 1d-15) exit
    enddo
    erf_H = ct * erf_H

  else

    erf_H = (0.d0, 0.d0)

  endif 

end function erf_H

! ---


