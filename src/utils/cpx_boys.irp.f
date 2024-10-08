
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
  rho_k        = (1.d0, 0.d0)
  crint_smallz = ct * rho_k / gamma(dble(n) + 1.5d0)

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
    stop 1
  endif

end

! ---

complex*16 function crint_2(n, rho)

  implicit none
  include 'constants.include.F'

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
  double precision, parameter :: tol = 1.0d-03

  double precision, parameter :: x0(1:50) = (/ 0.5d0,                &
                                               0.8660254037844386d0, &
                                               1.2331060371652351d0, &
                                               1.6005429364718398d0, &
                                               1.968141713517676d0,  &
                                               2.3358274254733784d0, &
                                               2.703565149324218d0,  &
                                               3.0713364398393708d0, &
                                               3.4391306442915526d0, &
                                               3.8069411832467615d0, &
                                               4.174763774674144d0,  &
                                               4.5425955121971775d0, &
                                               4.9104343539723185d0, &
                                               5.278278823606476d0,  &
                                               5.646127827158555d0,  &
                                               6.0139805368686865d0, &
                                               6.381836314862427d0,  &
                                               6.749694661668329d0,  &
                                               7.117555180620059d0,  &
                                               7.48541755270593d0,   &
                                               7.853281518456119d0,  &
                                               8.221146864672829d0,  &
                                               8.589013414557305d0,  &
                                               8.956881020260978d0,  &
                                               9.324749557193652d0,  &
                                               9.692618919623587d0,  &
                                               10.060489017239897d0, &
                                               10.428359772440366d0, &
                                               10.796231118172152d0, &
                                               11.164102996198254d0, &
                                               11.531975355694907d0, &
                                               11.899848152108484d0, &
                                               12.26772134621756d0,  &
                                               12.63559490335844d0,  &
                                               13.003468792781872d0, &
                                               13.371342987115684d0, &
                                               13.739217461913594d0, &
                                               14.107092195274445d0, &
                                               14.474967167519445d0, &
                                               14.842842360917235d0, &
                                               15.210717759448817d0, &
                                               15.578593348605757d0, &
                                               15.94646911521622d0,  &
                                               16.31434504729452d0,  &
                                               16.68222113391055d0,  &
                                               17.05009736507609d0,  &
                                               17.4179737316455d0,   &
                                               17.78585022522868d0,  &
                                               18.153726838114668d0, &
                                               18.521603563204277d0 /)

  double precision, parameter :: t(0:11) = (/ 0.20000000000000000d+01, &
                                              0.66666666666666663d+00, &
                                              0.40000000000000002d+00, &
                                              0.28571428571428570d+00, &
                                              0.22222222222222221d+00, &
                                              0.18181818181818182d+00, &
                                              0.15384615384615385d+00, &
                                              0.13333333333333333d+00, &
                                              0.11764705882352941d+00, &
                                              0.10526315789473684d+00, &
                                              0.95238095238095233d-01, &
                                              0.86956521739130432d-01 /)

  complex*16, parameter :: zz(1:10) = (/ ( 0.64304020652330500d+01, 0.18243694739308491d+02), &
                                         ( 0.64304020652330500d+01,-0.18243694739308491d+02), &
                                         (-0.12572081889410178d+01, 0.14121366415342502d+02), &
                                         (-0.12572081889410178d+01,-0.14121366415342502d+02), &
                                         (-0.54103079551670268d+01, 0.10457909575828442d+02), &
                                         (-0.54103079551670268d+01,-0.10457909575828442d+02), &
                                         (-0.78720025594983341d+01, 0.69309284623985663d+01), &
                                         (-0.78720025594983341d+01,-0.69309284623985663d+01), &
                                         (-0.92069621609035313d+01, 0.34559308619699376d+01), &
                                         (-0.92069621609035313d+01,-0.34559308619699376d+01) /)

  complex*16, parameter :: fact(1:10) = (/ ( 0.13249210991966042d-02, 0.91787356295447745d-03), &
                                           ( 0.13249210991966042d-02,-0.91787356295447745d-03), &
                                           ( 0.55545905103006735d-01,-0.35151540664451613d+01), &
                                           ( 0.55545905103006735d-01, 0.35151540664451613d+01), &
                                           (-0.11456407675096416d+03, 0.19213789620924834d+03), &
                                           (-0.11456407675096416d+03,-0.19213789620924834d+03), &
                                           ( 0.20915556220686653d+04,-0.15825742912360638d+04), &
                                           ( 0.20915556220686653d+04, 0.15825742912360638d+04), &
                                           (-0.94779394228935325d+04, 0.30814443710192086d+04), &
                                           (-0.94779394228935325d+04,-0.30814443710192086d+04) /)

  complex*16, parameter :: ww(1:10) = (/ (-0.83418049867878959d-08,-0.70958810331788253d-08), &
                                         (-0.83418050437598581d-08, 0.70958810084577824d-08), &
                                         ( 0.82436739552884774d-07,-0.27704117936134414d-06), &
                                         ( 0.82436739547688584d-07, 0.27704117938414886d-06), &
                                         ( 0.19838416382728666d-05, 0.78321058613942770d-06), &
                                         ( 0.19838416382681279d-05,-0.78321058613180811d-06), &
                                         (-0.47372729839268780d-05, 0.58076919074212929d-05), &
                                         (-0.47372729839287016d-05,-0.58076919074154416d-05), &
                                         (-0.68186014282131608d-05,-0.13515261354290787d-04), &
                                         (-0.68186014282138385d-05, 0.13515261354295612d-04) /)

  double precision, parameter :: rzz(1:1) = (/ -0.96321934290343840d+01 /)
  double precision, parameter :: rfact(1:1) = (/ 0.15247844519077540d+05 /)
  double precision, parameter :: rww(1:1) = (/ 0.18995875677635889d-04 /)

  integer, intent(in) :: n_max
  complex*16, intent(in) :: x
  complex*16, intent(out) :: vals(0:n_max)

  integer :: n, k
  double precision :: x0_nmax
  complex*16 :: y, yy, rtmp
  complex*16 :: p, q, tmp

  y = zexp(-x)
!  x0_nmax = x0(n_max)

!  if(abs(x) .ge. x0_nmax) then

!print*,'check'
    call zboysfun00(x, vals(0))
!print*, vals(0)
    yy = 0.5d0 * y
    do n = 1, n_max
      vals(n) = ((dble(n) - 0.5d0) * vals(n-1) - yy) / x 
!print*, n, x
!print*, vals(n)
    enddo

!  else
!
!    rtmp = (0.d0, 0.d0)
!    do k = 1, 10
!      if(abs(x + zz(k)) .ge. tol) then
!        rtmp = rtmp + ww(k) * (1.0d0 - fact(k)*y) / (x + zz(k))
!      else
!        q = x + zz(k)
!        p = 1.0d0 - 0.5d0 * q + q*q/6.0d0 - q*q*q/24.0d0 + q*q*q*q/120.0d0
!        rtmp = rtmp + ww(k) * p
!      endif
!    enddo
!
!    tmp = (0.d0, 0.d0)
!    do k = 1, 1
!      if(abs(x + rzz(k)) .ge. tol) then 
!        tmp = tmp + rww(k) * (1.0d0 - rfact(k)*y) / (x + rzz(k))
!      else
!        q = x + rzz(k)
!        p = 1.0d0 - 0.5d0 * q + q*q/6.0d0 - q*q*q/24.0d0 + q*q*q*q/120.0d0
!        tmp = tmp + rww(k) * p
!      endif
!    enddo
!
!    vals(n_max) = rtmp + tmp
!    print*, vals(n_max)
!    yy = 0.5d0 * y
!    do n = n_max-1, 0, -1
!      vals(n) = (x * vals(n+1) + yy) * t(n)
!    enddo
!
!  endif

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
  double precision, parameter :: tol = 1.0d-03

  double precision, parameter :: x0(1:50) = (/ 0.5d0,                &
                                               0.8660254037844386d0, &
                                               1.2331060371652351d0, &
                                               1.6005429364718398d0, &
                                               1.968141713517676d0,  &
                                               2.3358274254733784d0, &
                                               2.703565149324218d0,  &
                                               3.0713364398393708d0, &
                                               3.4391306442915526d0, &
                                               3.8069411832467615d0, &
                                               4.174763774674144d0,  &
                                               4.5425955121971775d0, &
                                               4.9104343539723185d0, &
                                               5.278278823606476d0,  &
                                               5.646127827158555d0,  &
                                               6.0139805368686865d0, &
                                               6.381836314862427d0,  &
                                               6.749694661668329d0,  &
                                               7.117555180620059d0,  &
                                               7.48541755270593d0,   &
                                               7.853281518456119d0,  &
                                               8.221146864672829d0,  &
                                               8.589013414557305d0,  &
                                               8.956881020260978d0,  &
                                               9.324749557193652d0,  &
                                               9.692618919623587d0,  &
                                               10.060489017239897d0, &
                                               10.428359772440366d0, &
                                               10.796231118172152d0, &
                                               11.164102996198254d0, &
                                               11.531975355694907d0, &
                                               11.899848152108484d0, &
                                               12.26772134621756d0,  &
                                               12.63559490335844d0,  &
                                               13.003468792781872d0, &
                                               13.371342987115684d0, &
                                               13.739217461913594d0, &
                                               14.107092195274445d0, &
                                               14.474967167519445d0, &
                                               14.842842360917235d0, &
                                               15.210717759448817d0, &
                                               15.578593348605757d0, &
                                               15.94646911521622d0,  &
                                               16.31434504729452d0,  &
                                               16.68222113391055d0,  &
                                               17.05009736507609d0,  &
                                               17.4179737316455d0,   &
                                               17.78585022522868d0,  &
                                               18.153726838114668d0, &
                                               18.521603563204277d0 /)

  double precision, parameter :: t(0:11) = (/ 0.20000000000000000D+01, &
                                              0.66666666666666663D+00, &
                                              0.40000000000000002D+00, &
                                              0.28571428571428570D+00, &
                                              0.22222222222222221D+00, &
                                              0.18181818181818182D+00, &
                                              0.15384615384615385D+00, &
                                              0.13333333333333333D+00, &
                                              0.11764705882352941D+00, &
                                              0.10526315789473684D+00, &
                                              0.95238095238095233D-01, &
                                              0.86956521739130432D-01 /)

  complex*16, parameter :: zz(1:10) = (/ ( 0.64304020652330500D+01, 0.18243694739308491D+02), &
                                         ( 0.64304020652330500D+01,-0.18243694739308491D+02), &
                                         (-0.12572081889410178D+01, 0.14121366415342502D+02), &
                                         (-0.12572081889410178D+01,-0.14121366415342502D+02), &
                                         (-0.54103079551670268D+01, 0.10457909575828442D+02), &
                                         (-0.54103079551670268D+01,-0.10457909575828442D+02), &
                                         (-0.78720025594983341D+01, 0.69309284623985663D+01), &
                                         (-0.78720025594983341D+01,-0.69309284623985663D+01), &
                                         (-0.92069621609035313D+01, 0.34559308619699376D+01), &
                                         (-0.92069621609035313D+01,-0.34559308619699376D+01) /)

  complex*16, parameter :: fact(1:10) = (/ ( 0.13249210991966042D-02, 0.91787356295447745D-03), &
                                           ( 0.13249210991966042D-02,-0.91787356295447745D-03), &
                                           ( 0.55545905103006735D-01,-0.35151540664451613D+01), &
                                           ( 0.55545905103006735D-01, 0.35151540664451613D+01), &
                                           (-0.11456407675096416D+03, 0.19213789620924834D+03), &
                                           (-0.11456407675096416D+03,-0.19213789620924834D+03), &
                                           ( 0.20915556220686653D+04,-0.15825742912360638D+04), &
                                           ( 0.20915556220686653D+04, 0.15825742912360638D+04), &
                                           (-0.94779394228935325D+04, 0.30814443710192086D+04), &
                                           (-0.94779394228935325D+04,-0.30814443710192086D+04) /)

  complex*16, parameter :: ww(1:10) = (/ (-0.83418049867878959D-08,-0.70958810331788253D-08), &
                                         (-0.83418050437598581D-08, 0.70958810084577824D-08), &
                                         ( 0.82436739552884774D-07,-0.27704117936134414D-06), &
                                         ( 0.82436739547688584D-07, 0.27704117938414886D-06), &
                                         ( 0.19838416382728666D-05, 0.78321058613942770D-06), &
                                         ( 0.19838416382681279D-05,-0.78321058613180811D-06), &
                                         (-0.47372729839268780D-05, 0.58076919074212929D-05), &
                                         (-0.47372729839287016D-05,-0.58076919074154416D-05), &
                                         (-0.68186014282131608D-05,-0.13515261354290787D-04), &
                                         (-0.68186014282138385D-05, 0.13515261354295612D-04) /)

  double precision, parameter :: rzz(1:1) = (/ -0.96321934290343840D+01 /)
  double precision, parameter :: rfact(1:1) = (/ 0.15247844519077540D+05 /)
  double precision, parameter :: rww(1:1) = (/ 0.18995875677635889D-04 /)

  integer, intent(in) :: n_max
  complex*16, intent(in) :: x
  complex*16, intent(out) :: vals(0:n_max)

  integer :: n, k
  double precision :: x0_nmax
  complex*16 :: y, rtmp
  complex*16 :: p, q, tmp

  y = exp(x)
!  x0_nmax = x0(n_max)

!  if(abs(x) .ge. x0_nmax) then
    call zboysfun00nrp(x, vals(0))
    do n = 1, n_max
      vals(n) = ((n - 0.5d0) * vals(n-1) - 0.5d0) / x 
    enddo
    return
!  endif
!
!  rtmp = (0.d0, 0.d0)
!  do k = 1, 10
!    if(abs(x + zz(k)) .ge. tol) then 
!      rtmp = rtmp + ww(k) * (y - fact(k)) / (x + zz(k))
!    else
!      q = x+zz(k)
!      p = 1.0d0 - 0.5d0*q + q*q/6.0d0 - q*q*q/24.0d0 + q*q*q*q/120.0d0
!      rtmp = rtmp + ww(k)*p
!    endif
!  enddo
!  tmp = (0.d0, 0.d0)
!  do k = 1, 1
!    if (abs(x + rzz(k)) .ge. tol) then 
!      tmp = tmp + rww(k)*(y-rfact(k))/(x + rzz(k))
!    else
!      q = x+rzz(k)
!      p = 1.0d0 - 0.5d0*q + q*q/6.0d0 - q*q*q/24.0d0 + q*q*q*q/120.0d0
!      tmp = tmp + rww(k) * p
!    endif
!  enddo
!  vals(n_max) = rtmp+tmp
!  do n = n_max-1, 0, -1
!    vals(n) = (x * vals(n+1) + 0.5d0) * t(n)
!  enddo
!  return

end

! ---

complex*16 function crint_sum_2(n_pt_out, rho, d1)

  implicit none

  integer,    intent(in) :: n_pt_out
  complex*16, intent(in) :: rho, d1(0:n_pt_out)

  integer                :: n, i

  complex*16, external   :: crint_2

  crint_sum_2 = (0.d0, 0.d0)
  do i = 0, n_pt_out, 2
    n = shiftr(i, 1)
    crint_sum_2 = crint_sum_2 + d1(i) * crint_2(n, rho)
  enddo

end

! ---

