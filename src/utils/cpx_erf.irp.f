
! ---

complex*16 function cpx_erf_1(x, y)

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

    cpx_erf_1 = (1.d0, 0.d0) * derf(x)
    return

  else

    erf_tmp1 = (1.d0, 0.d0) * derf(x)
    erf_tmp2 = erf_E(x, yabs) + erf_F(x, yabs)
    erf_tmp3 = zexp(-(0.d0, 2.d0) * x * yabs) * (erf_G(x, yabs) + erf_H(x, yabs))
    erf_tot  = erf_tmp1 + erf_tmp2 - erf_tmp3

  endif

  if(y .gt. 0.d0) then
    cpx_erf_1 = erf_tot
  else
    cpx_erf_1 = conjg(erf_tot)
  endif

end

! ---

complex*16 function erf_E(x, yabs)
 
  implicit none
  include 'constants.include.F'

  double precision, intent(in) :: x, yabs

  if((dabs(x).gt.6.d0) .or. (x==0.d0)) then
    erf_E = (0.d0, 0.d0)
    return
  endif

  if(dabs(x) .lt. 1.d-7) then

    erf_E = -inv_pi * (0.d0, 1.d0) * yabs

  else

    erf_E = 0.5d0 * inv_pi * dexp(-x*x) &
          * ((1.d0, 0.d0) - zexp(-(2.d0, 0.d0) * x * yabs)) / x

  endif

end

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

end

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
    tmp2  = dexp(idble * yabs - tmp1 - dlog(tmp1) - dlog_2pi) * (x - (0.d0, 1.d0)*tmp0)

    erf_G = erf_G + tmp2

  enddo

end

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


  if((dabs(x) .lt. 10d0) .and. (yabs .lt. 6.1d0)) then

    x2 = x * x
    ct = 0.5d0 * inv_pi

    erf_H = 0.d0
    do i = 1, Nmax

      idble = dble(i)
      tmp0  = 0.5d0 * idble
      tmp1  = tmp0 * tmp0 + x2
      tmp2  = dexp(-tmp1-idble*yabs) * (x + (0.d0, 1.d0)*tmp0) / tmp1
      erf_H = erf_H + tmp2

      tmp_mod = dsqrt(real(tmp2)*real(tmp2) + aimag(tmp2)*aimag(tmp2))
      if(tmp_mod .lt. 1d-15) exit
    enddo
    erf_H = ct * erf_H

  else

    erf_H = (0.d0, 0.d0)

  endif 

end

! ---

subroutine zboysfun00(z, val)

  BEGIN_DOC
  !
  ! Computes values of the Boys function for n=0 
  ! for a complex valued argument
  !
  ! Input: z  --- argument, complex*16, Real(z) >= 0
  ! Output: val  --- value of the Boys function n=0
  !
  ! Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
  ! https://doi.org/10.1063/5.0062444
  !
  END_DOC

  implicit none

  double precision, parameter :: asymcoef(1:7) = (/ -0.499999999999999799d0, &
                                                     0.249999999999993161d0, &
                                                    -0.374999999999766599d0, &
                                                     0.937499999992027020d0, &
                                                    -3.28124999972738868d0,  &
                                                     14.7656249906697030d0,  &
                                                    -81.2109371803307752d0 /)

  double precision, parameter :: taylcoef(0:10) = (/ 1.0d0,                    &
                                                    -0.333333333333333333d0,   &
                                                     0.1d0,                    &
                                                    -0.238095238095238095d-01, &
                                                     0.462962962962962963d-02, &
                                                    -0.757575757575757576d-03, &
                                                     0.106837606837606838d-03, & 
                                                    -0.132275132275132275d-04, &
                                                     1.458916900093370682d-06, &
                                                    -1.450385222315046877d-07, &
                                                     1.3122532963802805073d-08 /)

  double precision, parameter :: sqpio2 = 0.886226925452758014d0

  double precision, parameter :: pp(1:22) = (/ 0.001477878263796956477d0, &
                                               0.013317276413725817441d0, &
                                               0.037063591452052541530d0, &
                                               0.072752512422882761543d0, &
                                               0.120236941228785688896d0, &
                                               0.179574293958937717967d0, &
                                               0.253534046984087292596d0, &
                                               0.350388652780721927513d0, &
                                               0.482109575931276669313d0, &
                                               0.663028993158374107103d0, &
                                               0.911814736856590885929d0, &
                                               1.2539502287919293d0,      &
                                               1.7244634233573395d0,      &
                                               2.3715248262781863d0,      &
                                               3.2613796996078355d0,      &
                                               4.485130169059591d0,       &
                                               6.168062135122484d0,       &
                                               8.48247187231787d0,        &
                                               11.665305486296793d0,      &
                                               16.042417132288328d0,      &
                                               22.06192951814709d0,       &
                                               30.340112094708307d0 /)

  double precision, parameter :: ff(1:22) = (/ 0.0866431027201416556d0, &
                                               0.0857720608434394764d0, &
                                               0.0839350436829178814d0, &
                                               0.0809661970413229146d0, &
                                               0.0769089548492978618d0, &
                                               0.0731552078711821626d0, &
                                               0.0726950035163157228d0, &
                                               0.0752842556089304050d0, &
                                               0.0770943953645196145d0, &
                                               0.0754250625677530441d0, &
                                               0.0689686192650315305d0, &
                                               0.05744480422143023d0,   &
                                               0.04208199434694545d0,   &
                                               0.025838539448223282d0,  &
                                               0.012445024157255563d0,  &
                                               0.004292541592599837d0,  &
                                               0.0009354342987735969d0, &
                                               0.10840885466502504d-03, &
                                               5.271867966761674d-06,   &
                                               7.765974039750418d-08,   &
                                               2.2138172422680093d-10,  &
                                               6.594161760037707d-14 /)

  complex*16, intent(in)  :: z
  complex*16, intent(out) :: val

  integer    :: k
  complex*16 :: z1, y

  if(abs(z) .ge. 100.0d0) then

    ! large |z|
    z1 = (1.d0, 0.d0) / zsqrt(z)
    y  = (1.d0, 0.d0) / z
    val = asymcoef(7)
    do k = 6, 1, -1
      val = val * y + asymcoef(k)
    enddo
    val = zexp(-z) * val * y + z1 * sqpio2

  else if(abs(z) .le. 0.35d0) then

    ! small |z|
    val = taylcoef(10) * (1.d0, 0.d0)
    do k = 9, 0, -1
      val = val * z + taylcoef(k)
    enddo

  else

    ! intermediate |z|
    val = sqpio2 / zsqrt(z) - 0.5d0 * zexp(-z) * sum(ff(1:22)/(z+pp(1:22)))

  endif

  return
end

! ---

subroutine zboysfun00nrp(z, val)

  BEGIN_DOC
  !
  ! Computes values of the exp(z) F(0,z)
  ! (where F(0,z) is the Boys function)
  ! for a complex valued argument with Real(z)<=0
  !
  ! Input: z  --- argument, complex*16, !!! Real(z)<=0 !!!
  ! Output: val  --- value of the function !!! exp(z) F(0,z) !!!, where F(0,z) is the Boys function
  !
  ! Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
  ! https://doi.org/10.1063/5.0062444
  !
  END_DOC

  include 'constants.include.F'

  implicit none

  double precision, parameter :: asymcoef(1:7) = (/ -0.499999999999999799d0, &
                                                     0.249999999999993161d0, &
                                                    -0.374999999999766599d0, &
                                                     0.937499999992027020d0, &
                                                    -3.28124999972738868d0,  &
                                                     14.7656249906697030d0,  &
                                                    -81.2109371803307752d0 /)

  double precision, parameter :: taylcoef(0:10) = (/ 1.0d0,                    &
                                                    -0.333333333333333333d0,   &
                                                     0.1d0,                    &
                                                    -0.238095238095238095d-01, &
                                                     0.462962962962962963d-02, &
                                                    -0.757575757575757576d-03, &
                                                     0.106837606837606838d-03, & 
                                                    -0.132275132275132275d-04, &
                                                     1.458916900093370682d-06, &
                                                    -1.450385222315046877d-07, &
                                                     1.3122532963802805073d-08 /)

  double precision, parameter :: tol     = 1.0d-03
  double precision, parameter :: sqpio2  = 0.886226925452758014d0 ! sqrt(pi)/2
  double precision, parameter :: etmax   = 25.7903399171930621d0
  double precision, parameter :: etmax1  = 26.7903399171930621d0
  complex*16, parameter :: ima = (0.d0, 1.d0)

  double precision, parameter :: pp(1:16) = (/ 0.005299532504175031d0, &
                                               0.0277124884633837d0,   &
                                               0.06718439880608407d0,  &
                                               0.12229779582249845d0,  &
                                               0.19106187779867811d0,  &
                                               0.27099161117138637d0,  &
                                               0.35919822461037054d0,  &
                                               0.45249374508118123d0,  &
                                               0.5475062549188188d0,   &
                                               0.6408017753896295d0,   &
                                               0.7290083888286136d0,   &
                                               0.8089381222013219d0,   &
                                               0.8777022041775016d0,   &
                                               0.9328156011939159d0,   &
                                               0.9722875115366163d0,   &
                                               0.994700467495825d0 /)

  double precision, parameter :: ww(1:16) = (/ 0.013576229705876844d0, &
                                               0.03112676196932382d0,  &
                                               0.04757925584124612d0,  &
                                               0.062314485627766904d0, &
                                               0.07479799440828848d0,  &
                                               0.08457825969750153d0,  &
                                               0.09130170752246194d0,  &
                                               0.0947253052275344d0,   &
                                               0.0947253052275344d0,   &
                                               0.09130170752246194d0,  &
                                               0.08457825969750153d0,  &
                                               0.07479799440828848d0,  &
                                               0.062314485627766904d0, &
                                               0.04757925584124612d0,  &
                                               0.03112676196932382d0,  &
                                               0.013576229705876844d0 /)

  double precision, parameter :: qq(1:16) = (/ 0.0007243228510223928d0, &
                                               0.01980651726441906d0,   &
                                               0.11641097769229371d0,   &
                                               0.38573968881461146d0,   &
                                               0.9414671037609641d0,    &
                                               1.8939510935716377d0,    &
                                               3.3275564293459383d0,    &
                                               5.280587297262129d0,     &
                                               7.730992222360452d0,     &
                                               10.590207725831563d0,    &
                                               13.706359477128965d0,    &
                                               16.876705473663804d0,    &
                                               19.867876155236257d0,    &
                                               22.441333930203022d0,    &
                                               24.380717439613566d0,    &
                                               25.51771075067431d0 /)


  double precision, parameter :: qq1(1:16) = (/ 0.0007524078957852004d0,&
                                                0.020574499281252233d0, &
                                                0.12092472113522865d0,  &
                                                0.40069643967765295d0,  &
                                                0.9779717449089211d0,   &
                                                1.9673875468969015d0,   &
                                                3.4565797939091802d0,   &
                                                5.485337886599723d0,    &
                                                8.030755321535683d0,    &
                                                11.000834641174064d0,   &
                                                14.237812708111456d0,   &
                                                17.531086359214406d0,   &
                                                20.6382373144543d0,     &
                                                23.31147887603379d0,    &
                                                25.326060444703632d0,   &
                                                26.507139770710722d0 /)

  double precision, parameter :: uu(1:16) = (/ 0.9992759394074501d0,      &
                                               0.9803883431758104d0,      &
                                               0.8901093330366746d0,      &
                                               0.6799475005849274d0,      &
                                               0.3900551639790145d0,      &
                                               0.15047608763371934d0,     &
                                               0.0358806749968974d0,      &
                                               0.005089440900100864d0,    &
                                               0.00043900830706867264d0,  &
                                               0.000025161192619824898d0, &
                                               1.1153308427285078d-6,     &
                                               4.68317018372038d-8,       &
                                               2.3522908467181876d-9,     &
                                               1.7941242138648815d-10,    &
                                               2.5798173021885247d-11,    &
                                               8.27559122014575d-12 /)


  double precision, parameter :: uu1(1:16) = (/ 0.999247875092057d0,       &
                                                0.979635711599488d0,       &
                                                0.8861006617341018d0,      &
                                                0.6698533710831932d0,      &
                                                0.3760730980014839d0,      &
                                                0.13982165701683388d0,     &
                                                0.031537442321301304d0,    &
                                                0.004147133581658446d0,    &
                                                0.0003253024081883165d0,   &
                                                0.000016687766678889653d0, &
                                                6.555359391864376d-7,      &
                                                2.4341421258295026d-8,     &
                                                1.0887481200652014d-9,     &
                                                7.51542178140961d-11,      &
                                                1.002378402152542d-11,     &
                                                3.0767730761654096d-12 /)

  complex*16, intent(in)  :: z
  complex*16, intent(out) :: val

  integer    :: k
  complex*16 :: z1, zz, y, zsum, tmp, zt, q, p

  zz = zexp(z)

  if(abs(z) .ge. 100.0d0) then
    ! large |z|
    z1 = (1.d0, 0.d0) / zsqrt(z)
    y  = (1.d0, 0.d0) / z
    val = asymcoef(7)
    do k = 6, 1, -1
      val = val * y + asymcoef(k)
    enddo
    val = val * y + z1 * sqpio2 * zz
    return
  endif

  if(abs(z) .le. 0.35d0) then
    ! small |z|
    val = taylcoef(10) * (1.d0, 0.d0)
    do k = 9, 0, -1
      val = val * z + taylcoef(k)
    enddo
    val = val * zz
    return
  endif

  if(abs(etmax+z) .ge. 0.5d0) then
    ! intermediate |z|
    zsum = (0.d0, 0.d0)
    do k = 1, 16
      if(abs(z + qq(k)) .ge. tol) then 
        zsum = zsum + ww(k) * (zz - uu(k)) / (qq(k) + z)
      else
        q = z + qq(k)
        p = q * (0.041666666666666664d0*q * (q * (0.2d0*q - 1.d0) + 4.d0) - 0.5d0) + 1.d0
        zsum = zsum + ww(k) * p * zz
      endif
    enddo
    zt = ima * zsqrt(z / etmax)
    tmp = 0.5d0 * ima * log((1.0d0 - zt) / (1.0d0 + zt))
    val = dsqrt(etmax) * zsum * inv_sq_pi + zz * tmp / zsqrt(pi*z)
  else
    zsum = (0.d0, 0.d0)
    do k = 1, 16
      if(abs(z + qq1(k)) .ge. tol) then 
        zsum = zsum + ww(k) * (zz - uu1(k)) / (qq1(k) + z)
      else
        q = z + qq1(k)
        !p = 1.0d0 - 0.5d0*q + q*q/6.0d0 - q*q*q/24.0d0 + q*q*q*q/120.0d0
        p = q * (0.041666666666666664d0*q * (q * (0.2d0*q - 1.d0) + 4.d0) - 0.5d0) + 1.d0
        zsum = zsum + ww(k) * p * zz
      endif
    enddo
    zt = ima * zsqrt(z / etmax1)
    tmp = 0.5d0 * ima * log((1.0d0 - zt) / (1.0d0 + zt))
    val = dsqrt(etmax1) * zsum * inv_sq_pi + zz * tmp / zsqrt(pi*z)
  endif

  return
end

! ---

