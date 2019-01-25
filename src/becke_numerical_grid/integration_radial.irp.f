 double precision function knowles_function(alpha,m,x)
 implicit none
 BEGIN_DOC
! Function proposed by Knowles (JCP, 104, 1996) for distributing the radial points :
! the Log "m" function ( equation (7) in the paper )
 END_DOC
 double precision, intent(in) :: alpha,x
 integer, intent(in) :: m
!print*, x
 knowles_function = -alpha * dlog(1.d0-x**m)
 end

 double precision function derivative_knowles_function(alpha,m,x)
 implicit none
 BEGIN_DOC
! Derivative of the function proposed by Knowles (JCP, 104, 1996) for distributing the radial points
 END_DOC
 double precision, intent(in) :: alpha,x
 integer, intent(in) :: m
 double precision :: f
 f = x**(m-1)
 derivative_knowles_function = alpha * dble(m) * f / (1.d0 - x*f)
 end

 BEGIN_PROVIDER [double precision, alpha_knowles, (100)]
 implicit none
 integer :: i
 BEGIN_DOC
! Recommended values for the alpha parameters according to the paper of Knowles (JCP, 104, 1996)
! as a function of the nuclear charge
 END_DOC

 ! H-He
 alpha_knowles(1) = 5.d0
 alpha_knowles(2) = 5.d0

 ! Li-Be
 alpha_knowles(3) = 7.d0
 alpha_knowles(4) = 7.d0

 ! B-Ne
 do i = 5, 10
  alpha_knowles(i) = 5.d0
 enddo

 ! Na-Mg
 do i = 11, 12
  alpha_knowles(i) = 7.d0
 enddo

 ! Al-Ar
 do i = 13, 18
  alpha_knowles(i) = 5.d0
 enddo

 ! K-Ca
 do i = 19, 20
  alpha_knowles(i) = 7.d0
 enddo

 ! Sc-Zn
 do i = 21, 30
  alpha_knowles(i) = 5.d0
 enddo

 ! Ga-Kr
 do i = 31, 36
  alpha_knowles(i) = 7.d0
 enddo

 END_PROVIDER
