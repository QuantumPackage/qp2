double precision function derf_mu_x(mu,x)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: mu,x
  if(dabs(x).gt.1.d-6)then
   derf_mu_x = derf(mu * x)/x
  else
   derf_mu_x =  inv_sq_pi * 2.d0 * mu * (1.d0 - mu*mu*x*x/3.d0)                                                           
  endif                      
end    


double precision function binom_func(i,j)
  implicit none
  BEGIN_DOC
  !.. math                       ::
  !
  !  \frac{i!}{j!(i-j)!}
  !
  END_DOC
  integer,intent(in)             :: i,j
  double precision               :: logfact
  integer, save                  :: ifirst
  double precision, save         :: memo(0:15,0:15)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: memo
  integer                        :: k,l
  if (ifirst == 0) then
    ifirst = 1
    do k=0,15
      do l=0,15
        memo(k,l) = dexp( logfact(k)-logfact(l)-logfact(k-l) )
      enddo
    enddo
  endif
  if ( (i<=15).and.(j<=15) ) then
    binom_func = memo(i,j)
  else
    binom_func = dexp( logfact(i)-logfact(j)-logfact(i-j) )
  endif

  ! To avoid .999999 numbers
  binom_func = floor(binom_func + 0.5d0)

end


 BEGIN_PROVIDER [ double precision, binom, (0:40,0:40) ]
&BEGIN_PROVIDER [ double precision, binom_transp, (0:40,0:40) ]
  implicit none
  BEGIN_DOC
  ! Binomial coefficients
  END_DOC
  integer                        :: k,l
  double precision               :: logfact
  do k=0,40
    do l=0,40
      binom(k,l) = dexp( logfact(k)-logfact(l)-logfact(k-l) )
      binom_transp(l,k) = binom(k,l)
    enddo
  enddo
END_PROVIDER


 BEGIN_PROVIDER [ integer*8, binom_int, (0:40,0:40) ]
&BEGIN_PROVIDER [ integer*8, binom_int_transp, (0:40,0:40) ]
  implicit none
  BEGIN_DOC
  ! Binomial coefficients, as integers*8
  END_DOC
  integer                        :: k,l
  double precision               :: logfact
  do l=0,40
    do k=0,40
      binom_int(k,l) = int(binom(k,l)+0.1d0,8)
    enddo
  enddo
END_PROVIDER



double precision function fact(n)
  implicit none
  BEGIN_DOC
  ! n!
  END_DOC
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1

  if (n<=memomax) then
    if (n<2) then
      fact = 1.d0
    else
      fact = memo(n)
    endif
    return
  endif

  integer                        :: i
  memo(1) = 1.d0
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)*dble(i)
  enddo
  memomax = min(n,100)
  double precision :: logfact
  fact = dexp(logfact(n))
end function

double precision function logfact(n)
  implicit none
  BEGIN_DOC
  ! n!
  END_DOC
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1

  if (n<=memomax) then
    if (n<2) then
      logfact = 0.d0
    else
      logfact = memo(n)
    endif
    return
  endif

  integer                        :: i
  memo(1) = 0.d0
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)+dlog(dble(i))
  enddo
  memomax = min(n,100)
  logfact = memo(memomax)
  do i=101,n
    logfact += dlog(dble(i))
  enddo
end function

! ---

BEGIN_PROVIDER [ double precision, fact_inv, (128) ]
  implicit none
  BEGIN_DOC
  ! 1/n!
  END_DOC
  integer                        :: i
  double precision               :: fact
  do i=1,size(fact_inv)
    fact_inv(i) = 1.d0/fact(i)
  enddo
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, shiftfact_op5_inv, (128) ]

  BEGIN_DOC
  !
  ! 1 / Gamma(n + 0.5)
  !
  END_DOC

  implicit none
  integer          :: i
  double precision :: tmp

  do i = 1, size(shiftfact_op5_inv)
    !tmp = dgamma(dble(i) + 0.5d0)
    tmp = gamma(dble(i) + 0.5d0)
    shiftfact_op5_inv(i) = 1.d0 / tmp
  enddo

END_PROVIDER

! ---

double precision function dble_fact(n)
  implicit none
  integer :: n
  double precision :: dble_fact_even, dble_fact_odd

  dble_fact = 1.d0

  if(n.lt.0) return

  if(iand(n,1).eq.0)then
    dble_fact = dble_fact_even(n)
  else
    dble_fact= dble_fact_odd(n)
  endif

end function

double precision function dble_fact_even(n) result(fact2)
  implicit none
  BEGIN_DOC
  ! n!!
  END_DOC
  integer                        :: n,k
  double precision, save         :: memo(0:100)
  integer, save                  :: memomax = 0
  double precision               :: prod

  ASSERT (iand(n,1) /= 1)

!  prod=1.d0
!  do k=2,n,2
!   prod=prod*dfloat(k)
!  enddo
!  fact2=prod
!  return
!
  if (n <= memomax) then
    if (n < 2) then
      fact2 = 1.d0
    else
      fact2 = memo(n)
    endif
    return
  endif

  integer                        :: i
  memo(0)=1.d0
  memo(1)=1.d0
  do i=memomax+2,min(n,100),2
    memo(i) = memo(i-2)* dble(i)
  enddo
  memomax = min(n,100)
  fact2 = memo(memomax)

  if (n > 100) then
    double precision :: dble_logfact
    fact2 = dexp(dble_logfact(n))
  endif

end function

double precision function dble_fact_odd(n) result(fact2)
  implicit none
  BEGIN_DOC
  ! n!!
  END_DOC
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1

  ASSERT (iand(n,1) /= 0)
  if (n<=memomax) then
    if (n<3) then
      fact2 = 1.d0
    else
      fact2 = memo(n)
    endif
    return
  endif

  integer                        :: i
  memo(1) = 1.d0
  do i=memomax+2,min(n,99),2
    memo(i) = memo(i-2)* dble(i)
  enddo
  memomax = min(n,99)
  fact2 = memo(memomax)

  if (n > 99) then
    double precision :: dble_logfact
    fact2 = dexp(dble_logfact(n))
  endif

end function

double precision function dble_logfact(n) result(logfact2)
  implicit none
  BEGIN_DOC
  ! n!!
  END_DOC
  integer                        :: n
  integer :: k
  double precision :: prod
  prod=0.d0
  do k=2,n,2
   prod=prod+dlog(dfloat(k))
  enddo
  logfact2=prod
  return

end function

subroutine write_git_log(iunit)
  implicit none
  BEGIN_DOC
  ! Write the last git commit in file iunit.
  END_DOC
  integer, intent(in)            :: iunit
  write(iunit,*) '----------------'
  write(iunit,*) 'Last git commit:'
  BEGIN_SHELL [ /bin/bash ]
  git log -1 2>/dev/null | sed "s/'//g"| sed "s/^/    write(iunit,*) '/g" | sed "s/$/'/g" || echo "Unknown"
  END_SHELL
  write(iunit,*) '----------------'
end

BEGIN_PROVIDER [ double precision, inv_int, (128) ]
  implicit none
  BEGIN_DOC
  ! 1/i
  END_DOC
  integer                        :: i
  do i=1,128
    inv_int(i) = 1.d0/dble(i)
  enddo
END_PROVIDER

subroutine wall_time(t)
  implicit none
  BEGIN_DOC
  ! The equivalent of cpu_time, but for the wall time.
  END_DOC
  double precision, intent(out)  :: t
  integer*8                        :: c
  integer*8, save                  :: rate = 0
  if (rate == 0) then
    CALL SYSTEM_CLOCK(count_rate=rate)
  endif
  CALL SYSTEM_CLOCK(count=c)
  t = dble(c)/dble(rate)
end

BEGIN_PROVIDER [ integer, nproc ]
  implicit none
  BEGIN_DOC
  ! Number of current OpenMP threads
  END_DOC

  integer, external              :: omp_get_num_threads
  nproc = 1
  !$OMP PARALLEL
  !$OMP MASTER
  !$ nproc = omp_get_num_threads()
  !$OMP END MASTER
  !$OMP END PARALLEL
END_PROVIDER


double precision function u_dot_v(u,v,sze)
  implicit none
  BEGIN_DOC
  ! Compute <u|v>
  END_DOC
  integer, intent(in)            :: sze
  double precision, intent(in)   :: u(sze),v(sze)
  double precision, external     :: ddot

  !DIR$ FORCEINLINE
  u_dot_v = ddot(sze,u,1,v,1)

end

double precision function u_dot_u(u,sze)
  implicit none
  BEGIN_DOC
  ! Compute <u|u>
  END_DOC
  integer, intent(in)            :: sze
  double precision, intent(in)   :: u(sze)
  double precision, external     :: ddot

  !DIR$ FORCEINLINE
  u_dot_u = ddot(sze,u,1,u,1)

end

subroutine normalize(u,sze)
  implicit none
  BEGIN_DOC
  ! Normalizes vector u
  END_DOC
  integer, intent(in)            :: sze
  double precision, intent(inout):: u(sze)
  double precision               :: d
  double precision, external     :: dnrm2
  integer                        :: i

  !DIR$ FORCEINLINE
  d = dnrm2(sze,u,1)
  if (d /= 0.d0) then
    d = 1.d0/d
  endif
  if (d /= 1.d0) then
    !DIR$ FORCEINLINE
    call dscal(sze,d,u,1)
  endif
end

double precision function approx_dble(a,n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: a
  double precision :: f
  integer :: i

  if (a == 0.d0) then
    approx_dble = 0.d0
    return
  endif
  f = 1.d0
  do i=1,-int(dlog10(dabs(a)))+n
    f = f*.1d0
  enddo
  do i=1,int(dlog10(dabs(a)))-n
    f = f*10.d0
  enddo
  approx_dble = dnint(a/f)*f

end



subroutine lowercase(txt,n)
  implicit none
  BEGIN_DOC
! Transform to lower case
  END_DOC
  character*(*), intent(inout)   :: txt
  integer, intent(in)            :: n
  character( * ), PARAMETER      :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  character( * ), PARAMETER      :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  integer                        :: i, ic
  do i=1,n
    ic = index( UPPER_CASE, txt(i:i) )
    if (ic /= 0) then
      txt(i:i) = LOWER_CASE(ic:ic)
    endif
  enddo
end

subroutine v2_over_x(v,x,res)

  !BEGIN_DOC
  ! Two by two diagonalization to avoid the divergence in v^2/x when x goes to 0
  !END_DOC

  implicit none

  double precision, intent(in)  :: v, x
  double precision, intent(out) :: res

  double precision :: delta_E, tmp, val

  res = 0d0
  delta_E = x
  if (v == 0.d0) return

  val = 2d0 * v
  tmp = dsqrt(delta_E * delta_E + val * val)
  if (delta_E < 0.d0) then
      tmp = -tmp
  endif
  res = 0.5d0 * (tmp - delta_E)

end
