subroutine domain_cheb_weights(n,a,b,w)
  implicit none
  BEGIN_DOC
  ! Weights of the terms $t_q = v^T K^q s$, q=0,n-1, of the series
  ! $v^T (1-K)^{-1} s = \sum_q t_q$.
  !
  ! The first term, which is the second order, is kept with a weight 1.
  ! The rest, $v^T (1-K)^{-1} K s$, is approximated with
  ! $(1-\pi(x))/(1-x)$, where $\pi$ is the Chebyshev polynomial of the
  ! interval [a,b], normalized to 1 at x=1, of degree n-1.
  END_DOC
  integer, intent(in)           :: n
  double precision, intent(in)  :: a, b
  double precision, intent(out) :: w(0:n-1)
  double precision :: t0(0:n), t1(0:n), t2(0:n), al, be, s
  integer :: m, k, d

  w(0) = 1.d0
  m = n-1
  if (m < 1) return
  ! Monomial coefficients of T_m( al x + be )
  al = 2.d0/(b-a)
  be = -(a+b)/(b-a)
  t0 = 0.d0
  t1 = 0.d0
  t0(0) = 1.d0
  t1(0) = be
  t1(1) = al
  do d=2,m
    t2 = 0.d0
    do k=0,d-1
      t2(k)   = t2(k)   + 2.d0*be*t1(k)
      t2(k+1) = t2(k+1) + 2.d0*al*t1(k)
    enddo
    t2 = t2 - t0
    t0 = t1
    t1 = t2
  enddo
  if (m == 0) t1 = t0
  s = sum(t1(0:m))
  t1 = t1/s
  do k=1,m
    w(k) = sum(t1(k:m))
  enddo
end

subroutine domain_cheb_interval(np,nj,delta,a,b,ok,positive)
  implicit none
  BEGIN_DOC
  ! Estimation of the extreme eigenvalues of K from the moments
  ! $t_q = \Delta_{q+2}$, with the Ritz values of the Krylov space
  ! of $s$. The interval is enlarged.
  ! ok is |false| if the estimation is not possible. positive is
  ! |false| if the moments are not those of a positive measure: the
  ! spectrum of K can then be complex, or the moments are too noisy.
  END_DOC
  integer, intent(in)           :: np, nj
  double precision, intent(in)  :: delta(2:np,0:nj)
  double precision, intent(out) :: a, b
  logical, intent(out)          :: ok, positive
  double precision, allocatable :: h0(:,:), h1(:,:), e(:), work(:)
  integer :: k, i, j, info

  ok = .False.
  positive = .True.
  ! Moments 0 to 2k-1 are needed, and 0 to np-2 are available
  do k=(np-1)/2,2,-1
    allocate(h0(k,k),h1(k,k),e(k),work(10*k))
    do j=1,k
      do i=1,k
        h0(i,j) = -delta(i+j,0)
        h1(i,j) = -delta(i+j+1,0)
      enddo
    enddo
    call dsygv(1,'N','U',k,h1,k,h0,k,e,work,10*k,info)
    if (info == 0) then
      if ((e(1) > -5.d0).and.(e(k) < 1.d0)) then
        a = e(1) - 0.1d0
        b = min(e(k) + 0.1d0, 0.95d0)
        ok = a < b
      endif
    endif
    if (info /= 0) positive = .False.
    deallocate(h0,h1,e,work)
    if (ok) return
  enddo
end

double precision function domain_solve(e0,eg,nj,d)
  implicit none
  BEGIN_DOC
  ! Solution of $E = E_0 + \sum_m d_m (E-E_g)^m$, by Newton's method
  ! starting from $E_g$.
  END_DOC
  integer, intent(in)          :: nj
  double precision, intent(in) :: e0, eg, d(0:nj)
  double precision :: x, f, df
  integer :: i, m
  x = 0.d0
  do i=1,50
    f = e0 - eg - x
    df = -1.d0
    do m=0,nj
      f = f + d(m)*x**m
      if (m > 0) df = df + dble(m)*d(m)*x**(m-1)
    enddo
    f = f/df
    x = x - f
    if (dabs(f) < 1.d-14) exit
  enddo
  domain_solve = eg + x
end

subroutine domain_bw_energies(istate,np,nj,eg,delta,a,b,ebw,ecb)
  implicit none
  BEGIN_DOC
  ! Energies obtained with the orders 2 to p, for p=2,np: ebw(p) is
  ! the plain Brillouin-Wigner series and ecb(p) is the one accelerated
  ! by Chebyshev polynomials.
  END_DOC
  integer, intent(in)           :: istate, np, nj
  double precision, intent(in)  :: eg, a, b, delta(2:np,0:nj)
  double precision, intent(out) :: ebw(2:np), ecb(2:np)
  double precision :: w(0:np), d(0:nj)
  double precision, external :: domain_solve
  integer :: p, q
  do p=2,np
    d = 0.d0
    do q=2,p
      d(:) = d(:) + delta(q,:)
    enddo
    ebw(p) = domain_solve(domain_e0(istate),eg,nj,d)
    call domain_cheb_weights(p-1,a,b,w)
    d = 0.d0
    do q=2,p
      d(:) = d(:) + w(q-2)*delta(q,:)
    enddo
    ecb(p) = domain_solve(domain_e0(istate),eg,nj,d)
  enddo
end

subroutine domain_state_cheb(istate,eg,a,b,auto,positive,            &
      ebw,err_bw,ecb,err_cb)
  implicit none
  BEGIN_DOC
  ! Chebyshev-accelerated series for one state. The terms are computed
  ! at the energy eg = $E_0 + E_2$.
  END_DOC
  integer, intent(in)           :: istate
  double precision, intent(out) :: eg, a, b
  logical, intent(out)          :: auto, positive
  double precision, intent(out) :: ebw(2:n_pert), err_bw(2:n_pert)
  double precision, intent(out) :: ecb(2:n_pert), err_cb(2:n_pert)

  double precision, allocatable :: delta2(:), dr(:,:,:), dmean(:,:)
  double precision, allocatable :: dl(:,:), x(:,:), y(:,:), em(:)
  integer :: nrep, r, p, njet
  logical :: exact
  double precision :: f

  njet = domain_njet
  exact = exact_coupling.or.(n_walkers == 0)
  nrep = n_replicas
  if (n_walkers == 0) nrep = 1
  allocate(delta2(0:njet),dr(2:n_pert,0:njet,nrep))
  allocate(dmean(2:n_pert,0:njet),dl(2:n_pert,0:njet))
  allocate(x(2:n_pert,nrep),y(2:n_pert,nrep),em(2:n_pert))

  ! Reference energy: E_0 + E_2
  call domain_samples(istate,domain_e0(istate),nrep,exact,.True.,    &
      delta2,dr)
  eg = domain_e0(istate) + delta2(0)
  call domain_samples(istate,eg,nrep,exact,.False.,delta2,dr)

  dmean = sum(dr,3)/dble(nrep)
  a = domain_cheb_min
  b = domain_cheb_max
  auto = .False.
  positive = .True.
  if (domain_cheb_auto) then
    call domain_cheb_interval(n_pert,njet,dmean,a,b,auto,positive)
    if (.not.auto) then
      a = domain_cheb_min
      b = domain_cheb_max
    endif
  endif
  call domain_bw_energies(istate,n_pert,njet,eg,dmean,a,b,ebw,ecb)
  err_bw = 0.d0
  err_cb = 0.d0
  if (nrep > 1) then
    do r=1,nrep
      dl = (dmean*dble(nrep) - dr(:,:,r))/dble(nrep-1)
      call domain_bw_energies(istate,n_pert,njet,eg,dl,a,b,          &
          x(2,r),y(2,r))
    enddo
    f = dble(nrep-1)/dble(nrep)
    em = sum(x,2)/dble(nrep)
    do p=2,n_pert
      err_bw(p) = dsqrt(f*sum((x(p,:)-em(p))**2))
    enddo
    ebw = dble(nrep)*ebw - dble(nrep-1)*em
    em = sum(y,2)/dble(nrep)
    do p=2,n_pert
      err_cb(p) = dsqrt(f*sum((y(p,:)-em(p))**2))
    enddo
    ecb = dble(nrep)*ecb - dble(nrep-1)*em
  endif
  deallocate(delta2,dr,dmean,dl,x,y,em)
end
