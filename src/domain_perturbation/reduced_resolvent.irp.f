use bitmasks

 BEGIN_PROVIDER [ double precision, domain_phi0, (N_det,N_states) ]
&BEGIN_PROVIDER [ double precision, domain_e0, (N_states) ]
&BEGIN_PROVIDER [ double precision, domain_residual, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Normalized reference states $\Phi_n$, their electronic energies
 ! $E_n$, and the norms of the residuals $P H P \Phi_n - E_n \Phi_n$.
 END_DOC
 integer :: k
 double precision, allocatable :: hu(:,:)
 allocate(hu(N_det,N_states))
 do k=1,N_states
   domain_phi0(1:N_det,k) = psi_coef(1:N_det,k) /                    &
       dsqrt(sum(psi_coef(1:N_det,k)**2))
 enddo
 call domain_h_u(hu,domain_phi0,N_states)
 do k=1,N_states
   domain_e0(k) = sum(hu(:,k)*domain_phi0(:,k))
   domain_residual(k) = dsqrt(sum((hu(:,k)                           &
       -domain_e0(k)*domain_phi0(:,k))**2))
 enddo
 deallocate(hu)
END_PROVIDER

subroutine domain_h_u(v,u,n)
  implicit none
  BEGIN_DOC
  ! $v = P H P u$ in the reference space, determinants in the
  ! order of psi_det.
  END_DOC
  integer, intent(in) :: n
  double precision, intent(in)  :: u(N_det,n)
  double precision, intent(out) :: v(N_det,n)
  double precision, allocatable :: u2(:,:)
  if (N_det == 1) then
    v(1,:) = u(1,:) * ref_h_det1
    return
  endif
  allocate(u2(N_det,n))
  u2 = u
  call H_u_0_nstates_openmp(v,u2,n,N_det)
  deallocate(u2)
end

BEGIN_PROVIDER [ double precision, ref_h_det1 ]
 implicit none
 BEGIN_DOC
 ! Diagonal element of H of the first determinant of psi_det
 END_DOC
 double precision, external :: diag_H_mat_elem
 ref_h_det1 = diag_H_mat_elem(psi_det(1,1,1),N_int)
END_PROVIDER

BEGIN_PROVIDER [ double precision, domain_hdiag, (N_det) ]
 implicit none
 BEGIN_DOC
 ! Diagonal of H in the reference space
 END_DOC
 integer :: i
 double precision, external :: diag_H_mat_elem
 PROVIDE ref_bitmask_energy
 !$OMP PARALLEL DO PRIVATE(i)
 do i=1,N_det
   domain_hdiag(i) = diag_H_mat_elem(psi_det(1,1,i),N_int)
 enddo
 !$OMP END PARALLEL DO
END_PROVIDER

subroutine domain_r0(istate,e,b,z)
  implicit none
  BEGIN_DOC
  ! Applies the reduced resolvent of the reference space:
  ! $z = R_0(e) b$, $R_0(e) = \sum_{m \ne n} |\Phi_m \rangle
  ! (e - E_m)^{-1} \langle \Phi_m |$,
  ! by the minimal residual method with a diagonal preconditioner.
  ! Works for any state: the matrix may be indefinite.
  END_DOC
  integer, intent(in)           :: istate
  double precision, intent(in)  :: e, b(N_det)
  double precision, intent(out) :: z(N_det)
  double precision, allocatable :: r1(:), r2(:), y(:), v(:), w(:)
  double precision, allocatable :: w1(:), w2(:), prec(:), q(:)
  double precision :: beta1, beta, oldb, alfa, dbar, epsln, oldeps
  double precision :: delta, gbar, gam, phi, phibar, cs, sn
  double precision :: bnorm, s
  integer :: iter, i
  integer, parameter :: itmax = 1000

  z = 0.d0
  if (N_det == 1) return
  PROVIDE domain_minres_tol domain_minres_stat
  allocate(r1(N_det),r2(N_det),y(N_det),v(N_det),w(N_det))
  allocate(w1(N_det),w2(N_det),prec(N_det),q(N_det))
  do i=1,N_det
    prec(i) = 1.d0/max(dabs(domain_hdiag(i)-e),0.1d0)
  enddo
  q = b - domain_phi0(:,istate)*sum(domain_phi0(:,istate)*b)
  bnorm = dsqrt(sum(q*q))
  if (bnorm == 0.d0) then
    deallocate(r1,r2,y,v,w,w1,w2,prec,q)
    return
  endif

  ! Preconditioned MINRES (Paige and Saunders)
  r1 = q
  y = prec*r1
  beta1 = dsqrt(sum(r1*y))
  r2 = r1
  oldb = 0.d0
  beta = beta1
  dbar = 0.d0
  epsln = 0.d0
  phibar = beta1
  cs = -1.d0
  sn = 0.d0
  w = 0.d0
  w2 = 0.d0
  do iter=1,itmax
    s = 1.d0/beta
    v = s*y
    call domain_r0_matvec(istate,e,v,y)
    if (iter >= 2) y = y - (beta/oldb)*r1
    alfa = sum(v*y)
    y = y - (alfa/beta)*r2
    r1 = r2
    r2 = y
    y = prec*r2
    oldb = beta
    beta = dsqrt(max(sum(r2*y),0.d0))
    oldeps = epsln
    delta = cs*dbar + sn*alfa
    gbar = sn*dbar - cs*alfa
    epsln = sn*beta
    dbar = -cs*beta
    gam = max(dsqrt(gbar*gbar+beta*beta),1.d-300)
    cs = gbar/gam
    sn = beta/gam
    phi = cs*phibar
    phibar = sn*phibar
    w1 = w2
    w2 = w
    w = (v - oldeps*w1 - delta*w2)/gam
    z = z + phi*w
    if (dabs(phibar) < domain_minres_tol*beta1) exit
    if (beta < 1.d-300) exit
  enddo
  if (iter > itmax) then
    print *, 'Warning: domain_r0 not converged', phibar/beta1
  endif
  domain_minres_stat(1) = domain_minres_stat(1) + 1_8
  domain_minres_stat(2) = domain_minres_stat(2) + int(min(iter,itmax),8)
  domain_minres_stat(3) = max(domain_minres_stat(3),int(min(iter,itmax),8))
  z = z - domain_phi0(:,istate)*sum(domain_phi0(:,istate)*z)
  deallocate(r1,r2,y,v,w,w1,w2,prec,q)
end

subroutine domain_r0_matvec(istate,e,x,y)
  implicit none
  BEGIN_DOC
  ! $y = (e - P H P) x$, where the eigenvalue of $\Phi_n$ is
  ! replaced by -1.
  END_DOC
  integer, intent(in)           :: istate
  double precision, intent(in)  :: e, x(N_det)
  double precision, intent(out) :: y(N_det)
  call domain_h_u(y,x,1)
  y = e*x - y - domain_phi0(:,istate)                               &
      * ((e-domain_e0(istate)+1.d0)*sum(domain_phi0(:,istate)*x))
end

subroutine domain_r0_jets(istate,e,d,z,njet)
  implicit none
  BEGIN_DOC
  ! Taylor expansion of $z(E) = R_0(E) d(E)$ around the energy e:
  ! $z^{(m)} = R_0(e) (d^{(m)} - z^{(m-1)})$.
  END_DOC
  integer, intent(in) :: istate, njet
  double precision, intent(in)  :: e, d(0:njet,N_det)
  double precision, intent(out) :: z(0:njet,N_det)
  double precision, allocatable :: b(:), x(:)
  integer :: m
  allocate(b(N_det),x(N_det))
  do m=0,njet
    b(:) = d(m,:)
    if (m > 0) b(:) = b(:) - z(m-1,:)
    call domain_r0(istate,e,b,x)
    z(m,:) = x(:)
  enddo
  deallocate(b,x)
end

BEGIN_PROVIDER [ integer(8), domain_minres_stat, (3) ]
 implicit none
 BEGIN_DOC
 ! Statistics of the minimal residual method in domain_r0: number of
 ! solutions, total and maximum number of iterations. Reset for each
 ! state by the main program.
 END_DOC
 domain_minres_stat(:) = 0_8
END_PROVIDER
