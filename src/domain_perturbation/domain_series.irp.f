subroutine domain_rs_energies(delta,np,nj,e)
  implicit none
  BEGIN_DOC
  ! Rayleigh-Schrodinger coefficients $E_p$, p=2,np, from the
  ! Taylor coefficients delta(p,m) of the Brillouin-Wigner-like terms
  ! $\Delta_p(E) = \sum_m \Delta_{p,m} (E-E_0)^m$, by solving order
  ! by order
  ! $E(\lambda) = E_0 + \sum_p \lambda^p \Delta_p(E(\lambda))$.
  END_DOC
  integer, intent(in)           :: np, nj
  double precision, intent(in)  :: delta(2:np,0:nj)
  double precision, intent(out) :: e(2:np)
  ! pw(k,m) : coefficient of lambda**k in (E(lambda)-E0)**m
  double precision :: pw(0:np,0:nj)
  integer :: n, p, m, k

  e(:) = 0.d0
  do n=2,np
    ! Powers of the series known up to order n-1
    pw = 0.d0
    pw(0,0) = 1.d0
    do m=1,nj
      do k=2,n-1
        do p=0,n-1-k
          pw(p+k,m) = pw(p+k,m) + e(k)*pw(p,m-1)
        enddo
      enddo
    enddo
    e(n) = 0.d0
    do p=2,n
      do m=0,nj
        e(n) = e(n) + delta(p,m)*pw(n-p,m)
      enddo
    enddo
  enddo
end

subroutine domain_rpt_energies(delta,np,nj,e)
  implicit none
  BEGIN_DOC
  ! Renormalized corrections rPTp, p=2,np: the equation
  ! $E - E_0 = \sum_{q=2}^p \Delta_q(E)$ linearized around $E_0$
  ! (one Newton step),
  ! $rPT_p = \sum_q \Delta_{q,0} / (1 - \sum_q \Delta_{q,1})$.
  ! For p=2 with singletons this is the rPT2 of the CIPSI programs,
  ! as $\Delta_{2,1} = -\langle \psi_1 | \psi_1 \rangle$.
  ! Requires nj >= 1.
  END_DOC
  integer, intent(in)           :: np, nj
  double precision, intent(in)  :: delta(2:np,0:nj)
  double precision, intent(out) :: e(2:np)
  double precision :: s0, s1
  integer :: p
  s0 = 0.d0
  s1 = 0.d0
  do p=2,np
    s0 = s0 + delta(p,0)
    s1 = s1 + delta(p,1)
    e(p) = s0/(1.d0 - s1)
  enddo
end

subroutine domain_replica(istate,irep,e_ref,tv,exact,delta2,delta)
  use det_table_module
  implicit none
  BEGIN_DOC
  ! One evaluation of the Taylor coefficients delta(p,m) for a state.
  ! If n_walkers = 0 the result is exact and does not depend on irep.
  !
  ! $c_0 = v = Q H \Phi_0$, $c_{q+1} = K^T c_q$, and
  ! $\Delta_{q+2} = c_q^T s$ with $s = - G v$.
  ! K contains the reduced resolvent of the reference space on the rows
  ! of the reference space.
  !
  ! If exact, tv contains the exact v, and delta2 is the
  ! exact $\Delta_2$. Otherwise $c_0$ is a sample of v, and the elements
  ! of v needed by the scores are computed exactly on the fly.
  END_DOC
  integer, intent(in)           :: istate, irep
  double precision, intent(in)  :: e_ref
  type(det_table), intent(in)   :: tv
  logical, intent(in)           :: exact
  double precision, intent(in)  :: delta2(0:domain_njet)
  double precision, intent(out) :: delta(2:n_pert,0:domain_njet)
  type(det_table) :: c, cn, t1
  double precision, allocatable :: d0(:,:), z(:,:), zphi(:,:), score(:)
  double precision :: mw
  integer(8) :: seed
  integer :: q, njet
  logical :: spawn

  njet = domain_njet
  mw = dble(n_walkers)
  allocate(d0(0:njet,N_det),z(0:njet,N_det),zphi(0:njet,N_det))
  allocate(score(0:njet))
  delta = 0.d0
  z = 0.d0

  if (.not.exact) then
    zphi = 0.d0
    zphi(0,:) = domain_phi0(:,istate)
    call table_create(t1,2*N_int,njet)
    call domain_make_seed(istate,irep,1,seed)
    call domain_spawn_d0(zphi,njet,mw,seed,t1)
  endif

  do q=0,n_pert-2
    spawn = q < n_pert-2
    call domain_make_seed(istate,irep,10+q,seed)
    d0 = 0.d0
    call table_create(cn,2*N_int,njet)
    if (exact) then
      if (q == 0) then
        score = delta2
        if (spawn) then
          call domain_generation(istate,e_ref,tv,.False.,.False., &
              tv,score,.True.,mw,seed,cn,d0)
          score = delta2
        endif
      else
        call domain_generation(istate,e_ref,c,.True.,.False.,tv,   &
            score,spawn,mw,seed,cn,d0)
      endif
    else
      if (q == 0) then
        call domain_generation(istate,e_ref,t1,.True.,.True.,tv,   &
            score,spawn,mw,seed,cn,d0)
      else
        call domain_generation(istate,e_ref,c,.True.,.True.,tv,    &
            score,spawn,mw,seed,cn,d0)
      endif
    endif
    delta(q+2,:) = score(:)

    if (spawn) then
      ! The reference-space part of c_q is spawned after the
      ! application of the reduced resolvent
      if (q > 0) then
        call domain_make_seed(istate,irep,1000+q,seed)
        call domain_spawn_d0(z,njet,mw,seed,cn)
      endif
      ! Reference-space part of c_{q+1}
      if (q < n_pert-3) call domain_r0_jets(istate,e_ref,d0,z,njet)
    endif
    if (q > 0) call table_free(c)
    call table_create(c,2*N_int,njet)
    call move_alloc_table(cn,c)
    call table_free(cn)
  enddo
  call table_free(c)
  if (.not.exact) then
    call table_free(t1)
  endif
  deallocate(d0,z,zphi,score)
end

subroutine domain_make_seed(istate,irep,k,seed)
  implicit none
  BEGIN_DOC
  ! Seed of the sampling, different for each state, replica and step
  END_DOC
  integer, intent(in)     :: istate, irep, k
  integer(8), intent(out) :: seed
  integer(8) :: key(4), z
  integer    :: i
  key(1) = int(domain_seed,8)
  key(2) = int(istate,8)
  key(3) = int(irep,8)
  key(4) = int(k,8)
  z = 88172645463325252_8
  do i=1,4
    z = ieor(z,key(i)+1_8)
    z = z * (-7046029254386353131_8)
    z = ieor(z,shiftr(z,31))
    z = z * (-4658895280553007687_8)
    z = ieor(z,shiftr(z,29))
  enddo
  seed = z
end

subroutine domain_samples(istate,e_ref,nrep,exact,first,delta2,dr)
  use det_table_module
  implicit none
  BEGIN_DOC
  ! Taylor coefficients dr(p,m,r) of $\Delta_p(E)$ around e_ref, for
  ! all the replicas r. delta2 is $\Delta_2$, exact or sampled
  ! by the first replica.
  ! If first is |true|, only delta2 is computed.
  END_DOC
  integer, intent(in)           :: istate, nrep
  double precision, intent(in)  :: e_ref
  logical, intent(in)           :: exact, first
  double precision, intent(out) :: delta2(0:domain_njet)
  double precision, intent(out) :: dr(2:n_pert,0:domain_njet,nrep)
  type(det_table) :: tv, tdum, t1
  double precision, allocatable :: zphi(:,:), d0(:,:)
  integer(8) :: seed
  integer :: r, njet
  njet = domain_njet
  delta2 = 0.d0
  dr = 0.d0
  call table_create(tv,2*N_int,njet)
  allocate(zphi(0:njet,N_det),d0(0:njet,N_det))
  zphi = 0.d0
  zphi(0,:) = domain_phi0(:,istate)
  d0 = 0.d0
  call table_create(tdum,2*N_int,njet)
  if (exact) then
    call domain_spawn_d0(zphi,njet,0.d0,0_8,tv)
    if (.not.first) then
      print *, 'Number of external determinants :', table_size(tv)
    endif
  else if (first) then
    call table_create(t1,2*N_int,njet)
    call domain_make_seed(istate,0,1,seed)
    call domain_spawn_d0(zphi,njet,dble(n_walkers),seed,t1)
    call domain_generation(istate,e_ref,t1,.True.,.True.,tv,delta2,  &
        .False.,0.d0,0_8,tdum,d0)
    call table_free(t1)
  endif
  if (exact) then
    call domain_generation(istate,e_ref,tv,.True.,.False.,tv,        &
        delta2,.False.,0.d0,0_8,tdum,d0)
  endif
  call table_free(tdum)
  if (.not.first) then
    do r=1,nrep
      call domain_replica(istate,r,e_ref,tv,exact,delta2,dr(2,0,r))
    enddo
  endif
  call table_free(tv)
  deallocate(zphi,d0)
end

subroutine domain_state(istate,e_order,err_order,err_tot,e_rpt,err_rpt)
  implicit none
  BEGIN_DOC
  ! Rayleigh-Schrodinger series for one state. e_order(p) is $E_p$, and
  ! err_order(p) and err_tot(p) are the statistical errors on $E_p$ and
  ! on the sum of the orders 2 to p, estimated by a jackknife over
  ! the replicas. e_rpt(p) is the renormalized correction rPTp
  ! (domain_rpt_energies) and err_rpt(p) its error; they are zero if
  ! the Taylor expansions are of order 0 (n_pert < 4).
  END_DOC
  integer, intent(in)           :: istate
  double precision, intent(out) :: e_order(2:n_pert)
  double precision, intent(out) :: err_order(2:n_pert), err_tot(2:n_pert)
  double precision, intent(out) :: e_rpt(2:n_pert), err_rpt(2:n_pert)

  double precision, allocatable :: delta2(:)
  double precision, allocatable :: dr(:,:,:), dmean(:,:), dl(:,:)
  double precision, allocatable :: el(:,:), etot(:,:), em(:), erl(:,:)
  integer :: nrep, r, p, njet
  logical :: exact
  double precision :: f

  njet = domain_njet
  exact = exact_coupling.or.(n_walkers == 0)
  nrep = n_replicas
  if (n_walkers == 0) nrep = 1

  allocate(delta2(0:njet))
  allocate(dr(2:n_pert,0:njet,nrep),dmean(2:n_pert,0:njet))
  allocate(dl(2:n_pert,0:njet),el(2:n_pert,nrep),etot(2:n_pert,nrep))
  allocate(em(2:n_pert),erl(2:n_pert,nrep))
  call domain_samples(istate,domain_e0(istate),nrep,exact,.False.,   &
      delta2,dr)

  dmean = sum(dr,3)/dble(nrep)
  call domain_rs_energies(dmean,n_pert,njet,e_order)
  e_rpt   = 0.d0
  err_rpt = 0.d0
  if (njet >= 1) call domain_rpt_energies(dmean,n_pert,njet,e_rpt)
  err_order = 0.d0
  err_tot   = 0.d0
  if (nrep > 1) then
    do r=1,nrep
      dl = (dmean*dble(nrep) - dr(:,:,r))/dble(nrep-1)
      call domain_rs_energies(dl,n_pert,njet,el(2,r))
      etot(2,r) = el(2,r)
      do p=3,n_pert
        etot(p,r) = etot(p-1,r) + el(p,r)
      enddo
      if (njet >= 1) call domain_rpt_energies(dl,n_pert,njet,erl(2,r))
    enddo
    f = dble(nrep-1)/dble(nrep)
    em = sum(el,2)/dble(nrep)
    do p=2,n_pert
      err_order(p) = dsqrt(f*sum((el(p,:)-em(p))**2))
    enddo
    ! The products of averages of E_p are biased by 1/n_replicas:
    ! jackknife correction
    e_order = dble(nrep)*e_order - dble(nrep-1)*em
    em = sum(etot,2)/dble(nrep)
    do p=2,n_pert
      err_tot(p) = dsqrt(f*sum((etot(p,:)-em(p))**2))
    enddo
    if (njet >= 1) then
      em = sum(erl,2)/dble(nrep)
      do p=2,n_pert
        err_rpt(p) = dsqrt(f*sum((erl(p,:)-em(p))**2))
      enddo
      e_rpt = dble(nrep)*e_rpt - dble(nrep-1)*em
    endif
  endif
  deallocate(delta2,dr,dmean,dl,el,etot,em,erl)
end
