use bitmasks

subroutine domain_spawn_d0(z,njet,n_walk,seed,tout)
  use bitmasks
  use det_table_module
  implicit none
  BEGIN_DOC
  ! Adds to the table tout the vector $Q H z$, where z is a vector of
  ! the reference space. If n_walk = 0 the product is exact, otherwise
  ! it is sampled without bias with n_walk walkers on average.
  END_DOC
  integer, intent(in)            :: njet
  double precision, intent(in)   :: z(0:njet,N_det)
  double precision, intent(in)   :: n_walk
  integer(8), intent(in)         :: seed
  type(det_table), intent(inout) :: tout

  integer(bit_kind), allocatable :: list(:,:,:)
  double precision, allocatable  :: h(:), cum(:)
  double precision :: norm, a, u, x(0:njet), wt, fac
  integer, allocatable :: t_det(:), t_nd(:), t_id(:)
  double precision, allocatable :: t_a(:)
  integer          :: ntask, itask, m, nw
  integer, parameter :: chunk = 10000
  integer          :: i, k, n, nq, nd, it
  integer(8)       :: state
  integer          :: no(2), nv(2), npair, nlist
  integer          :: occ(N_int*bit_kind_size,2)
  integer          :: vir(N_int*bit_kind_size,2)
  integer, allocatable :: pij(:), ptype(:)
  double precision, allocatable :: pcum(:)
  integer(bit_kind) :: d2(N_int,2)
  double precision :: hij, pgen
  logical          :: hb
  integer, external :: domain_index_d0

  PROVIDE all_mo_integrals N_int generators_bitmask domain_jet_scale
  PROVIDE psi_det psi_det_sorted_bit psi_coef_sorted_bit_order
  PROVIDE n_connected_max domain_hij_threshold
  hb = domain_heat_bath.and.(n_walk > 0.d0)
  nlist = n_connected_max
  if (hb) then
    PROVIDE hb_cum hb_ab hb_sum hb_start hb_inv hb_orb hb_n_pair_max
    nlist = 1
  endif

  norm = 0.d0
  do i=1,N_det
    norm = norm + sum(dabs(z(:,i))*domain_jet_scale(0:njet))
  enddo
  if (norm == 0.d0) return

  ! List of tasks: the walkers of a determinant are split in chunks, for
  ! load balancing. The random numbers of a chunk depend only on the
  ! determinant and on the chunk.
  ntask = 0
  do m=1,2
    if (m == 2) then
      allocate(t_det(ntask),t_nd(ntask),t_id(ntask),t_a(ntask))
      ntask = 0
    endif
    do i=1,N_det
      wt = sum(dabs(z(:,i))*domain_jet_scale(0:njet))
      if (wt == 0.d0) cycle
      nd = 1
      a = 1.d0
      if (n_walk > 0.d0) then
        state = ieor(det_hash(psi_det(1,1,i),2*N_int),seed)
        u = det_rnd(state)
        u = det_rnd(state)
        a = n_walk*wt/norm
        nd = int(a+u)
        if (nd == 0) cycle
      endif
      k = 0
      do while (nd > 0)
        nw = nd
        if (n_walk > 0.d0) nw = min(nd,chunk)
        if ((n_walk > 0.d0).and.(.not.domain_heat_bath)) nw = nd
        k = k+1
        ntask = ntask+1
        if (m == 2) then
          t_det(ntask) = i
          t_nd(ntask) = nw
          t_id(ntask) = k
          t_a(ntask) = a
        endif
        nd = nd-nw
      enddo
    enddo
  enddo

  !$OMP PARALLEL DEFAULT(SHARED)                                     &
      !$OMP PRIVATE(list,h,cum,i,k,n,nq,nd,it,a,u,x,fac,state,itask)  &
      !$OMP PRIVATE(no,nv,npair,occ,vir,pij,ptype,pcum,d2,hij,pgen)
  allocate(list(N_int,2,nlist),h(nlist),cum(0:nlist))
  allocate(pij(hb_n_pair_max),ptype(hb_n_pair_max))
  allocate(pcum(0:hb_n_pair_max))
  !$OMP DO SCHEDULE(DYNAMIC,1)
  do itask=1,ntask
    i = t_det(itask)
    nd = t_nd(itask)
    a = t_a(itask)
    state = ieor(det_hash(psi_det(1,1,i),2*N_int),seed)
    state = ieor(state,int(t_id(itask),8)*(-7046029254386353131_8))
    state = det_mix(state)
    if (hb) then
      call hb_prepare(psi_det(1,1,i),no,occ,nv,vir,npair,pij,ptype,  &
          pcum)
      do it=1,nd
        call hb_draw(psi_det(1,1,i),no,occ,nv,vir,npair,pij,ptype,   &
            pcum,state,d2,hij,pgen)
        if (hij == 0.d0) cycle
        if (domain_index_d0(d2) > 0) cycle
        x(:) = z(:,i)*(hij/(a*pgen))
        call table_add(tout,d2,x)
      enddo
      cycle
    endif
    call domain_connected(psi_det(1,1,i),n,list,h)
    ! Keep the external determinants
    nq = 0
    cum(0) = 0.d0
    do k=1,n
      if (domain_index_d0(list(1,1,k)) > 0) cycle
      nq = nq+1
      list(:,:,nq) = list(:,:,k)
      h(nq) = h(k)
      cum(nq) = cum(nq-1)+dabs(h(nq))
    enddo
    if (nq == 0) cycle
    if (n_walk <= 0.d0) then
      do k=1,nq
        x(:) = z(:,i)*h(k)
        call table_add(tout,list(1,1,k),x)
      enddo
    else
      fac = cum(nq)/a
      do it=1,nd
        u = det_rnd(state)*cum(nq)
        call domain_search(cum,nq,u,k)
        x(:) = z(:,i)*dsign(fac,h(k))
        call table_add(tout,list(1,1,k),x)
      enddo
    endif
  enddo
  !$OMP END DO
  deallocate(list,h,cum,pij,ptype,pcum)
  !$OMP END PARALLEL
  deallocate(t_det,t_nd,t_id,t_a)
end

subroutine domain_search(cum,n,u,k)
  implicit none
  BEGIN_DOC
  ! Smallest k such that cum(k) > u
  END_DOC
  integer, intent(in) :: n
  double precision, intent(in) :: cum(0:n), u
  integer, intent(out) :: k
  integer :: lo, hi, mid
  lo = 1
  hi = n
  do while (lo < hi)
    mid = (lo+hi)/2
    if (cum(mid) > u) then
      hi = mid
    else
      lo = mid+1
    endif
  enddo
  k = lo
end

subroutine domain_generation(istate,e_ref,tin,do_score,otf,tv,score, &
      do_spawn,n_walk,seed,tout,d0)
  use bitmasks
  use det_table_module
  implicit none
  BEGIN_DOC
  ! Loops over the external determinants j of the table tin, with
  ! amplitudes $c_j(E)$, and computes
  ! $u_{jl}(E) = - c_j(E) \langle j | (P_j H P_j - E)^{-1} | l \rangle$
  ! for all l in the domain of j.
  !
  ! * if do_score: score = $\sum_{jl} u_{jl} v_l$. The elements of v are
  !   in the table tv, or computed on the fly if otf.
  !
  ! * if do_spawn: $\sum_{jl} u_{jl} H_{lk}$ for k outside the domain of j
  !   is added to tout if k is an external determinant, and to d0
  !   if k is in the reference space. If n_walk = 0 the sum is exact,
  !   otherwise it is sampled without bias, with n_walk walkers on
  !   average.
  !
  ! All the quantities are Taylor expansions around e_ref, up to the
  ! order njet of the table.
  END_DOC
  integer, intent(in)            :: istate
  double precision, intent(in)   :: e_ref
  type(det_table), intent(in)    :: tin, tv
  logical, intent(in)            :: do_score, do_spawn, otf
  double precision, intent(in)   :: n_walk
  integer(8), intent(in)         :: seed
  double precision, intent(out)  :: score(0:tin%njet)
  type(det_table), intent(inout) :: tout
  double precision, intent(inout):: d0(0:tin%njet,N_det)

  integer(bit_kind), allocatable :: list(:,:,:), dom(:,:,:)
  double precision, allocatable  :: h(:), cum(:), g(:,:), u(:,:)
  double precision, allocatable  :: imp(:,:), pl(:)
  integer, allocatable           :: nl(:)
  double precision :: sc(0:tin%njet), x(0:tin%njet), norm, a, r, hd
  double precision :: v, fac, dmin
  integer    :: is, l, k, m, mm, n, nd, nw, it, njet
  integer(8) :: ii, state, nmax_sh
  logical    :: skip
  double precision, external :: diag_H_mat_elem
  integer          :: no(2), nv(2), npair, nlist
  integer          :: occ(N_int*bit_kind_size,2)
  integer          :: vir(N_int*bit_kind_size,2)
  integer, allocatable :: pij(:), ptype(:)
  double precision, allocatable :: pcum(:)
  integer(bit_kind) :: d2(N_int,2)
  double precision :: hij, pgen, va(1)
  logical          :: hb

  PROVIDE all_mo_integrals N_int generators_bitmask domain_jet_scale
  PROVIDE psi_det psi_det_sorted_bit psi_coef_sorted_bit_order
  PROVIDE n_connected_max domain_hij_threshold n_det_domain
  PROVIDE ref_bitmask_energy mo_two_e_integrals_jj_anti
  PROVIDE domain_min_gap_seen domain_phi0
  hb = domain_heat_bath.and.(n_walk > 0.d0)
  nlist = n_connected_max
  if (hb) then
    PROVIDE hb_cum hb_ab hb_sum hb_start hb_inv hb_orb hb_n_pair_max
    if (n_det_domain == 1) nlist = 1
  endif
  if (.not.do_spawn.and.(n_det_domain == 1)) nlist = 1
  PROVIDE mo_two_e_integrals_jj mo_one_e_integrals

  njet = tin%njet
  score = 0.d0

  ! Importance of the determinants
  nmax_sh = 0_8
  do is=1,n_shard
    nmax_sh = max(nmax_sh,tin%sh(is)%n)
  enddo
  if (nmax_sh == 0_8) return
  allocate(imp(nmax_sh,n_shard))
  norm = 0.d0
  dmin = huge(1.d0)
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(is,ii,hd)                &
      !$OMP REDUCTION(+:norm) REDUCTION(min:dmin) SCHEDULE(DYNAMIC,4)
  do is=1,n_shard
    do ii=1,tin%sh(is)%n
      hd = diag_H_mat_elem(tin%sh(is)%key(1,ii),N_int) - e_ref
      dmin = min(dmin,dabs(hd))
      imp(ii,is) = sum(dabs(tin%sh(is)%w(0:njet,ii))                 &
          *domain_jet_scale(0:njet)) / max(dabs(hd),0.1d0)
      norm = norm + imp(ii,is)
    enddo
  enddo
  !$OMP END PARALLEL DO
  domain_min_gap_seen = min(domain_min_gap_seen,dmin)
  if (norm == 0.d0) then
    deallocate(imp)
    return
  endif

  !$OMP PARALLEL DEFAULT(SHARED)                                     &
      !$OMP PRIVATE(list,dom,h,cum,g,u,pl,nl,sc,x,a,r,hd,v,fac,is,l) &
      !$OMP PRIVATE(k,m,mm,n,nd,nw,it,ii,state,skip,no,nv,npair,occ)   &
      !$OMP PRIVATE(vir,pij,ptype,pcum,d2,hij,pgen,va)
  allocate(list(N_int,2,nlist),h(nlist))
  allocate(cum(0:nlist),dom(N_int,2,n_det_domain))
  allocate(pij(hb_n_pair_max),ptype(hb_n_pair_max))
  allocate(pcum(0:hb_n_pair_max))
  allocate(g(n_det_domain,0:njet),u(0:njet,n_det_domain))
  allocate(pl(0:n_det_domain),nl(n_det_domain))
  sc = 0.d0
  !$OMP DO SCHEDULE(DYNAMIC,1)
  do is=1,n_shard
    do ii=1,tin%sh(is)%n
      state = ieor(det_hash(tin%sh(is)%key(1,ii),2*N_int),seed)
      r = det_rnd(state)
      r = det_rnd(state)
      nw = 0
      a = 1.d0
      if (do_spawn) then
        nw = 1
        if (n_walk > 0.d0) then
          a = n_walk*imp(ii,is)/norm
          nw = int(a+r)
        endif
      endif
      if ((nw == 0).and.(.not.do_score)) cycle

      hd = diag_H_mat_elem(tin%sh(is)%key(1,ii),N_int)
      call domain_build(tin%sh(is)%key(1,ii),hd,nd,dom,list,h)
      call domain_resolvent(e_ref,nd,dom,njet,g,n_det_domain,hd)
      do l=1,nd
        do m=0,njet
          u(m,l) = 0.d0
          do mm=0,m
            u(m,l) = u(m,l) - tin%sh(is)%w(m-mm,ii)*g(l,mm)
          enddo
        enddo
      enddo

      if (do_score) then
        do l=1,nd
          if (otf) then
            call i_H_psi(dom(1,1,l),psi_det,domain_phi0(1,istate),   &
                N_int,N_det,N_det,1,va)
            v = va(1)
          else
            v = table_get(tv,dom(1,1,l))
          endif
          sc(:) = sc(:) + u(:,l)*v
        enddo
      endif
      if (nw == 0) cycle

      ! Number of walkers on each determinant of the domain
      if (n_walk > 0.d0) then
        pl(0) = 0.d0
        do l=1,nd
          pl(l) = pl(l-1) + sum(dabs(u(:,l))*domain_jet_scale(0:njet))
        enddo
        if (pl(nd) == 0.d0) cycle
        nl = 0
        do it=1,nw
          r = det_rnd(state)*pl(nd)
          call domain_search(pl,nd,r,l)
          nl(l) = nl(l)+1
        enddo
      endif

      do l=1,nd
        if (n_walk > 0.d0) then
          if (nl(l) == 0) cycle
        endif
        if (hb) then
          call hb_prepare(dom(1,1,l),no,occ,nv,vir,npair,pij,ptype,  &
              pcum)
          fac = pl(nd)/(a*(pl(l)-pl(l-1)))
          do it=1,nl(l)
            call hb_draw(dom(1,1,l),no,occ,nv,vir,npair,pij,ptype,   &
                pcum,state,d2,hij,pgen)
            if (hij == 0.d0) cycle
            call domain_in_list(d2,dom,nd,skip)
            if (skip) cycle
            x(:) = u(:,l)*(hij*fac/pgen)
            call domain_deposit(d2,x,njet,tout,d0)
          enddo
          cycle
        endif
        call domain_connected(dom(1,1,l),n,list,h)
        if (n_walk <= 0.d0) then
          do k=1,n
            call domain_in_list(list(1,1,k),dom,nd,skip)
            if (skip) cycle
            x(:) = u(:,l)*h(k)
            call domain_deposit(list(1,1,k),x,njet,tout,d0)
          enddo
        else
          cum(0) = 0.d0
          do k=1,n
            cum(k) = cum(k-1)+dabs(h(k))
          enddo
          if (n == 0) cycle
          fac = cum(n)*pl(nd)/(a*(pl(l)-pl(l-1)))
          do it=1,nl(l)
            r = det_rnd(state)*cum(n)
            call domain_search(cum,n,r,k)
            call domain_in_list(list(1,1,k),dom,nd,skip)
            if (skip) cycle
            x(:) = u(:,l)*dsign(fac,h(k))
            call domain_deposit(list(1,1,k),x,njet,tout,d0)
          enddo
        endif
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP CRITICAL
  score(:) = score(:) + sc(:)
  !$OMP END CRITICAL
  deallocate(list,h,cum,dom,g,u,pl,nl,pij,ptype,pcum)
  !$OMP END PARALLEL
  deallocate(imp)
end

subroutine domain_in_list(det,dom,nd,found)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! found is |true| if det is one of the nd determinants of dom.
  END_DOC
  integer, intent(in)           :: nd
  integer(bit_kind), intent(in) :: det(N_int,2), dom(N_int,2,nd)
  logical, intent(out)          :: found
  integer :: l
  found = .True.
  do l=1,nd
    if (all(det(:,:) == dom(:,:,l))) return
  enddo
  found = .False.
end

subroutine domain_deposit(det,x,njet,tout,d0)
  use bitmasks
  use det_table_module
  implicit none
  BEGIN_DOC
  ! Adds the contribution x of a walker on the determinant det.
  END_DOC
  integer, intent(in)            :: njet
  integer(bit_kind), intent(in)  :: det(N_int,2)
  double precision, intent(in)   :: x(0:njet)
  type(det_table), intent(inout) :: tout
  double precision, intent(inout):: d0(0:njet,N_det)
  integer, external :: domain_index_d0
  integer :: idx, m
  idx = domain_index_d0(det)
  if (idx > 0) then
    do m=0,njet
      !$OMP ATOMIC
      d0(m,idx) = d0(m,idx) + x(m)
    enddo
  else
    call table_add(tout,det,x)
  endif
end

BEGIN_PROVIDER [ integer, domain_njet ]
 implicit none
 BEGIN_DOC
 ! Order of the Taylor expansions in the energy. For the
 ! Rayleigh-Schrodinger series, as $E_1 = 0$, each power of $E - E_0$
 ! costs two orders of perturbation. With Chebyshev polynomials, the
 ! expansion is used only to solve the self-consistent equation
 ! in the vicinity of $E_0 + E_2$.
 END_DOC
 domain_njet = max(0,(n_pert-2)/2)
 if (domain_chebyshev) domain_njet = max(2,domain_njet)
END_PROVIDER

BEGIN_PROVIDER [ double precision, domain_jet_scale, (0:domain_njet) ]
 implicit none
 BEGIN_DOC
 ! Importance of each order of the Taylor expansions in the sampling.
 END_DOC
 integer :: m
 do m=0,domain_njet
   domain_jet_scale(m) = 0.1d0**m
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, domain_min_gap_seen ]
 implicit none
 BEGIN_DOC
 ! Smallest $|H_{jj}-E_n|$ met among the external determinants
 END_DOC
 domain_min_gap_seen = huge(1.d0)
END_PROVIDER
