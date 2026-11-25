use bitmasks

 BEGIN_PROVIDER [ integer, hb_n_orb ]
&BEGIN_PROVIDER [ integer, hb_orb, (mo_num) ]
 implicit none
 BEGIN_DOC
 ! Orbitals which can be holes or particles.
 END_DOC
 integer :: i, k
 hb_n_orb = 0
 hb_orb = 0
 do i=1,mo_num
   k = shiftr(i-1,bit_kind_shift)+1
   if (btest(ior(ior(generators_bitmask(k,1,s_hole),                 &
                     generators_bitmask(k,1,s_part)),                &
                 ior(generators_bitmask(k,2,s_hole),                 &
                     generators_bitmask(k,2,s_part))),               &
             iand(i-1,bit_kind_size-1))) then
     hb_n_orb = hb_n_orb+1
     hb_orb(hb_n_orb) = i
   endif
 enddo
END_PROVIDER

BEGIN_PROVIDER [ integer, hb_inv, (mo_num) ]
 implicit none
 BEGIN_DOC
 ! Index of an orbital in hb_orb
 END_DOC
 integer :: i
 hb_inv = 0
 do i=1,hb_n_orb
   hb_inv(hb_orb(i)) = i
 enddo
END_PROVIDER

 BEGIN_PROVIDER [ integer(8), hb_start, (0:hb_n_orb*hb_n_orb,2) ]
&BEGIN_PROVIDER [ integer(8), hb_size_max ]
 implicit none
 BEGIN_DOC
 ! Number of non-zero doubles ij -> ab, for all the pairs of orbitals.
 ! hb_start(ij,s) is the position of the last element of the pair.
 ! s=1 : opposite spins, i for the first electron. s=2 : same spin, i<j
 ! (pair index (i-1)*n+j), a<b.
 END_DOC
 integer :: i, j, a, b, ij, n
 integer(8) :: cnt(hb_n_orb*hb_n_orb,2)
 double precision :: x
 double precision, external :: get_two_e_integral
 n = hb_n_orb
 cnt = 0_8
 PROVIDE all_mo_integrals hb_orb
 !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,i,j,a,b,x)             &
     !$OMP SCHEDULE(DYNAMIC,8)
 do ij=1,n*n
   i = (ij-1)/n+1
   j = ij-(i-1)*n
   do a=1,n
     if (a == i) cycle
     do b=1,n
       if (b == j) cycle
       x = get_two_e_integral(hb_orb(i),hb_orb(j),hb_orb(a),         &
           hb_orb(b),mo_integrals_map)
       if (dabs(x) > domain_hij_threshold) cnt(ij,1) = cnt(ij,1)+1_8
     enddo
   enddo
   if (i >= j) cycle
   do a=1,n
     do b=a+1,n
       if ((a == i).or.(a == j).or.(b == i).or.(b == j)) cycle
       x = get_two_e_integral(hb_orb(i),hb_orb(j),hb_orb(a),         &
           hb_orb(b),mo_integrals_map)                               &
         - get_two_e_integral(hb_orb(i),hb_orb(j),hb_orb(b),         &
           hb_orb(a),mo_integrals_map)
       if (dabs(x) > domain_hij_threshold) cnt(ij,2) = cnt(ij,2)+1_8
     enddo
   enddo
 enddo
 !$OMP END PARALLEL DO
 do b=1,2
   hb_start(0,b) = 0_8
   do ij=1,n*n
     hb_start(ij,b) = hb_start(ij-1,b)+cnt(ij,b)
   enddo
 enddo
 hb_size_max = max(1_8,maxval(hb_start(n*n,:)))
END_PROVIDER

 BEGIN_PROVIDER [ double precision, hb_cum, (hb_size_max,2) ]
&BEGIN_PROVIDER [ integer, hb_ab, (hb_size_max,2) ]
&BEGIN_PROVIDER [ double precision, hb_sum, (hb_n_orb*hb_n_orb,2) ]
 implicit none
 BEGIN_DOC
 ! Cumulated absolute values of the integrals of the doubles ij -> ab,
 ! orbitals a and b ((a-1)*n+b), and their sum for each pair.
 END_DOC
 integer :: i, j, a, b, ij, n
 integer(8) :: k
 double precision :: x, s
 double precision, external :: get_two_e_integral
 n = hb_n_orb
 hb_sum = 0.d0
 !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,i,j,a,b,x,s,k)         &
     !$OMP SCHEDULE(DYNAMIC,8)
 do ij=1,n*n
   i = (ij-1)/n+1
   j = ij-(i-1)*n
   s = 0.d0
   k = hb_start(ij-1,1)
   do a=1,n
     if (a == i) cycle
     do b=1,n
       if (b == j) cycle
       x = get_two_e_integral(hb_orb(i),hb_orb(j),hb_orb(a),         &
           hb_orb(b),mo_integrals_map)
       if (dabs(x) > domain_hij_threshold) then
         k = k+1_8
         s = s+dabs(x)
         hb_cum(k,1) = s
         hb_ab(k,1) = (a-1)*n+b
       endif
     enddo
   enddo
   hb_sum(ij,1) = s
   if (i >= j) cycle
   s = 0.d0
   k = hb_start(ij-1,2)
   do a=1,n
     do b=a+1,n
       if ((a == i).or.(a == j).or.(b == i).or.(b == j)) cycle
       x = get_two_e_integral(hb_orb(i),hb_orb(j),hb_orb(a),         &
           hb_orb(b),mo_integrals_map)                               &
         - get_two_e_integral(hb_orb(i),hb_orb(j),hb_orb(b),         &
           hb_orb(a),mo_integrals_map)
       if (dabs(x) > domain_hij_threshold) then
         k = k+1_8
         s = s+dabs(x)
         hb_cum(k,2) = s
         hb_ab(k,2) = (a-1)*n+b
       endif
     enddo
   enddo
   hb_sum(ij,2) = s
 enddo
 !$OMP END PARALLEL DO
END_PROVIDER

BEGIN_PROVIDER [ integer, hb_n_pair_max ]
 implicit none
 BEGIN_DOC
 ! Maximum number of pairs of electrons
 END_DOC
 hb_n_pair_max = max(1,(elec_num*(elec_num-1))/2)
END_PROVIDER

subroutine hb_prepare(det,no,occ,nv,vir,npair,pij,ptype,pcum)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Occupied and empty orbitals of det, and cumulated weights of its
  ! pairs of electrons. ptype : 1 alpha-alpha, 2 beta-beta, 3 alpha-beta.
  END_DOC
  integer(bit_kind), intent(in) :: det(N_int,2)
  integer, intent(out) :: no(2), nv(2), npair
  integer, intent(out) :: occ(N_int*bit_kind_size,2)
  integer, intent(out) :: vir(N_int*bit_kind_size,2)
  integer, intent(out) :: pij(hb_n_pair_max), ptype(hb_n_pair_max)
  double precision, intent(out) :: pcum(0:hb_n_pair_max)
  integer(bit_kind) :: hm(N_int,2), pm(N_int,2)
  integer :: s, k, i, j, n, a, b, c
  n = hb_n_orb
  do s=1,2
    do k=1,N_int
      hm(k,s) = iand(det(k,s),generators_bitmask(k,s,s_hole))
      pm(k,s) = iand(not(det(k,s)),generators_bitmask(k,s,s_part))
    enddo
  enddo
  call bitstring_to_list_ab(hm,occ,no,N_int)
  call bitstring_to_list_ab(pm,vir,nv,N_int)
  npair = 0
  pcum(0) = 0.d0
  do s=1,2
    do i=1,no(s)
      a = hb_inv(occ(i,s))
      do j=i+1,no(s)
        b = hb_inv(occ(j,s))
        npair = npair+1
        c = (min(a,b)-1)*n+max(a,b)
        pij(npair) = c
        ptype(npair) = s
        pcum(npair) = pcum(npair-1)+hb_sum(c,2)
      enddo
    enddo
  enddo
  do i=1,no(1)
    a = hb_inv(occ(i,1))
    do j=1,no(2)
      b = hb_inv(occ(j,2))
      npair = npair+1
      c = (a-1)*n+b
      pij(npair) = c
      ptype(npair) = 3
      pcum(npair) = pcum(npair-1)+hb_sum(c,1)
    enddo
  enddo
end

subroutine hb_draw(det,no,occ,nv,vir,npair,pij,ptype,pcum,state,     &
      det2,hij,pgen)
  use bitmasks
  use det_table_module
  implicit none
  BEGIN_DOC
  ! Draws a single or double excitation det2 of det, with probability
  ! pgen. Doubles are drawn with a probability proportional to the
  ! absolute value of the integral. hij = 0 if the excitation is not
  ! possible.
  END_DOC
  integer(bit_kind), intent(in) :: det(N_int,2)
  integer, intent(in) :: no(2), nv(2), npair
  integer, intent(in) :: occ(N_int*bit_kind_size,2)
  integer, intent(in) :: vir(N_int*bit_kind_size,2)
  integer, intent(in) :: pij(hb_n_pair_max), ptype(hb_n_pair_max)
  double precision, intent(in)   :: pcum(0:hb_n_pair_max)
  integer(8), intent(inout)      :: state
  integer(bit_kind), intent(out) :: det2(N_int,2)
  double precision, intent(out)  :: hij, pgen
  integer :: ns(2), nsing, i, j, s, a, b, k, n, ip, t, ij
  integer(8) :: lo, hi, mid, first
  double precision :: r, ps, st
  logical, external :: hb_free

  hij = 0.d0
  pgen = 1.d0
  n = hb_n_orb
  ns(1) = no(1)*nv(1)
  ns(2) = no(2)*nv(2)
  nsing = ns(1)+ns(2)
  st = pcum(npair)
  ps = 0.1d0
  if (nsing == 0) ps = 0.d0
  if (st <= 0.d0) ps = 1.d0
  if (nsing == 0 .and. st <= 0.d0) return
  det2 = det

  if (det_rnd(state) < ps) then
    k = min(int(det_rnd(state)*nsing),nsing-1)
    s = 1
    if (k >= ns(1)) then
      s = 2
      k = k-ns(1)
    endif
    i = k/nv(s)+1
    a = k-(i-1)*nv(s)+1
    call domain_move(det2,occ(i,s),vir(a,s),s)
    call i_H_j_single_spin(det,det2,N_int,s,hij)
    pgen = ps/dble(nsing)
    return
  endif

  r = det_rnd(state)*st
  call domain_search(pcum,npair,r,ip)
  ij = pij(ip)
  t = 1
  if (ptype(ip) < 3) t = 2
  lo = hb_start(ij-1,t)+1_8
  hi = hb_start(ij,t)
  first = lo
  if (hi < lo) return
  r = det_rnd(state)*hb_sum(ij,t)
  do while (lo < hi)
    mid = (lo+hi)/2_8
    if (hb_cum(mid,t) > r) then
      hi = mid
    else
      lo = mid+1_8
    endif
  enddo
  pgen = hb_cum(lo,t)
  if (lo > first) pgen = pgen - hb_cum(lo-1_8,t)
  pgen = (1.d0-ps)*pgen/st
  if (pgen <= 0.d0) return
  a = (hb_ab(lo,t)-1)/n+1
  b = hb_ab(lo,t)-(a-1)*n
  i = (ij-1)/n+1
  j = ij-(i-1)*n
  a = hb_orb(a)
  b = hb_orb(b)
  i = hb_orb(i)
  j = hb_orb(j)
  if (ptype(ip) == 3) then
    if (.not.hb_free(det,a,1)) return
    if (.not.hb_free(det,b,2)) return
    call domain_move(det2,i,a,1)
    call domain_move(det2,j,b,2)
    call i_H_j_double_alpha_beta(det,det2,N_int,hij)
  else
    s = ptype(ip)
    if (.not.hb_free(det,a,s)) return
    if (.not.hb_free(det,b,s)) return
    call domain_move(det2,i,a,s)
    call domain_move(det2,j,b,s)
    call i_H_j_double_spin(det(1,s),det2(1,s),N_int,hij)
  endif
end

logical function hb_free(det,a,s)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! |true| if the orbital a can receive an electron of spin s
  END_DOC
  integer(bit_kind), intent(in) :: det(N_int,2)
  integer, intent(in) :: a, s
  integer :: k, l
  k = shiftr(a-1,bit_kind_shift)+1
  l = iand(a-1,bit_kind_size-1)
  hb_free = (.not.btest(det(k,s),l)).and.                            &
      btest(generators_bitmask(k,s,s_part),l)
end
