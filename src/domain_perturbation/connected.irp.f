use bitmasks

BEGIN_PROVIDER [ integer, n_connected_max ]
 implicit none
 BEGIN_DOC
 ! Upper bound for the number of singles and doubles of a
 ! determinant, with the orbitals of generators_bitmask.
 END_DOC
 integer :: nh, np, no(2), s
 integer :: occ(N_int*bit_kind_size,2)
 call bitstring_to_list_ab(generators_bitmask(1,1,s_hole),occ,no,N_int)
 integer :: nhole(2), npart(2)
 do s=1,2
   nhole(s) = min(no(s), elec_alpha_num)
   if (s == 2) nhole(s) = min(no(s), elec_beta_num)
 enddo
 call bitstring_to_list_ab(generators_bitmask(1,1,s_part),occ,no,N_int)
 do s=1,2
   npart(s) = no(s)
 enddo
 n_connected_max = 0
 do s=1,2
   nh = nhole(s)
   np = npart(s)
   n_connected_max += nh*np + ((nh*(nh-1))/2)*((np*(np-1))/2)
 enddo
 n_connected_max += nhole(1)*npart(1)*nhole(2)*npart(2)
 n_connected_max = max(n_connected_max,1)
END_PROVIDER

BEGIN_PROVIDER [ double precision, domain_hij_threshold ]
 implicit none
 BEGIN_DOC
 ! Matrix elements of H smaller than this value are neglected in
 ! the domain perturbation theory.
 END_DOC
 domain_hij_threshold = 1.d-12
END_PROVIDER

subroutine domain_connected(det,n,list,h)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! All the determinants connected to det by a single or double
  ! excitation allowed by generators_bitmask, with
  ! $|\langle l |H| k \rangle|$ above the threshold.
  ! The number of determinants is n.
  END_DOC
  integer(bit_kind), intent(in)  :: det(N_int,2)
  integer, intent(out)           :: n
  integer(bit_kind), intent(out) :: list(N_int,2,n_connected_max)
  double precision, intent(out)  :: h(n_connected_max)

  integer(bit_kind) :: hm(N_int,2), pm(N_int,2), d1(N_int,2), d2(N_int,2)
  integer :: ho(N_int*bit_kind_size,2), pa(N_int*bit_kind_size,2)
  integer :: nh(2), np(2)
  integer :: s, i, j, a, b, k
  double precision :: hij

  do s=1,2
    do k=1,N_int
      hm(k,s) = iand(det(k,s),generators_bitmask(k,s,s_hole))
      pm(k,s) = iand(not(det(k,s)),generators_bitmask(k,s,s_part))
    enddo
  enddo
  call bitstring_to_list_ab(hm,ho,nh,N_int)
  call bitstring_to_list_ab(pm,pa,np,N_int)

  n = 0
  ! Singles
  do s=1,2
    do i=1,nh(s)
      do a=1,np(s)
        d1 = det
        call domain_move(d1,ho(i,s),pa(a,s),s)
        call i_H_j_single_spin(det,d1,N_int,s,hij)
        if (dabs(hij) > domain_hij_threshold) then
          n = n+1
          list(:,:,n) = d1(:,:)
          h(n) = hij
        endif
      enddo
    enddo
  enddo

  ! Same-spin doubles
  do s=1,2
    do i=1,nh(s)
      do j=i+1,nh(s)
        do a=1,np(s)
          do b=a+1,np(s)
            d1 = det
            call domain_move(d1,ho(i,s),pa(a,s),s)
            call domain_move(d1,ho(j,s),pa(b,s),s)
            call i_H_j_double_spin(det(1,s),d1(1,s),N_int,hij)
            if (dabs(hij) > domain_hij_threshold) then
              n = n+1
              list(:,:,n) = d1(:,:)
              h(n) = hij
            endif
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Opposite-spin doubles
  do i=1,nh(1)
    do a=1,np(1)
      d1 = det
      call domain_move(d1,ho(i,1),pa(a,1),1)
      do j=1,nh(2)
        do b=1,np(2)
          d2 = d1
          call domain_move(d2,ho(j,2),pa(b,2),2)
          call i_H_j_double_alpha_beta(det,d2,N_int,hij)
          if (dabs(hij) > domain_hij_threshold) then
            n = n+1
            list(:,:,n) = d2(:,:)
            h(n) = hij
          endif
        enddo
      enddo
    enddo
  enddo
end

subroutine domain_move(det,ih,ip,s)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Moves an electron of spin s from orbital ih to orbital ip.
  END_DOC
  integer(bit_kind), intent(inout) :: det(N_int,2)
  integer, intent(in) :: ih, ip, s
  integer :: k
  k = shiftr(ih-1,bit_kind_shift)+1
  det(k,s) = ibclr(det(k,s),iand(ih-1,bit_kind_size-1))
  k = shiftr(ip-1,bit_kind_shift)+1
  det(k,s) = ibset(det(k,s),iand(ip-1,bit_kind_size-1))
end

integer function domain_index_d0(det)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Index in psi_det of a determinant, 0 if it is not in the
  ! reference space.
  END_DOC
  integer(bit_kind), intent(in) :: det(N_int,2)
  integer, external :: get_index_in_psi_det_sorted_bit
  integer :: k
  k = get_index_in_psi_det_sorted_bit(det,N_int)
  domain_index_d0 = 0
  if (k > 0) domain_index_d0 = psi_coef_sorted_bit_order(k)
end

subroutine domain_build(det,hdiag,nd,dom,list,h)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Domain of an external determinant: the determinant itself,
  ! followed by at most n_det_domain-1 determinants
  !
  ! * connected to it by H,
  ! * not in the reference space,
  ! * with a lower diagonal element of H.
  !
  ! The determinants with the lowest diagonal elements are kept.
  ! Ties (1.d-9) are broken by comparing the bit strings, so the
  ! domain does not depend on the order of the determinants.
  ! list and h are work arrays.
  END_DOC
  integer(bit_kind), intent(in)  :: det(N_int,2)
  double precision, intent(in)   :: hdiag
  integer, intent(out)           :: nd
  integer(bit_kind), intent(out) :: dom(N_int,2,n_det_domain)
  integer(bit_kind), intent(inout) :: list(N_int,2,n_connected_max)
  double precision, intent(inout)  :: h(n_connected_max)

  double precision, external :: diag_H_mat_elem_fock
  integer, external          :: domain_index_d0
  double precision :: fock(2,mo_num+1), x, val(n_det_domain)
  integer :: n, i, p
  logical, external :: domain_lower

  nd = 1
  dom(:,:,1) = det(:,:)
  val(1) = hdiag
  if (n_det_domain == 1) return

  call build_fock_tmp(fock,det,N_int)
  call domain_connected(det,n,list,h)
  do i=1,n
    if (dabs(h(i)) < 1.d-11) cycle
    x = diag_H_mat_elem_fock(det,list(1,1,i),fock,N_int)
    if (x > hdiag-1.d-9) cycle
    if (nd == n_det_domain) then
      if (.not.domain_lower(x,list(1,1,i),val(nd),dom(1,1,nd))) cycle
    endif
    if (domain_index_d0(list(1,1,i)) > 0) cycle
    if (nd < n_det_domain) nd = nd+1
    ! Insertion
    p = nd
    do while (p > 2)
      if (.not.domain_lower(x,list(1,1,i),val(p-1),dom(1,1,p-1))) exit
      dom(:,:,p) = dom(:,:,p-1)
      val(p) = val(p-1)
      p = p-1
    enddo
    dom(:,:,p) = list(:,:,i)
    val(p) = x
  enddo
end

logical function domain_lower(xa,da,xb,db)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Order used to choose the determinants of a domain: lowest
  ! diagonal element, then comparison of the bit strings.
  END_DOC
  double precision, intent(in)  :: xa, xb
  integer(bit_kind), intent(in) :: da(N_int,2), db(N_int,2)
  integer :: ii, ss
  domain_lower = .False.
  if (xa < xb-1.d-9) then
    domain_lower = .True.
  else if (xa <= xb+1.d-9) then
    do ss=1,2
      do ii=N_int,1,-1
        if (da(ii,ss) /= db(ii,ss)) then
          domain_lower = da(ii,ss) < db(ii,ss)
          return
        endif
      enddo
    enddo
  endif
end

subroutine domain_resolvent(e_ref,nd,dom,njet,g,ldg,hdiag,ok)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Row of the resolvent of the domain, for its first determinant:
  ! $g_l(E) = \langle 1 | (P H P - E)^{-1} | l \rangle$, expanded as
  ! $\sum_m g_l^{(m)} (E - E_{ref})^m$.
  ! $g^{(0)} = A^{-1} e_1$, $g^{(m)} = A^{-1} g^{(m-1)}$, with
  ! $A = P H P - E_{ref}$.
  ! ok is |false| if the domain matrix is singular (an external
  ! determinant degenerate with the state): the determinant is then
  ! skipped by the caller.
  END_DOC
  double precision, intent(in)  :: e_ref
  integer, intent(in)           :: nd, njet, ldg
  integer(bit_kind), intent(in) :: dom(N_int,2,nd)
  double precision, intent(in)  :: hdiag
  double precision, intent(out) :: g(ldg,0:njet)
  logical, intent(out)          :: ok
  double precision :: A(nd,nd), hij
  double precision, external :: diag_H_mat_elem
  integer :: i, j, m, ipiv(nd), info

  ok = .True.
  if (nd == 1) then
    if (dabs(hdiag-e_ref) < domain_singular_threshold) then
      ok = .False.
      return
    endif
    g(1,0) = 1.d0/(hdiag-e_ref)
    do m=1,njet
      g(1,m) = g(1,m-1)*g(1,0)
    enddo
    return
  endif
  do j=1,nd
    do i=1,j-1
      call i_H_j(dom(1,1,i),dom(1,1,j),N_int,hij)
      A(i,j) = hij
      A(j,i) = hij
    enddo
    A(j,j) = diag_H_mat_elem(dom(1,1,j),N_int) - e_ref
  enddo
  call dgetrf(nd,nd,A,nd,ipiv,info)
  if (info /= 0) then
    ok = .False.
    return
  endif
  do i=1,nd
    if (dabs(A(i,i)) < domain_singular_threshold) then
      ok = .False.
      return
    endif
  enddo
  g(1:nd,0) = 0.d0
  g(1,0) = 1.d0
  do m=0,njet
    if (m > 0) g(1:nd,m) = g(1:nd,m-1)
    call dgetrs('N',nd,1,A,nd,ipiv,g(1,m),ldg,info)
  enddo
end

BEGIN_PROVIDER [ double precision, domain_singular_threshold ]
 implicit none
 BEGIN_DOC
 ! External determinants whose domain matrix $P_j(H-E)P_j$ has a pivot
 ! smaller than this value (determinants degenerate with the state) are
 ! skipped, and counted in domain_n_singular.
 END_DOC
 domain_singular_threshold = 1.d-8
END_PROVIDER

BEGIN_PROVIDER [ integer(8), domain_n_singular ]
 implicit none
 BEGIN_DOC
 ! Number of external determinants skipped because their domain matrix
 ! is singular
 END_DOC
 domain_n_singular = 0_8
END_PROVIDER
