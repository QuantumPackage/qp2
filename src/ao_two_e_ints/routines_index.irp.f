
subroutine two_e_integrals_index(i,j,k,l,i1)
  use map_module
  implicit none
  BEGIN_DOC
! Gives a unique index for i,j,k,l using permtuation symmetry.
! i <-> k, j <-> l, and (i,k) <-> (j,l) for non-periodic systems
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind), intent(out) :: i1
  integer(key_kind)              :: p,q,r,s,i2
  p = min(i,k)
  r = max(i,k)
  p = p+shiftr(r*r-r,1)
  q = min(j,l)
  s = max(j,l)
  q = q+shiftr(s*s-s,1)
  i1 = min(p,q)
  i2 = max(p,q)
  i1 = i1+shiftr(i2*i2-i2,1)
end



subroutine two_e_integrals_index_reverse(i,j,k,l,i1)
  use map_module
  implicit none
  BEGIN_DOC
! Computes the 4 indices $i,j,k,l$ from a unique index $i_1$.
! For 2 indices $i,j$ and $i \le j$, we have
! $p = i(i-1)/2 + j$.
! The key point is that because $j < i$,
! $i(i-1)/2 < p \le i(i+1)/2$. So $i$ can be found by solving
! $i^2 - i - 2p=0$. One obtains $i=1 + \sqrt{1+8p}/2$
! and $j = p - i(i-1)/2$.
! This rule is applied 3 times. First for the symmetry of the
! pairs (i,k) and (j,l), and then for the symmetry within each pair.
  END_DOC
  integer, intent(out)           :: i(8),j(8),k(8),l(8)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind)              :: i2,i3
  i = 0
  i2   = ceiling(0.5d0*(dsqrt(dble(shiftl(i1,3)+1))-1.d0))
  l(1) = ceiling(0.5d0*(dsqrt(dble(shiftl(i2,3)+1))-1.d0))
  i3   = i1 - shiftr(i2*i2-i2,1)
  k(1) = ceiling(0.5d0*(dsqrt(dble(shiftl(i3,3)+1))-1.d0))
  j(1) = int(i2 - shiftr(l(1)*l(1)-l(1),1),4)
  i(1) = int(i3 - shiftr(k(1)*k(1)-k(1),1),4)

              !ijkl
  i(2) = i(1) !ilkj
  j(2) = l(1)
  k(2) = k(1)
  l(2) = j(1)

  i(3) = k(1) !kjil
  j(3) = j(1)
  k(3) = i(1)
  l(3) = l(1)

  i(4) = k(1) !klij
  j(4) = l(1)
  k(4) = i(1)
  l(4) = j(1)

  i(5) = j(1) !jilk
  j(5) = i(1)
  k(5) = l(1)
  l(5) = k(1)

  i(6) = j(1) !jkli
  j(6) = k(1)
  k(6) = l(1)
  l(6) = i(1)

  i(7) = l(1) !lijk
  j(7) = i(1)
  k(7) = j(1)
  l(7) = k(1)

  i(8) = l(1) !lkji
  j(8) = k(1)
  k(8) = j(1)
  l(8) = i(1)

  integer :: ii, jj
  do ii=2,8
    do jj=1,ii-1
      if ( (i(ii) == i(jj)).and. &
           (j(ii) == j(jj)).and. &
           (k(ii) == k(jj)).and. &
           (l(ii) == l(jj)) ) then
         i(ii) = 0
         exit
      endif
    enddo
  enddo
! This has been tested with up to 1000 AOs, and all the reverse indices are
! correct ! We can remove the test
!    do ii=1,8
!      if (i(ii) /= 0) then
!        call two_e_integrals_index(i(ii),j(ii),k(ii),l(ii),i2)
!        if (i1 /= i2) then
!          print *,  i1, i2
!          print *,  i(ii), j(ii), k(ii), l(ii)
!          stop 'two_e_integrals_index_reverse failed'
!        endif
!      endif
!    enddo


end



subroutine ao_idx2_sq(i,j,ij)
  implicit none
  integer, intent(in)  :: i,j
  integer, intent(out) :: ij
  if (i<j) then
    ij=(j-1)*(j-1)+2*i-mod(j+1,2)
  else if (i>j) then
    ij=(i-1)*(i-1)+2*j-mod(i,2)
  else
    ij=i*i
  endif
end

subroutine idx2_tri_int(i,j,ij)
  implicit none
  integer, intent(in)  :: i,j
  integer, intent(out) :: ij
  integer :: p,q
  p = max(i,j)
  q = min(i,j)
  ij = q+ishft(p*p-p,-1)
end

subroutine ao_idx2_tri_key(i,j,ij)
  use map_module
  implicit none
  integer, intent(in)  :: i,j
  integer(key_kind), intent(out) :: ij
  integer(key_kind) :: p,q
  p = max(i,j)
  q = min(i,j)
  ij = q+ishft(p*p-p,-1)
end

subroutine two_e_integrals_index_2fold(i,j,k,l,i1)
  use map_module
  implicit none
  integer, intent(in)            :: i,j,k,l
  integer(key_kind), intent(out) :: i1
  integer                        :: ik,jl

  call ao_idx2_sq(i,k,ik)
  call ao_idx2_sq(j,l,jl)
  call ao_idx2_tri_key(ik,jl,i1)
end

subroutine ao_idx2_sq_rev(i,k,ik)
  BEGIN_DOC
  ! reverse square compound index
  END_DOC
!  p = ceiling(dsqrt(dble(ik)))
!  q = ceiling(0.5d0*(dble(ik)-dble((p-1)*(p-1))))
!  if (mod(ik,2)==0) then
!    k=p
!    i=q
!  else
!    i=p
!    k=q
!  endif
  integer, intent(in)           :: ik
  integer, intent(out)          :: i,k
  integer                       :: pq(0:1),i1,i2
  pq(0) = ceiling(dsqrt(dble(ik)))
  pq(1) = ceiling(0.5d0*(dble(ik)-dble((pq(0)-1)*(pq(0)-1))))
  i1=mod(ik,2)
  i2=mod(ik+1,2)

  k=pq(i1)
  i=pq(i2)
end

subroutine ao_idx2_tri_rev_key(i,k,ik)
  use map_module
  BEGIN_DOC
  !return i<=k
  END_DOC
  integer(key_kind), intent(in) :: ik
  integer, intent(out)          :: i,k
  integer(key_kind) :: tmp_k
  k = ceiling(0.5d0*(dsqrt(8.d0*dble(ik)+1.d0)-1.d0))
  tmp_k = k
  i = int(ik - ishft(tmp_k*tmp_k-tmp_k,-1))
end

subroutine idx2_tri_rev_int(i,k,ik)
  BEGIN_DOC
  !return i<=k
  END_DOC
  integer, intent(in)           :: ik
  integer, intent(out)          :: i,k
  k = ceiling(0.5d0*(dsqrt(8.d0*dble(ik)+1.d0)-1.d0))
  i = int(ik - ishft(k*k-k,-1))
end

subroutine two_e_integrals_index_reverse_2fold(i,j,k,l,i1)
  use map_module
  implicit none
  integer, intent(out)           :: i(2),j(2),k(2),l(2)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind)              :: i0
  integer                        :: i2,i3
  i = 0
  call ao_idx2_tri_rev_key(i3,i2,i1)

  call ao_idx2_sq_rev(j(1),l(1),i2)
  call ao_idx2_sq_rev(i(1),k(1),i3)

              !ijkl
  i(2) = j(1) !jilk
  j(2) = i(1)
  k(2) = l(1)
  l(2) = k(1)

!  i(3) = k(1) !klij   complex conjugate
!  j(3) = l(1)
!  k(3) = i(1)
!  l(3) = j(1)
!
!  i(4) = l(1) !lkji   complex conjugate
!  j(4) = k(1)
!  k(4) = j(1)
!  l(4) = i(1)

  integer :: ii
  if ( (i(1)==i(2)).and. &
       (j(1)==j(2)).and. &
       (k(1)==k(2)).and. &
       (l(1)==l(2)) ) then
    i(2) = 0
  endif
! This has been tested with up to 1000 AOs, and all the reverse indices are
! correct ! We can remove the test
!  do ii=1,2
!    if (i(ii) /= 0) then
!      call two_e_integrals_index_2fold(i(ii),j(ii),k(ii),l(ii),i0)
!      if (i1 /= i0) then
!        print *,  i1, i0
!        print *,  i(ii), j(ii), k(ii), l(ii)
!        stop 'two_e_integrals_index_reverse_2fold failed'
!      endif
!    endif
!  enddo
end


