use map_module

!! AO Map
!! ======

BEGIN_PROVIDER [ type(map_type), ao_integrals_map ]
  implicit none
  BEGIN_DOC
  ! AO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call two_e_integrals_index(ao_num,ao_num,ao_num,ao_num,key_max)
  sze = key_max
  call map_init(ao_integrals_map,sze)
  print*,  'AO map initialized : ', sze
END_PROVIDER

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




 BEGIN_PROVIDER [ integer, ao_integrals_cache_min ]
&BEGIN_PROVIDER [ integer, ao_integrals_cache_max ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the AOs for which the integrals are in the cache
 END_DOC
 ao_integrals_cache_min = max(1,ao_num - 63)
 ao_integrals_cache_max = ao_num

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_integrals_cache, (0:64*64*64*64) ]
 implicit none
 BEGIN_DOC
 ! Cache of AO integrals for fast access
 END_DOC
 PROVIDE ao_two_e_integrals_in_map
 integer                        :: i,j,k,l,ii
 integer(key_kind)              :: idx, idx2
 real(integral_kind)            :: integral
 real(integral_kind)            :: tmp_re, tmp_im
 integer(key_kind)              :: idx_re,idx_im

  !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx,ii,integral)
  do l=ao_integrals_cache_min,ao_integrals_cache_max
    do k=ao_integrals_cache_min,ao_integrals_cache_max
      do j=ao_integrals_cache_min,ao_integrals_cache_max
        do i=ao_integrals_cache_min,ao_integrals_cache_max
          !DIR$ FORCEINLINE
          call two_e_integrals_index(i,j,k,l,idx)
          !DIR$ FORCEINLINE
          call map_get(ao_integrals_map,idx,integral)
          ii = l-ao_integrals_cache_min
          ii = ior( shiftl(ii,6), k-ao_integrals_cache_min)
          ii = ior( shiftl(ii,6), j-ao_integrals_cache_min)
          ii = ior( shiftl(ii,6), i-ao_integrals_cache_min)
          ao_integrals_cache(ii) = integral
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO
END_PROVIDER

! ---

double precision function get_ao_two_e_integral(i, j, k, l, map) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one AO bi-electronic integral from the AO map in PHYSICIST NOTATION
  !
  ! <1:k, 2:l |1:i, 2:j> 
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  real(integral_kind)            :: tmp
  logical, external              :: ao_two_e_integral_zero
  PROVIDE ao_two_e_integrals_in_map ao_integrals_cache ao_integrals_cache_min
  !DIR$ FORCEINLINE
  if (ao_two_e_integral_zero(i,j,k,l)) then
    tmp = 0.d0
  else
    ii = l-ao_integrals_cache_min
    ii = ior(ii, k-ao_integrals_cache_min)
    ii = ior(ii, j-ao_integrals_cache_min)
    ii = ior(ii, i-ao_integrals_cache_min)
    if (iand(ii, -64) /= 0) then
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,idx)
      !DIR$ FORCEINLINE
      call map_get(map,idx,tmp)
    else
      ii = l-ao_integrals_cache_min
      ii = ior( shiftl(ii,6), k-ao_integrals_cache_min)
      ii = ior( shiftl(ii,6), j-ao_integrals_cache_min)
      ii = ior( shiftl(ii,6), i-ao_integrals_cache_min)
      tmp = ao_integrals_cache(ii)
    endif
  endif
  result = tmp
end

BEGIN_PROVIDER [ complex*16, ao_integrals_cache_periodic, (0:64*64*64*64) ]
 implicit none
 BEGIN_DOC
 ! Cache of AO integrals for fast access
 END_DOC
 PROVIDE ao_two_e_integrals_in_map
 integer                        :: i,j,k,l,ii
 integer(key_kind)              :: idx1, idx2
 real(integral_kind)            :: tmp_re, tmp_im
 integer(key_kind)              :: idx_re,idx_im
 complex(integral_kind)         :: integral


 !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx1,idx2,tmp_re,tmp_im,idx_re,idx_im,ii,integral)
 do l=ao_integrals_cache_min,ao_integrals_cache_max
   do k=ao_integrals_cache_min,ao_integrals_cache_max
     do j=ao_integrals_cache_min,ao_integrals_cache_max
       do i=ao_integrals_cache_min,ao_integrals_cache_max
         !DIR$ FORCEINLINE
         call two_e_integrals_index_2fold(i,j,k,l,idx1)
         !DIR$ FORCEINLINE
         call two_e_integrals_index_2fold(k,l,i,j,idx2)
         idx_re = min(idx1,idx2)
         idx_im = max(idx1,idx2)
         !DIR$ FORCEINLINE
         call map_get(ao_integrals_map,idx_re,tmp_re)
         if (idx_re /= idx_im) then
           call map_get(ao_integrals_map,idx_im,tmp_im)
           if (idx1 < idx2) then
             integral = dcmplx(tmp_re,tmp_im)
           else
             integral = dcmplx(tmp_re,-tmp_im)
           endif
         else
           tmp_im = 0.d0
           integral = dcmplx(tmp_re,tmp_im)
         endif

         ii = l-ao_integrals_cache_min
         ii = ior( shiftl(ii,6), k-ao_integrals_cache_min)
         ii = ior( shiftl(ii,6), j-ao_integrals_cache_min)
         ii = ior( shiftl(ii,6), i-ao_integrals_cache_min)
         ao_integrals_cache_periodic(ii) = integral
       enddo
     enddo
   enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER


complex*16 function get_ao_two_e_integral_periodic(i,j,k,l,map) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one AO bi-electronic integral from the AO map
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx1,idx2
  real(integral_kind)            :: tmp_re, tmp_im
  integer(key_kind)              :: idx_re,idx_im
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  complex(integral_kind)         :: tmp
  PROVIDE ao_two_e_integrals_in_map ao_integrals_cache_periodic ao_integrals_cache_min
  !DIR$ FORCEINLINE
  logical, external              :: ao_two_e_integral_zero
  if (ao_two_e_integral_zero(i,j,k,l)) then
    tmp = (0.d0,0.d0)
  else
    ii = l-ao_integrals_cache_min
    ii = ior(ii, k-ao_integrals_cache_min)
    ii = ior(ii, j-ao_integrals_cache_min)
    ii = ior(ii, i-ao_integrals_cache_min)
    if (iand(ii, -64) /= 0) then
         !DIR$ FORCEINLINE
         call two_e_integrals_index_2fold(i,j,k,l,idx1)
         !DIR$ FORCEINLINE
         call two_e_integrals_index_2fold(k,l,i,j,idx2)
         idx_re = min(idx1,idx2)
         idx_im = max(idx1,idx2)
         !DIR$ FORCEINLINE
         call map_get(ao_integrals_map,idx_re,tmp_re)
         if (idx_re /= idx_im) then
           call map_get(ao_integrals_map,idx_im,tmp_im)
           if (idx1 < idx2) then
             tmp = dcmplx(tmp_re,tmp_im)
           else
             tmp = dcmplx(tmp_re,-tmp_im)
           endif
         else
           tmp_im = 0.d0
           tmp = dcmplx(tmp_re,tmp_im)
         endif
    else
      ii = l-ao_integrals_cache_min
      ii = ior( shiftl(ii,6), k-ao_integrals_cache_min)
      ii = ior( shiftl(ii,6), j-ao_integrals_cache_min)
      ii = ior( shiftl(ii,6), i-ao_integrals_cache_min)
      tmp = ao_integrals_cache_periodic(ii)
    endif
    result = tmp
  endif
end


subroutine get_ao_two_e_integrals(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All i are retrieved for j,k,l fixed.
  ! physicist convention : <ij|kl>
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)

  integer                        :: i
  integer(key_kind)              :: hash
  logical, external              :: ao_one_e_integral_zero
  PROVIDE ao_two_e_integrals_in_map ao_integrals_map 

  if (ao_one_e_integral_zero(j,l)) then
    out_val = 0.d0
    return
  endif

  double precision :: get_ao_two_e_integral
  do i=1,sze
    out_val(i) = get_ao_two_e_integral(i,j,k,l,ao_integrals_map)
  enddo

end


subroutine get_ao_two_e_integrals_periodic(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All i are retrieved for j,k,l fixed.
  ! physicist convention : <ij|kl>
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  complex(integral_kind), intent(out) :: out_val(sze)

  integer                        :: i
  integer(key_kind)              :: hash
  logical, external              :: ao_one_e_integral_zero
  PROVIDE ao_two_e_integrals_in_map ao_integrals_map

  if (ao_one_e_integral_zero(j,l)) then
    out_val = 0.d0
    return
  endif

  double precision :: get_ao_two_e_integral
  do i=1,sze
    out_val(i) = get_ao_two_e_integral(i,j,k,l,ao_integrals_map)
  enddo

end

subroutine get_ao_two_e_integrals_non_zero(j,k,l,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  integer, intent(out)           :: out_val_index(sze),non_zero_int

  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: tmp
  logical, external              :: ao_one_e_integral_zero
  logical, external              :: ao_two_e_integral_zero
  PROVIDE ao_two_e_integrals_in_map

  non_zero_int = 0
  if (ao_one_e_integral_zero(j,l)) then
    out_val = 0.d0
    return
  endif

  non_zero_int = 0
  do i=1,sze
    integer, external :: ao_l4
    double precision, external :: ao_two_e_integral
    !DIR$ FORCEINLINE
    if (ao_two_e_integral_zero(i,j,k,l)) then
      cycle
    endif
    call two_e_integrals_index(i,j,k,l,hash)
    call map_get(ao_integrals_map, hash,tmp)
    if (dabs(tmp) < ao_integrals_threshold) cycle
    non_zero_int = non_zero_int+1
    out_val_index(non_zero_int) = i
    out_val(non_zero_int) = tmp
  enddo

end


subroutine get_ao_two_e_integrals_non_zero_jl(j,l,thresh,sze_max,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  double precision, intent(in)   :: thresh
  integer, intent(in)            :: j,l, sze,sze_max
  real(integral_kind), intent(out) :: out_val(sze_max)
  integer, intent(out)           :: out_val_index(2,sze_max),non_zero_int

  integer                        :: i,k
  integer(key_kind)              :: hash
  double precision               :: tmp
  logical, external              :: ao_one_e_integral_zero
  logical, external              :: ao_two_e_integral_zero

  PROVIDE ao_two_e_integrals_in_map
  non_zero_int = 0
  if (ao_one_e_integral_zero(j,l)) then
    out_val = 0.d0
    return
  endif

  non_zero_int = 0
  do k = 1, sze
   do i = 1, sze
     integer, external :: ao_l4
     double precision, external :: ao_two_e_integral
     !DIR$ FORCEINLINE
     if (ao_two_e_integral_zero(i,j,k,l)) then
       cycle
     endif
     call two_e_integrals_index(i,j,k,l,hash)
     call map_get(ao_integrals_map, hash,tmp)
     if (dabs(tmp) < thresh ) cycle
     non_zero_int = non_zero_int+1
     out_val_index(1,non_zero_int) = i
     out_val_index(2,non_zero_int) = k
     out_val(non_zero_int) = tmp
   enddo
  enddo

end


subroutine get_ao_two_e_integrals_non_zero_jl_from_list(j,l,thresh,list,n_list,sze_max,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO two-electron integrals from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  double precision, intent(in)   :: thresh
  integer, intent(in)            :: sze_max
  integer, intent(in)            :: j,l, n_list,list(2,sze_max)
  real(integral_kind), intent(out) :: out_val(sze_max)
  integer, intent(out)           :: out_val_index(2,sze_max),non_zero_int

  integer                        :: i,k
  integer(key_kind)              :: hash
  double precision               :: tmp
  logical, external              :: ao_one_e_integral_zero
  logical, external              :: ao_two_e_integral_zero

  PROVIDE ao_two_e_integrals_in_map
  non_zero_int = 0
  if (ao_one_e_integral_zero(j,l)) then
    out_val = 0.d0
    return
  endif

  non_zero_int = 0
 integer :: kk
  do kk = 1, n_list
   k = list(1,kk)
   i = list(2,kk)
   integer, external :: ao_l4
   double precision, external :: ao_two_e_integral
   !DIR$ FORCEINLINE
   if (ao_two_e_integral_zero(i,j,k,l)) then
     cycle
   endif
   call two_e_integrals_index(i,j,k,l,hash)
   call map_get(ao_integrals_map, hash,tmp)
   if (dabs(tmp) < thresh ) cycle
   non_zero_int = non_zero_int+1
   out_val_index(1,non_zero_int) = i
   out_val_index(2,non_zero_int) = k
   out_val(non_zero_int) = tmp
  enddo

end




function get_ao_map_size()
  implicit none
  integer (map_size_kind) :: get_ao_map_size
  BEGIN_DOC
  ! Returns the number of elements in the AO map
  END_DOC
  get_ao_map_size = ao_integrals_map % n_elements
end

subroutine clear_ao_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the AO map
  END_DOC
  call map_deinit(ao_integrals_map)
  FREE ao_integrals_map
end


subroutine insert_into_ao_integrals_map(n_integrals,buffer_i, buffer_values)
  use map_module
  implicit none
  BEGIN_DOC
  ! Create new entry into AO map
  END_DOC

  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)

  call map_append(ao_integrals_map, buffer_i, buffer_values, n_integrals)
end


