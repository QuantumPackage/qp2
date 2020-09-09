use map_module


subroutine idx2_tri_int(i,j,ij)
  implicit none
  integer, intent(in)  :: i,j
  integer, intent(out) :: ij
  integer :: p,q
  p = max(i,j)
  q = min(i,j)
  ij = q+ishft(p*p-p,-1)
end

subroutine idx2_tri_key(i,j,ij)
  use map_module
  implicit none
  integer, intent(in)  :: i,j
  integer(key_kind), intent(out) :: ij
  integer(key_kind) :: p,q
  p = max(i,j)
  q = min(i,j)
  ij = q+ishft(p*p-p,-1)
end
subroutine two_e_integrals_index_complex(i,j,k,l,i1,p,q)
  use map_module
  implicit none
  BEGIN_DOC
! Gives a unique index for i,j,k,l using permtuation symmetry.
! i <-> k, j <-> l, and (i,k) <-> (j,l) 
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind), intent(out) :: i1
  integer(key_kind)              :: r,s,i2
  integer(key_kind),intent(out)   :: p,q
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



subroutine two_e_integrals_index_reverse_complex_1(i,j,k,l,i1)
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
! always returns first set such that i<=k, j<=l, ik<=jl
  END_DOC
  integer, intent(out)           :: i(4),j(4),k(4),l(4)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind)              :: i2,i3
  i = 0
  i2   = ceiling(0.5d0*(dsqrt(dble(shiftl(i1,3)+1))-1.d0))
  l(1) = ceiling(0.5d0*(dsqrt(dble(shiftl(i2,3)+1))-1.d0))
  i3   = i1 - shiftr(i2*i2-i2,1)
  k(1) = ceiling(0.5d0*(dsqrt(dble(shiftl(i3,3)+1))-1.d0))
  j(1) = int(i2 - shiftr(l(1)*l(1)-l(1),1),4)
  i(1) = int(i3 - shiftr(k(1)*k(1)-k(1),1),4)

              !ijkl a+ib
  i(2) = j(1) !jilk a+ib
  j(2) = i(1)
  k(2) = l(1)
  l(2) = k(1)

  i(3) = k(1) !klij a-ib
  j(3) = l(1)
  k(3) = i(1)
  l(3) = j(1)

  i(4) = l(1) !lkji a-ib
  j(4) = k(1)
  k(4) = j(1)
  l(4) = i(1)

  integer :: ii, jj
  do ii=2,4
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
end

subroutine two_e_integrals_index_reverse_complex_2(i,j,k,l,i1)
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
! always returns first set such that k<=i, j<=l, ik<=jl
  END_DOC
  integer, intent(out)           :: i(4),j(4),k(4),l(4)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind)              :: i2,i3
  i = 0
  i2   = ceiling(0.5d0*(dsqrt(dble(shiftl(i1,3)+1))-1.d0))
  l(1) = ceiling(0.5d0*(dsqrt(dble(shiftl(i2,3)+1))-1.d0))
  i3   = i1 - shiftr(i2*i2-i2,1)
  i(1) = ceiling(0.5d0*(dsqrt(dble(shiftl(i3,3)+1))-1.d0))
  j(1) = int(i2 - shiftr(l(1)*l(1)-l(1),1),4)
  k(1) = int(i3 - shiftr(i(1)*i(1)-i(1),1),4)

              !kjil a+ib
  i(2) = j(1) !jkli a+ib
  j(2) = i(1)
  k(2) = l(1)
  l(2) = k(1)

  i(3) = k(1) !ilkj a-ib
  j(3) = l(1)
  k(3) = i(1)
  l(3) = j(1)

  i(4) = l(1) !lijk a-ib
  j(4) = k(1)
  k(4) = j(1)
  l(4) = i(1)

  integer :: ii, jj
  do ii=2,4
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
end


BEGIN_PROVIDER [ complex*16, ao_integrals_cache_complex, (0:64*64*64*64) ]
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
  integer(key_kind)              :: p,q,r,s,ik,jl
  logical :: ilek, jlel, iklejl
  complex*16 :: get_ao_two_e_integral_complex_simple


 !$OMP PARALLEL DO PRIVATE (ilek,jlel,p,q,r,s, ik,jl,iklejl, &
 !$OMP                     i,j,k,l,idx1,idx2,tmp_re,tmp_im,idx_re,idx_im,ii,integral)
 do l=ao_integrals_cache_min,ao_integrals_cache_max
   do k=ao_integrals_cache_min,ao_integrals_cache_max
     do j=ao_integrals_cache_min,ao_integrals_cache_max
       do i=ao_integrals_cache_min,ao_integrals_cache_max
         !DIR$ FORCEINLINE
         integral = get_ao_two_e_integral_complex_simple(i,j,k,l,&
                    ao_integrals_map,ao_integrals_map_2)
         
         ii = l-ao_integrals_cache_min
         ii = ior( shiftl(ii,6), k-ao_integrals_cache_min)
         ii = ior( shiftl(ii,6), j-ao_integrals_cache_min)
         ii = ior( shiftl(ii,6), i-ao_integrals_cache_min)
         ao_integrals_cache_complex(ii) = integral
       enddo
     enddo
   enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER

subroutine ao_two_e_integral_complex_map_idx_sign(i,j,k,l,use_map1,idx,sign) 
  use map_module
  implicit none
  BEGIN_DOC
  ! get position of periodic AO integral <ij|kl> 
  ! use_map1: true if integral is in first ao map, false if integral is in second ao map
  ! idx: position of real part of integral in map (imag part is at idx+1)
  ! sign: sign of imaginary part
  !
  !
  !  for <ab|cd>, conditionals are [a<c, b<d, ac<bd]
  !  last two rows are real (ab==cd)
  ! +---------+---------+---------+---------+---------+---------+---------+---------+---------+
  ! | NEW     | <ij|kl> | <ji|lk> | <kl|ij> | <lk|ji> | <kj|il> | <jk|li> | <il|kj> | <li|jk> |
  ! +---------+---------+---------+---------+---------+---------+---------+---------+---------+
  ! |         |         m1        |         m1*       |         m2        |         m2*       |
  ! +---------+---------+---------+---------+---------+---------+---------+---------+---------+
  ! | <ij|kl> | TTT     | TTF     | FFT     | FFF     | FTT     | TFF     | TFT     | FTF     |
  ! | <ij|il> | 0TT     | T0F     | 0FT     | F0F     |         |         |         |         |
  ! | <ij|kj> | T0T     | 0TF     | F0T     | 0FF     |         |         |         |         |
  ! | <ii|jj> | TT0     |         | FF0     |         | FT0(r)  | TF0(r)  |         |         |
  ! +---------+---------+---------+---------+---------+---------+---------+---------+---------+
  ! | <ij|ij> |         |         |         |         | 00T(r)  | 00F(r)  |         |         |
  ! | <ii|ii> |         |         |         |         | 000     |         |         |         |
  ! +---------+---------+---------+---------+---------+---------+---------+---------+---------+
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind), intent(out) :: idx
  logical, intent(out)           :: use_map1
  double precision, intent(out)  :: sign
  integer(key_kind)              :: p,q,r,s,ik,jl,ij,kl
  !DIR$ FORCEINLINE
  call two_e_integrals_index_complex(i,j,k,l,idx,ik,jl)
  p = min(i,j)
  r = max(i,j)
  ij = p+shiftr(r*r-r,1)
  q = min(k,l)
  s = max(k,l)
  kl = q+shiftr(s*s-s,1)

  idx = 2*idx-1

  if (ij==kl) then !real, J -> map1, K -> map2
    sign=0.d0
    use_map1=.False.
  else
    if (ik.eq.jl) then
      if (i.lt.k) then   !TT0
        sign=1.d0
        use_map1=.True. 
      else               !FF0
        sign=-1.d0
        use_map1=.True.
      endif
    else if (i.eq.k) then
      if (j.lt.l) then   !0T* 
        sign=1.d0
        use_map1=.True.
      else               !0F*
        sign=-1.d0
        use_map1=.True.
      endif
    else if (j.eq.l) then
      if (i.lt.k) then
        sign=1.d0
        use_map1=.True.
      else
        sign=-1.d0
        use_map1=.True.
      endif
    else if ((i.lt.k).eqv.(j.lt.l)) then
      if (i.lt.k) then
        sign=1.d0
        use_map1=.True.
      else
        sign=-1.d0
        use_map1=.True.
      endif
    else
      if ((j.lt.l).eqv.(ik.lt.jl)) then
        sign=1.d0
        use_map1=.False.
      else
        sign=-1.d0
        use_map1=.False.
      endif
    endif
  endif
end

complex*16 function get_ao_two_e_integral_complex_simple(i,j,k,l,map,map2) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one AO bi-electronic integral from the AO map
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx1,idx2,idx
  real(integral_kind)            :: tmp_re, tmp_im
  integer(key_kind)              :: idx_re,idx_im
  type(map_type), intent(inout)  :: map,map2
  integer                        :: ii
  complex(integral_kind)         :: tmp
  integer(key_kind)              :: p,q,r,s,ik,jl
  logical :: ilek, jlel, iklejl,use_map1
  double precision :: sign
  ! a.le.c, b.le.d, tri(a,c).le.tri(b,d)
  PROVIDE ao_two_e_integrals_in_map
  call ao_two_e_integral_complex_map_idx_sign(i,j,k,l,use_map1,idx,sign)
  if (use_map1) then
    call map_get(map,idx,tmp_re)
    call map_get(map,idx+1,tmp_im)
    tmp_im *= sign
  else
    call map_get(map2,idx,tmp_re)
    if (sign/=0.d0) then
      call map_get(map2,idx+1,tmp_im)
      tmp_im *= sign
    else
      tmp_im=0.d0
    endif
  endif
  tmp = dcmplx(tmp_re,tmp_im)
  result = tmp
end


complex*16 function get_ao_two_e_integral_complex(i,j,k,l,map,map2) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one AO bi-electronic integral from the AO map
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx1,idx2
  real(integral_kind)            :: tmp_re, tmp_im
  integer(key_kind)              :: idx_re,idx_im
  type(map_type), intent(inout)  :: map,map2
  integer                        :: ii
  complex(integral_kind)         :: tmp
  complex(integral_kind)         :: get_ao_two_e_integral_complex_simple
  integer(key_kind)              :: p,q,r,s,ik,jl
  logical :: ilek, jlel, iklejl
  ! a.le.c, b.le.d, tri(a,c).le.tri(b,d)
  PROVIDE ao_two_e_integrals_in_map ao_integrals_cache_complex ao_integrals_cache_min
  !DIR$ FORCEINLINE
  !logical, external              :: ao_two_e_integral_zero
  !if (ao_two_e_integral_zero(i,j,k,l)) then
  !  tmp = (0.d0,0.d0)
  !else
  if (.True.) then
    ii = l-ao_integrals_cache_min
    ii = ior(ii, k-ao_integrals_cache_min)
    ii = ior(ii, j-ao_integrals_cache_min)
    ii = ior(ii, i-ao_integrals_cache_min)
    if (iand(ii, -64) /= 0) then
      tmp = get_ao_two_e_integral_complex_simple(i,j,k,l,map,map2)
    else
      ii = l-ao_integrals_cache_min
      ii = ior( shiftl(ii,6), k-ao_integrals_cache_min)
      ii = ior( shiftl(ii,6), j-ao_integrals_cache_min)
      ii = ior( shiftl(ii,6), i-ao_integrals_cache_min)
      tmp = ao_integrals_cache_complex(ii)
    endif
  endif
  result = tmp
end


subroutine get_ao_two_e_integrals_complex(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All i are retrieved for j,k,l fixed.
  ! physicist convention : <ij|kl> 
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  complex*16, intent(out) :: out_val(sze)

  integer                        :: i
  integer(key_kind)              :: hash
  !logical, external              :: ao_one_e_integral_zero
  PROVIDE ao_two_e_integrals_in_map ao_integrals_map

  !if (ao_one_e_integral_zero(j,l)) then
  !  out_val = (0.d0,0.d0)
  !  return
  !endif

  complex*16 :: get_ao_two_e_integral_complex
  do i=1,sze
    out_val(i) = get_ao_two_e_integral_complex(i,j,k,l,ao_integrals_map,ao_integrals_map_2)
  enddo

end

subroutine get_ao_two_e_integrals_non_zero_complex(j,k,l,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  complex(integral_kind), intent(out) :: out_val(sze)
  integer, intent(out)           :: out_val_index(sze),non_zero_int
  print*,'not implemented for periodic',irp_here
  stop -1
  !placeholder to keep compiler from complaining about out values not assigned
  out_val=0.d0
  out_val_index=0
  non_zero_int=0
!
!  integer                        :: i
!  integer(key_kind)              :: hash
!  double precision               :: thresh,tmp
!  if(is_complex) then
!    print*,'not implemented for periodic:',irp_here
!    stop -1
!  endif
!  PROVIDE ao_two_e_integrals_in_map
!  thresh = ao_integrals_threshold
!
!  non_zero_int = 0
!  if (ao_overlap_abs(j,l) < thresh) then
!    out_val = 0.d0
!    return
!  endif
!
!  non_zero_int = 0
!  do i=1,sze
!    integer, external :: ao_l4
!    double precision, external :: ao_two_e_integral
!    !DIR$ FORCEINLINE
!    if (ao_two_e_integral_schwartz(i,k)*ao_two_e_integral_schwartz(j,l) < thresh) then
!      cycle
!    endif
!    call two_e_integrals_index(i,j,k,l,hash)
!    call map_get(ao_integrals_map, hash,tmp)
!    if (dabs(tmp) < thresh ) cycle
!    non_zero_int = non_zero_int+1
!    out_val_index(non_zero_int) = i
!    out_val(non_zero_int) = tmp
!  enddo

end


subroutine get_ao_two_e_integrals_non_zero_jl_complex(j,l,thresh,sze_max,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  double precision, intent(in)   :: thresh
  integer, intent(in)            :: j,l, sze,sze_max
  complex(integral_kind), intent(out) :: out_val(sze_max)
  integer, intent(out)           :: out_val_index(2,sze_max),non_zero_int
  print*,'not implemented for periodic',irp_here
  stop -1
  !placeholder to keep compiler from complaining about out values not assigned
  out_val=0.d0
  out_val_index=0
  non_zero_int=0
!
!  integer                        :: i,k
!  integer(key_kind)              :: hash
!  double precision               :: tmp
!
!  if(is_complex) then
!    print*,'not implemented for periodic:',irp_here
!    stop -1
!  endif
!  PROVIDE ao_two_e_integrals_in_map
!  non_zero_int = 0
!  if (ao_overlap_abs(j,l) < thresh) then
!    out_val = 0.d0
!    return
!  endif
!
!  non_zero_int = 0
!  do k = 1, sze
!   do i = 1, sze
!     integer, external :: ao_l4
!     double precision, external :: ao_two_e_integral
!     !DIR$ FORCEINLINE
!     if (ao_two_e_integral_schwartz(i,k)*ao_two_e_integral_schwartz(j,l) < thresh) then
!       cycle
!     endif
!     call two_e_integrals_index(i,j,k,l,hash)
!     call map_get(ao_integrals_map, hash,tmp)
!     if (dabs(tmp) < thresh ) cycle
!     non_zero_int = non_zero_int+1
!     out_val_index(1,non_zero_int) = i
!     out_val_index(2,non_zero_int) = k
!     out_val(non_zero_int) = tmp
!   enddo
!  enddo

end


subroutine get_ao_two_e_integrals_non_zero_jl_from_list_complex(j,l,thresh,list,n_list,sze_max,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO two-electron integrals from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  double precision, intent(in)   :: thresh
  integer, intent(in)            :: sze_max
  integer, intent(in)            :: j,l, n_list,list(2,sze_max)
  complex(integral_kind), intent(out) :: out_val(sze_max)
  integer, intent(out)           :: out_val_index(2,sze_max),non_zero_int
  print*,'not implemented for periodic',irp_here
  stop -1
  !placeholder to keep compiler from complaining about out values not assigned
  out_val=0.d0
  out_val_index=0
  non_zero_int=0
!
!  integer                        :: i,k
!  integer(key_kind)              :: hash
!  double precision               :: tmp
!
!  if(is_complex) then
!    print*,'not implemented for periodic:',irp_here
!    stop -1
!  endif
!  PROVIDE ao_two_e_integrals_in_map
!  non_zero_int = 0
!  if (ao_overlap_abs(j,l) < thresh) then
!    out_val = 0.d0
!    return
!  endif
!
!  non_zero_int = 0
! integer :: kk
!  do kk = 1, n_list
!   k = list(1,kk)
!   i = list(2,kk)
!   integer, external :: ao_l4
!   double precision, external :: ao_two_e_integral
!   !DIR$ FORCEINLINE
!   if (ao_two_e_integral_schwartz(i,k)*ao_two_e_integral_schwartz(j,l) < thresh) then
!     cycle
!   endif
!   call two_e_integrals_index(i,j,k,l,hash)
!   call map_get(ao_integrals_map, hash,tmp)
!   if (dabs(tmp) < thresh ) cycle
!   non_zero_int = non_zero_int+1
!   out_val_index(1,non_zero_int) = i
!   out_val_index(2,non_zero_int) = k
!   out_val(non_zero_int) = tmp
!  enddo

end

subroutine insert_into_ao_integrals_map_2(n_integrals,buffer_i, buffer_values)
  use map_module
  implicit none
  BEGIN_DOC
  ! Create new entry into AO map
  END_DOC

  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)

  call map_append(ao_integrals_map_2, buffer_i, buffer_values, n_integrals)
end


