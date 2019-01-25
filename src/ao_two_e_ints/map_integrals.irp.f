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
  integer, intent(out)           :: i(8),j(8),k(8),l(8)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind)              :: i2,i3
  i = 0
  i2   = ceiling(0.5d0*(dsqrt(8.d0*dble(i1)+1.d0)-1.d0))
  l(1) = ceiling(0.5d0*(dsqrt(8.d0*dble(i2)+1.d0)-1.d0))
  i3   = i1 - shiftr(i2*i2-i2,1)
  k(1) = ceiling(0.5d0*(dsqrt(8.d0*dble(i3)+1.d0)-1.d0))
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
  do ii=1,8
    if (i(ii) /= 0) then
      call two_e_integrals_index(i(ii),j(ii),k(ii),l(ii),i2)
      if (i1 /= i2) then
        print *,  i1, i2
        print *,  i(ii), j(ii), k(ii), l(ii)
        stop 'two_e_integrals_index_reverse failed'
      endif
    endif
  enddo


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
 integer(key_kind)              :: idx
 real(integral_kind)            :: integral
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


double precision function get_ao_two_e_integral(i,j,k,l,map) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one AO bi-electronic integral from the AO map
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  real(integral_kind)            :: tmp
  PROVIDE ao_two_e_integrals_in_map ao_integrals_cache ao_integrals_cache_min
  !DIR$ FORCEINLINE
  if (ao_overlap_abs(i,k)*ao_overlap_abs(j,l) < ao_integrals_threshold ) then
    tmp = 0.d0
  else if (ao_two_e_integral_schwartz(i,k)*ao_two_e_integral_schwartz(j,l) < ao_integrals_threshold) then
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


subroutine get_ao_two_e_integrals(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All i are retrieved for j,k,l fixed.
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)

  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh
  PROVIDE ao_two_e_integrals_in_map ao_integrals_map
  thresh = ao_integrals_threshold

  if (ao_overlap_abs(j,l) < thresh) then
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
  double precision               :: thresh,tmp
  PROVIDE ao_two_e_integrals_in_map
  thresh = ao_integrals_threshold

  non_zero_int = 0
  if (ao_overlap_abs(j,l) < thresh) then
    out_val = 0.d0
    return
  endif

  non_zero_int = 0
  do i=1,sze
    integer, external :: ao_l4
    double precision, external :: ao_two_e_integral
    !DIR$ FORCEINLINE
    if (ao_two_e_integral_schwartz(i,k)*ao_two_e_integral_schwartz(j,l) < thresh) then
      cycle
    endif
    call two_e_integrals_index(i,j,k,l,hash)
    call map_get(ao_integrals_map, hash,tmp)
    if (dabs(tmp) < thresh ) cycle
    non_zero_int = non_zero_int+1
    out_val_index(non_zero_int) = i
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


subroutine dump_ao_integrals(filename)
  use map_module
  implicit none
  BEGIN_DOC
  ! Save to disk the |AO| integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer*8                      :: i,j, n
  if (.not.mpi_master) then
    return
  endif
  call ezfio_set_work_empty(.False.)
  open(unit=66,file=filename,FORM='unformatted')
  write(66) integral_kind, key_kind
  write(66) ao_integrals_map%sorted, ao_integrals_map%map_size,    &
      ao_integrals_map%n_elements
  do i=0_8,ao_integrals_map%map_size
    write(66) ao_integrals_map%map(i)%sorted, ao_integrals_map%map(i)%map_size,&
        ao_integrals_map%map(i)%n_elements
  enddo
  do i=0_8,ao_integrals_map%map_size
    key => ao_integrals_map%map(i)%key
    val => ao_integrals_map%map(i)%value
    n = ao_integrals_map%map(i)%n_elements
    write(66) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  close(66)

end


integer function load_ao_integrals(filename)
  implicit none
  BEGIN_DOC
  ! Read from disk the |AO| integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer*8                      :: i
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer                        :: iknd, kknd
  integer*8                      :: n, j
  load_ao_integrals = 1
  open(unit=66,file=filename,FORM='unformatted',STATUS='UNKNOWN')
  read(66,err=98,end=98) iknd, kknd
  if (iknd /= integral_kind) then
    print *,  'Wrong integrals kind in file :', iknd
    stop 1
  endif
  if (kknd /= key_kind) then
    print *,  'Wrong key kind in file :', kknd
    stop 1
  endif
  read(66,err=98,end=98) ao_integrals_map%sorted, ao_integrals_map%map_size,&
      ao_integrals_map%n_elements
  do i=0_8, ao_integrals_map%map_size
    read(66,err=99,end=99) ao_integrals_map%map(i)%sorted,          &
        ao_integrals_map%map(i)%map_size, ao_integrals_map%map(i)%n_elements
    call cache_map_reallocate(ao_integrals_map%map(i),ao_integrals_map%map(i)%map_size)
  enddo
  do i=0_8, ao_integrals_map%map_size
    key => ao_integrals_map%map(i)%key
    val => ao_integrals_map%map(i)%value
    n = ao_integrals_map%map(i)%n_elements
    read(66,err=99,end=99) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  call map_sort(ao_integrals_map)
  load_ao_integrals = 0
  return
  99 continue
  call map_deinit(ao_integrals_map)
  98 continue
  stop 'Problem reading ao_integrals_map file in work/'

end

