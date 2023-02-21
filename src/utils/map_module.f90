module map_module

! A map is an array of maps (cache_maps)
! A cache map is an array of keys and values sorted by keys
! A cache map has its own OpenMP lock
! To access a (key,value) pair in the map, the
! index of the cache_map in the map array is obtained
! by removing the first 15 bits of the key.
! The key in the cache_map is composed of the first
! 15 bits of the key. Therefore, it can be stored
! as integer*2 and is found by applying the map_mask
! to the initial key. The element are found in the
! cache_map using a binary search
!
! When using the map_update subroutine to build the map,
! the map_merge subroutine
! should be called before getting data from the map.

 use omp_lib

 integer, parameter             :: integral_kind = 8

 integer, parameter             :: cache_key_kind = 2
 integer, parameter             :: cache_map_size_kind = 4

 integer, parameter             :: key_kind = 8
 integer, parameter             :: map_size_kind = 8

 integer,   parameter           :: map_shift = 15
 integer*8, parameter           :: map_mask  = ibset(0_8,15)-1_8

 type cache_map_type
  real(integral_kind), pointer   :: value(:)
  integer(cache_key_kind), pointer :: key(:)
  logical                        :: sorted
  integer(cache_map_size_kind)   :: map_size
  integer(cache_map_size_kind)   :: n_elements
  integer(omp_lock_kind)         :: lock
 end type cache_map_type

 type map_type
  type(cache_map_type), allocatable :: map(:)
  real(integral_kind), pointer   :: consolidated_value(:)
  integer(cache_key_kind), pointer :: consolidated_key(:)
  integer*8, pointer             :: consolidated_idx(:)
  logical                        :: sorted
  logical                        :: consolidated
  integer(map_size_kind)         :: map_size
  integer(map_size_kind)         :: n_elements
  integer(omp_lock_kind)         :: lock
 end type map_type

end module map_module


double precision function map_mb(map)
  use map_module
  use omp_lib
  implicit none
  type (map_type), intent(in)    :: map
  integer(map_size_kind)         :: i

  map_mb = dble(8+map_size_kind+map_size_kind+omp_lock_kind+4)
  do i=0,map%map_size
    map_mb = map_mb + dble(map%map(i)%map_size*(cache_key_kind+integral_kind) +&
        8+8+4+cache_map_size_kind+cache_map_size_kind+omp_lock_kind)
  enddo
  map_mb = map_mb / (1024.d0*1024.d0)
end

subroutine cache_map_init(map,sze)
  use map_module
  implicit none
  type (cache_map_type), intent(inout) :: map
  integer(cache_map_size_kind)   :: sze
  call omp_set_lock(map%lock)
  map%n_elements = 0_8
  map%map_size = 0_8
  map%sorted = .True.
  NULLIFY(map%value, map%key)
  call cache_map_reallocate(map,sze)
  call omp_unset_lock(map%lock)
end

subroutine map_init(map,keymax)
  use map_module
  implicit none
  integer*8, intent(in)          :: keymax
  type (map_type), intent(inout) :: map
  integer(map_size_kind)         :: i
  integer(cache_map_size_kind)   :: sze
  integer                        :: err

  call omp_init_lock(map%lock)
  call omp_set_lock(map%lock)

  map%n_elements = 0_8
  map%map_size = shiftr(keymax,map_shift)
  map%consolidated = .False.

  allocate(map%map(0_8:map%map_size),stat=err)
  if (err /= 0) then
    print *,  'Unable to allocate map'
    stop 5
  endif
  sze = 2
  do i=0_8,map%map_size
    call omp_init_lock(map%map(i)%lock)
  enddo
  !$OMP PARALLEL DEFAULT(NONE) SHARED(map,sze) PRIVATE(i)
  !$OMP DO SCHEDULE(STATIC,512)
  do i=0_8,map%map_size
    call cache_map_init(map%map(i),sze)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  map%sorted = .True.

  call omp_unset_lock(map%lock)

end

subroutine cache_map_reallocate(map,sze)
  use map_module
  implicit none
  integer(cache_map_size_kind), intent(in) :: sze
  type (cache_map_type), intent(inout) :: map

  integer(cache_key_kind), pointer :: key_new(:)
  real(integral_kind), pointer   :: value_new(:)
  integer(map_size_kind)         :: i
  integer                        :: err
  !DIR$ ATTRIBUTES ALIGN : 64    :: key_new, value_new

  if (sze < map%n_elements) then
    print *,  'Unable to resize map : map too large'
    stop 3
  endif

  ! Resize keys
  allocate( key_new(sze), stat=err )
  if (err /= 0) then
    print *,  'Unable to allocate map', sze
    stop 1
  endif
  if (associated(map%key)) then
    do i=1_8,min(size(map%key),map%n_elements)
      key_new(i) = map%key(i)
    enddo
    deallocate(map%key)
  endif

  ! Resize values
  allocate( value_new(sze), stat=err )
  if (err /= 0) then
    print *,  'Unable to allocate map', sze
    stop 2
  endif
  if (associated(map%value)) then
    do i=1_8,min(size(map%key),map%n_elements)
      value_new(i) = map%value(i)
    enddo
    deallocate(map%value)
  endif

  ! Set new pointers
  map%key => key_new
  map%value => value_new
  map%map_size = sze

end


subroutine cache_map_deinit(map)
  use map_module
  implicit none
  type (cache_map_type), intent(inout) :: map

  integer                        :: err

  if (associated( map % value )) then
    deallocate( map % value, stat=err )
    if (err /= 0) then
      print *,  'Unable to deallocate map'
      stop 2
    endif
    NULLIFY(map%value)
  endif

  if (associated( map % key )) then
    deallocate( map % key, stat=err )
    if (err /= 0) then
      print *,  'Unable to deallocate map'
      stop 4
    endif
    NULLIFY(map%key)
  endif

  map%n_elements = 0_8
  map%map_size   = 0_8
  call omp_destroy_lock(map%lock)
end

subroutine map_deinit(map)
  use map_module
  implicit none
  type (map_type), intent(inout) :: map
  integer                        :: err
  integer(map_size_kind)         :: i

  if (allocated( map % map )) then
    do i=0_8, map%map_size
      call cache_map_deinit(map%map(i))
    enddo
    deallocate( map % map, stat=err )
    if (err /= 0) then
      print *,  'Unable to deallocate map'
      stop 6
    endif
  endif

  map%n_elements = 0_8
  map%map_size   = 0_8
  call omp_destroy_lock(map%lock)
end

subroutine cache_map_sort(map)
  use map_module
  implicit none
  type (cache_map_type), intent(inout) :: map
  integer(cache_map_size_kind), allocatable :: iorder(:)
  integer(cache_map_size_kind)   :: i
  !DIR$ ATTRIBUTES ALIGN : 64    :: iorder

  if (.not.map%sorted) then
    allocate(iorder(map%n_elements))
    do i=1,map%n_elements
      iorder(i) = i
    enddo
    if (cache_key_kind == 2) then
      call i2sort(map%key,iorder,map%n_elements,-1)
    else if (cache_key_kind == 4) then
      call isort(map%key,iorder,map%n_elements,-1)
    else if (cache_key_kind == 8) then
      call i8sort(map%key,iorder,map%n_elements,-1)
    endif
    if (integral_kind == 4) then
      call set_order(map%value,iorder,map%n_elements)
    else if (integral_kind == 8) then
      call dset_order(map%value,iorder,map%n_elements)
    endif
    deallocate(iorder)
    map%sorted = .True.
  endif

end

subroutine map_sort(map)
  use map_module
  implicit none
  type (map_type), intent(inout) :: map
  integer(map_size_kind)         :: i

  if (.not.map%sorted) then
    !$OMP PARALLEL DO SCHEDULE(static,1024) DEFAULT(SHARED) PRIVATE(i)
    do i=0_8,map%map_size
      call omp_set_lock(map%map(i)%lock)
      call cache_map_sort(map%map(i))
      call omp_unset_lock(map%map(i)%lock)
    enddo
    !$OMP END PARALLEL  DO
    map%sorted = .True.
  endif

end

subroutine cache_map_merge(map)
  use map_module
  implicit none
  type (cache_map_type), intent(inout) :: map
  integer(cache_key_kind)        :: prev_key
  integer(cache_map_size_kind)   :: i, j

  call cache_map_sort(map)
  prev_key = -1_8
  j=0
  do i=1,map%n_elements
    if (map%key(i) /= prev_key) then
      j = j+1
      map%value(j) = map%value(i)
      map%key(j) = map%key(i)
      prev_key = map%key(i)
    else
      map%value(j) = map%value(j)+map%value(i)
    endif
  enddo
  map%n_elements = j

end

subroutine cache_map_unique(map)
  use map_module
  implicit none
  type (cache_map_type), intent(inout) :: map
  integer(cache_key_kind)        :: prev_key
  integer(cache_map_size_kind)   :: i, j

  call cache_map_sort(map)
  prev_key = -1_8
  j=0
  do i=1,map%n_elements
    if (map%key(i) /= prev_key) then
      j = j+1
      map%value(j) = map%value(i)
      map%key(j) = map%key(i)
      prev_key = map%key(i)
    endif
  enddo
  map%n_elements = j

end

subroutine cache_map_shrink(map,thr)
  use map_module
  implicit none
  type (cache_map_type), intent(inout) :: map
  real(integral_kind)  , intent(in) :: thr
  integer(cache_map_size_kind)   :: i,j

  j=0
  do i=1,map%n_elements
    if (abs(map%value(i)) > thr) then
      j = j+1
      map%value(j) = map%value(i)
      map%key(j) = map%key(i)
    endif
  enddo
  map%n_elements = j

end

subroutine map_unique(map)
  use map_module
  implicit none
  type (map_type), intent(inout) :: map
  integer(map_size_kind)         :: i
  integer(map_size_kind)         :: icount

  icount = 0_8
  !$OMP PARALLEL DO SCHEDULE(dynamic,1000) DEFAULT(SHARED) PRIVATE(i)&
      !$OMP REDUCTION(+:icount)
  do i=0_8,map%map_size
    call omp_set_lock(map%map(i)%lock)
    call cache_map_unique(map%map(i))
    call omp_unset_lock(map%map(i)%lock)
    icount = icount + map%map(i)%n_elements
  enddo
  !$OMP END PARALLEL DO
  map%n_elements = icount

end

subroutine map_merge(map)
  use map_module
  implicit none
  type (map_type), intent(inout) :: map
  integer(map_size_kind)         :: i
  integer(map_size_kind)         :: icount

  icount = 0_8
  !$OMP PARALLEL DO SCHEDULE(dynamic,1000) DEFAULT(SHARED) PRIVATE(i)&
      !$OMP REDUCTION(+:icount)
  do i=0_8,map%map_size
    call omp_set_lock(map%map(i)%lock)
    call cache_map_merge(map%map(i))
    call omp_unset_lock(map%map(i)%lock)
    icount = icount + map%map(i)%n_elements
  enddo
  !$OMP END PARALLEL DO
  map%n_elements = icount

end

subroutine map_shrink(map,thr)
  use map_module
  implicit none
  type (map_type), intent(inout) :: map
  real(integral_kind), intent(in) :: thr
  integer(map_size_kind)         :: i
  integer(map_size_kind)         :: icount

  icount = 0_8
  !$OMP PARALLEL DO SCHEDULE(dynamic,1000) DEFAULT(SHARED) PRIVATE(i)&
      !$OMP REDUCTION(+:icount)
  do i=0_8,map%map_size
    call omp_set_lock(map%map(i)%lock)
    call cache_map_shrink(map%map(i),thr)
    call omp_unset_lock(map%map(i)%lock)
    icount = icount + map%map(i)%n_elements
  enddo
  !$OMP END PARALLEL DO
  map%n_elements = icount

end

subroutine map_update(map, key, value, sze, thr)
  use map_module
  implicit none
  type (map_type), intent(inout) :: map
  integer, intent(in)            :: sze
  integer(key_kind), intent(inout) :: key(sze)
  real(integral_kind), intent(inout) :: value(sze)
  real(integral_kind), intent(in) :: thr

  integer                        :: i
  integer(map_size_kind)         :: idx_cache, idx_cache_new
  integer(cache_map_size_kind)   :: idx
  integer                        :: sze2
  integer(cache_key_kind)        :: cache_key
  integer(map_size_kind)         :: n_elements_temp
  type (cache_map_type)          :: local_map
  logical                        :: map_sorted

  sze2 = sze
  map_sorted = .True.

  n_elements_temp = 0_8
  n_elements_temp = n_elements_temp + 1_8
  do while (sze2>0)
    i=1
    do while (i<=sze)
      if (key(i) /= 0_8) then
        idx_cache = shiftr(key(i),map_shift)
        if (omp_test_lock(map%map(idx_cache)%lock)) then
          local_map%key => map%map(idx_cache)%key
          local_map%value => map%map(idx_cache)%value
          local_map%sorted = map%map(idx_cache)%sorted
          local_map%map_size = map%map(idx_cache)%map_size
          local_map%n_elements = map%map(idx_cache)%n_elements
          do
          !DIR$ FORCEINLINE
          call search_key_big_interval(key(i),local_map%key, local_map%n_elements, idx, 1, local_map%n_elements)
          if (idx > 0_8) then
            local_map%value(idx) = local_map%value(idx) + value(i)
          else
            ! Assert that the map has a proper size
            if (local_map%n_elements == local_map%map_size) then
              call cache_map_merge(local_map)
              call cache_map_reallocate(local_map, local_map%n_elements + local_map%n_elements)
              call cache_map_shrink(local_map,thr)
            endif
            cache_key = int(iand(key(i),map_mask),2)
            local_map%n_elements = local_map%n_elements + 1
            local_map%value(local_map%n_elements) = value(i)
            local_map%key(local_map%n_elements) = cache_key
            local_map%sorted = .False.
            n_elements_temp = n_elements_temp + 1_8
          endif  ! idx > 0
          key(i) = 0_8
          i = i+1
          sze2 = sze2-1
          if (i>sze) then
            i=1
          endif
          if ( (shiftr(key(i),map_shift) /= idx_cache).or.(key(i)==0_8)) then
            exit
          endif
        enddo
        map%map(idx_cache)%key => local_map%key
        map%map(idx_cache)%value => local_map%value
        map%map(idx_cache)%sorted = local_map%sorted
        map%map(idx_cache)%n_elements = local_map%n_elements
        map%map(idx_cache)%map_size = local_map%map_size
        map_sorted = map_sorted .and. local_map%sorted
        call omp_unset_lock(map%map(idx_cache)%lock)
      endif  ! omp_test_lock
    else
      i=i+1
    endif  ! key = 0
  enddo  ! i
enddo  ! sze2 > 0
call omp_set_lock(map%lock)
map%n_elements = map%n_elements + n_elements_temp
map%sorted = map%sorted .and. map_sorted
call omp_unset_lock(map%lock)

end

subroutine map_append(map, key, value, sze)
  use map_module
  implicit none
  type (map_type), intent(inout) :: map
  integer, intent(in)            :: sze
  integer(key_kind), intent(inout) :: key(sze)
  real(integral_kind), intent(inout) :: value(sze)

  integer                        :: i
  integer(cache_map_size_kind)   :: n_elements
  integer(map_size_kind)         :: idx_cache
  integer(cache_key_kind)        :: cache_key

  do i=1,sze
    idx_cache = shiftr(key(i),map_shift)
    call omp_set_lock(map%map(idx_cache)%lock)
    n_elements = map%map(idx_cache)%n_elements + 1
    ! Assert that the map has a proper size
    if (n_elements == map%map(idx_cache)%map_size) then
      call cache_map_reallocate(map%map(idx_cache), n_elements+ shiftr(n_elements,1))
    endif
    cache_key = int(iand(key(i),map_mask),2)
    map%map(idx_cache)%value(n_elements) = value(i)
    map%map(idx_cache)%key(n_elements) = cache_key
    map%map(idx_cache)%n_elements = n_elements
    if (map%map(idx_cache)%sorted.and.n_elements > 1) then
      map%map(idx_cache)%sorted = (map%map(idx_cache)%key(n_elements-1) <= cache_key)
      map%sorted = map%sorted .and. map%map(idx_cache)%sorted
    endif
    call omp_unset_lock(map%map(idx_cache)%lock)
  enddo
  call omp_set_lock(map%lock)
  map%n_elements = map%n_elements + sze
  call omp_unset_lock(map%lock)

end

subroutine map_get(map, key, value)
  use map_module
  implicit none
  type (map_type), intent(inout) :: map
  integer(key_kind), intent(in)  :: key
  real(integral_kind), intent(out) :: value
  integer(map_size_kind)         :: idx_cache
  integer(cache_map_size_kind)   :: idx

  ! index in tha pointers array
  idx_cache = shiftr(key,map_shift)
  !DIR$ FORCEINLINE
  call cache_map_get_interval(map%map(idx_cache), key, value, 1, map%map(idx_cache)%n_elements,idx)
end

subroutine cache_map_get_interval(map, key, value, ibegin, iend, idx)
  use map_module
  implicit none
  type (cache_map_type), intent(inout) :: map
  integer(key_kind), intent(in)  :: key
  integer(cache_map_size_kind), intent(in) :: ibegin, iend
  real(integral_kind), intent(out) :: value
  integer(cache_map_size_kind), intent(inout) :: idx
  double precision, pointer :: v(:)
  integer :: i

  call search_key_big_interval(key,map%key, map%n_elements, idx, ibegin, iend)
  if (idx > 0) then
    value = map%value(idx)
  else
    value = 0._integral_kind
  endif
!  call search_key_value_big_interval(key, value, map%key, map%value, map%n_elements, idx, ibegin, iend)
end


subroutine map_get_many(map, key, value, sze)
  use map_module
  implicit none
  type (map_type), intent(inout) :: map
  integer, intent(in)            :: sze
  integer(key_kind), intent(in)  :: key(sze)
  real(integral_kind), intent(out) :: value(sze)
  integer                        :: i
  integer(map_size_kind)         :: idx_cache
  integer(cache_map_size_kind)   :: ibegin, iend
  integer(cache_map_size_kind), allocatable :: idx(:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: idx

  allocate(idx(sze))
  do i=1,sze
    idx_cache = shiftr(key(i),map_shift)
    iend = map%map(idx_cache)%n_elements
    !DIR$ FORCEINLINE
    call search_key_big_interval(key(i),map%map(idx_cache)%key, iend, idx(i), 1, iend)
  enddo
  do i=1,sze
    idx_cache = shiftr(key(i),map_shift)
    if (idx(i) > 0) then
      value(i) = map%map(idx_cache)%value(idx(i))
    else
      value(i) = 0.
    endif
  enddo
  deallocate(idx)
end

subroutine map_exists_many(map, key, sze)
  use map_module
  implicit none
  type (map_type), intent(inout) :: map
  integer, intent(in)            :: sze
  integer(key_kind), intent(inout) :: key(sze)
  integer                        :: i
  integer(map_size_kind)         :: idx_cache, idx_cache_prev
  integer(cache_map_size_kind)   :: ibegin, iend
  integer(cache_map_size_kind), allocatable :: idx(:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: idx

  idx_cache_prev = -1_map_size_kind
  allocate(idx(sze))
  do i=1,sze
    idx_cache = shiftr(key(i),map_shift)
    iend = map%map(idx_cache)%n_elements
    if (idx_cache == idx_cache_prev) then
      if ((idx(i-1) > 0_cache_map_size_kind).and.(idx(i-1) < iend)) then
        if ((key(i) == key(i-1)+1).and.(map%map(idx_cache)%key(idx(i-1))+1) == key(i)) then
          idx(i) = idx(i-1)+1
          cycle
        endif
      endif
    endif
    !DIR$ FORCEINLINE
    call search_key_big_interval(key(i),map%map(idx_cache)%key, iend, idx(i), 1, iend)
    idx_cache_prev = idx_cache
  enddo
  do i=1,sze
    idx_cache = shiftr(key(i),map_shift)
    if (idx(i) <= 0) then
      key(i) = 0_key_kind
    endif
  enddo
  deallocate(idx)
end

subroutine search_key_big(key,X,sze,idx)
  use map_module
  implicit none
  integer(cache_map_size_kind), intent(in) :: sze
  integer(key_kind)           , intent(in) :: key
  integer(cache_key_kind)     , intent(in) :: X(sze)
  integer(cache_map_size_kind), intent(out) :: idx

  call  search_key_big_interval(key,X,sze,idx,1,sze)
end


subroutine search_key_big_interval(key,X,sze,idx,ibegin_in,iend_in)
  use map_module
  implicit none
  integer(cache_map_size_kind), intent(in) :: sze
  integer(key_kind)           , intent(in) :: key
  integer(cache_key_kind)     , intent(in) :: X(sze)
  integer(cache_map_size_kind), intent(in) :: ibegin_in, iend_in
  integer(cache_map_size_kind), intent(out) :: idx

  integer(cache_map_size_kind)   :: istep, ibegin, iend, i
  integer(cache_key_kind)        :: cache_key

  if (sze /= 0) then
    continue
  else
    idx = -1
    return
  endif
  cache_key = int(iand(key,map_mask),2)
  ibegin = min(ibegin_in,sze)
  iend   = min(iend_in,sze)
  if ((cache_key > X(ibegin)) .and. (cache_key < X(iend))) then

    istep = shiftr(iend-ibegin,1)
    idx = ibegin + istep
    do while (istep > 4)
      idx = ibegin + istep
      ! TODO : Cache misses
      if (cache_key < X(idx)) then
        iend = idx
        istep = shiftr(idx-ibegin,1)
        idx = ibegin + istep
        if (cache_key < X(idx)) then
          iend = idx
          istep = shiftr(idx-ibegin,1)
          cycle
        else if (cache_key > X(idx)) then
          ibegin = idx
          istep = shiftr(iend-idx,1)
          cycle
        else
          return
        endif
      else if (cache_key > X(idx)) then
        ibegin = idx
        istep = shiftr(iend-idx,1)
        idx = idx + istep
        if (cache_key < X(idx)) then
          iend = idx
          istep = shiftr(idx-ibegin,1)
          cycle
        else if (cache_key > X(idx)) then
          ibegin = idx
          istep = shiftr(iend-idx,1)
          cycle
        else
          return
        endif
      else
        return
      endif
    enddo
    idx = ibegin
    if (min(iend_in,sze) > ibegin+4) then
      iend = ibegin+4
      !DIR$ LOOP COUNT MAX(4)
      do while (cache_key > X(idx))
        idx = idx+1
      end do
    else
      !DIR$ LOOP COUNT MAX(4)
      do while (cache_key > X(idx))
        idx = idx+1
        if (idx == iend) then
          exit
        endif
      end do
    endif
    if (cache_key /= X(idx)) then
      idx = 1-idx
    endif
    return

  else

    if (cache_key < X(ibegin)) then
      idx = -ibegin
      return
    endif
    if (cache_key > X(iend)) then
      idx = -iend
      return
    endif
    if (cache_key == X(ibegin)) then
      idx = ibegin
      return
    endif
    if (cache_key == X(iend)) then
      idx = iend
      return
    endif
  endif

end

subroutine search_key_value_big_interval(key,value,X,Y,sze,idx,ibegin_in,iend_in)
  use map_module
  implicit none
  integer(cache_map_size_kind), intent(in) :: sze
  integer(key_kind)           , intent(in) :: key
  real(integral_kind)         , intent(out) :: value
  integer(cache_key_kind)     , intent(in) :: X(sze)
  real(integral_kind)         , intent(in) :: Y(sze)
  integer(cache_map_size_kind), intent(in) :: ibegin_in, iend_in
  integer(cache_map_size_kind), intent(out) :: idx

  integer(cache_map_size_kind)   :: istep, ibegin, iend, i
  integer(cache_key_kind)        :: cache_key

  if (sze /= 0) then
    continue
  else
    idx = -1
    value = 0.d0
    return
  endif
  cache_key = int(iand(key,map_mask),2)
  ibegin = min(ibegin_in,sze)
  iend   = min(iend_in,sze)
  if ((cache_key > X(ibegin)) .and. (cache_key < X(iend))) then

    istep = shiftr(iend+ibegin,1)
    idx = ibegin + istep
    do while (istep > 4)
      idx = ibegin + istep
      ! TODO : Cache misses
      if (cache_key < X(idx)) then
        iend = idx
        istep = shiftr(idx-ibegin,1)
        idx = ibegin + istep
        if (cache_key < X(idx)) then
          iend = idx
          istep = shiftr(idx-ibegin,1)
          cycle
        else if (cache_key > X(idx)) then
          ibegin = idx
          istep = shiftr(iend-idx,1)
          cycle
        else
          value = Y(idx)
          return
        endif
      else if (cache_key > X(idx)) then
        ibegin = idx
        istep = shiftr(iend-idx,1)
        idx = idx + istep
        if (cache_key < X(idx)) then
          iend = idx
          istep = shiftr(idx-ibegin,1)
          cycle
        else if (cache_key > X(idx)) then
          ibegin = idx
          istep = shiftr(iend-idx,1)
          cycle
        else
          value = Y(idx)
          return
        endif
      else
        value = Y(idx)
        return
      endif
    enddo
    idx = ibegin
    if (min(iend_in,sze) > ibegin+4) then
      iend = ibegin+4
      !DIR$ LOOP COUNT MAX(4)
      do while (cache_key > X(idx))
        idx = idx+1
      end do
    else
      !DIR$ LOOP COUNT MAX(4)
      do while (cache_key > X(idx))
        idx = idx+1
        if (idx == iend) then
          exit
        endif
      end do
    endif
    if (cache_key /= X(idx)) then
      idx = 1-idx
      value = 0.d0
    else
      value = Y(idx)
    endif
    return

  else

    if (cache_key < X(ibegin)) then
      idx = -ibegin
      value = 0.d0
      return
    endif
    if (cache_key > X(iend)) then
      idx = -iend
      value = 0.d0
      return
    endif
    if (cache_key == X(ibegin)) then
      idx = ibegin
      value = Y(idx)
      return
    endif
    if (cache_key == X(iend)) then
      idx = iend
      value = Y(idx)
      return
    endif
  endif

end


subroutine get_cache_map_n_elements_max(map,n_elements_max)
  use map_module
  implicit none
  ! Returns the size of the largest cache_map
  type (map_type), intent(in)    :: map
  integer(cache_map_size_kind), intent(out) :: n_elements_max
  integer(map_size_kind)         :: i
  n_elements_max = 0_cache_map_size_kind
  do i=0_8,map%map_size
    n_elements_max = max(n_elements_max, map%map(i)%n_elements)
  enddo
end



subroutine get_cache_map(map,map_idx,keys,values,n_elements)
  use map_module
  implicit none
  type (map_type), intent(in)    :: map
  integer(map_size_kind), intent(in) :: map_idx
  integer(cache_map_size_kind), intent(inout) :: n_elements
  integer(key_kind), intent(out) :: keys(n_elements)
  double precision, intent(out)  :: values(n_elements)
  integer(cache_map_size_kind)   :: i
  integer(key_kind)              :: shift

  shift = shiftl(map_idx,map_shift)

  n_elements = map%map(map_idx)%n_elements
  do i=1,n_elements
    keys(i)   = map%map(map_idx)%key(i)   + shift
    values(i) = map%map(map_idx)%value(i)
  enddo

end

