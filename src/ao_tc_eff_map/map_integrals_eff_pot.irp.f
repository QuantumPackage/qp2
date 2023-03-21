use map_module

!! AO Map
!! ======

BEGIN_PROVIDER [ type(map_type), ao_tc_sym_two_e_pot_map ]
  implicit none
  BEGIN_DOC
  ! |AO| integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call two_e_integrals_index(ao_num,ao_num,ao_num,ao_num,key_max)
  sze = key_max
  call map_init(ao_tc_sym_two_e_pot_map,sze)
  print*,  'ao_tc_sym_two_e_pot_map map initialized : ', sze
END_PROVIDER

 BEGIN_PROVIDER [ integer, ao_tc_sym_two_e_pot_cache_min ]
&BEGIN_PROVIDER [ integer, ao_tc_sym_two_e_pot_cache_max ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the AOs for which the integrals are in the cache
 END_DOC
 ao_tc_sym_two_e_pot_cache_min = max(1,ao_num - 63)
 ao_tc_sym_two_e_pot_cache_max = ao_num

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ao_tc_sym_two_e_pot_cache, (0:64*64*64*64) ]

  use map_module
  implicit none

  BEGIN_DOC
  ! Cache of |AO| integrals for fast access
  END_DOC

  integer                        :: i,j,k,l,ii
  integer(key_kind)              :: idx
  real(integral_kind)            :: integral

  PROVIDE ao_tc_sym_two_e_pot_in_map

 !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx,ii,integral)
  do l = ao_tc_sym_two_e_pot_cache_min, ao_tc_sym_two_e_pot_cache_max
    do k = ao_tc_sym_two_e_pot_cache_min, ao_tc_sym_two_e_pot_cache_max
      do j = ao_tc_sym_two_e_pot_cache_min, ao_tc_sym_two_e_pot_cache_max
        do i = ao_tc_sym_two_e_pot_cache_min, ao_tc_sym_two_e_pot_cache_max
          !DIR$ FORCEINLINE
          call two_e_integrals_index(i, j, k, l, idx)
          !DIR$ FORCEINLINE
          call map_get(ao_tc_sym_two_e_pot_map, idx, integral)
          ii = l-ao_tc_sym_two_e_pot_cache_min
          ii = ior( ishft(ii,6), k-ao_tc_sym_two_e_pot_cache_min)
          ii = ior( ishft(ii,6), j-ao_tc_sym_two_e_pot_cache_min)
          ii = ior( ishft(ii,6), i-ao_tc_sym_two_e_pot_cache_min)
          ao_tc_sym_two_e_pot_cache(ii) = integral
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER

! ---

subroutine insert_into_ao_tc_sym_two_e_pot_map(n_integrals, buffer_i, buffer_values)

  use map_module
  implicit none

  BEGIN_DOC
  ! Create new entry into |AO| map
  END_DOC

  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)

  call map_append(ao_tc_sym_two_e_pot_map, buffer_i, buffer_values, n_integrals)

end

! ---

double precision function get_ao_tc_sym_two_e_pot(i, j, k, l, map) result(result)

  use map_module

  implicit none

  BEGIN_DOC
  ! Gets one |AO| two-electron integral from the |AO| map
  END_DOC

  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  real(integral_kind)            :: tmp
  logical, external              :: ao_two_e_integral_zero

  PROVIDE ao_tc_sym_two_e_pot_in_map ao_tc_sym_two_e_pot_cache ao_tc_sym_two_e_pot_cache_min

  !DIR$ FORCEINLINE
!  if (ao_two_e_integral_zero(i,j,k,l)) then
  if (.False.) then
    tmp = 0.d0
  !else if (ao_two_e_integral_erf_schwartz(i,k)*ao_two_e_integral_erf_schwartz(j,l) < ao_integrals_threshold) then
  !  tmp = 0.d0
  else
    ii = l-ao_tc_sym_two_e_pot_cache_min
    ii = ior(ii, k-ao_tc_sym_two_e_pot_cache_min)
    ii = ior(ii, j-ao_tc_sym_two_e_pot_cache_min)
    ii = ior(ii, i-ao_tc_sym_two_e_pot_cache_min)
    if (iand(ii, -64) /= 0) then
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i, j, k, l, idx)
      !DIR$ FORCEINLINE
      call map_get(map, idx, tmp)
      tmp = tmp
    else
      ii = l-ao_tc_sym_two_e_pot_cache_min
      ii = ior( ishft(ii,6), k-ao_tc_sym_two_e_pot_cache_min)
      ii = ior( ishft(ii,6), j-ao_tc_sym_two_e_pot_cache_min)
      ii = ior( ishft(ii,6), i-ao_tc_sym_two_e_pot_cache_min)
      tmp = ao_tc_sym_two_e_pot_cache(ii)
    endif
  endif

  result = tmp

end

! ---

subroutine get_many_ao_tc_sym_two_e_pot(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple |AO| two-electron integral from the |AO| map .
  ! All i are retrieved for j,k,l fixed.
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)

  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh
! logical, external              :: ao_one_e_integral_zero
  PROVIDE ao_tc_sym_two_e_pot_in_map ao_tc_sym_two_e_pot_map
  thresh = ao_integrals_threshold

! if (ao_one_e_integral_zero(j,l)) then
  if (.False.) then
    out_val = 0.d0
    return
  endif

  double precision :: get_ao_tc_sym_two_e_pot
  do i=1,sze
    out_val(i) = get_ao_tc_sym_two_e_pot(i,j,k,l,ao_tc_sym_two_e_pot_map)
  enddo

end

subroutine get_many_ao_tc_sym_two_e_pot_non_zero(j,k,l,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple |AO| two-electron integrals from the |AO| map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  integer, intent(out)           :: out_val_index(sze),non_zero_int

  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh,tmp
! logical, external              :: ao_one_e_integral_zero
  PROVIDE ao_tc_sym_two_e_pot_in_map
  thresh = ao_integrals_threshold

  non_zero_int = 0
! if (ao_one_e_integral_zero(j,l)) then
  if (.False.) then
    out_val = 0.d0
    return
  endif

  non_zero_int = 0
  do i=1,sze
    integer, external :: ao_l4
    double precision, external :: ao_two_e_integral_eff_pot
    !DIR$ FORCEINLINE
    !if (ao_two_e_integral_erf_schwartz(i,k)*ao_two_e_integral_erf_schwartz(j,l) < thresh) then
    !  cycle
    !endif
    call two_e_integrals_index(i,j,k,l,hash)
    call map_get(ao_tc_sym_two_e_pot_map, hash,tmp)
    if (dabs(tmp) < thresh ) cycle
    non_zero_int = non_zero_int+1
    out_val_index(non_zero_int) = i
    out_val(non_zero_int) = tmp
  enddo

end


function get_ao_tc_sym_two_e_pot_map_size()
  implicit none
  integer (map_size_kind) :: get_ao_tc_sym_two_e_pot_map_size
  BEGIN_DOC
  ! Returns the number of elements in the |AO| map
  END_DOC
  get_ao_tc_sym_two_e_pot_map_size = ao_tc_sym_two_e_pot_map % n_elements
end

subroutine clear_ao_tc_sym_two_e_pot_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the |AO| map
  END_DOC
  call map_deinit(ao_tc_sym_two_e_pot_map)
  FREE ao_tc_sym_two_e_pot_map
end



subroutine dump_ao_tc_sym_two_e_pot(filename)
  use map_module
  implicit none
  BEGIN_DOC
  ! Save to disk the |AO| eff_pot integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer*8                      :: i,j, n
  call ezfio_set_work_empty(.False.)
  open(unit=66,file=filename,FORM='unformatted')
  write(66) integral_kind, key_kind
  write(66) ao_tc_sym_two_e_pot_map%sorted, ao_tc_sym_two_e_pot_map%map_size,    &
      ao_tc_sym_two_e_pot_map%n_elements
  do i=0_8,ao_tc_sym_two_e_pot_map%map_size
    write(66) ao_tc_sym_two_e_pot_map%map(i)%sorted, ao_tc_sym_two_e_pot_map%map(i)%map_size,&
        ao_tc_sym_two_e_pot_map%map(i)%n_elements
  enddo
  do i=0_8,ao_tc_sym_two_e_pot_map%map_size
    key => ao_tc_sym_two_e_pot_map%map(i)%key
    val => ao_tc_sym_two_e_pot_map%map(i)%value
    n = ao_tc_sym_two_e_pot_map%map(i)%n_elements
    write(66) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  close(66)

end



integer function load_ao_tc_sym_two_e_pot(filename)
  implicit none
  BEGIN_DOC
  ! Read from disk the |AO| eff_pot integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer*8                      :: i
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer                        :: iknd, kknd
  integer*8                      :: n, j
  load_ao_tc_sym_two_e_pot = 1
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
  read(66,err=98,end=98) ao_tc_sym_two_e_pot_map%sorted, ao_tc_sym_two_e_pot_map%map_size,&
      ao_tc_sym_two_e_pot_map%n_elements
  do i=0_8, ao_tc_sym_two_e_pot_map%map_size
    read(66,err=99,end=99) ao_tc_sym_two_e_pot_map%map(i)%sorted,          &
        ao_tc_sym_two_e_pot_map%map(i)%map_size, ao_tc_sym_two_e_pot_map%map(i)%n_elements
    call cache_map_reallocate(ao_tc_sym_two_e_pot_map%map(i),ao_tc_sym_two_e_pot_map%map(i)%map_size)
  enddo
  do i=0_8, ao_tc_sym_two_e_pot_map%map_size
    key => ao_tc_sym_two_e_pot_map%map(i)%key
    val => ao_tc_sym_two_e_pot_map%map(i)%value
    n = ao_tc_sym_two_e_pot_map%map(i)%n_elements
    read(66,err=99,end=99) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  call map_sort(ao_tc_sym_two_e_pot_map)
  load_ao_tc_sym_two_e_pot = 0
  return
  99 continue
  call map_deinit(ao_tc_sym_two_e_pot_map)
  98 continue
  stop 'Problem reading ao_tc_sym_two_e_pot_map file in work/'

end




