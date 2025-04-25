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
    out_val(1:sze) = 0.d0
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


