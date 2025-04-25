
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

