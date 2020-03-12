use map_module

subroutine insert_into_mo_integrals_map_2(n_integrals,               &
      buffer_i, buffer_values, thr)
  use map_module
  implicit none

  BEGIN_DOC
  ! Create new entry into MO map, or accumulate in an existing entry
  END_DOC

  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  real(integral_kind), intent(in)    :: thr
  call map_update(mo_integrals_map_2, buffer_i, buffer_values, n_integrals, thr)
end

BEGIN_PROVIDER [ complex*16, mo_integrals_cache_complex, (0_8:128_8*128_8*128_8*128_8) ]
 implicit none
 BEGIN_DOC
 ! Cache of MO integrals for fast access
 END_DOC
 PROVIDE mo_two_e_integrals_in_map
 integer*8                      :: i,j,k,l
 integer*4                      :: i4,j4,k4,l4
 integer*8                      :: ii
 integer(key_kind)              :: idx
 complex(integral_kind)            :: integral
 complex*16 :: get_two_e_integral_complex_simple
 FREE ao_integrals_cache
 !$OMP PARALLEL DO PRIVATE (i,j,k,l,i4,j4,k4,l4,idx,ii,integral)
 do l=mo_integrals_cache_min_8,mo_integrals_cache_max_8
   l4 = int(l,4)
   do k=mo_integrals_cache_min_8,mo_integrals_cache_max_8
     k4 = int(k,4)
     do j=mo_integrals_cache_min_8,mo_integrals_cache_max_8
       j4 = int(j,4)
       do i=mo_integrals_cache_min_8,mo_integrals_cache_max_8
         i4 = int(i,4)
         !DIR$ FORCEINLINE
         integral = get_two_e_integral_complex_simple(i4,j4,k4,l4,&
                    mo_integrals_map,mo_integrals_map_2)
         ii = l-mo_integrals_cache_min_8
         ii = ior( shiftl(ii,7), k-mo_integrals_cache_min_8)
         ii = ior( shiftl(ii,7), j-mo_integrals_cache_min_8)
         ii = ior( shiftl(ii,7), i-mo_integrals_cache_min_8)
         mo_integrals_cache_complex(ii) = integral
       enddo
     enddo
   enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER


complex*16 function get_two_e_integral_complex_simple(i,j,k,l,map,map2) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one MO bi-electronic integral from the MO map
  ! reuse ao map/idx/sign function
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  real(integral_kind)            :: tmp_re, tmp_im
  type(map_type), intent(inout)  :: map,map2
  complex(integral_kind)         :: tmp
  logical                        :: use_map1
  double precision               :: sign
  PROVIDE mo_two_e_integrals_in_map
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

complex*16 function get_two_e_integral_complex(i,j,k,l,map,map2)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  ! TODO: finish this
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  integer                        :: ii
  integer*8                      :: ii_8
  type(map_type), intent(inout)  :: map,map2
  complex(integral_kind)         :: tmp
  complex(integral_kind)         :: get_two_e_integral_complex_simple
  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache_complex
  ii = l-mo_integrals_cache_min
  ii = ior(ii, k-mo_integrals_cache_min)
  ii = ior(ii, j-mo_integrals_cache_min)
  ii = ior(ii, i-mo_integrals_cache_min)
  if (iand(ii, -128) /= 0) then
    tmp = get_two_e_integral_complex_simple(i,j,k,l,map,map2)
  else
    ii_8 = int(l,8)-mo_integrals_cache_min_8
    ii_8 = ior( shiftl(ii_8,7), int(k,8)-mo_integrals_cache_min_8)
    ii_8 = ior( shiftl(ii_8,7), int(j,8)-mo_integrals_cache_min_8)
    ii_8 = ior( shiftl(ii_8,7), int(i,8)-mo_integrals_cache_min_8)
    tmp = mo_integrals_cache_complex(ii_8)
  endif
  get_two_e_integral_complex = tmp
end

complex*16 function mo_two_e_integral_complex(i,j,k,l)
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  complex*16                     :: get_two_e_integral_complex
  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache_complex
  PROVIDE mo_two_e_integrals_in_map
  !DIR$ FORCEINLINE
  mo_two_e_integral_complex = get_two_e_integral_complex(i,j,k,l,mo_integrals_map,mo_integrals_map_2)
  return
end

subroutine get_mo_two_e_integrals_complex(j,k,l,sze,out_val,map,map2)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  complex*16, intent(out)        :: out_val(sze)
  type(map_type), intent(inout)  :: map,map2
  integer                        :: i
  complex*16, external           :: get_two_e_integral_complex_simple

  integer                        :: ii, ii0
  integer*8                      :: ii_8, ii0_8
  complex(integral_kind)         :: tmp
  integer(key_kind)              :: i1, idx
  integer(key_kind)              :: p,q,r,s,i2
  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache_complex

!DEBUG
!  do i=1,sze
!    out_val(i) = get_two_e_integral_complex(i,j,k,l,map,map2)
!  enddo
!  return
!DEBUG

  ii0 = l-mo_integrals_cache_min
  ii0 = ior(ii0, k-mo_integrals_cache_min)
  ii0 = ior(ii0, j-mo_integrals_cache_min)

  ii0_8 = int(l,8)-mo_integrals_cache_min_8
  ii0_8 = ior( shiftl(ii0_8,7), int(k,8)-mo_integrals_cache_min_8)
  ii0_8 = ior( shiftl(ii0_8,7), int(j,8)-mo_integrals_cache_min_8)

  do i=1,sze
    ii = ior(ii0, i-mo_integrals_cache_min)
    if (iand(ii, -128) == 0) then
      ii_8 = ior( shiftl(ii0_8,7), int(i,8)-mo_integrals_cache_min_8)
      out_val(i) = mo_integrals_cache_complex(ii_8)
    else
      out_val(i) = get_two_e_integral_complex_simple(i,j,k,l,map,map2)
    endif
  enddo
end

!subroutine get_mo_two_e_integrals_ij_complex(k,l,sze,out_array,map)
!  use map_module
!  implicit none
!  BEGIN_DOC
!  ! Returns multiple integrals <ij|kl> in the MO basis, all
!  ! i(1)j(2) 1/r12 k(1)l(2)
!  ! i, j for k,l fixed.
!  END_DOC
!  integer, intent(in)            :: k,l, sze
!  double precision, intent(out)  :: out_array(sze,sze)
!  type(map_type), intent(inout)  :: map
!  integer                        :: i,j,kk,ll,m
!  integer(key_kind),allocatable  :: hash(:)
!  integer  ,allocatable          :: pairs(:,:), iorder(:)
!  real(integral_kind), allocatable :: tmp_val(:)
!
!  PROVIDE mo_two_e_integrals_in_map
!  allocate (hash(sze*sze), pairs(2,sze*sze),iorder(sze*sze), &
!  tmp_val(sze*sze))
!
!  kk=0
!  out_array = 0.d0
!  do j=1,sze
!   do i=1,sze
!    kk += 1
!    !DIR$ FORCEINLINE
!    call two_e_integrals_index(i,j,k,l,hash(kk))
!    pairs(1,kk) = i
!    pairs(2,kk) = j
!    iorder(kk) = kk
!   enddo
!  enddo
!
!  logical :: integral_is_in_map
!  if (key_kind == 8) then
!    call i8radix_sort(hash,iorder,kk,-1)
!  else if (key_kind == 4) then
!    call iradix_sort(hash,iorder,kk,-1)
!  else if (key_kind == 2) then
!    call i2radix_sort(hash,iorder,kk,-1)
!  endif
!
!  call map_get_many(mo_integrals_map, hash, tmp_val, kk)
!
!  do ll=1,kk
!    m = iorder(ll)
!    i=pairs(1,m)
!    j=pairs(2,m)
!    out_array(i,j) = tmp_val(ll)
!  enddo
!
!  deallocate(pairs,hash,iorder,tmp_val)
!end

!subroutine get_mo_two_e_integrals_i1j1_complex(k,l,sze,out_array,map)
!  use map_module
!  implicit none
!  BEGIN_DOC
!  ! Returns multiple integrals <ik|jl> in the MO basis, all
!  ! i(1)j(1) 1/r12 k(2)l(2)
!  ! i, j for k,l fixed.
!  END_DOC
!  integer, intent(in)            :: k,l, sze
!  double precision, intent(out)  :: out_array(sze,sze)
!  type(map_type), intent(inout)  :: map
!  integer                        :: i,j,kk,ll,m
!  integer(key_kind),allocatable  :: hash(:)
!  integer  ,allocatable          :: pairs(:,:), iorder(:)
!  real(integral_kind), allocatable :: tmp_val(:)
!
!  PROVIDE mo_two_e_integrals_in_map
!  allocate (hash(sze*sze), pairs(2,sze*sze),iorder(sze*sze), &
!  tmp_val(sze*sze))
!
!  kk=0
!  out_array = 0.d0
!  do j=1,sze
!   do i=1,sze
!    kk += 1
!    !DIR$ FORCEINLINE
!    call two_e_integrals_index(i,k,j,l,hash(kk))
!    pairs(1,kk) = i
!    pairs(2,kk) = j
!    iorder(kk) = kk
!   enddo
!  enddo
!
!  logical :: integral_is_in_map
!  if (key_kind == 8) then
!    call i8radix_sort(hash,iorder,kk,-1)
!  else if (key_kind == 4) then
!    call iradix_sort(hash,iorder,kk,-1)
!  else if (key_kind == 2) then
!    call i2radix_sort(hash,iorder,kk,-1)
!  endif
!
!  call map_get_many(mo_integrals_map, hash, tmp_val, kk)
!
!  do ll=1,kk
!    m = iorder(ll)
!    i=pairs(1,m)
!    j=pairs(2,m)
!    out_array(i,j) = tmp_val(ll)
!  enddo
!
!  deallocate(pairs,hash,iorder,tmp_val)
!end

subroutine get_mo_two_e_integrals_coulomb_ii_complex(k,l,sze,out_val,map,map2)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|li>
  ! k(1)i(2) 1/r12 l(1)i(2) :: out_val(i1)
  ! for k,l fixed.
  ! real and in map2 if k==l
  ! complex and in map1 otherwise
  ! take conjugate if k>l
  ! TODO: determine best way to structure code 
  !       to account for single/double integral_kind, real/complex, and +/- imag part
  END_DOC
  integer, intent(in)            :: k,l, sze
  complex*16, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map,map2
  integer                        :: i
  integer(key_kind)              :: hash(sze),hash_re(sze),hash_im(sze)
  real(integral_kind)            :: tmp_re(sze),tmp_im(sze)
  double precision               :: out_re(sze),out_im(sze)
  double precision :: sign
  PROVIDE mo_two_e_integrals_in_map

  if (k.eq.l) then ! real, call other function
    call get_mo_two_e_integrals_coulomb_ijij_complex(k,sze,out_re,map2)
    do i=1,sze
      out_val(i) = dcmplx(out_re(i),0.d0)
    enddo
  else ! complex
    if (k.gt.l) then
      sign = -1.d0
    else
      sign = 1.d0
    endif

    do i=1,sze
      !DIR$ FORCEINLINE
      call two_e_integrals_index(k,i,l,i,hash(i))
      !hash_im(i) = hash(i)*2
      hash_im(i) = shiftl(hash(i),1)
      hash_re(i) = hash_im(i)-1
    enddo

    if (integral_kind == 8) then
      call map_get_many(map, hash_re, out_re, sze)
      call map_get_many(map, hash_im, out_im, sze)
      do i=1,sze
        out_val(i) = dcmplx(out_re(i),sign*out_im(i))
      enddo
    else
      call map_get_many(map, hash_re, tmp_re, sze)
      call map_get_many(map, hash_im, tmp_im, sze)
      ! Conversion to double complex
      do i=1,sze
        out_val(i) = dcmplx(tmp_re(i),sign*tmp_im(i))
      enddo
    endif
  endif
end

subroutine get_mo_two_e_integrals_coulomb_ijij_complex(j,sze,out_val,map2)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|ij>
  ! i*(1)j*(2) 1/r12 i(1)j(2) :: out_val(i)
  ! for j fixed.
  ! always in map2, always real
  END_DOC
  integer, intent(in)            :: j, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map2
  integer                        :: i
  integer(key_kind)              :: hash(sze),hash_re(sze)
  real(integral_kind)            :: tmp_re(sze)
  PROVIDE mo_two_e_integrals_in_map

  do i=1,sze
    !DIR$ FORCEINLINE
    call two_e_integrals_index(i,j,i,j,hash(i))
    !hash_re(i) = hash(i)*2 - 1
    hash_re(i) = shiftl(hash(i),1) - 1
  enddo

  if (integral_kind == 8) then
    call map_get_many(map2, hash_re, out_val, sze)
  else
    call map_get_many(map2, hash_re, tmp_re, sze)
    ! Conversion to double complex
    do i=1,sze
      out_val(i) = dble(tmp_re(i))
    enddo
  endif
end

subroutine get_mo_two_e_integrals_exch_ii_complex(k,l,sze,out_val,map,map2)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|il>
  ! k*(1)i*(2) 1/r12 i(1)l(2) :: out_val(i1)
  ! for k,l fixed.
  ! 
  ! if k<l, then:
  ! i < k        map2 +
  ! k <= i <= l  map1 +
  ! l < i        map2 -
  !
  ! if l<k, then same maps as above, but take complex conjugate
  END_DOC
  integer, intent(in)            :: k,l, sze
  complex*16, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map,map2
  integer                        :: i,klmin,klmax
  integer(key_kind)              :: hash(sze),hash_re(sze),hash_im(sze)
  real(integral_kind)            :: tmp_re(sze),tmp_im(sze)
  double precision               :: out_re(sze),out_im(sze)
  double precision :: sign,sign2(sze)

  PROVIDE mo_two_e_integrals_in_map

  
  if (k.eq.l) then ! real, call other function
    call get_mo_two_e_integrals_exch_ijji_complex(k,sze,out_re,map,map2)
    do i=1,sze
      out_val(i) = dcmplx(out_re(i),0.d0)
    enddo
  else ! complex
    if (k.gt.l) then
      sign = -1.d0
    else
      sign = 1.d0
    endif
    klmin = min(k,l) ! put these in conditional above?
    klmax = max(k,l)
    sign2(1:klmax) = 1.d0
    sign2(klmax+1:sze) = -1.d0
    do i=1,sze
      !DIR$ FORCEINLINE
      call two_e_integrals_index(k,i,i,l,hash(i))
      !hash_im(i) = hash(i)*2
      hash_im(i) = shiftl(hash(i),1)
      hash_re(i) = hash_im(i)-1
    enddo

    if (integral_kind == 8) then
      call map_get_many(map2, hash_re(1:klmin-1),   out_re(1:klmin-1),   klmin-1)
      call map_get_many(map2, hash_im(1:klmin-1),   out_im(1:klmin-1),   klmin-1)
      call map_get_many(map,  hash_re(klmin:klmax), out_re(klmin:klmax), klmax-klmin+1)
      call map_get_many(map,  hash_im(klmin:klmax), out_im(klmin:klmax), klmax-klmin+1)
      if (klmax.lt.size) then
        call map_get_many(map2, hash_re(klmax+1:sze), out_re(klmax+1:sze), sze-klmax)
        call map_get_many(map2, hash_im(klmax+1:sze), out_im(klmax+1:sze), sze-klmax)
      endif
      do i=1,sze
        out_val(i) = dcmplx(out_re(i),sign*sign2(i)*out_im(i))
      enddo
    else
      call map_get_many(map2, hash_re(1:klmin-1),   tmp_re(1:klmin-1),   klmin-1)
      call map_get_many(map2, hash_im(1:klmin-1),   tmp_im(1:klmin-1),   klmin-1)
      call map_get_many(map,  hash_re(klmin:klmax), tmp_re(klmin:klmax), klmax-klmin+1)
      call map_get_many(map,  hash_im(klmin:klmax), tmp_im(klmin:klmax), klmax-klmin+1)
      if (klmax.lt.size) then
        call map_get_many(map2, hash_re(klmax+1:sze), tmp_re(klmax+1:sze), sze-klmax)
        call map_get_many(map2, hash_im(klmax+1:sze), tmp_im(klmax+1:sze), sze-klmax)
      endif
      ! Conversion to double complex
      do i=1,sze
        out_val(i) = dcmplx(tmp_re(i),sign*sign2(i)*tmp_im(i))
      enddo
    endif
  endif
end

subroutine get_mo_two_e_integrals_exch_ijji_complex(j,sze,out_val,map,map2)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|ji>
  ! i*(1)j*(2) 1/r12 j(1)i(2) :: out_val(i)
  ! for j fixed.
  ! always real, always in map2 
  END_DOC
  integer, intent(in)            :: j, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map,map2
  integer                        :: i
  integer(key_kind)              :: hash(sze),hash_re(sze)
  real(integral_kind)            :: tmp_val(sze)
  PROVIDE mo_two_e_integrals_in_map

  do i=1,sze
    !DIR$ FORCEINLINE
    call two_e_integrals_index(i,j,j,i,hash(i))
    !hash_re(i) = 2*hash(i) - 1
    hash_re(i) = shiftl(hash(i),1) - 1
  enddo

  if (integral_kind == 8) then
    call map_get_many(map2, hash_re, out_val, sze)
  else
    call map_get_many(map2, hash_re, tmp_val, sze)
    ! Conversion to double precision
    do i=1,sze
      out_val(i) = dble(tmp_val(i))
    enddo
  endif
end

