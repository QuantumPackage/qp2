use map_module

BEGIN_PROVIDER [ logical, all_mo_integrals ]
  implicit none
  BEGIN_DOC
! Used to provide everything needed before using MO integrals
! PROVIDE all_mo_integrals
  END_DOC
  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache mo_two_e_integrals_jj_exchange mo_two_e_integrals_jj_anti mo_two_e_integrals_jj big_array_exchange_integrals big_array_coulomb_integrals mo_one_e_integrals
END_PROVIDER

BEGIN_PROVIDER [ logical, mo_cholesky_double ]
  implicit none
  BEGIN_DOC
! If true, use double precision to compute integrals from cholesky vectors
  END_DOC
  mo_cholesky_double = .True.
END_PROVIDER


!! MO Map
!! ======

BEGIN_PROVIDER [ type(map_type), mo_integrals_map ]
  implicit none
  BEGIN_DOC
  ! MO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call two_e_integrals_index(mo_num,mo_num,mo_num,mo_num,key_max)
  sze = key_max
  call map_init(mo_integrals_map,sze)
  print*, 'MO map initialized: ', sze
END_PROVIDER

subroutine insert_into_mo_integrals_map(n_integrals,                 &
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
  call map_update(mo_integrals_map, buffer_i, buffer_values, n_integrals, thr)
end

 BEGIN_PROVIDER [ integer, mo_integrals_cache_min ]
&BEGIN_PROVIDER [ integer, mo_integrals_cache_max ]
&BEGIN_PROVIDER [ integer, mo_integrals_cache_size ]
&BEGIN_PROVIDER [ integer*8, mo_integrals_cache_size_8 ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the MOs for which the integrals are in the cache
 END_DOC

 mo_integrals_cache_size  = shiftl(1,mo_integrals_cache_shift)
 mo_integrals_cache_size_8  = shiftl(1_8, mo_integrals_cache_shift*4)


 mo_integrals_cache_min = max(1,elec_alpha_num - (mo_integrals_cache_size/2 - 1) )
 mo_integrals_cache_max = min(mo_num, mo_integrals_cache_min + mo_integrals_cache_size - 1)
 print *, 'MO integrals cache: (', mo_integrals_cache_min, ', ', mo_integrals_cache_max, '), ', &
          shiftr(mo_integrals_cache_size_8, 17), 'MiB'

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_integrals_cache, (0_8:mo_integrals_cache_size_8) ]
 implicit none
 BEGIN_DOC
 ! Cache of MO integrals for fast access
 END_DOC
 PROVIDE mo_two_e_integrals_in_map
 integer                        :: i,j,k,l
 integer*8                      :: ii
 integer(key_kind)              :: idx
 real(integral_kind)            :: integral
 FREE ao_integrals_cache

 if (do_mo_cholesky) then

   call set_multiple_levels_omp(.False.)

   !$OMP PARALLEL DO PRIVATE(k,l,ii) SCHEDULE(dynamic)
   do l=mo_integrals_cache_min,mo_integrals_cache_max
     do k=mo_integrals_cache_min,mo_integrals_cache_max
         ii = int(l-mo_integrals_cache_min,8)
         ii = ior( shiftl(ii,mo_integrals_cache_shift), int(k-mo_integrals_cache_min,8))
         ii = shiftl(ii,2*mo_integrals_cache_shift)
         call dgemm('T','N', mo_integrals_cache_max-mo_integrals_cache_min+1, &
                             mo_integrals_cache_max-mo_integrals_cache_min+1, &
           cholesky_mo_num, 1.d0, &
           cholesky_mo_transp(1,mo_integrals_cache_min,k), cholesky_mo_num, &
           cholesky_mo_transp(1,mo_integrals_cache_min,l), cholesky_mo_num, 0.d0, &
           mo_integrals_cache(ii), mo_integrals_cache_size)
     enddo
   enddo
   !$OMP END PARALLEL DO

 else
   !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx,ii,integral) SCHEDULE(dynamic)
   do l=mo_integrals_cache_min,mo_integrals_cache_max
     do k=mo_integrals_cache_min,mo_integrals_cache_max
       do j=mo_integrals_cache_min,mo_integrals_cache_max
         do i=mo_integrals_cache_min,mo_integrals_cache_max
           !DIR$ FORCEINLINE
           call two_e_integrals_index(i,j,k,l,idx)
           !DIR$ FORCEINLINE
           call map_get(mo_integrals_map,idx,integral)
           ii = int(l-mo_integrals_cache_min,8)
           ii = ior( shiftl(ii,mo_integrals_cache_shift), int(k-mo_integrals_cache_min,8))
           ii = ior( shiftl(ii,mo_integrals_cache_shift), int(j-mo_integrals_cache_min,8))
           ii = ior( shiftl(ii,mo_integrals_cache_shift), int(i-mo_integrals_cache_min,8))
           mo_integrals_cache(ii) = integral
         enddo
       enddo
     enddo
   enddo
   !$OMP END PARALLEL DO
 endif

END_PROVIDER


double precision function get_two_e_integral_cache(i,j,k,l)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis taken from the cache
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer*8                      :: ii

  ii = int(l-mo_integrals_cache_min,8)
  ii = ior( shiftl(ii,mo_integrals_cache_shift), int(k-mo_integrals_cache_min,8))
  ii = ior( shiftl(ii,mo_integrals_cache_shift), int(j-mo_integrals_cache_min,8))
  ii = ior( shiftl(ii,mo_integrals_cache_shift), int(i-mo_integrals_cache_min,8))
  get_two_e_integral_cache = mo_integrals_cache(ii)

end


double precision function get_two_e_integral(i,j,k,l,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  integer                        :: ii, kk
  type(map_type), intent(inout)  :: map
  real(integral_kind)            :: tmp

  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache do_mo_cholesky

  if (use_banned_excitation) then
    if (banned_excitation(i,k)) then
      get_two_e_integral = 0.d0
      return
    endif
    if (banned_excitation(j,l)) then
      get_two_e_integral = 0.d0
      return
    endif
  endif


  ii = l-mo_integrals_cache_min
  ii = ior(ii, k-mo_integrals_cache_min)
  ii = ior(ii, j-mo_integrals_cache_min)
  ii = ior(ii, i-mo_integrals_cache_min)

  if (iand(ii, -mo_integrals_cache_size) == 0) then

    double precision, external :: get_two_e_integral_cache
    get_two_e_integral = get_two_e_integral_cache(i,j,k,l)

  else

    ! Integral is not in the cache

    if  (do_mo_cholesky) then

      double precision, external :: ddot
      real, external :: sdot
      integer :: isplit
      if (mo_cholesky_double) then
        get_two_e_integral = ddot(cholesky_mo_num, cholesky_mo_transp(1,i,k), 1, cholesky_mo_transp(1,j,l), 1)
      else
        get_two_e_integral = 0.d0
        do isplit=1,4
          get_two_e_integral = get_two_e_integral + &
                               sdot(cholesky_mo_num_split(isplit+1) - cholesky_mo_num_split(isplit), &
                                    cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),i,k), 1, &
                                    cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),j,l), 1)
        enddo
      endif

    else

      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,idx)
      !DIR$ FORCEINLINE
      call map_get(map,idx,tmp)
      get_two_e_integral = dble(tmp)
    endif

  endif
end


subroutine get_mo_two_e_integrals(j,k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i

  integer                        :: ii
  real(integral_kind)            :: tmp
  integer(key_kind)              :: i1, idx
  integer(key_kind)              :: p,q,r,s,i2
  real, allocatable :: out_val_sp(:)
  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache cholesky_mo_transp cholesky_mo_transp_sp

  if (banned_excitation(j,l)) then
      out_val(1:sze) = 0.d0
      return
  endif
!
  ii = l-mo_integrals_cache_min
  ii = ior(ii, k-mo_integrals_cache_min)
  ii = ior(ii, j-mo_integrals_cache_min)

  if (do_mo_cholesky.and. .not.mo_cholesky_double) then
    allocate(out_val_sp(sze))
  endif

  if (iand(ii, -mo_integrals_cache_size) == 0) then
    ! Some integrals are in the cache

    if (mo_integrals_cache_min > 1) then

      if (do_mo_cholesky) then

        !TODO: bottleneck here
        if (mo_cholesky_double) then
          call dgemv('T', cholesky_mo_num, mo_integrals_cache_min-1, 1.d0, &
            cholesky_mo_transp(1,1,k), cholesky_mo_num, &
            cholesky_mo_transp(1,j,l), 1, 0.d0, &
            out_val, 1)
        else
          integer :: isplit
          call sgemv('T', cholesky_mo_num_split(isplit+1) - cholesky_mo_num_split(isplit), &
              mo_integrals_cache_min-1, 1., &
              cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),1,k), cholesky_mo_num, &
              cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),j,l), 1, 0., &
              out_val_sp, 1)
            out_val(1:mo_integrals_cache_min-1) += out_val_sp(1:mo_integrals_cache_min-1)
          out_val(1:mo_integrals_cache_min-1) = out_val_sp(1:mo_integrals_cache_min-1)
          do isplit=2,4
            call sgemv('T', cholesky_mo_num_split(isplit+1) - cholesky_mo_num_split(isplit), &
              mo_integrals_cache_min-1, 1., &
              cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),1,k), cholesky_mo_num, &
              cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),j,l), 1, 0., &
              out_val_sp, 1)
            out_val(1:mo_integrals_cache_min-1) += out_val_sp(1:mo_integrals_cache_min-1)
          enddo
        endif

      else

        q = min(j,l)
        s = max(j,l)
        q = q+shiftr(s*s-s,1)

        do i=1,mo_integrals_cache_min-1
          if (banned_excitation(i,k)) then
            out_val(i) = 0.d0
            cycle
          endif
          p = min(i,k)
          r = max(i,k)
          p = p+shiftr(r*r-r,1)
          i1 = min(p,q)
          i2 = max(p,q)
          idx = i1+shiftr(i2*i2-i2,1)
          !DIR$ FORCEINLINE
          call map_get(map,idx,tmp)
          out_val(i) = dble(tmp)
        enddo

      endif

    endif

    call get_mo_two_e_integrals_cache(j,k,l,sze,out_val)

    if (mo_integrals_cache_max < mo_num) then

      if (do_mo_cholesky) then

        !TODO: bottleneck here
        if (mo_cholesky_double) then
          call dgemv('T', cholesky_mo_num, mo_num-mo_integrals_cache_max, 1.d0, &
             cholesky_mo_transp(1,mo_integrals_cache_max+1,k), cholesky_mo_num, &
             cholesky_mo_transp(1,j,l), 1, 0.d0, &
             out_val(mo_integrals_cache_max+1), 1)
        else
          out_val = 0.d0
          do isplit=1,4
            call sgemv('T', cholesky_mo_num_split(isplit+1) - cholesky_mo_num_split(isplit), &
              mo_num-mo_integrals_cache_max, 1., &
              cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),mo_integrals_cache_max+1,k), cholesky_mo_num, &
              cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),j,l), 1, 0., &
              out_val_sp(mo_integrals_cache_max+1), 1)
            out_val(mo_integrals_cache_max+1:sze) += out_val_sp(mo_integrals_cache_max+1:sze)
          enddo
        endif

      else

        q = min(j,l)
        s = max(j,l)
        q = q+shiftr(s*s-s,1)

        do i=mo_integrals_cache_max+1,mo_num
          if (banned_excitation(i,k)) then
            out_val(i) = 0.d0
            cycle
          endif
          p = min(i,k)
          r = max(i,k)
          p = p+shiftr(r*r-r,1)
          i1 = min(p,q)
          i2 = max(p,q)
          idx = i1+shiftr(i2*i2-i2,1)
          !DIR$ FORCEINLINE
          call map_get(map,idx,tmp)
          out_val(i) = dble(tmp)
        enddo

      endif

    endif

  else

    if (do_mo_cholesky) then

      !TODO: bottleneck here
      if (mo_cholesky_double) then
          call dgemv('T', cholesky_mo_num, sze, 1.d0, &
               cholesky_mo_transp(1,1,k), cholesky_mo_num, &
               cholesky_mo_transp(1,j,l), 1, 0.d0, &
               out_val, 1)
      else
          out_val = 0.d0
          do isplit=1,4
            call sgemv('T', cholesky_mo_num_split(isplit+1) - cholesky_mo_num_split(isplit), &
              sze, 1., &
              cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),1,k), cholesky_mo_num, &
              cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),j,l), 1, 0., &
              out_val_sp, 1)
            out_val(1:sze) += out_val_sp(1:sze)
          enddo
      endif

    else

      q = min(j,l)
      s = max(j,l)
      q = q+shiftr(s*s-s,1)

      do i=1,sze
        if (banned_excitation(i,k)) cycle
        p = min(i,k)
        r = max(i,k)
        p = p+shiftr(r*r-r,1)
        i1 = min(p,q)
        i2 = max(p,q)
        idx = i1+shiftr(i2*i2-i2,1)
        !DIR$ FORCEINLINE
        call map_get(map,idx,tmp)
        out_val(i) = dble(tmp)
      enddo

    endif

  endif

end

double precision function mo_two_e_integral(i,j,k,l)
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  double precision               :: get_two_e_integral
  PROVIDE all_mo_integrals
  !DIR$ FORCEINLINE
  mo_two_e_integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
  return
end


subroutine get_mo_two_e_integrals_cache(j,k,l,sze,out_val)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i for j,k,l fixed, all integrals from the cache
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  double precision, intent(out)  :: out_val(sze)
  integer*8                      :: ii

  ii = int(l-mo_integrals_cache_min,8)
  ii = ior( shiftl(ii, mo_integrals_cache_shift), int(k-mo_integrals_cache_min,8))
  ii = ior( shiftl(ii, mo_integrals_cache_shift), int(j-mo_integrals_cache_min,8))
  ii = shiftl(ii, mo_integrals_cache_shift)
  out_val(mo_integrals_cache_min:mo_integrals_cache_max) = &
       mo_integrals_cache(ii:ii+int(mo_integrals_cache_max-mo_integrals_cache_min,8))

end

subroutine get_mo_two_e_integrals_ij(k,l,sze,out_array,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i(1)j(2) 1/r12 k(1)l(2)
  ! i, j for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_array(sze,sze)
  type(map_type), intent(inout)  :: map
  integer                        :: j

  if ( (mo_integrals_cache_min>1).or.(mo_integrals_cache_max<mo_num) ) then

    if (do_mo_cholesky) then

      if (mo_cholesky_double) then
          call dgemm('T', 'N', mo_num, mo_num, cholesky_mo_num, 1.d0, &
             cholesky_mo_transp(1,1,k), cholesky_mo_num, &
             cholesky_mo_transp(1,1,l), cholesky_mo_num, 0.d0, &
             out_array, sze)
      else
          integer :: isplit
          double precision, allocatable :: out_array_sp(:,:)
          allocate(out_array_sp(sze,sze))
          call sgemm('T', 'N', mo_num, mo_num, &
             cholesky_mo_num_split(isplit+1) - cholesky_mo_num_split(isplit), 1.d0, &
             cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),1,k), cholesky_mo_num, &
             cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),1,l), cholesky_mo_num, 0.d0, &
             out_array_sp, sze)
          out_array(1:sze,1:sze) = out_array_sp(1:sze,1:sze)
          do isplit=2,4
            call sgemm('T', 'N', mo_num, mo_num, &
             cholesky_mo_num_split(isplit+1) - cholesky_mo_num_split(isplit), 1.d0, &
             cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),1,k), cholesky_mo_num, &
             cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),1,l), cholesky_mo_num, 0.d0, &
             out_array_sp, sze)
            out_array(1:sze,1:sze) = out_array(1:sze,1:sze) + out_array_sp(1:sze,1:sze)
         enddo
         deallocate(out_array_sp)
      endif

    else

      do j=1,sze
        call get_mo_two_e_integrals(j,k,l,sze,out_array(1,j),map)
      enddo

    endif

  else

    double precision, external :: get_two_e_integral_cache
    do j=1,sze
      call get_mo_two_e_integrals_cache(j,k,l,sze,out_array(1,j))
    enddo

  endif

end

subroutine get_mo_two_e_integrals_i1j1(k,l,sze,out_array,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ik|jl> in the MO basis, all
  ! i(1)j(1) 1/r12 k(2)l(2)
  ! i, j for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_array(sze,sze)
  type(map_type), intent(inout)  :: map
  integer                        :: j
  PROVIDE mo_two_e_integrals_in_map

  if ( (mo_integrals_cache_min>1).or.(mo_integrals_cache_max<mo_num) ) then

    if (do_mo_cholesky) then

      call dgemv('T', cholesky_mo_num, mo_num*mo_num, 1.d0, &
         cholesky_mo_transp(1,1,1), cholesky_mo_num, &
         cholesky_mo_transp(1,k,l), 1, 0.d0, &
         out_array, 1)

    else

      do j=1,sze
        call get_mo_two_e_integrals(k,j,l,sze,out_array(1,j),map)
      enddo

    endif

  else

    double precision, external :: get_two_e_integral_cache
    do j=1,sze
      call get_mo_two_e_integrals_cache(k,j,l,sze,out_array(1,j))
    enddo

  endif

end


subroutine get_mo_two_e_integrals_coulomb_ii(k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|li>
  ! k(1)i(2) 1/r12 l(1)i(2) :: out_val(i1)
  ! for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  double precision, external     :: get_two_e_integral
  PROVIDE mo_two_e_integrals_in_map

  if ( (mo_integrals_cache_min>1).or.(mo_integrals_cache_max<mo_num) ) then

    if (do_mo_cholesky) then

      call dgemv('T', cholesky_mo_num, mo_num, 1.d0, &
         cholesky_mo_transp(1,1,1), cholesky_mo_num*(mo_num+1), &
         cholesky_mo_transp(1,k,l), 1, 0.d0, &
         out_val, 1)

    else

      do i=1,sze
        out_val(i) = get_two_e_integral(i,k,i,l,map)
      enddo

    endif

  else

    double precision, external :: get_two_e_integral_cache
    do i=1,sze
      out_val(i) = get_two_e_integral_cache(i,k,i,l)
    enddo

  endif


end

subroutine get_mo_two_e_integrals_exch_ii(k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|il>
  ! k(1)i(2) 1/r12 i(1)l(2) :: out_val(i1)
  ! for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  double precision, external     :: get_two_e_integral
  PROVIDE mo_two_e_integrals_in_map mo_cholesky_double

  if ( (mo_integrals_cache_min>1).or.(mo_integrals_cache_max<mo_num) ) then

    if (do_mo_cholesky) then

      if ( (k>=mo_integrals_cache_min).and.(k<=mo_integrals_cache_max).and. &
           (l>=mo_integrals_cache_min).and.(l<=mo_integrals_cache_max) ) then

        double precision, external :: ddot
        real, external :: sdot
        integer :: kk

        if (mo_cholesky_double) then

          do i=1,mo_integrals_cache_min-1
            out_val(i) = ddot(cholesky_mo_num, cholesky_mo_transp(1,i,k), 1, &
                                               cholesky_mo_transp(1,i,l), 1)
          enddo

          do i=mo_integrals_cache_min,mo_integrals_cache_max
            out_val(i) = get_two_e_integral_cache(i,i,k,l)
          enddo

          do i=mo_integrals_cache_max, sze
            out_val(i) = ddot(cholesky_mo_num, cholesky_mo_transp(1,i,k), 1, &
                                               cholesky_mo_transp(1,i,l), 1)
          enddo

        else

          integer :: isplit
          do i=1,mo_integrals_cache_min-1
             out_val(i) = 0.d0
             do isplit=1,4
               out_val(i) += sdot(cholesky_mo_num_split(isplit+1) - cholesky_mo_num_split(isplit), &
                                          cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),i,k), 1, &
                                          cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),i,l), 1)
             enddo
          enddo

          do i=mo_integrals_cache_min,mo_integrals_cache_max
            out_val(i) = get_two_e_integral_cache(i,i,k,l)
          enddo

          do i=mo_integrals_cache_max, sze
             out_val(i) = 0.d0
             do isplit=1,4
               out_val(i) += sdot(cholesky_mo_num_split(isplit+1) - cholesky_mo_num_split(isplit), &
                                          cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),i,k), 1, &
                                          cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),i,l), 1)
             enddo
          enddo

        endif

      else

        if (mo_cholesky_double) then
          do i=1,sze
            out_val(i) = ddot(cholesky_mo_num, cholesky_mo_transp(1,i,k), 1, &
                                               cholesky_mo_transp(1,i,l), 1)
          enddo
        else
          do i=1,sze
             out_val(i) = 0.d0
             do isplit=1,4
               out_val(i) += sdot(cholesky_mo_num_split(isplit+1) - cholesky_mo_num_split(isplit), &
                                          cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),i,k), 1, &
                                          cholesky_mo_transp_sp(cholesky_mo_num_split(isplit),i,l), 1)
             enddo
          enddo
        endif

      endif

    else

      do i=1,sze
        out_val(i) = get_two_e_integral(i,i,k,l,map)
      enddo

    endif

  else

    double precision, external :: get_two_e_integral_cache
    do i=1,sze
      out_val(i) = get_two_e_integral_cache(i,i,k,l)
    enddo

  endif


end

 BEGIN_PROVIDER [ logical, banned_excitation, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ logical, use_banned_excitation  ]
 implicit none
 use map_module
 BEGIN_DOC
 ! If true, the excitation is banned in the selection. Useful with local MOs.
 END_DOC
 banned_excitation = .False.
 use_banned_excitation = .False.

 integer :: i,j, icount
 integer(key_kind)              :: idx
 double precision :: tmp

!icount = 1 ! Avoid division by zero
!do j=1,mo_num
!  do i=1,j-1
!   call two_e_integrals_index(i,j,j,i,idx)
!   !DIR$ FORCEINLINE
!   call map_get(mo_integrals_map,idx,tmp)
!   banned_excitation(i,j) = dabs(tmp) < 1.d-14
!   banned_excitation(j,i) = banned_excitation(i,j)
!   if (banned_excitation(i,j)) icount = icount+2
! enddo
!enddo
!use_banned_excitation =  (mo_num*mo_num) / icount <= 100  !1%
!if (use_banned_excitation) then
!  print *, 'Using sparsity of exchange integrals'
!endif

END_PROVIDER



integer*8 function get_mo_map_size()
  implicit none
  BEGIN_DOC
  ! Return the number of elements in the MO map
  END_DOC
  get_mo_map_size = mo_integrals_map % n_elements
end


subroutine dump_mo_integrals(filename)
  use map_module
  implicit none
  BEGIN_DOC
  ! Save to disk the |MO| integrals
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
  write(66) mo_integrals_map%sorted, mo_integrals_map%map_size,    &
      mo_integrals_map%n_elements
  do i=0_8,mo_integrals_map%map_size
    write(66) mo_integrals_map%map(i)%sorted, mo_integrals_map%map(i)%map_size,&
        mo_integrals_map%map(i)%n_elements
  enddo
  do i=0_8,mo_integrals_map%map_size
    key => mo_integrals_map%map(i)%key
    val => mo_integrals_map%map(i)%value
    n = mo_integrals_map%map(i)%n_elements
    write(66) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  close(66)

end


integer function load_mo_integrals(filename)
  implicit none
  BEGIN_DOC
  ! Read from disk the |MO| integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer*8                      :: i
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer                        :: iknd, kknd
  integer*8                      :: n, j
  load_mo_integrals = 1
  open(unit=66,file=filename,FORM='unformatted',STATUS='UNKNOWN')
  call lock_io()
  read(66,err=98,end=98) iknd, kknd
  if (iknd /= integral_kind) then
    print *,  'Wrong integrals kind in file :', iknd
    stop 1
  endif
  if (kknd /= key_kind) then
    print *,  'Wrong key kind in file :', kknd
    stop 1
  endif
  read(66,err=98,end=98) mo_integrals_map%sorted, mo_integrals_map%map_size,&
      mo_integrals_map%n_elements
  do i=0_8, mo_integrals_map%map_size
    read(66,err=99,end=99) mo_integrals_map%map(i)%sorted,          &
        mo_integrals_map%map(i)%map_size, mo_integrals_map%map(i)%n_elements
    call cache_map_reallocate(mo_integrals_map%map(i),mo_integrals_map%map(i)%map_size)
  enddo
  do i=0_8, mo_integrals_map%map_size
    key => mo_integrals_map%map(i)%key
    val => mo_integrals_map%map(i)%value
    n = mo_integrals_map%map(i)%n_elements
    read(66,err=99,end=99) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  call unlock_io()
  call map_sort(mo_integrals_map)
  load_mo_integrals = 0
  return
  99 continue
  call map_deinit(mo_integrals_map)
  98 continue
  stop 'Problem reading mo_integrals_map file in work/'

end

