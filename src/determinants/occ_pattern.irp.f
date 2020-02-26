use bitmasks
subroutine occ_pattern_of_det(d,o,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Transforms a determinant to an occupation pattern
  !
  ! occ(:,1) : Single occupations
  !
  ! occ(:,2) : Double occupations
  !
  END_DOC
  integer          ,intent(in)   :: Nint
  integer(bit_kind),intent(in)   :: d(Nint,2)
  integer(bit_kind),intent(out)  :: o(Nint,2)

  integer                        :: k

  do k=1,Nint
    o(k,1) = ieor(d(k,1),d(k,2))
    o(k,2) = iand(d(k,1),d(k,2))
  enddo
end


subroutine occ_pattern_to_dets_size(o,sze,n_alpha,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
!  Number of possible determinants for a given occ_pattern
  END_DOC
  integer          ,intent(in)   :: Nint, n_alpha
  integer(bit_kind),intent(in)   :: o(Nint,2)
  integer, intent(out)           :: sze
  integer                        :: amax,bmax,k
  double precision, external     :: binom_func

  bmax = 0
  amax = n_alpha
  do k=1,Nint
    bmax += popcnt( o(k,1) )
    amax -= popcnt( o(k,2) )
  enddo
  if (binom_int(bmax, amax) > huge(1)) then
    print *,  irp_here, ': Too many determinants to generate'
    stop 1
  endif
  sze = int(binom_int(bmax, amax),4)
end


subroutine occ_pattern_to_dets(o,d,sze,n_alpha,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Generate all possible determinants for a given occ_pattern
  ! 
  ! Input :
  !    o   : occupation pattern : (doubly occupied, singly occupied)
  !    sze : Number of produced determinants, computed by `occ_pattern_to_dets_size`
  !    n_alpha : Number of $\alpha$ electrons
  !    Nint    : N_int
  !
  ! Output:
  !    d : determinants 
  !
  END_DOC
  integer          ,intent(in)   :: Nint
  integer          ,intent(in)   :: n_alpha        ! Number of alpha electrons
  integer         ,intent(inout) :: sze            ! Dimension of the output dets
  integer(bit_kind),intent(in)   :: o(Nint,2)      ! Occ patters
  integer(bit_kind),intent(out)  :: d(Nint,2,sze)  ! Output determinants

  integer                        :: i, k, n, ispin, ispin2

  ! Extract list of singly occupied MOs as (int,pos) pairs
  ! ------------------------------------------------------

  integer           :: iint(2*n_alpha), ipos(2*n_alpha)
  integer(bit_kind) :: v, t, tt, diff, v_prev
  integer           :: n_alpha_in_single

  n=0
  n_alpha_in_single = n_alpha
  do i=1,Nint
    v = o(i,1)
    do while(v /= 0_bit_kind)
      n = n+1
      iint(n) = i
      ipos(n) = trailz(v)
      v = iand(v,v-1)
    enddo
    n_alpha_in_single = n_alpha_in_single - popcnt( o(i,2) )
  enddo

  v = shiftl(1,n_alpha_in_single) - 1

  ! Initialize first determinant
  d(:,1,1) = o(:,2)
  d(:,2,1) = o(:,2)

  do k=1,n_alpha_in_single
    d(iint(k),1,1) = ibset( d(iint(k),1,1), ipos(k) )
  enddo

  do k=n_alpha_in_single+1,n
    d(iint(k),2,1) = ibset( d(iint(k),2,1), ipos(k) )
  enddo

  sze = int(binom_int(n,n_alpha_in_single),4)

  if ( (shiftl(n_alpha_in_single,1) == n).and.n>0 ) then

    ! Time reversal symmetry
    d(:,1,2) = d(:,2,1)
    d(:,2,2) = d(:,1,1)

    do i=3,sze,2
      ! Generate next permutation with Anderson's algorithm
      v_prev = v
      t = ior(v,v-1)
      tt = t+1
      v = ior(tt, shiftr( and(not(t),tt) - 1, trailz(v)+1) )

      ! Find what has changed between v_prev and v
      diff = ieor(v,v_prev)

      ! Initialize with previous determinant
      d(:,1,i) = d(:,1,i-2)
      d(:,2,i) = d(:,2,i-2)

      ! Swap bits only where they have changed from v_prev to v
      do while (diff /= 0_bit_kind)
        k = trailz(diff)+1
        if (btest(v,k-1)) then
          d(iint(k),1,i) = ibset( d(iint(k),1,i), ipos(k) )
          d(iint(k),2,i) = ibclr( d(iint(k),2,i), ipos(k) )
        else
          d(iint(k),1,i) = ibclr( d(iint(k),1,i), ipos(k) )
          d(iint(k),2,i) = ibset( d(iint(k),2,i), ipos(k) )
        endif
        diff = iand(diff,diff-1_bit_kind)
      enddo

      ! Time reversal symmetry
      d(:,1,i+1) = d(:,2,i)
      d(:,2,i+1) = d(:,1,i)

    enddo

  else

    do i=2,sze
      ! Generate next permutation with Anderson's algorithm
      v_prev = v
      t = ior(v,v-1)
      tt = t+1
      v = ior(tt, shiftr( and(not(t),tt) - 1, trailz(v)+1) )

      ! Find what has changed between v_prev and v
      diff = ieor(v,v_prev)

      ! Initialize with previous determinant
      d(:,1,i) = d(:,1,i-1)
      d(:,2,i) = d(:,2,i-1)

      ! Swap bits only where they have changed from v_prev to v
      do while (diff /= 0_bit_kind)
        k = trailz(diff)+1
        if (btest(v,k-1)) then
          d(iint(k),1,i) = ibset( d(iint(k),1,i), ipos(k) )
          d(iint(k),2,i) = ibclr( d(iint(k),2,i), ipos(k) )
        else
          d(iint(k),1,i) = ibclr( d(iint(k),1,i), ipos(k) )
          d(iint(k),2,i) = ibset( d(iint(k),2,i), ipos(k) )
        endif
        diff = iand(diff,diff-1_bit_kind)
      enddo

    enddo

  endif

end


 BEGIN_PROVIDER [ integer(bit_kind), psi_occ_pattern, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_occ_pattern ]
 implicit none
 BEGIN_DOC
  ! Array of the occ_patterns present in the wave function.
  !
  ! psi_occ_pattern(:,1,j) = j-th occ_pattern of the wave function : represents all the single occupations
  !
  ! psi_occ_pattern(:,2,j) = j-th occ_pattern of the wave function : represents all the double occupations
  !
  ! The occ patterns are sorted by :c:func:`occ_pattern_search_key`
 END_DOC
 integer :: i,j,k

 ! create
 do i = 1, N_det
  do k = 1, N_int
   psi_occ_pattern(k,1,i) = ieor(psi_det(k,1,i),psi_det(k,2,i))
   psi_occ_pattern(k,2,i) = iand(psi_det(k,1,i),psi_det(k,2,i))
  enddo
 enddo

 ! Sort
 integer, allocatable           :: iorder(:)
 integer*8, allocatable         :: bit_tmp(:)
 integer*8, external            :: occ_pattern_search_key
 integer(bit_kind), allocatable :: tmp_array(:,:,:)
 logical,allocatable            :: duplicate(:)
 logical :: dup


 allocate ( iorder(N_det), duplicate(N_det), bit_tmp(N_det), tmp_array(N_int,2,N_det) )

 do i=1,N_det
   iorder(i) = i
   bit_tmp(i) = occ_pattern_search_key(psi_occ_pattern(1,1,i),N_int)
 enddo

 call i8sort(bit_tmp,iorder,N_det)


 !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j,k,dup)

 !$OMP DO
 do i=1,N_det
  do k=1,N_int
    tmp_array(k,1,i) = psi_occ_pattern(k,1,iorder(i))
    tmp_array(k,2,i) = psi_occ_pattern(k,2,iorder(i))
  enddo
  duplicate(i) = .False.
 enddo
 !$OMP END DO

 ! Find duplicates
 !$OMP DO
 do i=1,N_det-1
  if (duplicate(i)) then
    cycle
  endif
  j = i+1
  do while (bit_tmp(j)==bit_tmp(i))
    if (duplicate(j)) then
      j+=1
      if (j>N_det) then
        exit
      endif
      cycle
    endif
    dup = .True.
    do k=1,N_int
      if ( (tmp_array(k,1,i) /= tmp_array(k,1,j)) &
      .or. (tmp_array(k,2,i) /= tmp_array(k,2,j)) ) then
         dup = .False.
         exit
      endif
    enddo
    if (dup) then
      duplicate(j) = .True.
    endif
    j+=1
    if (j>N_det) then
      exit
    endif
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 ! Copy filtered result
 N_occ_pattern=0
 do i=1,N_det
  if (duplicate(i)) then
    cycle
  endif
  N_occ_pattern += 1
  do k=1,N_int
    psi_occ_pattern(k,1,N_occ_pattern) = tmp_array(k,1,i)
    psi_occ_pattern(k,2,N_occ_pattern) = tmp_array(k,2,i)
  enddo
 enddo

!- Check
!  print *,  'Checking for duplicates in occ pattern'
!  do i=1,N_occ_pattern
!   do j=i+1,N_occ_pattern
!     duplicate(1) = .True.
!     do k=1,N_int
!       if (psi_occ_pattern(k,1,i) /= psi_occ_pattern(k,1,j)) then
!         duplicate(1) = .False.
!         exit
!       endif
!       if (psi_occ_pattern(k,2,i) /= psi_occ_pattern(k,2,j)) then
!         duplicate(1) = .False.
!         exit
!       endif
!     enddo
!     if (duplicate(1)) then
!       call debug_det(psi_occ_pattern(1,1,i),N_int)
!       call debug_det(psi_occ_pattern(1,1,j),N_int)
!       stop 'DUPLICATE'
!     endif
!   enddo
!  enddo
!  print *,  'No duplicates'
!-
 deallocate(iorder,duplicate,bit_tmp,tmp_array)

END_PROVIDER

BEGIN_PROVIDER [ integer, det_to_occ_pattern, (N_det) ]
 implicit none
 BEGIN_DOC
 ! Returns the index of the occupation pattern for each determinant
 END_DOC
 integer :: i,j,k,r,l
 integer*8 :: key
 integer(bit_kind) :: occ(N_int,2)
 logical :: found
 integer*8, allocatable :: bit_tmp(:)
 integer*8, external            :: occ_pattern_search_key

 allocate(bit_tmp(N_occ_pattern))
 do i=1,N_occ_pattern
   bit_tmp(i) = occ_pattern_search_key(psi_occ_pattern(1,1,i),N_int)
 enddo

 !$OMP PARALLEL DO DEFAULT(SHARED) &
 !$OMP PRIVATE(i,k,j,r,l,key,found,occ)
 do i=1,N_det
    do k = 1, N_int
      occ(k,1) = ieor(psi_det(k,1,i),psi_det(k,2,i))
      occ(k,2) = iand(psi_det(k,1,i),psi_det(k,2,i))
    enddo

    key = occ_pattern_search_key(occ,N_int)

    ! TODO: Binary search
    l = 1
    r = N_occ_pattern
!    do while(r-l > 32)
!      j = shiftr(r+l,1)
!      if (bit_tmp(j) < key) then
!        l = j
!      else
!        r = j
!      endif
!    enddo
    do j=l,r
      found = .True.
      do k=1,N_int
        if ( (occ(k,1) /= psi_occ_pattern(k,1,j)) &
        .or. (occ(k,2) /= psi_occ_pattern(k,2,j)) ) then
          found = .False.
          exit
        endif
      enddo
      if (found) then
        det_to_occ_pattern(i) = j
        exit
      endif
    enddo

    if (.not.found) then
      print *,  '3 bug in ',  irp_here
      stop -1
    endif
 enddo
 !$OMP END PARALLEL DO
 deallocate(bit_tmp)
END_PROVIDER


BEGIN_PROVIDER [ double precision, psi_occ_pattern_Hii, (N_occ_pattern) ]
 implicit none
 BEGIN_DOC
 ! $\langle I|H|I \rangle$ where $|I\rangle$ is an occupation pattern.
 ! This is the minimum $H_{ii}$, where the $|i\rangle$ are the
 ! determinants of $|I\rangle$.
 END_DOC
 integer :: j, i

 psi_occ_pattern_Hii(:) = huge(1.d0)
 do i=1,N_det
  j = det_to_occ_pattern(i)
  psi_occ_pattern_Hii(j) = min(psi_occ_pattern_Hii(j), psi_det_Hii(i))
 enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, weight_occ_pattern, (N_occ_pattern,N_states) ]
 implicit none
 BEGIN_DOC
 ! Weight of the occupation patterns in the wave function
 END_DOC
 integer :: i,j,k
 weight_occ_pattern = 0.d0
 if (is_complex) then
 do i=1,N_det
  j = det_to_occ_pattern(i)
  do k=1,N_states
    weight_occ_pattern(j,k) += cdabs(psi_coef_complex(i,k) * psi_coef_complex(i,k))
  enddo
 enddo
 else
 do i=1,N_det
  j = det_to_occ_pattern(i)
  do k=1,N_states
    weight_occ_pattern(j,k) += psi_coef(i,k) * psi_coef(i,k)
  enddo
 enddo
 endif
END_PROVIDER

BEGIN_PROVIDER [ double precision, weight_occ_pattern_average, (N_occ_pattern) ]
 implicit none
 BEGIN_DOC
 ! State-average weight of the occupation patterns in the wave function
 END_DOC
 integer :: i,j,k
 weight_occ_pattern_average(:) = 0.d0
 if (is_complex) then
  do i=1,N_det
   j = det_to_occ_pattern(i)
   do k=1,N_states
     weight_occ_pattern_average(j) += cdabs(psi_coef_complex(i,k) * psi_coef_complex(i,k)) * state_average_weight(k)
   enddo
  enddo
 else
 do i=1,N_det
  j = det_to_occ_pattern(i)
  do k=1,N_states
    weight_occ_pattern_average(j) += psi_coef(i,k) * psi_coef(i,k) * state_average_weight(k)
  enddo
 enddo
 endif
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_occ_pattern_sorted, (N_int,2,N_occ_pattern) ]
&BEGIN_PROVIDER [ double precision, weight_occ_pattern_average_sorted, (N_occ_pattern) ]
&BEGIN_PROVIDER [ integer, psi_occ_pattern_sorted_order, (N_occ_pattern) ]
&BEGIN_PROVIDER [ integer, psi_occ_pattern_sorted_order_reverse, (N_occ_pattern) ]
 implicit none
 BEGIN_DOC
 ! Occupation patterns sorted by weight 
 END_DOC
 integer                        :: i,j,k
 integer, allocatable           :: iorder(:)
 allocate ( iorder(N_occ_pattern) )
 do i=1,N_occ_pattern
   weight_occ_pattern_average_sorted(i) = -weight_occ_pattern_average(i)
   iorder(i) = i
 enddo
 call dsort(weight_occ_pattern_average_sorted,iorder,N_occ_pattern)
 do i=1,N_occ_pattern
   do j=1,N_int
     psi_occ_pattern_sorted(j,1,i) = psi_occ_pattern(j,1,iorder(i)) 
     psi_occ_pattern_sorted(j,2,i) = psi_occ_pattern(j,2,iorder(i)) 
   enddo
   psi_occ_pattern_sorted_order(iorder(i)) = i
   psi_occ_pattern_sorted_order_reverse(i) = iorder(i)
   weight_occ_pattern_average_sorted(i) = -weight_occ_pattern_average_sorted(i)
 enddo

 deallocate(iorder)

END_PROVIDER


subroutine make_s2_eigenfunction
  implicit none
  integer                        :: i,j,k
  integer                        :: smax, s
  integer(bit_kind), allocatable :: d(:,:,:), det_buffer(:,:,:)
  integer                        :: N_det_new, ithread, omp_get_thread_num
  integer, parameter             :: bufsze = 1000
  logical, external              :: is_in_wavefunction
  logical                        :: update

  update=.False.
  call write_int(6,N_occ_pattern,'Number of occupation patterns')

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP  SHARED(N_occ_pattern, psi_occ_pattern, elec_alpha_num,N_int,update) &
  !$OMP  PRIVATE(s,ithread, d, det_buffer, smax, N_det_new,i,j,k)
  N_det_new = 0
  call occ_pattern_to_dets_size(psi_occ_pattern(1,1,1),s,elec_alpha_num,N_int)
  allocate (d(N_int,2,s+64), det_buffer(N_int,2,bufsze) )
  smax = s
  ithread=0
  !$ ithread = omp_get_thread_num()
  !$OMP DO SCHEDULE (dynamic,1000)
  do i=1,N_occ_pattern
    call occ_pattern_to_dets_size(psi_occ_pattern(1,1,i),s,elec_alpha_num,N_int)
    s += 1
    if (s > smax) then
      deallocate(d)
      allocate ( d(N_int,2,s+64) )
      smax = s
    endif
    call occ_pattern_to_dets(psi_occ_pattern(1,1,i),d,s,elec_alpha_num,N_int)
    do j=1,s
      if ( is_in_wavefunction(d(1,1,j), N_int) ) then
        cycle
      endif
      update = .true.
      N_det_new += 1
      det_buffer(:,:,N_det_new) = d(:,:,j)
      if (N_det_new == bufsze) then
        call fill_h_apply_buffer_no_selection(bufsze,det_buffer,N_int,ithread)
        N_det_new = 0
      endif
    enddo
  enddo
  !$OMP END DO NOWAIT

  if (N_det_new > 0) then
    call fill_H_apply_buffer_no_selection(N_det_new,det_buffer,N_int,ithread)
  endif
  !$OMP BARRIER
  deallocate(d,det_buffer)
  !$OMP END PARALLEL

  if (update) then
    call copy_h_apply_buffer_to_wf
    if (is_complex) then
      TOUCH N_det psi_coef_complex psi_det psi_occ_pattern N_occ_pattern
    else
    TOUCH N_det psi_coef psi_det psi_occ_pattern N_occ_pattern
    endif
  endif
  call write_time(6)

end



