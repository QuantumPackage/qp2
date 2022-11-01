use bitmasks

BEGIN_PROVIDER [ integer, spin_multiplicity ]
 implicit none
 BEGIN_DOC
 ! n_alpha - n_beta + 1
 END_DOC
 spin_multiplicity = elec_alpha_num - elec_beta_num + 1
END_PROVIDER

subroutine configuration_of_det(d,o,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Transforms a determinant to a configuration
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


subroutine configuration_to_dets_size(o,sze,n_alpha,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
!  Number of possible determinants for a given configuration
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
    print *, bmax, amax
    print *,  irp_here, ': Too many determinants to generate'
    stop 1
  endif
  sze = int(binom_int(bmax, amax),4)
end


subroutine configuration_to_dets(o,d,sze,n_alpha,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Generate all possible determinants for a given configuration
  !
  ! Input :
  !    o   : configuration : (doubly occupied, singly occupied)
  !    sze : Number of produced determinants, computed by `configuration_to_dets_size`
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
  integer(bit_kind),intent(in)   :: o(Nint,2)      ! Configurations
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

subroutine configuration_to_dets_tree_addressing(o,d,sze,n_alpha,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Generate all possible determinants for a given configuration
  !
  ! This function preserves the tree addressing i.e.
  ! the time-reversal determinants are at the opposite ends
  ! and not one after the other as in the parent function.
  !
  ! Input :
  !    o   : configuration : (doubly occupied, singly occupied)
  !    sze : Number of produced determinants, computed by `configuration_to_dets_size`
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
  integer(bit_kind),intent(in)   :: o(Nint,2)      ! Configurations
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
    d(:,1,sze) = d(:,2,1)
    d(:,2,sze) = d(:,1,1)

    do i=2,sze/2,1
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

      ! Time reversal symmetry
      d(:,1,sze-i+1) = d(:,2,i)
      d(:,2,sze-i+1) = d(:,1,i)

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


 BEGIN_PROVIDER [ integer(bit_kind), psi_configuration, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_configuration ]
 implicit none
 BEGIN_DOC
  ! Array of the configurations present in the wave function.
  !
  ! psi_configuration(:,1,j) = j-th configuration of the wave function : represents all the single occupations
  !
  ! psi_configuration(:,2,j) = j-th configuration of the wave function : represents all the double occupations
  !
  ! The occ patterns are sorted by :c:func:`configuration_search_key`
 END_DOC
 integer :: i,j,k

 ! create
 do i = 1, N_det
  do k = 1, N_int
   psi_configuration(k,1,i) = ieor(psi_det(k,1,i),psi_det(k,2,i))
   psi_configuration(k,2,i) = iand(psi_det(k,1,i),psi_det(k,2,i))
  enddo
 enddo

 ! Sort
 integer, allocatable           :: iorder(:)
 integer*8, allocatable         :: bit_tmp(:)
 integer*8, external            :: configuration_search_key
 integer(bit_kind), allocatable :: tmp_array(:,:,:)
 logical,allocatable            :: duplicate(:)
 logical :: dup


 allocate ( iorder(N_det), duplicate(N_det), bit_tmp(N_det), tmp_array(N_int,2,N_det) )

 do i=1,N_det
   iorder(i) = i
   bit_tmp(i) = configuration_search_key(psi_configuration(1,1,i),N_int)
 enddo

 call i8sort(bit_tmp,iorder,N_det)


 !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j,k,dup)

 !$OMP DO
 do i=1,N_det
  do k=1,N_int
    tmp_array(k,1,i) = psi_configuration(k,1,iorder(i))
    tmp_array(k,2,i) = psi_configuration(k,2,iorder(i))
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
      dup = dup .and. (tmp_array(k,1,i) == tmp_array(k,1,j)) &
                .and. (tmp_array(k,2,i) == tmp_array(k,2,j))
    enddo
    if (dup) then
      duplicate(j) = .True.
    endif
    j = j+1
    if (j>N_det) then
      exit
    endif
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 ! Copy filtered result
 N_configuration=0
 do i=1,N_det
  if (duplicate(i)) then
    cycle
  endif
  N_configuration += 1
  do k=1,N_int
    psi_configuration(k,1,N_configuration) = tmp_array(k,1,i)
    psi_configuration(k,2,N_configuration) = tmp_array(k,2,i)
  enddo
 enddo

!- Check
!  print *,  'Checking for duplicates in configuration'
!  do i=1,N_configuration
!   do j=i+1,N_configuration
!     duplicate(1) = .True.
!     do k=1,N_int
!       if (psi_configuration(k,1,i) /= psi_configuration(k,1,j)) then
!         duplicate(1) = .False.
!         exit
!       endif
!       if (psi_configuration(k,2,i) /= psi_configuration(k,2,j)) then
!         duplicate(1) = .False.
!         exit
!       endif
!     enddo
!     if (duplicate(1)) then
!       call debug_det(psi_configuration(1,1,i),N_int)
!       call debug_det(psi_configuration(1,1,j),N_int)
!       stop 'DUPLICATE'
!     endif
!   enddo
!  enddo
!  print *,  'No duplicates'
!-
 deallocate(iorder,duplicate,bit_tmp,tmp_array)

END_PROVIDER

 BEGIN_PROVIDER [ integer, cfg_seniority_index, (0:elec_num+2) ]
&BEGIN_PROVIDER [ integer, cfg_nsomo_max ]
&BEGIN_PROVIDER [ integer, cfg_nsomo_min ]
  implicit none
  BEGIN_DOC
 ! Returns the index in psi_configuration of the first cfg with
 ! the requested seniority
 !
 ! cfg_nsomo_max : Max number of SOMO in the current wave function
 END_DOC
 integer :: i, k, s, sold, soldmin
 cfg_seniority_index(:) = -1
 sold = -1
 soldmin = 2000
 cfg_nsomo_max = 0
 do i=1,N_configuration
   s = 0
   do k=1,N_int
     if (psi_configuration(k,1,i) == 0_bit_kind) cycle
     s = s + popcnt(psi_configuration(k,1,i))
   enddo
   if (s /= sold) then
     sold = s
     cfg_seniority_index(s) = i
     cfg_nsomo_max = s
   endif
   if (soldmin .GT. s ) then
     soldmin = s
     cfg_nsomo_min = s
   endif
 enddo
END_PROVIDER

BEGIN_PROVIDER [ integer, det_to_configuration, (N_det) ]
 implicit none
 BEGIN_DOC
 ! Returns the index of the configuration for each determinant
 END_DOC
 integer :: i,j,k,r,l
 integer*8 :: key, key2
 integer(bit_kind) :: occ(N_int,2)
 logical :: found
 integer*8, allocatable :: bit_tmp(:)
 integer*8, external            :: configuration_search_key

 allocate(bit_tmp(0:N_configuration))
 bit_tmp(0) = 0
 do i=1,N_configuration
   bit_tmp(i) = configuration_search_key(psi_configuration(1,1,i),N_int)
 enddo

 !$OMP PARALLEL DO DEFAULT(SHARED) &
 !$OMP PRIVATE(i,k,j,r,l,key,found,occ)
 do i=1,N_det
    do k = 1, N_int
      occ(k,1) = ieor(psi_det(k,1,i),psi_det(k,2,i))
      occ(k,2) = iand(psi_det(k,1,i),psi_det(k,2,i))
    enddo

    key = configuration_search_key(occ,N_int)

    ! Binary search
    l = 0
    r = N_configuration+1
    j = shiftr(r-l,1)
    do while (j>=1)
      j = j+l
      if (bit_tmp(j) == key) then
        do while (bit_tmp(j) == bit_tmp(j-1))
          j = j-1
        enddo
        do while (bit_tmp(j) == key)
          found = .True.
          do k=1,N_int
            found = found .and. (psi_configuration(k,1,j) == occ(k,1)) &
                          .and. (psi_configuration(k,2,j) == occ(k,2))
          enddo
          if (found) then
            det_to_configuration(i) = j
            exit
          endif
          j = j+1
        enddo
        if (found) exit
      else if (bit_tmp(j) > key) then
        r = j
      else
        l = j
      endif
      j = shiftr(r-l,1)
    enddo

 enddo
 !$OMP END PARALLEL DO
 deallocate(bit_tmp)
END_PROVIDER


BEGIN_PROVIDER [ double precision, psi_configuration_Hii, (N_configuration) ]
 implicit none
 BEGIN_DOC
 ! $\langle I|H|I \rangle$ where $|I\rangle$ is a configuration.
 ! This is the minimum $H_{ii}$, where the $|i\rangle$ are the
 ! determinants of $|I\rangle$.
 END_DOC
 integer :: j, i

 psi_configuration_Hii(:) = huge(1.d0)
 do i=1,N_det
  j = det_to_configuration(i)
  psi_configuration_Hii(j) = min(psi_configuration_Hii(j), psi_det_Hii(i))
 enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, weight_configuration, (N_configuration,N_states) ]
 implicit none
 BEGIN_DOC
 ! Weight of the configurations in the wave function
 END_DOC
 integer :: i,j,k
 weight_configuration = 0.d0
 do i=1,N_det
  j = det_to_configuration(i)
  do k=1,N_states
    weight_configuration(j,k) += psi_coef(i,k) * psi_coef(i,k)
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, weight_configuration_average, (N_configuration) ]
 implicit none
 BEGIN_DOC
 ! State-average weight of the configurations in the wave function
 END_DOC
 integer :: i,j,k
 weight_configuration_average(:) = 0.d0
 do i=1,N_det
  j = det_to_configuration(i)
  do k=1,N_states
    weight_configuration_average(j) += psi_coef(i,k) * psi_coef(i,k) * state_average_weight(k)
  enddo
 enddo
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_configuration_sorted, (N_int,2,N_configuration) ]
&BEGIN_PROVIDER [ double precision, weight_configuration_average_sorted, (N_configuration) ]
&BEGIN_PROVIDER [ integer, psi_configuration_sorted_order, (N_configuration) ]
&BEGIN_PROVIDER [ integer, psi_configuration_sorted_order_reverse, (N_configuration) ]
 implicit none
 BEGIN_DOC
 ! Configurations sorted by weight
 END_DOC
 integer                        :: i,j,k
 integer, allocatable           :: iorder(:)
 allocate ( iorder(N_configuration) )
 do i=1,N_configuration
   weight_configuration_average_sorted(i) = -weight_configuration_average(i)
   iorder(i) = i
 enddo
 call dsort(weight_configuration_average_sorted,iorder,N_configuration)
 do i=1,N_configuration
   do j=1,N_int
     psi_configuration_sorted(j,1,i) = psi_configuration(j,1,iorder(i))
     psi_configuration_sorted(j,2,i) = psi_configuration(j,2,iorder(i))
   enddo
   psi_configuration_sorted_order(iorder(i)) = i
   psi_configuration_sorted_order_reverse(i) = iorder(i)
   weight_configuration_average_sorted(i) = -weight_configuration_average_sorted(i)
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
  call write_int(6,N_configuration,'Number of configurations')

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP  SHARED(N_configuration, psi_configuration, elec_alpha_num,N_int,update) &
  !$OMP  PRIVATE(s,ithread, d, det_buffer, smax, N_det_new,i,j,k)
  N_det_new = 0
  call configuration_to_dets_size(psi_configuration(1,1,1),s,elec_alpha_num,N_int)
  allocate (d(N_int,2,s+64), det_buffer(N_int,2,bufsze) )
  smax = s
  ithread=0
  !$ ithread = omp_get_thread_num()
  !$OMP DO SCHEDULE (dynamic,1000)
  do i=1,N_configuration
    call configuration_to_dets_size(psi_configuration(1,1,i),s,elec_alpha_num,N_int)
    s += 1
    if (s > smax) then
      deallocate(d)
      allocate ( d(N_int,2,s+64) )
      smax = s
    endif
    call configuration_to_dets(psi_configuration(1,1,i),d,s,elec_alpha_num,N_int)
    do j=1,s
      if ( is_in_wavefunction(d(1,1,j), N_int) ) then
        cycle
      endif
      update = .true.
      N_det_new += 1
      det_buffer(:,:,N_det_new) = d(:,:,j)
      if (N_det_new == bufsze) then
        call fill_H_apply_buffer_no_selection(bufsze,det_buffer,N_int,ithread)
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
    call copy_H_apply_buffer_to_wf
    TOUCH N_det psi_coef psi_det psi_configuration N_configuration
  endif
  call write_time(6)

end



BEGIN_PROVIDER [ integer, dominant_cfg, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Configuration of the determinants with the largest weight, for each state
 END_DOC
 integer :: k
 dominant_cfg(1) = det_to_configuration(dominant_det(1))
 if (N_det < N_states) then
   dominant_cfg(:) = dominant_cfg(1)
 else
   do k=1,N_states
     dominant_cfg(k) = det_to_configuration(dominant_det(k))
   enddo
 endif
END_PROVIDER


BEGIN_PROVIDER [ integer, N_dominant_dets_of_cfgs ]
 implicit none
 BEGIN_DOC
 ! Number of determinants in all dominant determinants
 END_DOC
 integer                        :: k, sze

 N_dominant_dets_of_cfgs = 0
 do k=1,N_states
   call configuration_to_dets_size(  &
          psi_configuration(1,1,dominant_cfg(k)), &
          sze, elec_alpha_num, N_int)
   N_dominant_dets_of_cfgs += sze
 enddo
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), dominant_dets_of_cfgs, (N_int,2,N_dominant_dets_of_cfgs) ]
 implicit none
 BEGIN_DOC
 ! Configuration of the determinants with the largest weight, for each state
 END_DOC
 integer :: i,k,sze
 i=1
 do k=1,N_states
   sze = N_dominant_dets_of_cfgs
   call configuration_to_dets( &
          psi_configuration(1,1,dominant_cfg(k)), &
          dominant_dets_of_cfgs(1,1,i), &
          sze,elec_alpha_num,N_int)
   i += sze
 enddo
END_PROVIDER

subroutine binary_search_cfg(cfgInp,addcfg,bit_tmp)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Documentation for binary_search
  !
  ! Does a binary search to find
  ! the address of a configuration in a list of
  ! configurations.
  END_DOC
  integer(bit_kind), intent(in)  :: cfgInp(N_int,2)
  integer          , intent(out) :: addcfg
  integer*8,          intent(in) :: bit_tmp(0:N_configuration+1)

  logical                        :: found
  integer                        :: l, r, j, k
  integer*8                      :: key

  integer*8, external            :: configuration_search_key

  key = configuration_search_key(cfgInp,N_int)

  ! Binary search
  l = 0
  r = N_configuration+1
IRP_IF WITHOUT_SHIFTRL
  j = ishft(r-l,-1)
IRP_ELSE
  j = shiftr(r-l,1)
IRP_ENDIF
  do while (j>=1)
    j = j+l
    if (bit_tmp(j) == key) then
      ! Find 1st element which matches the key
      if (j > 1) then
        do while (j>1 .and. bit_tmp(j-1) == key)
          j = j-1
        enddo
      endif
      ! Find correct element matching the key
      do while (bit_tmp(j) == key)
        found = .True.
        do k=1,N_int
          found = found .and. (psi_configuration(k,1,j) == cfgInp(k,1))&
                        .and. (psi_configuration(k,2,j) == cfgInp(k,2))
        enddo
        if (found) then
          addcfg = j
          return
        endif
        j = j+1
      enddo
      addcfg = -1
      return
    else if (bit_tmp(j) > key) then
      r = j
    else
      l = j
    endif
IRP_IF WITHOUT_SHIFTRL
    j = ishft(r-l,-1)
IRP_ELSE
    j = shiftr(r-l,1)
IRP_ENDIF
  enddo

  addcfg = -1
  return

end subroutine

!subroutine binary_search_cfg(cfgInp,addcfg)
!  use bitmasks
!  implicit none
!  BEGIN_DOC
!  ! Documentation for binary_search
!  ! 
!  ! Does a binary search to find 
!  ! the address of a configuration in a list of
!  ! configurations.
!  END_DOC
!  integer(bit_kind), intent(in)  :: cfgInp(N_int,2)
!  integer          , intent(out) :: addcfg
!  integer :: i,j,k,r,l
!  integer*8 :: key, key2
!  logical :: found
!  !integer*8, allocatable :: bit_tmp(:)
!  !integer*8, external            :: configuration_search_key
!
!  !allocate(bit_tmp(0:N_configuration))
!  !bit_tmp(0) = 0
!  do i=1,N_configuration
!    !bit_tmp(i) = configuration_search_key(psi_configuration(1,1,i),N_int)
!    found = .True.
!    do k=1,N_int
!      found = found .and. (psi_configuration(k,1,i) == cfgInp(k,1)) &
!                    .and. (psi_configuration(k,2,i) == cfgInp(k,2))
!    enddo
!    if (found) then
!      addcfg = i
!      exit
!    endif
!  enddo
!
!end subroutine
!
 BEGIN_PROVIDER [ integer, psi_configuration_to_psi_det, (2,N_configuration) ]
&BEGIN_PROVIDER [ integer, psi_configuration_n_det, (N_configuration) ]
&BEGIN_PROVIDER [ integer, psi_configuration_to_psi_det_data, (N_det) ]

 implicit none
 BEGIN_DOC
 ! psi_configuration_to_psi_det_data(k) -> i : i is the index of the
 ! determinant in psi_det.
 !
 ! psi_configuration_to_psi_det(1:2,k) gives the first and last index of the
 ! determinants of configuration k in array psi_configuration_to_psi_det_data.
 END_DOC

 integer :: i, k, iorder
 integer, allocatable :: confs(:)
 allocate (confs(N_det))

 do i=1,N_det
   psi_configuration_to_psi_det_data(i) = i
   confs(i) = det_to_configuration(i)
 enddo

 call isort(confs, psi_configuration_to_psi_det_data, N_det)
 k=1
 psi_configuration_to_psi_det(1,1) = 1
 do i=2,N_det
   if (confs(i) /= confs(i-1)) then
     psi_configuration_to_psi_det(2,k) = i-1
     k = k+1
     psi_configuration_to_psi_det(1,k) = i
   endif
 enddo
 psi_configuration_to_psi_det(2,k) = N_det


 ! Reorder determinants according to generation 
 ! --------------------------------------------

 integer(bit_kind), allocatable :: dets(:,:,:)
 integer                        :: nmax, sze, degree, istart, iend, j
 integer, allocatable           :: old_order(:)


 nmax = 1000
 allocate(dets(N_int,2,nmax), old_order(nmax))

 do k=1,N_configuration
   istart = psi_configuration_to_psi_det(1,k)
   iend   = psi_configuration_to_psi_det(2,k)

   if (iend-istart+1 > nmax) then
      nmax = iend-istart+1
      deallocate(dets)
      allocate(dets(N_int,2,nmax))
   endif

   sze = nmax
   call configuration_to_dets_tree_addressing(                       &
       psi_configuration(1,1,k),                                     &
       dets, sze, elec_alpha_num, N_int)

   if (sze /= iend-istart+1) then
      print *, 'bug in ', irp_here
      stop -1
   endif

   do i=1,sze
     old_order(i) = psi_configuration_to_psi_det_data(i-1+istart)
   enddo

   do i=1,sze
     do j=1,sze

       if (old_order(j) == 0) cycle

       call get_excitation_degree(dets(1,1,i),                       &
           psi_det(1, 1, old_order(j)), degree, N_int)

       if (degree == 0) then
         psi_configuration_to_psi_det_data(i-1+istart) = old_order(j)
         old_order(j) = 0
         exit
       endif

     enddo
   enddo

 enddo

 deallocate(dets, old_order)
 integer :: ndet_conf
 do i = 1, N_configuration
  ndet_conf = psi_configuration_to_psi_det(2,i) - psi_configuration_to_psi_det(1,i) + 1
  psi_configuration_n_det(i) = ndet_conf
 enddo

END_PROVIDER


BEGIN_PROVIDER [ integer, n_elec_alpha_for_psi_configuration, (N_configuration)]
 implicit none
 integer :: i,j,k,l
 integer(bit_kind) :: det_tmp(N_int,2),det_alpha(N_int)
 n_elec_alpha_for_psi_configuration = 0
 do i = 1, N_configuration
  j = psi_configuration_to_psi_det(2,i) 
  det_tmp(:,:) = psi_det(:,:,j)
  k = 0
  do l = 1, N_int
   det_alpha(N_int) = iand(det_tmp(l,1),psi_configuration(l,1,i))
   k += popcnt(det_alpha(l))
  enddo
  n_elec_alpha_for_psi_configuration(i) = k
 enddo

END_PROVIDER 
