
subroutine double_allowed_mo_kpts(h1,h2,p1,p2,is_allowed)
  implicit none
  integer, intent(in) :: h1,h2,p1,p2
  logical, intent(out) :: is_allowed
  integer :: kh1,kh2,kp1,kp2

  kh1 = (h1-1)/mo_num_per_kpt+1
  kh2 = (h2-1)/mo_num_per_kpt+1
  kp1 = (p1-1)/mo_num_per_kpt+1
  kp2 = (p2-1)/mo_num_per_kpt+1
  call double_allowed_kpts(kh1,kh2,kp1,kp2,is_allowed)
end subroutine

subroutine add_integrals_to_map_complex(mask_ijkl)
  use map_module
  use bitmasks
  implicit none

  BEGIN_DOC
  ! Adds integrals to tha MO map according to some bitmask
  END_DOC

  integer(bit_kind), intent(in)  :: mask_ijkl(N_int,4)

  integer                        :: i,j,k,l
  integer                        :: i0,j0,k0,l0
  double precision               :: c, cpu_1, cpu_2, wall_1, wall_2, wall_0

  integer, allocatable           :: list_ijkl(:,:)
  integer                        :: n_i, n_j, n_k, n_l
  integer, allocatable           :: two_e_tmp_0_idx(:)
  real(integral_kind), allocatable :: two_e_tmp_0(:,:)
  double precision, allocatable  :: two_e_tmp_1(:)
  double precision, allocatable  :: two_e_tmp_2(:,:)
  double precision, allocatable  :: two_e_tmp_3(:,:,:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: two_e_tmp_1, two_e_tmp_2, two_e_tmp_3

  integer                        :: n_integrals
  integer                        :: size_buffer
  integer(key_kind),allocatable  :: buffer_i(:)
  real(integral_kind),allocatable :: buffer_value(:)
  double precision, external     :: map_mb

  integer                        :: i1,j1,k1,l1, ii1, kmax, thread_num
  integer                        :: i2,i3,i4
  double precision,parameter     :: thr_coef = 1.d-10

  print*,'not implemented for complex',irp_here
  stop -1
!  PROVIDE ao_two_e_integrals_in_map  mo_coef
!
!  !Get list of MOs for i,j,k and l
!  !-------------------------------
!
!  allocate(list_ijkl(mo_num,4))
!  call bitstring_to_list( mask_ijkl(1,1), list_ijkl(1,1), n_i, N_int )
!  call bitstring_to_list( mask_ijkl(1,2), list_ijkl(1,2), n_j, N_int )
!  call bitstring_to_list( mask_ijkl(1,3), list_ijkl(1,3), n_k, N_int )
!  call bitstring_to_list( mask_ijkl(1,4), list_ijkl(1,4), n_l, N_int )
!  j = 0
!  do i = 1, N_int
!    j += popcnt(mask_ijkl(i,1))
!  enddo
!  if(j==0)then
!    return
!  endif
!
!  j = 0
!  do i = 1, N_int
!    j += popcnt(mask_ijkl(i,2))
!  enddo
!  if(j==0)then
!    return
!  endif
!
!  j = 0
!  do i = 1, N_int
!    j += popcnt(mask_ijkl(i,3))
!  enddo
!  if(j==0)then
!    return
!  endif
!
!  j = 0
!  do i = 1, N_int
!    j += popcnt(mask_ijkl(i,4))
!  enddo
!  if(j==0)then
!    return
!  endif
!
!  size_buffer = min(ao_num*ao_num*ao_num,16000000)
!  print*, 'Buffers : ', 8.*(mo_num*(n_j)*(n_k+1) + mo_num+&
!      ao_num+ao_num*ao_num+ size_buffer*3)/(1024*1024), 'MB / core'
!
!  double precision               :: accu_bis
!  accu_bis = 0.d0
!  call wall_time(wall_1)
!
!  !$OMP PARALLEL PRIVATE(l1,k1,j1,i1,i2,i3,i4,i,j,k,l,c, ii1,kmax,   &
!      !$OMP  two_e_tmp_0_idx, two_e_tmp_0, two_e_tmp_1,two_e_tmp_2,two_e_tmp_3,&
!      !$OMP  buffer_i,buffer_value,n_integrals,wall_2,i0,j0,k0,l0,   &
!      !$OMP  wall_0,thread_num,accu_bis)                             &
!      !$OMP  DEFAULT(NONE)                                           &
!      !$OMP  SHARED(size_buffer,ao_num,mo_num,n_i,n_j,n_k,n_l,   &
!      !$OMP  mo_coef_transp,                                         &
!      !$OMP  mo_coef_transp_is_built, list_ijkl,                     &
!      !$OMP  mo_coef_is_built, wall_1,                               &
!      !$OMP  mo_coef,mo_integrals_threshold,mo_integrals_map)
!  n_integrals = 0
!  wall_0 = wall_1
!  allocate(two_e_tmp_3(mo_num, n_j, n_k),                 &
!      two_e_tmp_1(mo_num),                                &
!      two_e_tmp_0(ao_num,ao_num),                                   &
!      two_e_tmp_0_idx(ao_num),                                      &
!      two_e_tmp_2(mo_num, n_j),                           &
!      buffer_i(size_buffer),                                         &
!      buffer_value(size_buffer) )
!
!  thread_num = 0
!  !$  thread_num = omp_get_thread_num()
!  !$OMP DO SCHEDULE(guided)
!  do l1 = 1,ao_num
!    two_e_tmp_3 = 0.d0
!    do k1 = 1,ao_num
!      two_e_tmp_2 = 0.d0
!      do j1 = 1,ao_num
!        call get_ao_two_e_integrals(j1,k1,l1,ao_num,two_e_tmp_0(1,j1))
!        ! call compute_ao_two_e_integrals(j1,k1,l1,ao_num,two_e_tmp_0(1,j1))
!      enddo
!      do j1 = 1,ao_num
!        kmax = 0
!        do i1 = 1,ao_num
!          c = two_e_tmp_0(i1,j1)
!          if (c == 0.d0) then
!            cycle
!          endif
!          kmax += 1
!          two_e_tmp_0(kmax,j1) = c
!          two_e_tmp_0_idx(kmax) = i1
!        enddo
!
!        if (kmax==0) then
!          cycle
!        endif
!
!        two_e_tmp_1 = 0.d0
!        ii1=1
!        do ii1 = 1,kmax-4,4
!          i1 = two_e_tmp_0_idx(ii1)
!          i2 = two_e_tmp_0_idx(ii1+1)
!          i3 = two_e_tmp_0_idx(ii1+2)
!          i4 = two_e_tmp_0_idx(ii1+3)
!          do i = list_ijkl(1,1), list_ijkl(n_i,1)
!            two_e_tmp_1(i)  =  two_e_tmp_1(i) +                    &
!                mo_coef_transp(i,i1) * two_e_tmp_0(ii1,j1) +        &
!                mo_coef_transp(i,i2) * two_e_tmp_0(ii1+1,j1) +      &
!                mo_coef_transp(i,i3) * two_e_tmp_0(ii1+2,j1) +      &
!                mo_coef_transp(i,i4) * two_e_tmp_0(ii1+3,j1)
!          enddo ! i
!        enddo  ! ii1
!
!        i2 = ii1
!        do ii1 = i2,kmax
!          i1 = two_e_tmp_0_idx(ii1)
!          do i = list_ijkl(1,1), list_ijkl(n_i,1)
!            two_e_tmp_1(i) = two_e_tmp_1(i) + mo_coef_transp(i,i1) * two_e_tmp_0(ii1,j1)
!          enddo ! i
!        enddo  ! ii1
!        c = 0.d0
!
!        do i = list_ijkl(1,1), list_ijkl(n_i,1)
!          c = max(c,abs(two_e_tmp_1(i)))
!          if (c>mo_integrals_threshold) exit
!        enddo
!        if ( c < mo_integrals_threshold ) then
!          cycle
!        endif
!
!        do j0 = 1, n_j
!          j = list_ijkl(j0,2)
!          c = mo_coef_transp(j,j1)
!          if (abs(c) < thr_coef) then
!            cycle
!          endif
!          do i = list_ijkl(1,1), list_ijkl(n_i,1)
!            two_e_tmp_2(i,j0)  = two_e_tmp_2(i,j0) + c * two_e_tmp_1(i)
!          enddo ! i
!        enddo  ! j
!      enddo !j1
!      if ( maxval(abs(two_e_tmp_2)) < mo_integrals_threshold ) then
!        cycle
!      endif
!
!
!      do k0 = 1, n_k
!        k = list_ijkl(k0,3)
!        c = mo_coef_transp(k,k1)
!        if (abs(c) < thr_coef) then
!          cycle
!        endif
!
!        do j0 = 1, n_j
!          j = list_ijkl(j0,2)
!          do i = list_ijkl(1,1), k
!            two_e_tmp_3(i,j0,k0) = two_e_tmp_3(i,j0,k0) + c* two_e_tmp_2(i,j0)
!          enddo!i
!        enddo !j
!
!      enddo  !k
!    enddo   !k1
!
!
!
!    do l0 = 1,n_l
!      l = list_ijkl(l0,4)
!      c = mo_coef_transp(l,l1)
!      if (abs(c) < thr_coef) then
!        cycle
!      endif
!      j1 = shiftr((l*l-l),1)
!      do j0 = 1, n_j
!        j = list_ijkl(j0,2)
!        if (j > l)  then
!          exit
!        endif
!        j1 += 1
!        do k0 = 1, n_k
!          k = list_ijkl(k0,3)
!          i1 = shiftr((k*k-k),1)
!          if (i1<=j1) then
!            continue
!          else
!            exit
!          endif
!          two_e_tmp_1 = 0.d0
!          do i0 = 1, n_i
!            i = list_ijkl(i0,1)
!            if (i>k) then
!              exit
!            endif
!            two_e_tmp_1(i) = c*two_e_tmp_3(i,j0,k0)
!            !           i1+=1
!          enddo
!
!          do i0 = 1, n_i
!            i = list_ijkl(i0,1)
!            if(i> min(k,j1-i1+list_ijkl(1,1)-1))then
!              exit
!            endif
!            if (abs(two_e_tmp_1(i)) < mo_integrals_threshold) then
!              cycle
!            endif
!            n_integrals += 1
!            buffer_value(n_integrals) = two_e_tmp_1(i)
!            !DIR$ FORCEINLINE
!            call mo_two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
!            if (n_integrals == size_buffer) then
!              call insert_into_mo_integrals_map(n_integrals,buffer_i,buffer_value,&
!                  real(mo_integrals_threshold,integral_kind))
!              n_integrals = 0
!            endif
!          enddo
!        enddo
!      enddo
!    enddo
!
!    call wall_time(wall_2)
!    if (thread_num == 0) then
!      if (wall_2 - wall_0 > 1.d0) then
!        wall_0 = wall_2
!        print*, 100.*float(l1)/float(ao_num), '% in ',               &
!            wall_2-wall_1, 's', map_mb(mo_integrals_map) ,'MB'
!      endif
!    endif
!  enddo
!  !$OMP END DO NOWAIT
!  deallocate (two_e_tmp_1,two_e_tmp_2,two_e_tmp_3)
!
!  integer                        :: index_needed
!
!  call insert_into_mo_integrals_map(n_integrals,buffer_i,buffer_value,&
!      real(mo_integrals_threshold,integral_kind))
!  deallocate(buffer_i, buffer_value)
!  !$OMP END PARALLEL
!  call map_merge(mo_integrals_map)
!
!  call wall_time(wall_2)
!  call cpu_time(cpu_2)
!  integer*8                      :: get_mo_map_size, mo_map_size
!  mo_map_size = get_mo_map_size()
!
!  deallocate(list_ijkl)


end


subroutine add_integrals_to_map_three_indices_complex(mask_ijk)
  use map_module
  use bitmasks
  implicit none

  BEGIN_DOC
  ! Adds integrals to tha MO map according to some bitmask
  END_DOC

  integer(bit_kind), intent(in)  :: mask_ijk(N_int,3)

  integer                        :: i,j,k,l
  integer                        :: i0,j0,k0,l0
  double precision               :: c, cpu_1, cpu_2, wall_1, wall_2, wall_0

  integer, allocatable           :: list_ijkl(:,:)
  integer                        :: n_i, n_j, n_k
  integer                        :: m
  integer, allocatable           :: two_e_tmp_0_idx(:)
  real(integral_kind), allocatable :: two_e_tmp_0(:,:)
  double precision, allocatable  :: two_e_tmp_1(:)
  double precision, allocatable  :: two_e_tmp_2(:,:)
  double precision, allocatable  :: two_e_tmp_3(:,:,:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: two_e_tmp_1, two_e_tmp_2, two_e_tmp_3

  integer                        :: n_integrals
  integer                        :: size_buffer
  integer(key_kind),allocatable  :: buffer_i(:)
  real(integral_kind),allocatable :: buffer_value(:)
  double precision               :: map_mb

  integer                        :: i1,j1,k1,l1, ii1, kmax, thread_num
  integer                        :: i2,i3,i4
  double precision,parameter     :: thr_coef = 1.d-10

  print*,'not implemented for complex',irp_here
  stop -1
!  PROVIDE ao_two_e_integrals_in_map  mo_coef
!
!  !Get list of MOs for i,j,k and l
!  !-------------------------------
!
!  allocate(list_ijkl(mo_num,4))
!  call bitstring_to_list( mask_ijk(1,1), list_ijkl(1,1), n_i, N_int )
!  call bitstring_to_list( mask_ijk(1,2), list_ijkl(1,2), n_j, N_int )
!  call bitstring_to_list( mask_ijk(1,3), list_ijkl(1,3), n_k, N_int )
!  j = 0
!  do i = 1, N_int
!    j += popcnt(mask_ijk(i,1))
!  enddo
!  if(j==0)then
!    return
!  endif
!
!  j = 0
!  do i = 1, N_int
!    j += popcnt(mask_ijk(i,2))
!  enddo
!  if(j==0)then
!    return
!  endif
!
!  j = 0
!  do i = 1, N_int
!    j += popcnt(mask_ijk(i,3))
!  enddo
!  if(j==0)then
!    return
!  endif
!
!  size_buffer = min(ao_num*ao_num*ao_num,16000000)
!  print*, 'Providing the molecular integrals '
!  print*, 'Buffers : ', 8.*(mo_num*(n_j)*(n_k+1) + mo_num+&
!      ao_num+ao_num*ao_num+ size_buffer*3)/(1024*1024), 'MB / core'
!
!  call wall_time(wall_1)
!  call cpu_time(cpu_1)
!  double precision               :: accu_bis
!  accu_bis = 0.d0
!  !$OMP PARALLEL PRIVATE(m,l1,k1,j1,i1,i2,i3,i4,i,j,k,l,c, ii1,kmax, &
!      !$OMP  two_e_tmp_0_idx, two_e_tmp_0, two_e_tmp_1,two_e_tmp_2,two_e_tmp_3,&
!      !$OMP  buffer_i,buffer_value,n_integrals,wall_2,i0,j0,k0,l0,   &
!      !$OMP  wall_0,thread_num,accu_bis)                             &
!      !$OMP  DEFAULT(NONE)                                           &
!      !$OMP  SHARED(size_buffer,ao_num,mo_num,n_i,n_j,n_k,       &
!      !$OMP  mo_coef_transp,                                         &
!      !$OMP  mo_coef_transp_is_built, list_ijkl,                     &
!      !$OMP  mo_coef_is_built, wall_1,                               &
!      !$OMP  mo_coef,mo_integrals_threshold,mo_integrals_map)
!  n_integrals = 0
!  wall_0 = wall_1
!  allocate(two_e_tmp_3(mo_num, n_j, n_k),                 &
!      two_e_tmp_1(mo_num),                                &
!      two_e_tmp_0(ao_num,ao_num),                             &
!      two_e_tmp_0_idx(ao_num),                                &
!      two_e_tmp_2(mo_num, n_j),                           &
!      buffer_i(size_buffer),                                   &
!      buffer_value(size_buffer) )
!
!  thread_num = 0
!  !$  thread_num = omp_get_thread_num()
!  !$OMP DO SCHEDULE(guided)
!  do l1 = 1,ao_num
!    two_e_tmp_3 = 0.d0
!    do k1 = 1,ao_num
!      two_e_tmp_2 = 0.d0
!      do j1 = 1,ao_num
!        call get_ao_two_e_integrals(j1,k1,l1,ao_num,two_e_tmp_0(1,j1))
!      enddo
!      do j1 = 1,ao_num
!        kmax = 0
!        do i1 = 1,ao_num
!          c = two_e_tmp_0(i1,j1)
!          if (c == 0.d0) then
!            cycle
!          endif
!          kmax += 1
!          two_e_tmp_0(kmax,j1) = c
!          two_e_tmp_0_idx(kmax) = i1
!        enddo
!
!        if (kmax==0) then
!          cycle
!        endif
!
!        two_e_tmp_1 = 0.d0
!        ii1=1
!        do ii1 = 1,kmax-4,4
!          i1 = two_e_tmp_0_idx(ii1)
!          i2 = two_e_tmp_0_idx(ii1+1)
!          i3 = two_e_tmp_0_idx(ii1+2)
!          i4 = two_e_tmp_0_idx(ii1+3)
!          do i = list_ijkl(1,1), list_ijkl(n_i,1)
!            two_e_tmp_1(i)  =  two_e_tmp_1(i) +                    &
!                mo_coef_transp(i,i1) * two_e_tmp_0(ii1,j1) +        &
!                mo_coef_transp(i,i2) * two_e_tmp_0(ii1+1,j1) +      &
!                mo_coef_transp(i,i3) * two_e_tmp_0(ii1+2,j1) +      &
!                mo_coef_transp(i,i4) * two_e_tmp_0(ii1+3,j1)
!          enddo ! i
!        enddo  ! ii1
!
!        i2 = ii1
!        do ii1 = i2,kmax
!          i1 = two_e_tmp_0_idx(ii1)
!          do i = list_ijkl(1,1), list_ijkl(n_i,1)
!            two_e_tmp_1(i) = two_e_tmp_1(i) + mo_coef_transp(i,i1) * two_e_tmp_0(ii1,j1)
!          enddo ! i
!        enddo  ! ii1
!        c = 0.d0
!
!        do i = list_ijkl(1,1), list_ijkl(n_i,1)
!          c = max(c,abs(two_e_tmp_1(i)))
!          if (c>mo_integrals_threshold) exit
!        enddo
!        if ( c < mo_integrals_threshold ) then
!          cycle
!        endif
!
!        do j0 = 1, n_j
!          j = list_ijkl(j0,2)
!          c = mo_coef_transp(j,j1)
!          if (abs(c) < thr_coef) then
!            cycle
!          endif
!          do i = list_ijkl(1,1), list_ijkl(n_i,1)
!            two_e_tmp_2(i,j0)  = two_e_tmp_2(i,j0) + c * two_e_tmp_1(i)
!          enddo ! i
!        enddo  ! j
!      enddo !j1
!      if ( maxval(abs(two_e_tmp_2)) < mo_integrals_threshold ) then
!        cycle
!      endif
!
!
!      do k0 = 1, n_k
!        k = list_ijkl(k0,3)
!        c = mo_coef_transp(k,k1)
!        if (abs(c) < thr_coef) then
!          cycle
!        endif
!
!        do j0 = 1, n_j
!          j = list_ijkl(j0,2)
!          do i = list_ijkl(1,1), k
!            two_e_tmp_3(i,j0,k0) = two_e_tmp_3(i,j0,k0) + c* two_e_tmp_2(i,j0)
!          enddo!i
!        enddo !j
!
!      enddo  !k
!    enddo   !k1
!
!
!
!    do l0 = 1,n_j
!      l = list_ijkl(l0,2)
!      c = mo_coef_transp(l,l1)
!      if (abs(c) < thr_coef) then
!        cycle
!      endif
!      do k0 = 1, n_k
!        k = list_ijkl(k0,3)
!        i1 = shiftr((k*k-k),1)
!        two_e_tmp_1 = 0.d0
!        j0 = l0
!        j = list_ijkl(j0,2)
!        do i0 = 1, n_i
!          i = list_ijkl(i0,1)
!          if (i>k) then
!            exit
!          endif
!          two_e_tmp_1(i) = c*two_e_tmp_3(i,j0,k0)
!        enddo
!
!        do i0 = 1, n_i
!          i = list_ijkl(i0,1)
!          if (i>k) then !min(k,j1-i1)
!            exit
!          endif
!          if (abs(two_e_tmp_1(i)) < mo_integrals_threshold) then
!            cycle
!          endif
!          n_integrals += 1
!          buffer_value(n_integrals) = two_e_tmp_1(i)
!          if(i==k .and. j==l .and. i.ne.j)then
!            buffer_value(n_integrals) = buffer_value(n_integrals) *0.5d0
!          endif
!          !DIR$ FORCEINLINE
!          call mo_two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
!          if (n_integrals == size_buffer) then
!            call insert_into_mo_integrals_map(n_integrals,buffer_i,buffer_value,&
!                real(mo_integrals_threshold,integral_kind))
!            n_integrals = 0
!          endif
!        enddo
!      enddo
!    enddo
!
!    do l0 = 1,n_j
!      l = list_ijkl(l0,2)
!      c = mo_coef_transp(l,l1)
!      if (abs(c) < thr_coef) then
!        cycle
!      endif
!      do k0 = 1, n_k
!        k = list_ijkl(k0,3)
!        i1 = shiftr((k*k-k),1)
!        two_e_tmp_1 = 0.d0
!        j0 = k0
!        j = list_ijkl(k0,2)
!        i0 = l0
!        i = list_ijkl(i0,2)
!        if (k==l) then
!          cycle
!        endif
!        two_e_tmp_1(i) = c*two_e_tmp_3(i,j0,k0)
!
!        n_integrals += 1
!        buffer_value(n_integrals) = two_e_tmp_1(i)
!        !DIR$ FORCEINLINE
!        call mo_two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
!        if (n_integrals == size_buffer) then
!          call insert_into_mo_integrals_map(n_integrals,buffer_i,buffer_value,&
!              real(mo_integrals_threshold,integral_kind))
!          n_integrals = 0
!        endif
!      enddo
!    enddo
!
!    call wall_time(wall_2)
!    if (thread_num == 0) then
!      if (wall_2 - wall_0 > 1.d0) then
!        wall_0 = wall_2
!        print*, 100.*float(l1)/float(ao_num), '% in ',               &
!            wall_2-wall_1, 's', map_mb(mo_integrals_map) ,'MB'
!      endif
!    endif
!  enddo
!  !$OMP END DO NOWAIT
!  deallocate (two_e_tmp_1,two_e_tmp_2,two_e_tmp_3)
!
!  integer                        :: index_needed
!
!  call insert_into_mo_integrals_map(n_integrals,buffer_i,buffer_value,&
!      real(mo_integrals_threshold,integral_kind))
!  deallocate(buffer_i, buffer_value)
!  !$OMP END PARALLEL
!  call map_merge(mo_integrals_map)
!
!  call wall_time(wall_2)
!  call cpu_time(cpu_2)
!  integer*8                      :: get_mo_map_size, mo_map_size
!  mo_map_size = get_mo_map_size()
!
!  deallocate(list_ijkl)
!
!
!  print*,'Molecular integrals provided:'
!  print*,' Size of MO map           ', map_mb(mo_integrals_map) ,'MB'
!  print*,' Number of MO integrals: ',  mo_map_size
!  print*,' cpu  time :',cpu_2 - cpu_1, 's'
!  print*,' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1), ')'

end


subroutine add_integrals_to_map_no_exit_34_complex(mask_ijkl)
  use map_module
  use bitmasks
  implicit none

  BEGIN_DOC
  ! Adds integrals to tha MO map according to some bitmask
  END_DOC

  integer(bit_kind), intent(in)  :: mask_ijkl(N_int,4)

  integer                        :: i,j,k,l
  integer                        :: i0,j0,k0,l0
  double precision               :: c, cpu_1, cpu_2, wall_1, wall_2, wall_0

  integer, allocatable           :: list_ijkl(:,:)
  integer                        :: n_i, n_j, n_k, n_l
  integer, allocatable           :: two_e_tmp_0_idx(:)
  real(integral_kind), allocatable :: two_e_tmp_0(:,:)
  double precision, allocatable  :: two_e_tmp_1(:)
  double precision, allocatable  :: two_e_tmp_2(:,:)
  double precision, allocatable  :: two_e_tmp_3(:,:,:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: two_e_tmp_1, two_e_tmp_2, two_e_tmp_3

  integer                        :: n_integrals
  integer                        :: size_buffer
  integer(key_kind),allocatable  :: buffer_i(:)
  real(integral_kind),allocatable :: buffer_value(:)
  double precision               :: map_mb

  integer                        :: i1,j1,k1,l1, ii1, kmax, thread_num
  integer                        :: i2,i3,i4
  double precision,parameter     :: thr_coef = 1.d-10

  print*,'not implemented for complex',irp_here
  stop -1
!  PROVIDE ao_two_e_integrals_in_map  mo_coef
!
!  !Get list of MOs for i,j,k and l
!  !-------------------------------
!
!  allocate(list_ijkl(mo_num,4))
!  call bitstring_to_list( mask_ijkl(1,1), list_ijkl(1,1), n_i, N_int )
!  call bitstring_to_list( mask_ijkl(1,2), list_ijkl(1,2), n_j, N_int )
!  call bitstring_to_list( mask_ijkl(1,3), list_ijkl(1,3), n_k, N_int )
!  call bitstring_to_list( mask_ijkl(1,4), list_ijkl(1,4), n_l, N_int )
!
!  size_buffer = min(ao_num*ao_num*ao_num,16000000)
!  print*, 'Providing the molecular integrals '
!  print*, 'Buffers : ', 8.*(mo_num*(n_j)*(n_k+1) + mo_num+&
!      ao_num+ao_num*ao_num+ size_buffer*3)/(1024*1024), 'MB / core'
!
!  call wall_time(wall_1)
!  call cpu_time(cpu_1)
!
!  !$OMP PARALLEL PRIVATE(l1,k1,j1,i1,i2,i3,i4,i,j,k,l,c, ii1,kmax,   &
!      !$OMP  two_e_tmp_0_idx, two_e_tmp_0, two_e_tmp_1,two_e_tmp_2,two_e_tmp_3,&
!      !$OMP  buffer_i,buffer_value,n_integrals,wall_2,i0,j0,k0,l0,   &
!      !$OMP  wall_0,thread_num)                                      &
!      !$OMP  DEFAULT(NONE)                                           &
!      !$OMP  SHARED(size_buffer,ao_num,mo_num,n_i,n_j,n_k,n_l,   &
!      !$OMP  mo_coef_transp,                                         &
!      !$OMP  mo_coef_transp_is_built, list_ijkl,                     &
!      !$OMP  mo_coef_is_built, wall_1,                               &
!      !$OMP  mo_coef,mo_integrals_threshold,mo_integrals_map)
!  n_integrals = 0
!  wall_0 = wall_1
!  allocate(two_e_tmp_3(mo_num, n_j, n_k),                 &
!      two_e_tmp_1(mo_num),                                &
!      two_e_tmp_0(ao_num,ao_num),                                   &
!      two_e_tmp_0_idx(ao_num),                                      &
!      two_e_tmp_2(mo_num, n_j),                           &
!      buffer_i(size_buffer),                                         &
!      buffer_value(size_buffer) )
!
!  thread_num = 0
!  !$  thread_num = omp_get_thread_num()
!  !$OMP DO SCHEDULE(guided)
!  do l1 = 1,ao_num
!    !IRP_IF COARRAY
!    !    if (mod(l1-this_image(),num_images()) /= 0 ) then
!    !      cycle
!    !    endif
!    !IRP_ENDIF
!    two_e_tmp_3 = 0.d0
!    do k1 = 1,ao_num
!      two_e_tmp_2 = 0.d0
!      do j1 = 1,ao_num
!        call get_ao_two_e_integrals(j1,k1,l1,ao_num,two_e_tmp_0(1,j1))
!        ! call compute_ao_two_e_integrals(j1,k1,l1,ao_num,two_e_tmp_0(1,j1))
!      enddo
!      do j1 = 1,ao_num
!        kmax = 0
!        do i1 = 1,ao_num
!          c = two_e_tmp_0(i1,j1)
!          if (c == 0.d0) then
!            cycle
!          endif
!          kmax += 1
!          two_e_tmp_0(kmax,j1) = c
!          two_e_tmp_0_idx(kmax) = i1
!        enddo
!
!        if (kmax==0) then
!          cycle
!        endif
!
!        two_e_tmp_1 = 0.d0
!        ii1=1
!        do ii1 = 1,kmax-4,4
!          i1 = two_e_tmp_0_idx(ii1)
!          i2 = two_e_tmp_0_idx(ii1+1)
!          i3 = two_e_tmp_0_idx(ii1+2)
!          i4 = two_e_tmp_0_idx(ii1+3)
!          do i = list_ijkl(1,1), list_ijkl(n_i,1)
!            two_e_tmp_1(i)  =  two_e_tmp_1(i) +                    &
!                mo_coef_transp(i,i1) * two_e_tmp_0(ii1,j1) +        &
!                mo_coef_transp(i,i2) * two_e_tmp_0(ii1+1,j1) +      &
!                mo_coef_transp(i,i3) * two_e_tmp_0(ii1+2,j1) +      &
!                mo_coef_transp(i,i4) * two_e_tmp_0(ii1+3,j1)
!          enddo ! i
!        enddo  ! ii1
!
!        i2 = ii1
!        do ii1 = i2,kmax
!          i1 = two_e_tmp_0_idx(ii1)
!          do i = list_ijkl(1,1), list_ijkl(n_i,1)
!            two_e_tmp_1(i) = two_e_tmp_1(i) + mo_coef_transp(i,i1) * two_e_tmp_0(ii1,j1)
!          enddo ! i
!        enddo  ! ii1
!        c = 0.d0
!
!        do i = list_ijkl(1,1), list_ijkl(n_i,1)
!          c = max(c,abs(two_e_tmp_1(i)))
!          if (c>mo_integrals_threshold) exit
!        enddo
!        if ( c < mo_integrals_threshold ) then
!          cycle
!        endif
!
!        do j0 = 1, n_j
!          j = list_ijkl(j0,2)
!          c = mo_coef_transp(j,j1)
!          if (abs(c) < thr_coef) then
!            cycle
!          endif
!          do i = list_ijkl(1,1), list_ijkl(n_i,1)
!            two_e_tmp_2(i,j0)  = two_e_tmp_2(i,j0) + c * two_e_tmp_1(i)
!          enddo ! i
!        enddo  ! j
!      enddo !j1
!      if ( maxval(abs(two_e_tmp_2)) < mo_integrals_threshold ) then
!        cycle
!      endif
!
!
!      do k0 = 1, n_k
!        k = list_ijkl(k0,3)
!        c = mo_coef_transp(k,k1)
!        if (abs(c) < thr_coef) then
!          cycle
!        endif
!
!        do j0 = 1, n_j
!          j = list_ijkl(j0,2)
!          do i = list_ijkl(1,1), k
!            two_e_tmp_3(i,j0,k0) = two_e_tmp_3(i,j0,k0) + c* two_e_tmp_2(i,j0)
!          enddo!i
!        enddo !j
!
!      enddo  !k
!    enddo   !k1
!
!
!
!    do l0 = 1,n_l
!      l = list_ijkl(l0,4)
!      c = mo_coef_transp(l,l1)
!      if (abs(c) < thr_coef) then
!        cycle
!      endif
!      j1 = shiftr((l*l-l),1)
!      do j0 = 1, n_j
!        j = list_ijkl(j0,2)
!        if (j > l)  then
!          exit
!        endif
!        j1 += 1
!        do k0 = 1, n_k
!          k = list_ijkl(k0,3)
!          i1 = shiftr((k*k-k),1)
!          two_e_tmp_1 = 0.d0
!          do i0 = 1, n_i
!            i = list_ijkl(i0,1)
!            if (i>k) then
!              exit
!            endif
!            two_e_tmp_1(i) = c*two_e_tmp_3(i,j0,k0)
!          enddo
!
!          do i0 = 1, n_i
!            i = list_ijkl(i0,1)
!            if(i> k)then
!              exit
!            endif
!
!            if (abs(two_e_tmp_1(i)) < mo_integrals_threshold) then
!              cycle
!            endif
!            n_integrals += 1
!            buffer_value(n_integrals) = two_e_tmp_1(i)
!            !DIR$ FORCEINLINE
!            call mo_two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
!            if (n_integrals == size_buffer) then
!              call insert_into_mo_integrals_map(n_integrals,buffer_i,buffer_value,&
!                  real(mo_integrals_threshold,integral_kind))
!              n_integrals = 0
!            endif
!          enddo
!        enddo
!      enddo
!    enddo
!
!    call wall_time(wall_2)
!    if (thread_num == 0) then
!      if (wall_2 - wall_0 > 1.d0) then
!        wall_0 = wall_2
!        print*, 100.*float(l1)/float(ao_num), '% in ',               &
!            wall_2-wall_1, 's', map_mb(mo_integrals_map) ,'MB'
!      endif
!    endif
!  enddo
!  !$OMP END DO NOWAIT
!  deallocate (two_e_tmp_1,two_e_tmp_2,two_e_tmp_3)
!
!  call insert_into_mo_integrals_map(n_integrals,buffer_i,buffer_value,&
!      real(mo_integrals_threshold,integral_kind))
!  deallocate(buffer_i, buffer_value)
!  !$OMP END PARALLEL
!  !IRP_IF COARRAY
!  !  print*, 'Communicating the map'
!  !  call communicate_mo_integrals()
!  !IRP_ENDIF
!  call map_merge(mo_integrals_map)
!
!  call wall_time(wall_2)
!  call cpu_time(cpu_2)
!  integer*8                      :: get_mo_map_size, mo_map_size
!  mo_map_size = get_mo_map_size()
!
!  deallocate(list_ijkl)
!
!
!  print*,'Molecular integrals provided:'
!  print*,' Size of MO map           ', map_mb(mo_integrals_map) ,'MB'
!  print*,' Number of MO integrals: ',  mo_map_size
!  print*,' cpu  time :',cpu_2 - cpu_1, 's'
!  print*,' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1), ')'


end
