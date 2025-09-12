! 1,2-index integrals are always taken from:
! - mo_two_e_integrals_jj_exchange
! - mo_two_e_integrals_jj_anti
! - mo_two_e_integrals_jj
! 
! 3-index integrals are always taken from:
! - big_array_exchange_integrals
! - big_array_coulomb_integrals
!
! If (do_mo_cholesky):
! - Integrals with four 4 active orbitals are stored in the cache map,
!   all other integrals are used from cholesky vectors
! - 1,2,3-index arrays are built from cholesky vectors
! Else:
! - All integrals are stored in the map or cache map
! - 1,2,3-index arrays are built from the map
!

subroutine mo_two_e_integrals_index(i,j,k,l,i1)
  use map_module
  implicit none
  BEGIN_DOC
  ! Computes an unique index for i,j,k,l integrals
  END_DOC
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


BEGIN_PROVIDER [ logical, mo_two_e_integrals_in_map ]
  use map_module
  implicit none
  BEGIN_DOC
  ! If True, the map of MO two-electron integrals is provided
  END_DOC
  integer(bit_kind)              :: mask_ijkl(N_int,4)
  integer(bit_kind)              :: mask_ijk(N_int,3)
  double precision               :: cpu_1, cpu_2, wall_1, wall_2

  PROVIDE mo_class

  mo_two_e_integrals_in_map = .True.
  if (read_mo_two_e_integrals) then
    print*,'Reading the MO integrals'
    call map_load_from_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_map)
    print*, 'MO integrals provided'
    return
  endif

  call wall_time(wall_1)
  call cpu_time(cpu_1)

  if (do_mo_cholesky) then
     PROVIDE cholesky_mo_transp
  else
    if (do_ao_cholesky) then
      call add_integrals_to_map_cholesky
    else if (dble(ao_num)**4 * 32.d-9 < dble(qp_max_mem)) then
      call four_idx_dgemm
    else
      call add_integrals_to_map(full_ijkl_bitmask_4)
    endif

    integer*8                      :: get_mo_map_size, mo_map_size
    mo_map_size = get_mo_map_size()
    double precision, external     :: map_mb
    print*,'Molecular integrals provided:'
    print*,' Size of MO map           ', map_mb(mo_integrals_map) ,'MB'
    print*,' Number of MO integrals: ',  mo_map_size
  endif

  call wall_time(wall_2)
  call cpu_time(cpu_2)

  print*,' cpu  time :',cpu_2 - cpu_1, 's'
  print*,' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1), ')'

  if (write_mo_two_e_integrals.and.mpi_master) then
    call ezfio_set_work_empty(.False.)
    call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_map)
    call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
  endif

END_PROVIDER

subroutine four_idx_dgemm
  implicit none
  integer :: p,q,r,s,i,j,k,l
  double precision, allocatable :: a1(:,:,:,:)
  double precision, allocatable :: a2(:,:,:,:)

  PROVIDE ao_two_e_integrals_in_map mo_coef

  print *, ''
  print *, 'DGEMM-based AO->MO Transformation'
  print *, '---------------------------------'
  print *, ''

  if (ao_num > 1289) then
    print *,  irp_here, ': Integer overflow in ao_num**3. Set do_ao_cholesky=.True.'
  endif

  allocate (a1(ao_num,ao_num,ao_num,ao_num))

  print *, 'Getting AOs'
  call get_ao_two_e_integrals(1,1,1,ao_num,a1(1,1,1,1))
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(q,r,s) COLLAPSE(2)
  do s=1,ao_num
    do r=1,ao_num
      do q=1,ao_num
        a1(:,q,r,s) = 0.d0
        call get_ao_two_e_integrals(q,r,s,ao_num,a1(1,q,r,s))
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO


  print *, '1st transformation'
  ! 1st transformation
  allocate (a2(ao_num,ao_num,ao_num,mo_num))
  call dgemm('T','N', (ao_num*ao_num*ao_num), mo_num, ao_num, 1.d0, a1, ao_num, mo_coef, ao_num, 0.d0, a2, (ao_num*ao_num*ao_num))

  ! 2nd transformation
  print *, '2nd transformation'
  deallocate (a1)
  allocate (a1(ao_num,ao_num,mo_num,mo_num))
  call dgemm('T','N', (ao_num*ao_num*mo_num), mo_num, ao_num, 1.d0, a2, ao_num, mo_coef, ao_num, 0.d0, a1, (ao_num*ao_num*mo_num))

  ! 3rd transformation
  print *, '3rd transformation'
  deallocate (a2)
  allocate (a2(ao_num,mo_num,mo_num,mo_num))
  call dgemm('T','N', (ao_num*mo_num*mo_num), mo_num, ao_num, 1.d0, a1, ao_num, mo_coef, ao_num, 0.d0, a2, (ao_num*mo_num*mo_num))

  ! 4th transformation
  print *, '4th transformation'
  deallocate (a1)
  allocate (a1(mo_num,mo_num,mo_num,mo_num))
  call dgemm('T','N', (mo_num*mo_num*mo_num), mo_num, ao_num, 1.d0, a2, ao_num, mo_coef, ao_num, 0.d0, a1, (mo_num*mo_num*mo_num))

  deallocate (a2)

  integer :: n_integrals, size_buffer
  integer(key_kind)  , allocatable :: buffer_i(:)
  real(integral_kind), allocatable :: buffer_value(:)
  size_buffer = min(ao_num*ao_num*ao_num,16000000)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l,buffer_value,buffer_i,n_integrals)
  allocate ( buffer_i(size_buffer), buffer_value(size_buffer) )

  n_integrals = 0
  !$OMP DO
  do l=1,mo_num
    do k=1,mo_num
      do j=1,l
        do i=1,k
            if (abs(a1(i,j,k,l)) < mo_integrals_threshold) then
              cycle
            endif
            n_integrals += 1
            buffer_value(n_integrals) = a1(i,j,k,l)
            !DIR$ FORCEINLINE
            call mo_two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
            if (n_integrals == size_buffer) then
              call map_append(mo_integrals_map, buffer_i, buffer_value, n_integrals)
              n_integrals = 0
            endif
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO

  call map_append(mo_integrals_map, buffer_i, buffer_value, n_integrals)

  deallocate(buffer_i, buffer_value)
  !$OMP END PARALLEL

  deallocate (a1)

  call map_sort(mo_integrals_map)
  call map_unique(mo_integrals_map)

end subroutine

subroutine add_integrals_to_map(mask_ijkl)
  use bitmasks
  implicit none

  BEGIN_DOC
  ! Adds integrals to the MO map according to some bitmask
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

  PROVIDE ao_two_e_integrals_in_map  mo_coef


  print *, ''
  print *, 'Sparse AO->MO Transformation'
  print *, '----------------------------'
  print *, ''

  !Get list of MOs for i,j,k and l
  !-------------------------------

  allocate(list_ijkl(mo_num,4))
  call bitstring_to_list( mask_ijkl(1,1), list_ijkl(1,1), n_i, N_int )
  call bitstring_to_list( mask_ijkl(1,2), list_ijkl(1,2), n_j, N_int )
  call bitstring_to_list( mask_ijkl(1,3), list_ijkl(1,3), n_k, N_int )
  call bitstring_to_list( mask_ijkl(1,4), list_ijkl(1,4), n_l, N_int )
  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,1))
  enddo
  if(j==0)then
    return
  endif

  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,2))
  enddo
  if(j==0)then
    return
  endif

  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,3))
  enddo
  if(j==0)then
    return
  endif

  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,4))
  enddo
  if(j==0)then
    return
  endif

  call wall_time(wall_1)

  size_buffer = min(ao_num*ao_num,8000000)
  print*, 'Buffers : ', 8.*(mo_num*(n_j)*(n_k+1) + mo_num+&
      ao_num+ao_num*ao_num+ size_buffer*3)/(1024*1024), 'MB / core'

  !$OMP PARALLEL PRIVATE(l1,k1,j1,i1,i2,i3,i4,i,j,k,l,c, ii1,kmax,   &
      !$OMP  two_e_tmp_0_idx, two_e_tmp_0, two_e_tmp_1,two_e_tmp_2,two_e_tmp_3,&
      !$OMP  buffer_i,buffer_value,n_integrals,wall_2,i0,j0,k0,l0,   &
      !$OMP  wall_0,thread_num)                             &
      !$OMP  DEFAULT(NONE)                                           &
      !$OMP  SHARED(size_buffer,ao_num,mo_num,n_i,n_j,n_k,n_l,   &
      !$OMP  mo_coef_transp,                                         &
      !$OMP  mo_coef_transp_is_built, list_ijkl,                     &
      !$OMP  mo_coef_is_built, wall_1,                               &
      !$OMP  mo_coef,mo_integrals_threshold,mo_integrals_map)

  thread_num = 0
  !$  thread_num = omp_get_thread_num()

  n_integrals = 0
  wall_0 = wall_1
  allocate(two_e_tmp_3(mo_num, n_j, n_k),                 &
      two_e_tmp_1(mo_num),                                &
      two_e_tmp_0(ao_num,ao_num),                                   &
      two_e_tmp_0_idx(ao_num),                                      &
      two_e_tmp_2(mo_num, n_j),                           &
      buffer_i(size_buffer),                                         &
      buffer_value(size_buffer) )

  !$OMP DO SCHEDULE(guided)
  do l1 = 1,ao_num
    two_e_tmp_3 = 0.d0
    do k1 = 1,ao_num
      two_e_tmp_2 = 0.d0
      do j1 = 1,ao_num
        call get_ao_two_e_integrals(j1,k1,l1,ao_num,two_e_tmp_0(1,j1))
      enddo
      do j1 = 1,ao_num
        kmax = 0
        do i1 = 1,ao_num
          c = two_e_tmp_0(i1,j1)
          if (c == 0.d0) then
            cycle
          endif
          kmax += 1
          two_e_tmp_0(kmax,j1) = c
          two_e_tmp_0_idx(kmax) = i1
        enddo

        if (kmax==0) then
          cycle
        endif

        two_e_tmp_1 = 0.d0
        ii1=1
        do ii1 = 1,kmax-4,4
          i1 = two_e_tmp_0_idx(ii1)
          i2 = two_e_tmp_0_idx(ii1+1)
          i3 = two_e_tmp_0_idx(ii1+2)
          i4 = two_e_tmp_0_idx(ii1+3)
          do i = list_ijkl(1,1), list_ijkl(n_i,1)
            two_e_tmp_1(i)  =  two_e_tmp_1(i) +                    &
                mo_coef_transp(i,i1) * two_e_tmp_0(ii1,j1) +        &
                mo_coef_transp(i,i2) * two_e_tmp_0(ii1+1,j1) +      &
                mo_coef_transp(i,i3) * two_e_tmp_0(ii1+2,j1) +      &
                mo_coef_transp(i,i4) * two_e_tmp_0(ii1+3,j1)
          enddo ! i
        enddo  ! ii1

        i2 = ii1
        do ii1 = i2,kmax
          i1 = two_e_tmp_0_idx(ii1)
          do i = list_ijkl(1,1), list_ijkl(n_i,1)
            two_e_tmp_1(i) = two_e_tmp_1(i) + mo_coef_transp(i,i1) * two_e_tmp_0(ii1,j1)
          enddo ! i
        enddo  ! ii1
        c = 0.d0

        do i = list_ijkl(1,1), list_ijkl(n_i,1)
          c = max(c,abs(two_e_tmp_1(i)))
          if (c>mo_integrals_threshold) exit
        enddo
        if ( c < mo_integrals_threshold ) then
          cycle
        endif

        do j0 = 1, n_j
          j = list_ijkl(j0,2)
          c = mo_coef_transp(j,j1)
          if (abs(c) < thr_coef) then
            cycle
          endif
          do i = list_ijkl(1,1), list_ijkl(n_i,1)
            two_e_tmp_2(i,j0)  = two_e_tmp_2(i,j0) + c * two_e_tmp_1(i)
          enddo ! i
        enddo  ! j
      enddo !j1
      if ( maxval(abs(two_e_tmp_2)) < mo_integrals_threshold ) then
        cycle
      endif


      do k0 = 1, n_k
        k = list_ijkl(k0,3)
        c = mo_coef_transp(k,k1)
        if (abs(c) < thr_coef) then
          cycle
        endif

        do j0 = 1, n_j
          j = list_ijkl(j0,2)
          do i = list_ijkl(1,1), k
            two_e_tmp_3(i,j0,k0) = two_e_tmp_3(i,j0,k0) + c* two_e_tmp_2(i,j0)
          enddo!i
        enddo !j

      enddo  !k
    enddo   !k1



    do l0 = 1,n_l
      l = list_ijkl(l0,4)
      c = mo_coef_transp(l,l1)
      if (abs(c) < thr_coef) then
        cycle
      endif
      j1 = shiftr((l*l-l),1)
      do j0 = 1, n_j
        j = list_ijkl(j0,2)
        if (j > l)  then
          exit
        endif
        j1 += 1
        do k0 = 1, n_k
          k = list_ijkl(k0,3)
          i1 = shiftr((k*k-k),1)
          if (i1<=j1) then
            continue
          else
            exit
          endif
          two_e_tmp_1 = 0.d0
          do i0 = 1, n_i
            i = list_ijkl(i0,1)
            if (i>k) then
              exit
            endif
            two_e_tmp_1(i) = c*two_e_tmp_3(i,j0,k0)
            !           i1+=1
          enddo

          do i0 = 1, n_i
            i = list_ijkl(i0,1)
            if(i> min(k,j1-i1+list_ijkl(1,1)-1))then
              exit
            endif
            if (abs(two_e_tmp_1(i)) < mo_integrals_threshold) then
              cycle
            endif
            n_integrals += 1
            buffer_value(n_integrals) = two_e_tmp_1(i)
            !DIR$ FORCEINLINE
            call mo_two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
            if (n_integrals == size_buffer) then
              call insert_into_mo_integrals_map(n_integrals,buffer_i,buffer_value,&
                  real(mo_integrals_threshold,integral_kind))
              n_integrals = 0
            endif
          enddo
        enddo
      enddo
    enddo

    call wall_time(wall_2)
    if (thread_num == 0) then
      if (wall_2 - wall_0 > 1.d0) then
        wall_0 = wall_2
        print*, 100.*float(l1)/float(ao_num), '% in ',               &
            wall_2-wall_1, 's', map_mb(mo_integrals_map) ,'MB'
      endif
    endif
  enddo
  !$OMP END DO NOWAIT
  deallocate (two_e_tmp_1,two_e_tmp_2,two_e_tmp_3)

  if (n_integrals > 0) then
    call insert_into_mo_integrals_map(n_integrals,buffer_i,buffer_value,&
        real(mo_integrals_threshold,integral_kind))
  endif
  deallocate(buffer_i, buffer_value)
  !$OMP END PARALLEL
  call map_merge(mo_integrals_map)

  deallocate(list_ijkl)


end




subroutine add_integrals_to_map_cholesky
  use bitmasks
  implicit none

  BEGIN_DOC
  ! Adds integrals to the MO map using Cholesky vectors
  END_DOC

  integer :: i,j,k,l,m
  integer :: size_buffer, n_integrals
  size_buffer = min(mo_num*mo_num*mo_num,16000000)

  double precision, allocatable :: Vtmp(:,:,:)
  integer(key_kind)  , allocatable :: buffer_i(:)
  real(integral_kind), allocatable :: buffer_value(:)

  PROVIDE cholesky_mo_transp
  call set_multiple_levels_omp(.False.)

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(i,j,k,l,n_integrals,buffer_value, buffer_i, Vtmp)
  allocate (buffer_i(size_buffer), buffer_value(size_buffer))
  allocate (Vtmp(mo_num,mo_num,mo_num))
  n_integrals = 0

  !$OMP DO SCHEDULE(dynamic)
  do l=1,mo_num
    call dgemm('T','N',mo_num*mo_num,mo_num,cholesky_mo_num,1.d0, &
       cholesky_mo_transp, cholesky_mo_num, &
       cholesky_mo_transp(1,1,l), cholesky_mo_num, 0.d0, &
       Vtmp, mo_num*mo_num)

    do k=1,l
      do j=1,mo_num
        do i=1,j
          if (dabs(Vtmp(i,j,k)) > mo_integrals_threshold) then
            n_integrals = n_integrals + 1
            buffer_value(n_integrals) = Vtmp(i,j,k)
            !DIR$ FORCEINLINE
            call mo_two_e_integrals_index(i,k,j,l,buffer_i(n_integrals))
            if (n_integrals == size_buffer) then
              call map_append(mo_integrals_map, buffer_i, buffer_value, n_integrals)
              n_integrals = 0
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  if (n_integrals > 0) then
    call map_append(mo_integrals_map, buffer_i, buffer_value, n_integrals)
  endif
  deallocate(buffer_i, buffer_value, Vtmp)
  !$OMP BARRIER
  !$OMP END PARALLEL

  call map_sort(mo_integrals_map)
  call map_unique(mo_integrals_map)

end




 BEGIN_PROVIDER [ double precision, mo_two_e_integrals_jj, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, mo_two_e_integrals_jj_exchange, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, mo_two_e_integrals_jj_anti, (mo_num,mo_num) ]
  implicit none
  BEGIN_DOC
  ! mo_two_e_integrals_jj(i,j) = J_ij
  ! mo_two_e_integrals_jj_exchange(i,j) = K_ij
  ! mo_two_e_integrals_jj_anti(i,j) = J_ij - K_ij
  END_DOC

  integer                        :: i,j,k
  double precision               :: get_two_e_integral


  if (do_mo_cholesky) then
    double precision, allocatable :: buffer(:,:)
    allocate (buffer(cholesky_mo_num,mo_num))
    do k=1,cholesky_mo_num
      do i=1,mo_num
        buffer(k,i) = cholesky_mo_transp(k,i,i)
      enddo
    enddo
    call dgemm('T','N',mo_num,mo_num,cholesky_mo_num,1.d0, &
      buffer, cholesky_mo_num, buffer, cholesky_mo_num, 0.d0, mo_two_e_integrals_jj, mo_num)
    deallocate(buffer)

    do j=1,mo_num
      do i=1,mo_num
          mo_two_e_integrals_jj_exchange(i,j) = 0.d0
          do k=1,cholesky_mo_num
            mo_two_e_integrals_jj_exchange(i,j) = mo_two_e_integrals_jj_exchange(i,j) + &
                          cholesky_mo_transp(k,i,j)*cholesky_mo_transp(k,j,i)
          enddo
      enddo
    enddo

  else

    do j=1,mo_num
      do i=1,mo_num
        mo_two_e_integrals_jj(i,j) = get_two_e_integral(i,j,i,j,mo_integrals_map)
        mo_two_e_integrals_jj_exchange(i,j) = get_two_e_integral(i,j,j,i,mo_integrals_map)
      enddo
    enddo

  endif

  do j=1,mo_num
    do i=1,mo_num
        mo_two_e_integrals_jj_anti(i,j) = mo_two_e_integrals_jj(i,j) - mo_two_e_integrals_jj_exchange(i,j)
    enddo
  enddo

END_PROVIDER


subroutine clear_mo_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the MO map
  END_DOC
  call map_deinit(mo_integrals_map)
  FREE mo_integrals_map mo_two_e_integrals_jj mo_two_e_integrals_jj_anti
  FREE mo_two_e_integrals_jj_exchange mo_two_e_integrals_in_map
end

