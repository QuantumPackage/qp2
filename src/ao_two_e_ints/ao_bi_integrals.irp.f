subroutine ao_two_e_integrals_index(i,j,k,l,i1)
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

BEGIN_TEMPLATE 

BEGIN_PROVIDER [ logical, ao_two_e_integrals$_erf_in_map ]
  use map_module
  implicit none
  BEGIN_DOC
  ! If True, the map of AO two-electron integrals$_erf is provided
  END_DOC
  integer(bit_kind)              :: mask_ijkl(N_int,4)
  integer(bit_kind)              :: mask_ijk(N_int,3)
  double precision               :: cpu_1, cpu_2, wall_1, wall_2


  ao_two_e_integrals$_erf_in_map = .True.
  if (read_ao_two_e_integrals$_erf) then
    print*,'Reading the AO integrals$_erf'
    call map_load_from_disk(trim(ezfio_filename)//'/work/ao_ints$_erf',ao_integrals$_erf_map)
    print*, 'AO integrals$_erf provided'
    return
  endif

  call wall_time(wall_1)
  call cpu_time(cpu_1)

  if (do_ao_cholesky) then
     PROVIDE cholesky_ao_transp
  else
    call add_integrals$_erf_to_map_cholesky_ao
    integer*8                      :: get_ao_map_size, ao_map_size
    ao_map_size = get_ao_map_size()
    double precision, external     :: map_mb
    print*,'Atomic integrals$_erf provided:'
    print*,' Size of AO map           ', map_mb(ao_integrals$_erf_map) ,'MB'
    print*,' Number of AO integrals$_erf: ',  ao_map_size
  endif

  call wall_time(wall_2)
  call cpu_time(cpu_2)

  print*,' cpu  time :',cpu_2 - cpu_1, 's'
  print*,' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1), ')'

  if (write_ao_two_e_integrals$_erf.and.mpi_master) then
    call ezfio_set_work_empty(.False.)
    call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints$_erf',ao_integrals$_erf_map)
    call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals$_erf('Read')
  endif

END_PROVIDER

subroutine add_integrals$_erf_to_map_cholesky_ao
  use bitmasks
  implicit none

  BEGIN_DOC
  ! Adds integrals$_erf to the AO map using Cholesky vectors
  END_DOC

  integer :: i,j,k,l,m
  integer :: size_buffer, n_integrals$_erf
  size_buffer = min(ao_num*ao_num*ao_num,16000000)

  double precision, allocatable :: Vtmp(:,:,:)
  integer(key_kind)  , allocatable :: buffer_i(:)
  real(integral_kind), allocatable :: buffer_value(:)

  PROVIDE cholesky_ao_transp
  call set_multiple_levels_omp(.False.)

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(i,j,k,l,n_integrals$_erf,buffer_value, buffer_i, Vtmp)
  allocate (buffer_i(size_buffer), buffer_value(size_buffer))
  allocate (Vtmp(ao_num,ao_num,ao_num))
  n_integrals$_erf = 0

  !$OMP DO SCHEDULE(dynamic)
  do l=1,ao_num
    call dgemm('T','N',ao_num*ao_num,ao_num,cholesky_ao_num,1.d0, &
       cholesky_ao_transp, cholesky_ao_num, &
       cholesky_ao_transp(1,1,l), cholesky_ao_num, 0.d0, &
       Vtmp, ao_num*ao_num)

    do k=1,l
      do j=1,ao_num
        do i=1,j
          if (dabs(Vtmp(i,j,k)) > ao_integrals_threshold) then
            n_integrals$_erf = n_integrals$_erf + 1
            buffer_value(n_integrals$_erf) = Vtmp(i,j,k)
            !DIR$ FORCEINLINE
            call ao_two_e_integrals$_erf_index(i,k,j,l,buffer_i(n_integrals$_erf))
            if (n_integrals$_erf == size_buffer) then
              call map_append(ao_integrals$_erf_map, buffer_i, buffer_value, n_integrals$_erf)
              n_integrals$_erf = 0
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  if (n_integrals$_erf > 0) then
    call map_append(ao_integrals$_erf_map, buffer_i, buffer_value, n_integrals$_erf)
  endif
  deallocate(buffer_i, buffer_value, Vtmp)
  !$OMP BARRIER
  !$OMP END PARALLEL

  call map_sort(ao_integrals$_erf_map)
  call map_unique(ao_integrals$_erf_map)

end

subroutine clear_ao_map
  implicit none
  BEGIN_DOC
  ! Frees the meaory of the AO map
  END_DOC
  call map_deinit(ao_integrals$_erf_map)
  FREE ao_integrals$_erf_map 
  FREE ao_two_e_integrals$_erf_in_map
end

SUBST [ _erf ]
  
;;
_erf;;
_cgtos;;
  
END_TEMPLATE     
