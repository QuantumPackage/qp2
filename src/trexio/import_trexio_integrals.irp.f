program import_integrals_ao
  use trexio
  implicit none
  integer(trexio_t)              :: f ! TREXIO file handle
  integer(trexio_exit_code)      :: rc
  PROVIDE mo_num

  f = trexio_open(trexio_filename, 'r', TREXIO_AUTO, rc)
  if (f == 0_8) then
    print *, 'Unable to open TREXIO file for reading'
    print *, 'rc = ', rc
    stop -1
  endif

  call run(f)
  rc = trexio_close(f)
  call trexio_assert(rc, TREXIO_SUCCESS)
end

subroutine run(f)
  use trexio
  use map_module
  implicit none
  BEGIN_DOC
! Program to import integrals from TREXIO
  END_DOC

  integer(trexio_t), intent(in)  :: f ! TREXIO file handle
  integer(trexio_exit_code)      :: rc

  integer ::i,j,k,l
  integer(8) :: m, n_integrals
  double precision :: integral

  integer(key_kind), allocatable   :: buffer_i(:)
  real(integral_kind), allocatable :: buffer_values(:)


  double precision, allocatable :: A(:,:)
  double precision, allocatable :: V(:)
  integer         , allocatable :: Vi(:,:)
  double precision              :: s

  if (trexio_has_nucleus_repulsion(f) == TREXIO_SUCCESS) then
    rc = trexio_read_nucleus_repulsion(f, s)
    if (rc /= TREXIO_SUCCESS) then
      print *, irp_here, rc
      print *, 'Error reading nuclear repulsion'
      call trexio_assert(rc, TREXIO_SUCCESS)
      stop -1
    endif
    call ezfio_set_nuclei_nuclear_repulsion(s)
    call ezfio_set_nuclei_io_nuclear_repulsion('Read')
  endif

  ! AO integrals
  ! ------------

  allocate(A(ao_num, ao_num))


  if (trexio_has_ao_1e_int_overlap(f) == TREXIO_SUCCESS) then
    rc = trexio_read_ao_1e_int_overlap(f, A)
    if (rc /= TREXIO_SUCCESS) then
      print *, irp_here
      print *, 'Error reading AO overlap'
      call trexio_assert(rc, TREXIO_SUCCESS)
      stop -1
    endif
    call ezfio_set_ao_one_e_ints_ao_integrals_overlap(A)
    call ezfio_set_ao_one_e_ints_io_ao_integrals_overlap('Read')
  endif

  if (trexio_has_ao_1e_int_kinetic(f) == TREXIO_SUCCESS) then
    rc = trexio_read_ao_1e_int_kinetic(f, A)
    if (rc /= TREXIO_SUCCESS) then
      print *, irp_here
      print *, 'Error reading AO kinetic integrals'
      call trexio_assert(rc, TREXIO_SUCCESS)
      stop -1
    endif
    call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(A)
    call ezfio_set_ao_one_e_ints_io_ao_integrals_kinetic('Read')
  endif

!  if (trexio_has_ao_1e_int_ecp(f) == TREXIO_SUCCESS) then
!    rc = trexio_read_ao_1e_int_ecp(f, A)
!    if (rc /= TREXIO_SUCCESS) then
!      print *, irp_here
!      print *, 'Error reading AO ECP local integrals'
!      call trexio_assert(rc, TREXIO_SUCCESS)
!      stop -1
!    endif
!    call ezfio_set_ao_one_e_ints_ao_integrals_pseudo(A)
!    call ezfio_set_ao_one_e_ints_io_ao_integrals_pseudo('Read')
!  endif

  if (trexio_has_ao_1e_int_potential_n_e(f) == TREXIO_SUCCESS) then
    rc = trexio_read_ao_1e_int_potential_n_e(f, A)
    if (rc /= TREXIO_SUCCESS) then
      print *, irp_here
      print *, 'Error reading AO potential N-e integrals'
      call trexio_assert(rc, TREXIO_SUCCESS)
      stop -1
    endif
    call ezfio_set_ao_one_e_ints_ao_integrals_n_e(A)
    call ezfio_set_ao_one_e_ints_io_ao_integrals_n_e('Read')
  endif

  deallocate(A)

  ! AO 2e integrals
  ! ---------------

  rc = trexio_has_ao_2e_int(f)
  PROVIDE ao_num
  if (rc /= TREXIO_HAS_NOT) then
      PROVIDE ao_integrals_map

      integer*4 :: BUFSIZE
      BUFSIZE=ao_num**2
      allocate(buffer_i(BUFSIZE), buffer_values(BUFSIZE))
      allocate(Vi(4,BUFSIZE), V(BUFSIZE))

      integer*8 :: offset, icount

      offset = 0_8
      icount = BUFSIZE
      rc = TREXIO_SUCCESS
      do while (icount == size(V))
        rc = trexio_read_ao_2e_int_eri(f, offset, icount, Vi, V)
        do m=1,icount
          i = Vi(1,m)
          j = Vi(2,m)
          k = Vi(3,m)
          l = Vi(4,m)
          integral = V(m)
          call two_e_integrals_index(i, j, k, l, buffer_i(m) )
          buffer_values(m) = integral
        enddo
        call insert_into_ao_integrals_map(int(icount,4),buffer_i,buffer_values)
        offset = offset + icount
        if (rc /= TREXIO_SUCCESS) then
            exit
        endif
      end do
      n_integrals = offset

      call map_sort(ao_integrals_map)
      call map_unique(ao_integrals_map)

      call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
      call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read')

      deallocate(buffer_i, buffer_values, Vi, V)
      print *, 'AO integrals read from TREXIO file'
  else
      print *, 'AO integrals not found in TREXIO file'
  endif

  ! MO integrals
  ! ------------

  allocate(A(mo_num, mo_num))
  if (trexio_has_mo_1e_int_core_hamiltonian(f) == TREXIO_SUCCESS) then
    rc = trexio_read_mo_1e_int_core_hamiltonian(f, A)
    if (rc /= TREXIO_SUCCESS) then
      print *, irp_here
      print *, 'Error reading MO 1e integrals'
      call trexio_assert(rc, TREXIO_SUCCESS)
      stop -1
    endif
    call ezfio_set_mo_one_e_ints_mo_one_e_integrals(A)
    call ezfio_set_mo_one_e_ints_io_mo_one_e_integrals('Read')
  endif
  deallocate(A)

  ! MO 2e integrals
  ! ---------------

  rc = trexio_has_mo_2e_int(f)
  if (rc /= TREXIO_HAS_NOT) then

      BUFSIZE=mo_num**2
      allocate(buffer_i(BUFSIZE), buffer_values(BUFSIZE))
      allocate(Vi(4,BUFSIZE), V(BUFSIZE))


      offset = 0_8
      icount = BUFSIZE
      rc = TREXIO_SUCCESS
      do while (icount == size(V))
        rc = trexio_read_mo_2e_int_eri(f, offset, icount, Vi, V)
        do m=1,icount
          i = Vi(1,m)
          j = Vi(2,m)
          k = Vi(3,m)
          l = Vi(4,m)
          integral = V(m)
          call two_e_integrals_index(i, j, k, l, buffer_i(m) )
          buffer_values(m) = integral
        enddo
        call map_append(mo_integrals_map, buffer_i, buffer_values, int(icount,4))
        offset = offset + icount
        if (rc /= TREXIO_SUCCESS) then
            exit
        endif
      end do
      n_integrals = offset

      call map_sort(mo_integrals_map)
      call map_unique(mo_integrals_map)

      call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_map)
      call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
      deallocate(buffer_i, buffer_values, Vi, V)
      print *, 'MO integrals read from TREXIO file'
  else
      print *, 'MO integrals not found in TREXIO file'
  endif

end
