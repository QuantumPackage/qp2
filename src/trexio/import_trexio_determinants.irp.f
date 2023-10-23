program import_determinants_ao
  call run
end

subroutine run
  use trexio
  use map_module
  implicit none
  BEGIN_DOC
! Program to import determinants from TREXIO
  END_DOC

  integer(trexio_t)              :: f ! TREXIO file handle
  integer(trexio_exit_code)      :: rc

  integer :: m

  double precision, allocatable :: coef_buffer(:,:)
  integer*8       , allocatable :: det_buffer(:,:,:)

  f = trexio_open(trexio_filename, 'r', TREXIO_AUTO, rc)
  if (f == 0_8) then
    print *, 'Unable to open TREXIO file for reading'
    print *, 'rc = ', rc
    stop -1
  endif



  ! Determinants
  ! ------------

  integer :: nint, nstates
  integer :: bufsize

  rc = trexio_read_state_num(f, nstates)
  call trexio_assert(rc, TREXIO_SUCCESS)

!  rc = trexio_read_determinant_int64_num(f, nint)
!  call trexio_assert(rc, TREXIO_SUCCESS)
  nint = N_int
  if (nint /= N_int) then
     stop 'Problem with N_int'
  endif

  integer*8 :: offset, icount

  rc = trexio_read_determinant_num(f, bufsize)
  call trexio_assert(rc, TREXIO_SUCCESS)
  print *, 'N_det = ', bufsize

  allocate ( det_buffer(nint, 2, bufsize), coef_buffer(bufsize, n_states) )


  offset = 0_8
  icount = bufsize

  rc = trexio_read_determinant_list(f, offset, icount, det_buffer)
  call trexio_assert(rc, TREXIO_SUCCESS)
  if (icount /= bufsize) then
      print *, 'error: bufsize /= N_det: ', bufsize, icount
      stop -1
  endif

  do m=1,nstates
    rc = trexio_set_state(f, m-1)
    call trexio_assert(rc, TREXIO_SUCCESS)
    rc = trexio_read_determinant_coefficient(f, offset, icount, coef_buffer(1,m))
    call trexio_assert(rc, TREXIO_SUCCESS)
    if (icount /= bufsize) then
        print *, 'error: bufsize /= N_det for state', m, ':', icount, bufsize
        stop -1
    endif
  enddo

  call save_wavefunction_general(bufsize,nstates,det_buffer,size(coef_buffer,1),coef_buffer)


end
