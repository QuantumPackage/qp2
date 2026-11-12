
! ---

 BEGIN_PROVIDER [double precision, j1e_expo, (j1e_size, nucl_num)]
&BEGIN_PROVIDER [double precision, j1e_coef, (j1e_size, nucl_num)]

  BEGIN_DOC
  !
  ! parameters of the 1e-Jastrow
  !
  END_DOC

  implicit none
  logical :: exists
  integer :: i, j
  integer :: ierr

  PROVIDE ezfio_filename

  ! ---

  if (mpi_master) then
    call ezfio_has_jastrow_j1e_expo(exists)
  endif



  if (exists) then
    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: j1e_expo ] <<<<< ..'
      call ezfio_get_jastrow_j1e_expo(j1e_expo)
    endif
  else

    j1e_expo = 1.d0
    call ezfio_set_jastrow_j1e_expo(j1e_expo)
  endif

  ! ---

  if (mpi_master) then
    call ezfio_has_jastrow_j1e_coef(exists)
  endif



  if (exists) then
    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: j1e_coef ] <<<<< ..'
      call ezfio_get_jastrow_j1e_coef(j1e_coef)
    endif
  else

    j1e_coef = 0.d0
    call ezfio_set_jastrow_j1e_coef(j1e_coef)
  endif

  ! ---

  print *, ' parameters of the 1e-Jastrow'
  do i = 1, nucl_num
    print*, ' for Z = ', nucl_charge(i)
    do j = 1, j1e_size
      write(*,'(I4, 2x, 2(E15.7, 2X))') j, j1e_coef(j,i), j1e_expo(j,i)
    enddo
  enddo

END_PROVIDER

! ---

