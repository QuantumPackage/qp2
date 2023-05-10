
! ---

BEGIN_PROVIDER [double precision, mo_coef_aux, (ao_num,mo_num)]

  implicit none
  integer                        :: i, j
  logical                        :: exists
  double precision, allocatable  :: buffer(:,:)

  PROVIDE ezfio_filename

  if (mpi_master) then
    ! Coefs
    call ezfio_has_mo_basis_mo_coef_aux(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST(exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_coef_aux with MPI'
    endif
  IRP_ENDIF

  if (exists) then
    if (mpi_master) then
      call ezfio_get_mo_basis_mo_coef_aux(mo_coef_aux)
      write(*,*) 'Read  mo_coef_aux'
    endif
    IRP_IF MPI
      call MPI_BCAST(mo_coef_aux, mo_num*ao_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read mo_coef_aux with MPI'
      endif
    IRP_ENDIF
  else
    ! Orthonormalized AO basis
    do i = 1, mo_num
      do j = 1, ao_num
        mo_coef_aux(j,i) = ao_ortho_canonical_coef(j,i)
      enddo
    enddo
  endif

END_PROVIDER

