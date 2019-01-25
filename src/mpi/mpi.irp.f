BEGIN_PROVIDER [ logical, mpi_initialized ]
 implicit none
 BEGIN_DOC
 ! Always true. Initialized MPI
 END_DOC
 IRP_IF MPI
  include 'mpif.h'
  integer                        :: ierr
  call mpi_init(ierr)
  if (ierr /= MPI_SUCCESS) then
    print *,  'ierr = ', ierr
    stop 'Unable to initialize MPI'
  endif
 IRP_ENDIF
 mpi_initialized = .True.
END_PROVIDER


 BEGIN_PROVIDER [ integer, mpi_rank ]
&BEGIN_PROVIDER [ integer, mpi_size ]
 implicit none
 BEGIN_DOC
 ! Rank of MPI process and number of MPI processes
 END_DOC
 IRP_IF MPI
  include 'mpif.h'
  PROVIDE mpi_initialized
  integer                        :: ierr

  call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_rank, ierr)
  if (ierr /= MPI_SUCCESS) then
    print *,  'ierr = ', ierr
    stop 'Unable to get MPI rank'
  endif

  call MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_size, ierr)
  if (ierr /= MPI_SUCCESS) then
    print *,  'ierr = ', ierr
    stop 'Unable to get MPI size'
  endif

 IRP_ELSE
  mpi_rank = 0
  mpi_size = 1
 IRP_ENDIF
 ASSERT (mpi_rank >= 0)
 ASSERT (mpi_rank < mpi_size)

END_PROVIDER


BEGIN_PROVIDER [ logical, mpi_master ]
 implicit none
 BEGIN_DOC
 ! If true, rank is zero
 END_DOC
 mpi_master = (mpi_rank == 0)
 if (mpi_master.and.(mpi_size > 1)) then
    print *,  'MPI size: ', mpi_size
 endif

END_PROVIDER

BEGIN_TEMPLATE

subroutine broadcast_chunks_$double(A, LDA)
  implicit none
  integer*8, intent(in)             :: LDA
  $type, intent(inout) :: A(LDA)
  BEGIN_DOC
! Broadcast with chunks of ~2GB
  END_DOC
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: i, sze, ierr
    do i=1,LDA,200000000/$8
      sze = min(LDA-i+1, 200000000/$8)
      call MPI_BCAST (A(i), sze, MPI_$DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        print *,  irp_here//': Unable to broadcast chunks $double ', i
        stop -1
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    enddo
  IRP_ENDIF
end

SUBST [ double, type, 8, DOUBLE_PRECISION ]
double    ; double precision  ; 8 ; DOUBLE_PRECISION ;;
integer   ; integer           ; 4  ; INTEGER ;;
integer8  ; integer*8         ; 8  ; INTEGER8 ;;

END_TEMPLATE


subroutine mpi_print(string)
  implicit none
  BEGIN_DOC
! Print string to stdout if the MPI rank is zero.
  END_DOC
  character*(*)                  :: string
  if (mpi_master) then
    print *, string
  endif
end
