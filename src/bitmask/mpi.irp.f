BEGIN_PROVIDER [ integer, mpi_bit_kind ]
 use bitmasks
 implicit none
 BEGIN_DOC
 ! MPI bit kind type
 END_DOC
 IRP_IF MPI
  include 'mpif.h'
  if (bit_kind == 4) then
    mpi_bit_kind = MPI_INTEGER4
  else if (bit_kind == 8) then
    mpi_bit_kind = MPI_INTEGER8
  else
    stop 'Wrong bit kind in mpi_bit_kind'
  endif
 IRP_ELSE
  mpi_bit_kind = -1
 IRP_ENDIF
END_PROVIDER

subroutine broadcast_chunks_bit_kind(A, LDA)
  use bitmasks
  implicit none
  integer*8, intent(in)             :: LDA
  integer(bit_kind), intent(inout) :: A(LDA)
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
    do i=1,LDA,200000000/bit_kind_size
      sze = min(LDA-i+1, 200000000/bit_kind_size)
      call MPI_BCAST (A(i), sze, MPI_BIT_KIND, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        print *,  irp_here//': Unable to broadcast chunks bit_kind', i
        stop -1
      endif
    enddo
  IRP_ENDIF
end



