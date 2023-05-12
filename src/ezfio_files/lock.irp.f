use omp_lib

BEGIN_PROVIDER [ integer(omp_lock_kind), file_lock ]
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! OpenMP Lock for I/O
  END_DOC
  call omp_init_lock(file_lock)
END_PROVIDER

! These functions need to be called because internal read and write are not thread safe.
subroutine lock_io()
  implicit none
  BEGIN_DOC
! Needs to be called because before doing I/O because internal read and write
! are not thread safe.
  END_DOC
  call omp_set_lock(file_lock)
end subroutine lock_io()

subroutine unlock_io()
  implicit none
  BEGIN_DOC
! Needs to be called because afterdoing I/O because internal read and write
! are not thread safe.
  END_DOC
  call omp_unset_lock(file_lock)
end subroutine lock_io()
