use omp_lib

BEGIN_PROVIDER [ integer(omp_lock_kind), file_lock ]
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! OpenMP Lock for I/O
  END_DOC
  call omp_init_lock(file_lock)
END_PROVIDER


