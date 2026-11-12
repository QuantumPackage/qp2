BEGIN_PROVIDER [ logical, mpi_initialized ]
 implicit none
 BEGIN_DOC
 ! Always true. Initialized MPI
 END_DOC
 mpi_initialized = .True.
END_PROVIDER


 BEGIN_PROVIDER [ integer, mpi_rank ]
&BEGIN_PROVIDER [ integer, mpi_size ]
 implicit none
 BEGIN_DOC
 ! Rank of MPI process and number of MPI processes
 END_DOC
  mpi_rank = 0
  mpi_size = 1
 ASSERT (mpi_rank >= 0)
 ASSERT (mpi_rank < mpi_size)

END_PROVIDER


BEGIN_PROVIDER [ logical, mpi_master ]
 implicit none
 BEGIN_DOC
 ! If true, rank is zero
 END_DOC
 mpi_master = (mpi_rank == 0)
END_PROVIDER
