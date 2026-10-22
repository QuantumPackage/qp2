===
mpi
===

Initialization and basic routines for the parallelization with |MPI|.

The provider :c:data:`mpi_initialized` starts |MPI|. The providers
:c:data:`mpi_rank`, :c:data:`mpi_size` and :c:data:`mpi_master` give
the rank of the process, the number of processes, and the boolean that
is |true| only on the master rank.

:c:func:`mpi_print` prints a string only on the master rank. The
subroutines :c:func:`broadcast_chunks_double`,
:c:func:`broadcast_chunks_integer` and
:c:func:`broadcast_chunks_integer8` broadcast an array from the master
rank in chunks of about 2 GB.
