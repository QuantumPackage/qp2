===
mpi
===

Single-node placeholders that used to wrap |MPI|.

Multi-node parallelism (|MPI| and the distributed task server) was removed
in 2026-11; the code is now single-node, with in-process OpenMP
parallelism inside each executable. The providers
:c:data:`mpi_rank`, :c:data:`mpi_size` and :c:data:`mpi_master`
survive as constants (``0``, ``1`` and ``.True.``) because a large part
of the code base guards master-only I/O and file writes with
``mpi_master``. :c:data:`mpi_initialized` is always ``.True.``.
