BEGIN_PROVIDER [ integer, n_iter  ]
  implicit none
  BEGIN_DOC
! number of iterations
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then

      double precision :: zeros(N_states,100)
      integer :: izeros(100)
      zeros = 0.d0
      izeros = 0
      call ezfio_set_iterations_n_iter(0)
      call ezfio_set_iterations_energy_iterations(zeros)
      call ezfio_set_iterations_pt2_iterations(zeros)
      call ezfio_set_iterations_n_det_iterations(izeros)
      n_iter = 1
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( n_iter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read n_iter with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

