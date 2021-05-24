BEGIN_PROVIDER [ double precision, shell_normalization_factor , (shell_num) ]
  implicit none
  BEGIN_DOC
  ! Number of primitives per |AO|
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    if (size(shell_normalization_factor) == 0) return

    call ezfio_has_basis_shell_normalization_factor(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: shell_normalization_factor ] <<<<< ..'
      call ezfio_get_basis_shell_normalization_factor(shell_normalization_factor)
    else


      double precision               :: norm,overlap_x,overlap_y,overlap_z,C_A(3), c
      integer                        :: l, powA(3), nz
      integer                        :: i,j,k
      nz=100
      C_A(1) = 0.d0
      C_A(2) = 0.d0
      C_A(3) = 0.d0

      do i=1,shell_num
        powA(1) = shell_ang_mom(i)
        powA(2) = 0
        powA(3) = 0

        ! Normalization of the contracted basis functions
        norm = 0.d0
        do j=shell_prim_index(i), shell_prim_index(i)+shell_prim_num(i)-1
          do k=shell_prim_index(i), shell_prim_index(i)+shell_prim_num(i)-1
            call overlap_gaussian_xyz(C_A,C_A,shell_prim_expo(j),shell_prim_expo(k),&
                powA,powA,overlap_x,overlap_y,overlap_z,c,nz)
            norm = norm+c*shell_prim_coef(j)*shell_prim_coef(k)
          enddo
        enddo
        shell_normalization_factor(i) = dsqrt(norm)
      enddo


    endif
  endif
  IRP_IF MPI_DEBUG
  print *,  irp_here, mpi_rank
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
  include 'mpif.h'
  integer                        :: ierr
  call MPI_BCAST( shell_normalization_factor, (shell_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) then
    stop 'Unable to read shell_normalization_factor with MPI'
  endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
