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

        norm = 0.d0
        do k=shell_prim_index(i),shell_prim_index(i)+shell_prim_num(i)-1
          do j=shell_prim_index(i),shell_prim_index(i)+shell_prim_num(i)-1
            call overlap_gaussian_xyz(C_A,C_A,prim_expo(j),prim_expo(k), &
              powA,powA,overlap_x,overlap_y,overlap_z,c,nz)
            norm = norm+c*prim_coef(j)*prim_coef(k)
          enddo
        enddo
        shell_normalization_factor(i) = 1.d0/dsqrt(norm)
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


BEGIN_PROVIDER [ double precision, prim_normalization_factor , (prim_num) ]
  implicit none
  BEGIN_DOC
  ! Number of primitives per |AO|
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    if (size(prim_normalization_factor) == 0) return

    call ezfio_has_basis_prim_normalization_factor(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: prim_normalization_factor ] <<<<< ..'
      call ezfio_get_basis_prim_normalization_factor(prim_normalization_factor)
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

        do k=shell_prim_index(i),shell_prim_index(i)+shell_prim_num(i)-1
            call overlap_gaussian_xyz(C_A,C_A,prim_expo(k),prim_expo(k), &
              powA,powA,overlap_x,overlap_y,overlap_z,norm,nz)
            prim_normalization_factor(k) = 1.d0/dsqrt(norm)
        enddo
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
  call MPI_BCAST( prim_normalization_factor, (prim_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) then
    stop 'Unable to read prim_normalization_factor with MPI'
  endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
