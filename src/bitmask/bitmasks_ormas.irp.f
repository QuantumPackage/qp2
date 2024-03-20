use bitmasks

BEGIN_PROVIDER [integer, ormas_mstart, (ormas_n_space) ]
  implicit none
  BEGIN_DOC
! first orbital idx in each active space
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_bitmask_ormas_mstart(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: ormas_mstart ] <<<<< ..'
      call ezfio_get_bitmask_ormas_mstart(ormas_mstart)
      ASSERT (ormas_mstart(1).eq.1)
    else if (ormas_n_space.eq.1) then
      ormas_mstart = 1
    else
      print *, 'bitmask/ormas_mstart not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( ormas_mstart, ormas_n_space, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read ormas_mstart with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)
  

END_PROVIDER

BEGIN_PROVIDER [integer, ormas_min_e, (ormas_n_space) ]
  implicit none
  BEGIN_DOC
! min nelec in each active space
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_bitmask_ormas_min_e(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: ormas_min_e ] <<<<< ..'
      call ezfio_get_bitmask_ormas_min_e(ormas_min_e)
    else if (ormas_n_space.eq.1) then
      ormas_min_e = 0
    else
      print *, 'bitmask/ormas_min_e not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( ormas_min_e, ormas_n_space, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read ormas_min_e with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [integer, ormas_max_e, (ormas_n_space) ]
  implicit none
  BEGIN_DOC
! max nelec in each active space
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_bitmask_ormas_max_e(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: ormas_max_e ] <<<<< ..'
      call ezfio_get_bitmask_ormas_max_e(ormas_max_e)
    else if (ormas_n_space.eq.1) then
      ormas_max_e = elec_num
    else
      print *, 'bitmask/ormas_max_e not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( ormas_max_e, ormas_n_space, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read ormas_max_e with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)

END_PROVIDER

 BEGIN_PROVIDER [ integer, ormas_n_orb, (ormas_n_space) ]
&BEGIN_PROVIDER [ integer, ormas_max_n_orb ]
  implicit none
  BEGIN_DOC
  ! number of orbitals in each ormas space
  END_DOC
  integer :: i
  ormas_n_orb = 0
  ormas_n_orb(ormas_n_space) = mo_num + 1 - ormas_mstart(ormas_n_space) 
  do i = ormas_n_space-1, 1, -1
    ormas_n_orb(i) = ormas_mstart(i+1) - ormas_mstart(i)
    ASSERT (ormas_n_orb(i).ge.1)
  enddo
  ormas_max_n_orb = maxval(ormas_n_orb)
END_PROVIDER

BEGIN_PROVIDER [ integer, ormas_list_orb, (ormas_max_n_orb, ormas_n_space) ]
  implicit none
  BEGIN_DOC
  ! list of orbitals in each ormas space
  END_DOC
  integer :: i,j,k
  ormas_list_orb = 0
  i = 1
  do j = 1, ormas_n_space
    do k = 1, ormas_n_orb(j)
      ormas_list_orb(k,j) = i
      i += 1
    enddo
  enddo
END_PROVIDER
  
BEGIN_PROVIDER [ integer(bit_kind), ormas_bitmask, (N_int, ormas_n_space) ]
  implicit none
  BEGIN_DOC
  ! bitmask for each ormas space
  END_DOC
  integer :: j
  ormas_bitmask = 0_bit_kind
  do j = 1, ormas_n_space
    call list_to_bitstring(ormas_bitmask(1,j), ormas_list_orb(:,j), ormas_n_orb(j), N_int)
  enddo
END_PROVIDER

subroutine ormas_occ(key_in, occupancies)
  implicit none
  BEGIN_DOC
  ! number of electrons in each ormas space
  END_DOC
  integer(bit_kind), intent(in) :: key_in(N_int,2)
  integer, intent(out) :: occupancies(ormas_n_space)
  integer :: i,ispin,ispace

  occupancies = 0
  ! TODO: get start/end of each space within N_int
  do ispace=1,ormas_n_space
    do ispin=1,2
      do i=1,N_int
        occupancies(ispace) += popcnt(iand(ormas_bitmask(i,ispace),key_in(i,ispin)))
      enddo
    enddo
  enddo
end

logical function det_allowed_ormas(key_in)
  implicit none
  BEGIN_DOC
  ! return true if det has allowable ormas occupations
  END_DOC
  integer(bit_kind), intent(in) :: key_in(N_int,2)
  integer :: i,ispin,ispace,occ

  det_allowed_ormas = .True.
  if (ormas_n_space.eq.1) return
  det_allowed_ormas = .False.
  ! TODO: get start/end of each space within N_int
  do ispace=1,ormas_n_space
    occ = 0
    do ispin=1,2
      do i=1,N_int
        occ += popcnt(iand(ormas_bitmask(i,ispace),key_in(i,ispin)))
      enddo
    enddo
    if ((occ.lt.ormas_min_e(ispace)).or.(occ.gt.ormas_max_e(ispace))) return
  enddo
  det_allowed_ormas = .True.
end
  
