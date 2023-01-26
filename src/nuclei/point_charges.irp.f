! ---


BEGIN_PROVIDER [ integer, n_pts_charge  ]
  implicit none
  BEGIN_DOC
! Number of point charges to be added to the potential
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_nuclei_n_pts_charge(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: n_pts_charge ] <<<<< ..'
      call ezfio_get_nuclei_n_pts_charge(n_pts_charge)
    else
      print *, 'nuclei/n_pts_charge not found in EZFIO file'
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
    call MPI_BCAST( n_pts_charge, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read n_pts_charge with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ double precision, pts_charge_z, (n_pts_charge) ]

  BEGIN_DOC
  ! Charge associated to each point charge. 
  END_DOC

  implicit none
  logical :: exists

  PROVIDE ezfio_filename

  if (mpi_master) then
    call ezfio_has_nuclei_pts_charge_z(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST(pts_charge_z, (n_pts_charge), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read pts_charge_z with MPI'
    endif
  IRP_ENDIF

  if (exists) then

    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: pts_charge_z ] <<<<< ..'
      call ezfio_get_nuclei_pts_charge_z(pts_charge_z)
      IRP_IF MPI
        call MPI_BCAST(pts_charge_z, (n_pts_charge), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read pts_charge_z with MPI'
        endif
      IRP_ENDIF
    endif

  else
 
    integer :: i
    do i = 1, n_pts_charge
      pts_charge_z(i) = 0.d0
    enddo

  endif
 print*,'Point charges '
 do i = 1, n_pts_charge
  print*,'i,pts_charge_z(i)',i,pts_charge_z(i)
 enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, pts_charge_coord, (n_pts_charge,3) ]

  BEGIN_DOC
  ! Coordinates of each point charge. 
  END_DOC

  implicit none
  logical :: exists

  PROVIDE ezfio_filename

  if (mpi_master) then
    call ezfio_has_nuclei_pts_charge_coord(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST(pts_charge_coord, (n_pts_charge), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read pts_charge_coord with MPI'
    endif
  IRP_ENDIF

  if (exists) then

    if (mpi_master) then
     double precision, allocatable  :: buffer(:,:)
     allocate (buffer(n_pts_charge,3))
      write(6,'(A)') '.. >>>>> [ IO READ: pts_charge_coord ] <<<<< ..'
      call ezfio_get_nuclei_pts_charge_coord(buffer)
      integer :: i,j
      do i=1,3
        do j=1,n_pts_charge
          pts_charge_coord(j,i) = buffer(j,i)
        enddo
      enddo
     deallocate(buffer)
      IRP_IF MPI
        call MPI_BCAST(pts_charge_coord, (n_pts_charge), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read pts_charge_coord with MPI'
        endif
      IRP_ENDIF
    endif

  else
 
    do i = 1, n_pts_charge
      pts_charge_coord(i,:) = 0.d0
    enddo

  endif
 print*,'Coordinates for the point charges '
 do i = 1, n_pts_charge
  write(*,'(I3,X,3(F16.8,X))') i,pts_charge_coord(i,1:3)
 enddo

END_PROVIDER

! ---
BEGIN_PROVIDER [ double precision, pt_chrg_repulsion]
 implicit none
 BEGIN_DOC
 ! repulsion between the point charges 
 END_DOC
 integer :: i,j
 double precision               :: Z_A, z_B,A_center(3), B_center(3), dist
 pt_chrg_repulsion = 0.d0
 do  i = 1, n_pts_charge
   Z_A = pts_charge_z(i)
   A_center(1:3) = pts_charge_coord(i,1:3)
   do  j = i+1, n_pts_charge
     Z_B = pts_charge_z(j)
     B_center(1:3) = pts_charge_coord(j,1:3)
     dist = (A_center(1)-B_center(1))**2 + (A_center(2)-B_center(2))**2 + (A_center(3)-B_center(3))**2
     dist = dsqrt(dist)
     pt_chrg_repulsion += Z_A*Z_B/dist
   enddo
 enddo
 print*,'Repulsion of point charges '
 print*,'pt_chrg_repulsion = ',pt_chrg_repulsion
END_PROVIDER 

BEGIN_PROVIDER [ double precision, pt_chrg_nuclei_repulsion]
 implicit none
 BEGIN_DOC
 ! repulsion between the point charges and the nuclei
 END_DOC
 integer :: i,j
 double precision               :: Z_A, z_B,A_center(3), B_center(3), dist
 pt_chrg_nuclei_repulsion = 0.d0
 do  i = 1, n_pts_charge
   Z_A = pts_charge_z(i)
   A_center(1:3) = pts_charge_coord(i,1:3)
   do  j = 1, nucl_num
     Z_B = nucl_charge(j)
     B_center(1:3) = nucl_coord(j,1:3)
     dist = (A_center(1)-B_center(1))**2 + (A_center(2)-B_center(2))**2 + (A_center(3)-B_center(3))**2
     dist = dsqrt(dist)
     pt_chrg_nuclei_repulsion += Z_A*Z_B/dist
   enddo
 enddo
 print*,'Repulsion between point charges and nuclei'
 print*,'pt_chrg_nuclei_repulsion = ',pt_chrg_nuclei_repulsion
END_PROVIDER 

