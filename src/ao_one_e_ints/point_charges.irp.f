
! ---


BEGIN_PROVIDER [ integer, n_pts_charge  ]
  implicit none
  BEGIN_DOC
! Number of point charges to be added to the potential
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_ao_one_e_ints_n_pts_charge(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: n_pts_charge ] <<<<< ..'
      call ezfio_get_ao_one_e_ints_n_pts_charge(n_pts_charge)
    else
      print *, 'ao_one_e_ints/n_pts_charge not found in EZFIO file'
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
    call ezfio_has_ao_one_e_ints_pts_charge_z(exists)
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
      call ezfio_get_ao_one_e_ints_pts_charge_z(pts_charge_z)
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
    call ezfio_has_ao_one_e_ints_pts_charge_coord(exists)
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
      call ezfio_get_ao_one_e_ints_pts_charge_coord(buffer)
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

BEGIN_PROVIDER [ double precision, ao_integrals_pt_chrg, (ao_num,ao_num)]

  BEGIN_DOC
  !  Point charge-electron interaction, in the |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
  !
  !  These integrals also contain the pseudopotential integrals.
  END_DOC

  implicit none
  integer          :: num_A, num_B, power_A(3), power_B(3)
  integer          :: i, j, k, l, n_pt_in, m
  double precision :: alpha, beta
  double precision :: A_center(3),B_center(3),C_center(3)
  double precision :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult

  ao_integrals_pt_chrg = 0.d0

!  if (read_ao_integrals_pt_chrg) then
!
!    call ezfio_get_ao_one_e_ints_ao_integrals_pt_chrg(ao_integrals_pt_chrg)
!    print *,  'AO N-e integrals read from disk'
!
!  else

!    if(use_cosgtos) then
!      !print *, " use_cosgtos for ao_integrals_pt_chrg ?", use_cosgtos
!
!      do j = 1, ao_num
!        do i = 1, ao_num
!          ao_integrals_pt_chrg(i,j) = ao_integrals_pt_chrg_cosgtos(i,j)
!        enddo
!      enddo
!
!    else

      !$OMP PARALLEL                                                   &
          !$OMP DEFAULT (NONE)                                         &
          !$OMP PRIVATE (i,j,k,l,m,alpha,beta,A_center,B_center,C_center,power_A,power_B,&
          !$OMP          num_A,num_B,Z,c,c1,n_pt_in)                      &
          !$OMP SHARED (ao_num,ao_prim_num,ao_expo_ordered_transp,ao_power,ao_nucl,pts_charge_coord,ao_coef_normalized_ordered_transp,nucl_coord,&
          !$OMP         n_pt_max_integrals,ao_integrals_pt_chrg,n_pts_charge,pts_charge_z)

      n_pt_in = n_pt_max_integrals

      !$OMP DO SCHEDULE (dynamic)

      do j = 1, ao_num
        num_A = ao_nucl(j)
        power_A(1:3)= ao_power(j,1:3)
        A_center(1:3) = nucl_coord(num_A,1:3)

        do i = 1, ao_num

          num_B = ao_nucl(i)
          power_B(1:3)= ao_power(i,1:3)
          B_center(1:3) = nucl_coord(num_B,1:3)

          do l=1,ao_prim_num(j)
            alpha = ao_expo_ordered_transp(l,j)

            do m=1,ao_prim_num(i)
              beta = ao_expo_ordered_transp(m,i)

              double precision               :: c, c1
              c = 0.d0

              do  k = 1, n_pts_charge
                double precision               :: Z
                Z = pts_charge_z(k)

                C_center(1:3) = pts_charge_coord(k,1:3)

                c1 = NAI_pol_mult( A_center, B_center, power_A, power_B &
                                 , alpha, beta, C_center, n_pt_in )

                c = c + Z * c1

              enddo
              ao_integrals_pt_chrg(i,j) = ao_integrals_pt_chrg(i,j)  &
                  + ao_coef_normalized_ordered_transp(l,j)             &
                  * ao_coef_normalized_ordered_transp(m,i) * c
            enddo
          enddo
        enddo
      enddo

    !$OMP END DO
    !$OMP END PARALLEL

!    endif


!    IF(do_pseudo) THEN
!       ao_integrals_pt_chrg += ao_pseudo_integrals
!    ENDIF

!  endif


!  if (write_ao_integrals_pt_chrg) then
!    call ezfio_set_ao_one_e_ints_ao_integrals_pt_chrg(ao_integrals_pt_chrg)
!    print *,  'AO N-e integrals written to disk'
!  endif

END_PROVIDER
