BEGIN_PROVIDER [ double precision, extra_nucl_coord,  (extra_nucl_num,3) ]
   implicit none

   BEGIN_DOC
   ! Nuclear coordinates in the format (:, {x,y,z})
   END_DOC
   PROVIDE ezfio_filename extra_nucl_label extra_nucl_charge

   if (mpi_master) then
     double precision, allocatable  :: buffer(:,:)
     extra_nucl_coord = 0.d0
     allocate (buffer(extra_nucl_num,3))
     buffer = 0.d0
     logical                        :: has
     call ezfio_has_extra_nuclei_extra_nucl_coord(has)
     if (.not.has) then
       print *, irp_here
       stop 1
     endif
     call ezfio_get_extra_nuclei_extra_nucl_coord(buffer)
     integer                        :: i,j

     do i=1,3
       do j=1,extra_nucl_num
         extra_nucl_coord(j,i) = buffer(j,i)
       enddo
     enddo
     deallocate(buffer)

     character*(64), parameter      :: f = '(A16, 4(1X,F12.6))'
     character*(64), parameter      :: ft= '(A16, 4(1X,A12  ))'
     double precision, parameter    :: a0= 0.529177249d0

     call write_time(6)
     write(6,'(A)') ''
     write(6,'(A)') 'Extra Nuclear Coordinates (Angstroms)'
     write(6,'(A)') '====================================='
     write(6,'(A)') ''
     write(6,ft)                                         &
         '================','============','============','============','============'
     write(6,*)                                          &
         '     Atom          Charge          X            Y            Z '
     write(6,ft)                                         &
         '================','============','============','============','============'

     do i=1,extra_nucl_num
       write(6,f) extra_nucl_label(i), extra_nucl_charge(i),         &
           extra_nucl_coord(i,1)*a0,                                       &
           extra_nucl_coord(i,2)*a0,                                       &
           extra_nucl_coord(i,3)*a0
     enddo
     write(6,ft)                                         &
         '================','============','============','============','============'
     write(6,'(A)') ''

     if (extra_nucl_num > 1) then
       double precision               :: dist_min, x, y, z
       dist_min = huge(1.d0)
       do i=1,extra_nucl_num
         do j=i+1,extra_nucl_num
           x = extra_nucl_coord(i,1)-extra_nucl_coord(j,1)
           y = extra_nucl_coord(i,2)-extra_nucl_coord(j,2)
           z = extra_nucl_coord(i,3)-extra_nucl_coord(j,3)
           dist_min = min(x*x + y*y + z*z, dist_min)
         enddo
       enddo
       write(6,'(A,F12.4,A)') 'Minimal interatomic distance found: ', &
          dsqrt(dist_min)*a0,' Angstrom'
     endif

   endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
   IRP_IF MPI
     include 'mpif.h'
     integer                        :: ierr
     call MPI_BCAST( extra_nucl_coord, 3*extra_nucl_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     if (ierr /= MPI_SUCCESS) then
       stop 'Unable to read nucl_coord with MPI'
     endif
   IRP_ENDIF

END_PROVIDER


BEGIN_PROVIDER [ double precision, extra_nucl_coord_transp, (3,extra_nucl_num) ]
   implicit none
   BEGIN_DOC
   ! Transposed array of extra_nucl_coord
   END_DOC
   integer                        :: i, k
   extra_nucl_coord_transp = 0.d0

   do i=1,extra_nucl_num
     extra_nucl_coord_transp(1,i) = extra_nucl_coord(i,1)
     extra_nucl_coord_transp(2,i) = extra_nucl_coord(i,2)
     extra_nucl_coord_transp(3,i) = extra_nucl_coord(i,3)
   enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, extra_center_of_mass, (3) ]
  implicit none
  BEGIN_DOC
  ! Center of mass of the molecule
  END_DOC
  integer                        :: i,j
  double precision               :: s
  extra_center_of_mass(:) = 0.d0
  s = 0.d0
  do i=1,extra_nucl_num
    do j=1,3
      extra_center_of_mass(j) += extra_nucl_coord(i,j)* element_mass(int(extra_nucl_charge(i)))
    enddo
    s += element_mass(int(extra_nucl_charge(i)))
  enddo
  s = 1.d0/s
  extra_center_of_mass(:) = extra_center_of_mass(:)*s
END_PROVIDER


BEGIN_PROVIDER [ double precision, extra_nucl_dist, (extra_nucl_num,extra_nucl_num)]
 implicit none
 integer :: i,j
 double precision :: x,y,z
 do i = 1, extra_nucl_num
  do j = 1, extra_nucl_num
   x = extra_nucl_coord(i,1)-extra_nucl_coord(j,1)
   y = extra_nucl_coord(i,2)-extra_nucl_coord(j,2)
   z = extra_nucl_coord(i,3)-extra_nucl_coord(j,3)
   extra_nucl_dist(j,i) = x*x+y*y+z*z
   extra_nucl_dist(j,i) = dsqrt(extra_nucl_dist(j,i))
  enddo
 enddo

END_PROVIDER 
