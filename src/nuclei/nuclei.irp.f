BEGIN_PROVIDER [ double precision, nucl_coord,  (nucl_num,3) ]
   implicit none

   BEGIN_DOC
   ! Nuclear coordinates in the format (:, {x,y,z})
   END_DOC
   PROVIDE ezfio_filename nucl_label nucl_charge

   if (mpi_master) then
     double precision, allocatable  :: buffer(:,:)
     nucl_coord = 0.d0
     allocate (buffer(nucl_num,3))
     buffer = 0.d0
     logical                        :: has
     call ezfio_has_nuclei_nucl_coord(has)
     if (.not.has) then
       print *, irp_here
       stop 1
     endif
     call ezfio_get_nuclei_nucl_coord(buffer)
     integer                        :: i,j

     do i=1,3
       do j=1,nucl_num
         nucl_coord(j,i) = buffer(j,i)
       enddo
     enddo
     deallocate(buffer)

     character*(64), parameter      :: f = '(A16, 4(1X,F12.6))'
     character*(64), parameter      :: ft= '(A16, 4(1X,A12  ))'
     double precision, parameter    :: a0= 0.529177249d0

     call write_time(6)
     write(6,'(A)') ''
     write(6,'(A)') 'Nuclear Coordinates (Angstroms)'
     write(6,'(A)') '==============================='
     write(6,'(A)') ''
     write(6,ft)                                         &
         '================','============','============','============','============'
     write(6,*)                                          &
         '     Atom          Charge          X            Y            Z '
     write(6,ft)                                         &
         '================','============','============','============','============'

     do i=1,nucl_num
       write(6,f) nucl_label(i), nucl_charge(i),         &
           nucl_coord(i,1)*a0,                                       &
           nucl_coord(i,2)*a0,                                       &
           nucl_coord(i,3)*a0
     enddo
     write(6,ft)                                         &
         '================','============','============','============','============'
     write(6,'(A)') ''

     if (nucl_num > 1) then
       double precision               :: dist_min, x, y, z
       dist_min = huge(1.d0)
       do i=1,nucl_num
         do j=i+1,nucl_num
           x = nucl_coord(i,1)-nucl_coord(j,1)
           y = nucl_coord(i,2)-nucl_coord(j,2)
           z = nucl_coord(i,3)-nucl_coord(j,3)
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
     call MPI_BCAST( nucl_coord, 3*nucl_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     if (ierr /= MPI_SUCCESS) then
       stop 'Unable to read nucl_coord with MPI'
     endif
   IRP_ENDIF

END_PROVIDER


BEGIN_PROVIDER [ double precision, nucl_coord_transp, (3,nucl_num) ]
   implicit none
   BEGIN_DOC
   ! Transposed array of nucl_coord
   END_DOC
   integer                        :: i, k
   nucl_coord_transp = 0.d0

   do i=1,nucl_num
     nucl_coord_transp(1,i) = nucl_coord(i,1)
     nucl_coord_transp(2,i) = nucl_coord(i,2)
     nucl_coord_transp(3,i) = nucl_coord(i,3)
   enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, nucl_dist_inv, (nucl_num,nucl_num) ]
  implicit none
  BEGIN_DOC
  ! Inverse of the distance between nucleus I and nucleus J
  END_DOC

  integer                        :: ie1, ie2, l

  do ie1 = 1, nucl_num
    do ie2 = 1, nucl_num
      if(ie1 /= ie2) then
        nucl_dist_inv(ie2,ie1) = 1.d0/nucl_dist(ie2,ie1)
      else
        nucl_dist_inv(ie2,ie1) = 0.d0
      endif
    enddo
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, nucl_dist_2, (nucl_num,nucl_num) ]
&BEGIN_PROVIDER [ double precision, nucl_dist_vec_x, (nucl_num,nucl_num) ]
&BEGIN_PROVIDER [ double precision, nucl_dist_vec_y, (nucl_num,nucl_num) ]
&BEGIN_PROVIDER [ double precision, nucl_dist_vec_z, (nucl_num,nucl_num) ]
&BEGIN_PROVIDER [ double precision, nucl_dist, (nucl_num,nucl_num) ]
   implicit none
   BEGIN_DOC
   ! nucl_dist     : Nucleus-nucleus distances

   ! nucl_dist_2   : Nucleus-nucleus distances squared

   ! nucl_dist_vec : Nucleus-nucleus distances vectors
   END_DOC

   integer                        :: ie1, ie2, l

   do ie2 = 1,nucl_num
     do ie1 = 1,nucl_num
       nucl_dist_vec_x(ie1,ie2) = nucl_coord(ie1,1) - nucl_coord(ie2,1)
       nucl_dist_vec_y(ie1,ie2) = nucl_coord(ie1,2) - nucl_coord(ie2,2)
       nucl_dist_vec_z(ie1,ie2) = nucl_coord(ie1,3) - nucl_coord(ie2,3)
     enddo
     do ie1 = 1,nucl_num
       nucl_dist_2(ie1,ie2) = nucl_dist_vec_x(ie1,ie2)*nucl_dist_vec_x(ie1,ie2) +&
           nucl_dist_vec_y(ie1,ie2)*nucl_dist_vec_y(ie1,ie2) +       &
           nucl_dist_vec_z(ie1,ie2)*nucl_dist_vec_z(ie1,ie2)
       nucl_dist(ie1,ie2) = sqrt(nucl_dist_2(ie1,ie2))
     enddo
   enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, nuclear_repulsion ]
   implicit none
   BEGIN_DOC
   ! Nuclear repulsion energy
   END_DOC

   PROVIDE mpi_master nucl_coord nucl_charge nucl_num
   if (io_nuclear_repulsion == 'Read') then
     logical                        :: has

     if (mpi_master) then
       call ezfio_has_nuclei_nuclear_repulsion(has)
       if (has) then
         call ezfio_get_nuclei_nuclear_repulsion(nuclear_repulsion)
       else
         print *, 'nuclei/nuclear_repulsion not found in EZFIO file'
         stop 1
       endif
       print*, 'Read nuclear_repulsion'
     endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
     IRP_IF MPI
      include 'mpif.h'
      integer                        :: ierr
      call MPI_BCAST( nuclear_repulsion, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read nuclear_repulsion with MPI'
      endif
     IRP_ENDIF


   else

     integer                        :: k,l
     double precision               :: Z12, r2, x(3)
     nuclear_repulsion = 0.d0
     do l = 1, nucl_num
       do  k = 1, nucl_num
         if(k == l) then
           cycle
         endif
         Z12 = nucl_charge(k)*nucl_charge(l)
         x(1) = nucl_coord(k,1) - nucl_coord(l,1)
         x(2) = nucl_coord(k,2) - nucl_coord(l,2)
         x(3) = nucl_coord(k,3) - nucl_coord(l,3)
         r2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
         nuclear_repulsion += Z12/dsqrt(r2)
       enddo
     enddo
     nuclear_repulsion *= 0.5d0
     if(point_charges)then
      nuclear_repulsion += pt_chrg_nuclei_repulsion + pt_chrg_repulsion
     endif
   end if

   call write_time(6)
   call write_double(6,nuclear_repulsion,'Nuclear repulsion energy')

   if (io_nuclear_repulsion == 'Write') then
     if (mpi_master) then
       call ezfio_set_nuclei_nuclear_repulsion(nuclear_repulsion)
     endif
   endif


END_PROVIDER

 BEGIN_PROVIDER [ character*(4), element_name, (0:127)]
&BEGIN_PROVIDER [ double precision, element_mass, (0:127) ]
   implicit none
   BEGIN_DOC
   ! Array of the name of element, sorted by nuclear charge (integer)
   END_DOC
   integer                        :: iunit, i
   integer, external              :: getUnitAndOpen
   character*(1024)               :: filename
   if (mpi_master) then
     call getenv('QP_ROOT',filename)
     filename = trim(filename)//'/data/list_element.txt'
     iunit = getUnitAndOpen(filename,'r')
     element_mass(:) = 0.d0
     do i=0,127
       write(element_name(i),'(I4)') i
     enddo
     character*(80)                 :: buffer, dummy
     do
     read(iunit,'(A80)',end=10) buffer
     read(buffer,*) i ! First read i
     read(buffer,*) i, element_name(i), dummy, element_mass(i)
   enddo
   10 continue
   close(10)
 endif

 IRP_IF MPI_DEBUG
   print *,  irp_here, mpi_rank
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
 IRP_ENDIF
 IRP_IF MPI
  include 'mpif.h'
  integer                        :: ierr
  call MPI_BCAST( element_name, 128*4, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) then
    stop 'Unable to read element_name with MPI'
  endif
  call MPI_BCAST( element_mass, 128, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) then
    stop 'Unable to read element_name with MPI'
  endif
 IRP_ENDIF

END_PROVIDER

BEGIN_PROVIDER [ double precision, center_of_mass, (3) ]
  implicit none
  BEGIN_DOC
  ! Center of mass of the molecule
  END_DOC
  integer                        :: i,j
  double precision               :: s
  center_of_mass(:) = 0.d0
  s = 0.d0
  do i=1,nucl_num
    do j=1,3
      center_of_mass(j) += nucl_coord(i,j)* element_mass(int(nucl_charge(i)))
    enddo
    s += element_mass(int(nucl_charge(i)))
  enddo
  s = 1.d0/s
  center_of_mass(:) = center_of_mass(:)*s
END_PROVIDER

