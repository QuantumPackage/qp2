program import_kconserv

  PROVIDE ezfio_filename
  call run
end

subroutine run
  use map_module
  implicit none
  BEGIN_DOC
  ! read kconserv in physicists' notation order <ij|kl>
  ! if kconserv(i,j,k)=l, then <ij|kl> is allowed by symmetry
  ! NOTE: pyscf stores this internally in the order of chemists' notation (ik|jl)
  END_DOC
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  integer, allocatable :: A(:,:,:)

  allocate(A(kpt_num,kpt_num,kpt_num))
  
  A = 0
  iunit = getunitandopen('kconserv','r')
  do 
    read (iunit,*,end=10) i,j,k,l
    A(i,j,k) = l
  enddo
  10 continue
  close(iunit)
  call ezfio_set_nuclei_kconserv(A)
  call ezfio_set_nuclei_io_kconserv("Read")
  deallocate(A)

end
