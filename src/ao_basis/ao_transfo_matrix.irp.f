BEGIN_PROVIDER [ integer, ao_num  ]
  implicit none
  BEGIN_DOC
! Number of |AOs|
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    call ezfio_has_ao_basis_ao_num(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: ao_num ] <<<<< ..'
      call ezfio_get_ao_basis_ao_num(ao_num)
    else
      print *, 'Using the transformation matrix'
      if(.not.ao_cartesian)then
       ao_num = ao_sphe_num
      else
       ao_num = ao_cart_num 
      endif
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( ao_num, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read ao_num with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_cart_to_ao_basis_mat, ( ao_num,ao_cart_num)]
 implicit none
 if(.not.ao_cartesian)then
  ao_cart_to_ao_basis_mat(1:ao_num,1:ao_cart_num) = ao_cart_to_sphe_coef_inv(1:ao_sphe_num,1:ao_cart_num)
 else
  ao_cart_to_ao_basis_mat=0.d0
  integer :: i
  do i = 1, ao_num
   ao_cart_to_ao_basis_mat(i,i) = 1.d0
  enddo
 endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_cart_to_ao_basis_mat_transp, (ao_cart_num,ao_num)]
 implicit none
 integer :: i,j
 do i = 1, ao_num
  do j = 1, ao_cart_num
   ao_cart_to_ao_basis_mat_transp(j,i) = ao_cart_to_ao_basis_mat(i,j)
  enddo
 enddo
END_PROVIDER 

subroutine ao_cart_to_ao_basis(A_cart, LDA_cart, A_ao_basis, LDA_ao_basis)
 implicit none
 BEGIN_DOC
 ! transform a matrix "A_cart" expressed in the cartesian AO basis 
 ! 
 ! to the "AO" basis with the transformation matrix "ao_cart_to_ao_basis_mat" 
 !
 ! A_cart = ao_cart_to_ao_basis_mat A_cart ao_cart_to_ao_basis_mat^T
 END_DOC
 double precision, intent(in)  :: A_cart(LDA_cart, ao_cart_num)
 double precision, intent(out) :: A_ao_basis(LDA_ao_basis, ao_num)
 integer, intent(in) :: LDA_cart, LDA_ao_basis
 double precision, allocatable :: tmp(:,:)
 allocate (tmp(ao_num,ao_cart_num))

 call dgemm('N','N',ao_num,ao_cart_num,ao_cart_num, 1.d0, &
   ao_cart_to_ao_basis_mat,size(ao_cart_to_ao_basis_mat,1), &
   A_cart,size(A_cart,1), 0.d0, &
   tmp, size(tmp,1))

 call dgemm('N','T',ao_num,ao_num,ao_cart_num, 1.d0, &
   tmp, size(tmp,1), &
   ao_cart_to_ao_basis_mat,size(ao_cart_to_ao_basis_mat,1), 0.d0, &
   A_ao_basis,size(A_ao_basis,1))

 deallocate(tmp)

end

subroutine ao_cart_to_ao_basis_vec(V_cart, V_ao_basis)
 implicit none
 BEGIN_DOC
 ! transform a VECTOR "V_cart" expressed in the cartesian AO basis 
 ! 
 ! to the "AO" basis with the transformation matrix "ao_cart_to_ao_basis_mat" 
 !
 ! V_cart = ao_cart_to_ao_basis_mat V_cart 
 END_DOC
 double precision, intent(in)  :: V_cart(ao_cart_num)
 double precision, intent(out) :: V_ao_basis(ao_num)

 call dgemv('N',ao_num,ao_cart_num, 1.d0, &
   ao_cart_to_ao_basis_mat,size(ao_cart_to_ao_basis_mat,1), &
   V_cart,1, 0.d0, V_ao_basis, 1)

end

