BEGIN_PROVIDER [ double precision, mo_coef_imag, (ao_num,mo_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on |AO| basis set
  !
  ! mo_coef_imag(i,j) = coefficient of the i-th |AO| on the jth |MO|
  !
  ! mo_label : Label characterizing the |MOs| (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j
  double precision, allocatable  :: buffer(:,:)
  logical                        :: exists
  PROVIDE ezfio_filename


  if (mpi_master) then
    ! Coefs
    call ezfio_has_mo_basis_mo_coef_imag(exists)
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST(exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_coef_imag with MPI'
    endif
  IRP_ENDIF

  if (exists) then
    if (mpi_master) then
      call ezfio_get_mo_basis_mo_coef_imag(mo_coef_imag)
      write(*,*) 'Read  mo_coef_imag'
    endif
    IRP_IF MPI
      call MPI_BCAST( mo_coef_imag, mo_num*ao_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read mo_coef_imag with MPI'
      endif
    IRP_ENDIF
  else
    ! Orthonormalized AO basis
    do i=1,mo_num
      do j=1,ao_num
        mo_coef_imag(j,i) = 0.d0
      enddo
    enddo
  endif
END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_complex, (ao_num,mo_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on |AO| basis set
  !
  ! mo_coef_complex(i,j) = coefficient of the i-th |AO| on the jth |MO|
  !
  ! mo_label : Label characterizing the |MOs| (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j
  PROVIDE ezfio_filename

  provide mo_coef mo_coef_imag
  
  do i=1,mo_num
    do j=1,ao_num
      mo_coef_complex(j,i) = dcmplx(mo_coef(j,i),mo_coef_imag(j,i))
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_in_ao_ortho_basis_complex, (ao_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! |MO| coefficients in orthogonalized |AO| basis
 !
 ! $C^{-1}.C_{mo}$
 END_DOC
 call zgemm('N','N',ao_num,mo_num,ao_num,(1.d0,0.d0),                   &
     ao_ortho_canonical_coef_inv_complex, size(ao_ortho_canonical_coef_inv_complex,1),&
     mo_coef_complex, size(mo_coef_complex,1), (0.d0,0.d0),                                 &
     mo_coef_in_ao_ortho_basis_complex, size(mo_coef_in_ao_ortho_basis_complex,1))

END_PROVIDER

 BEGIN_PROVIDER [ complex*16, mo_coef_transp_complex, (mo_num,ao_num) ]
&BEGIN_PROVIDER [ complex*16, mo_coef_transp_complex_conjg, (mo_num,ao_num) ]
  implicit none
  BEGIN_DOC
  ! |MO| coefficients on |AO| basis set
  END_DOC
  integer                        :: i, j

  do j=1,ao_num
    do i=1,mo_num
      mo_coef_transp_complex(i,j) = mo_coef_complex(j,i)
      mo_coef_transp_complex_conjg(i,j) = dconjg(mo_coef_complex(j,i))
    enddo
  enddo

END_PROVIDER

subroutine ao_to_mo_complex(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the AO basis to the MO basis
  ! where A is complex in the AO basis
  !
  ! Ct.A_ao.C
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)   :: A_ao(LDA_ao,ao_num)
  complex*16, intent(out)  :: A_mo(LDA_mo,mo_num)
  complex*16, allocatable  :: T(:,:)
  
  allocate ( T(ao_num,mo_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  
  call zgemm('N','N', ao_num, mo_num, ao_num,                    &
      (1.d0,0.d0), A_ao,LDA_ao,                                      &
      mo_coef_complex, size(mo_coef_complex,1),                                      &
      (0.d0,0.d0), T, size(T,1))
  
  call zgemm('C','N', mo_num, mo_num, ao_num,                &
      (1.d0,0.d0), mo_coef_complex,size(mo_coef_complex,1),                          &
      T, ao_num,                                                     &   
      (0.d0,0.d0), A_mo, size(A_mo,1))
  
  deallocate(T)
end


subroutine ao_ortho_cano_to_ao_complex(A_ao,LDA_ao,A,LDA)
  implicit none
  BEGIN_DOC
  ! Transform A from the |AO| basis to the orthogonal |AO| basis
  !
  ! $C^{-1}.A_{ao}.C^{\dagger-1}$
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA
  complex*16, intent(in)   :: A_ao(LDA_ao,*)
  complex*16, intent(out)  :: A(LDA,*)
  complex*16, allocatable  :: T(:,:)

  allocate ( T(ao_num,ao_num) )

  call zgemm('C','N', ao_num, ao_num, ao_num,                        &
      (1.d0,0.d0),                                                     &
      ao_ortho_canonical_coef_inv_complex, size(ao_ortho_canonical_coef_inv_complex,1),&
      A_ao,size(A_ao,1),                                             &
      (0.d0,0.d0), T, size(T,1))

  call zgemm('N','N', ao_num, ao_num, ao_num, (1.d0,0.d0),             &
      T, size(T,1),                                                  &
      ao_ortho_canonical_coef_inv_complex,size(ao_ortho_canonical_coef_inv_complex,1),&
      (0.d0,0.d0), A, size(A,1))

  deallocate(T)
end

