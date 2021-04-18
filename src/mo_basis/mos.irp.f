BEGIN_PROVIDER [ integer, mo_num ]
  implicit none
  BEGIN_DOC
  ! Number of MOs
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    call ezfio_has_mo_basis_mo_num(has)
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer                        :: ierr
    call MPI_BCAST( has, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_num with MPI'
    endif
  IRP_ENDIF
  if (.not.has) then
    mo_num = ao_ortho_canonical_num
  else
    if (mpi_master) then
      call ezfio_get_mo_basis_mo_num(mo_num)
    endif
    IRP_IF MPI
      call MPI_BCAST( mo_num, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read mo_num with MPI'
      endif
    IRP_ENDIF
  endif
  call write_int(6,mo_num,'mo_num')
  ASSERT (mo_num > 0)

END_PROVIDER


BEGIN_PROVIDER [ double precision, mo_coef, (ao_num,mo_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on |AO| basis set
  !
  ! mo_coef(i,j) = coefficient of the i-th |AO| on the jth |MO|
  !
  ! mo_label : Label characterizing the |MOs| (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j
  double precision, allocatable  :: buffer(:,:)
  logical                        :: exists
  PROVIDE ezfio_filename


  if (mpi_master) then
    ! Coefs
    call ezfio_has_mo_basis_mo_coef(exists)
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
      stop 'Unable to read mo_coef with MPI'
    endif
  IRP_ENDIF

  if (exists) then
    if (mpi_master) then
      call ezfio_get_mo_basis_mo_coef(mo_coef)
      write(*,*) 'Read  mo_coef'
    endif
    IRP_IF MPI
      call MPI_BCAST( mo_coef, mo_num*ao_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read mo_coef with MPI'
      endif
    IRP_ENDIF
  else
    ! Orthonormalized AO basis
    do i=1,mo_num
      do j=1,ao_num
        mo_coef(j,i) = ao_ortho_canonical_coef(j,i)
      enddo
    enddo
  endif
END_PROVIDER

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

BEGIN_PROVIDER [ double precision, mo_coef_in_ao_ortho_basis, (ao_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! |MO| coefficients in orthogonalized |AO| basis
 !
 ! $C^{-1}.C_{mo}$
 END_DOC
 call dgemm('N','N',ao_num,mo_num,ao_num,1.d0,                   &
     ao_ortho_canonical_coef_inv, size(ao_ortho_canonical_coef_inv,1),&
     mo_coef, size(mo_coef,1), 0.d0,                                 &
     mo_coef_in_ao_ortho_basis, size(mo_coef_in_ao_ortho_basis,1))

END_PROVIDER

BEGIN_PROVIDER [ character*(64), mo_label ]
  implicit none
  BEGIN_DOC
  ! |MO| coefficients on |AO| basis set
  !
  ! mo_coef(i,j) = coefficient of the i-th |AO| on the j-th |MO|
  !
  ! mo_label : Label characterizing the |MOs| (local, canonical, natural, etc)
  END_DOC

  logical                        :: exists
  PROVIDE ezfio_filename
  if (mpi_master) then
    call ezfio_has_mo_basis_mo_label(exists)
    if (exists) then
      call ezfio_get_mo_basis_mo_label(mo_label)
      mo_label = trim(mo_label)
    else
      mo_label = 'no_label'
    endif
    write(*,*) '* mo_label          ', trim(mo_label)
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( mo_label, 64, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_label with MPI'
    endif
  IRP_ENDIF

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_coef_transp, (mo_num,ao_num) ]
  implicit none
  BEGIN_DOC
  ! |MO| coefficients on |AO| basis set
  END_DOC
  integer                        :: i, j

  do j=1,ao_num
    do i=1,mo_num
      mo_coef_transp(i,j) = mo_coef(j,i)
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, mo_occ, (mo_num) ]
  implicit none
  BEGIN_DOC
  ! |MO| occupation numbers
  END_DOC
  PROVIDE ezfio_filename elec_beta_num elec_alpha_num
  if (mpi_master) then
    logical :: exists
    call ezfio_has_mo_basis_mo_occ(exists)
    if (exists) then
      call ezfio_get_mo_basis_mo_occ(mo_occ)
    else
      mo_occ = 0.d0
      integer :: i
      do i=1,elec_beta_num
        mo_occ(i) = 2.d0
      enddo
      do i=elec_beta_num+1,elec_alpha_num
        mo_occ(i) = 1.d0
      enddo
    endif
    write(*,*) 'Read mo_occ'
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( mo_occ, mo_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_occ with MPI'
    endif
  IRP_ENDIF

END_PROVIDER


subroutine ao_to_mo(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the |AO| basis to the |MO| basis
  !
  ! $C^\dagger.A_{ao}.C$
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_ao(LDA_ao,ao_num)
  double precision, intent(out)  :: A_mo(LDA_mo,mo_num)
  double precision, allocatable  :: T(:,:)

  allocate ( T(ao_num,mo_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call dgemm('N','N', ao_num, mo_num, ao_num,                    &
      1.d0, A_ao,LDA_ao,                                             &
      mo_coef, size(mo_coef,1),                                      &
      0.d0, T, size(T,1))

  call dgemm('T','N', mo_num, mo_num, ao_num,                &
      1.d0, mo_coef,size(mo_coef,1),                                 &
      T, ao_num,                                                     &
      0.d0, A_mo, size(A_mo,1))

  call restore_symmetry(mo_num,mo_num,A_mo,size(A_mo,1),1.d-12)
  deallocate(T)
end


subroutine mix_mo_jk(j,k)
  implicit none
  integer, intent(in)            :: j,k
  integer                        :: i,i_plus,i_minus
  BEGIN_DOC
  ! Rotates the j-th |MO| with the k-th |MO| to give two new |MOs| that are
  !
  ! * $+ = \frac{1}{\sqrt{2}} ( | j\rangle +  | k\rangle)$
  !
  ! * $- = \frac{1}{\sqrt{2}} ( | j\rangle -  | k\rangle)$
  !
  ! by convention, the '+' |MO| is in the lowest  index (min(j,k))
  ! by convention, the '-' |MO| is in the highest index (max(j,k))
  END_DOC
  double precision               :: array_tmp(ao_num,2),dsqrt_2
  if(j==k)then
    print*,'You want to mix two orbitals that are the same !'
    print*,'It does not make sense ... '
    print*,'Stopping ...'
    stop
  endif
  array_tmp = 0.d0
  dsqrt_2 = 1.d0/dsqrt(2.d0)
  do i = 1, ao_num
    array_tmp(i,1) = dsqrt_2 * (mo_coef(i,j) + mo_coef(i,k))
    array_tmp(i,2) = dsqrt_2 * (mo_coef(i,j) - mo_coef(i,k))
  enddo
  i_plus = min(j,k)
  i_minus = max(j,k)
  do i = 1, ao_num
    mo_coef(i,i_plus) = array_tmp(i,1)
    mo_coef(i,i_minus) = array_tmp(i,2)
  enddo

end

subroutine ao_ortho_cano_to_ao(A_ao,LDA_ao,A,LDA)
  implicit none
  BEGIN_DOC
  ! Transform A from the |AO| basis to the orthogonal |AO| basis
  !
  ! $C^{-1}.A_{ao}.C^{\dagger-1}$
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA
  double precision, intent(in)   :: A_ao(LDA_ao,*)
  double precision, intent(out)  :: A(LDA,*)
  double precision, allocatable  :: T(:,:)

  allocate ( T(ao_num,ao_num) )

  call dgemm('T','N', ao_num, ao_num, ao_num,                        &
      1.d0,                                                          &
      ao_ortho_canonical_coef_inv, size(ao_ortho_canonical_coef_inv,1),&
      A_ao,size(A_ao,1),                                             &
      0.d0, T, size(T,1))

  call dgemm('N','N', ao_num, ao_num, ao_num, 1.d0,                  &
      T, size(T,1),                                                  &
      ao_ortho_canonical_coef_inv,size(ao_ortho_canonical_coef_inv,1),&
      0.d0, A, size(A,1))

  deallocate(T)
end

