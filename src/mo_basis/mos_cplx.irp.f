BEGIN_PROVIDER [ integer, mo_num_per_kpt ]
 implicit none
 BEGIN_DOC
 ! number of mos per kpt.
 END_DOC
 mo_num_per_kpt = mo_num/kpt_num
END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_complex, (ao_num,mo_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on |AO| basis set
  !
  ! mo_coef_imag(i,j) = coefficient of the i-th |AO| on the jth |MO|
  !
  ! mo_label : Label characterizing the |MOs| (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j
  logical                        :: exists
  PROVIDE ezfio_filename

  if (mpi_master) then
    ! Coefs
    call ezfio_has_mo_basis_mo_coef_complex(exists)
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
      stop 'Unable to read mo_coef_complex with MPI'
    endif
  IRP_ENDIF
  
  if (exists) then
    if (mpi_master) then
      call ezfio_get_mo_basis_mo_coef_complex(mo_coef_complex)
      write(*,*) 'Read  mo_coef_complex'
    endif
    IRP_IF MPI
      call MPI_BCAST( mo_coef_complex, mo_num*ao_num, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read mo_coef_complex with MPI'
      endif
    IRP_ENDIF
  else
    ! Orthonormalized AO basis
    do i=1,mo_num
      do j=1,ao_num
        mo_coef_complex(j,i) = ao_ortho_canonical_coef_complex(j,i)
      enddo
    enddo
  endif
END_PROVIDER


BEGIN_PROVIDER [ complex*16, mo_coef_in_ao_ortho_basis_complex, (ao_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! |MO| coefficients in orthogonalized |AO| basis
 !
 ! $C^{-1}.C_{mo}$
 END_DOC
 call zgemm('N','N',ao_num,mo_num,ao_num,(1.d0,0.d0),                   &
     ao_ortho_cano_coef_inv_cplx, size(ao_ortho_cano_coef_inv_cplx,1),&
     mo_coef_complex, size(mo_coef_complex,1), (0.d0,0.d0),                                 &
     mo_coef_in_ao_ortho_basis_complex, size(mo_coef_in_ao_ortho_basis_complex,1))

END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_complex_kpts, (ao_num_per_kpt, mo_num_per_kpt, kpt_num) ]
  implicit none
  BEGIN_DOC
  ! nonzero blocks of |MO| coefficients
  !
  END_DOC
  integer :: i,j,k, mo_shft, ao_shft
  mo_coef_complex_kpts = (0.d0,0.d0)

  do k=1,kpt_num
    mo_shft = (k-1)*mo_num_per_kpt
    ao_shft = (k-1)*ao_num_per_kpt
    do i=1,mo_num_per_kpt
      do j=1,ao_num_per_kpt
        mo_coef_complex_kpts(j,i,k) = mo_coef_complex(j+ao_shft,i+mo_shft)
      enddo
    enddo
  enddo

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
  ! C^\dagger.A_ao.C
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

subroutine ao_to_mo_noconjg_complex(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the AO basis to the MO basis
  ! where A is complex in the AO basis
  !
  ! C^T.A_ao.C
  ! needed for 4idx tranform in four_idx_novvvv
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
  
  call zgemm('T','N', mo_num, mo_num, ao_num,                &
      (1.d0,0.d0), mo_coef_complex,size(mo_coef_complex,1),                          &
      T, ao_num,                                                     &   
      (0.d0,0.d0), A_mo, size(A_mo,1))
  
  deallocate(T)
end


subroutine ao_ortho_cano_to_ao_cplx(A_ao,LDA_ao,A,LDA)
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
      ao_ortho_cano_coef_inv_cplx, size(ao_ortho_cano_coef_inv_cplx,1),&
      A_ao,size(A_ao,1),                                             &
      (0.d0,0.d0), T, size(T,1))

  call zgemm('N','N', ao_num, ao_num, ao_num, (1.d0,0.d0),             &
      T, size(T,1),                                                  &
      ao_ortho_cano_coef_inv_cplx,size(ao_ortho_cano_coef_inv_cplx,1),&
      (0.d0,0.d0), A, size(A,1))

  deallocate(T)
end

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!


BEGIN_PROVIDER [ complex*16, mo_coef_kpts, (ao_num_per_kpt, mo_num_per_kpt, kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on |AO| basis set
  !
  ! mo_coef_kpts(i,j,k) = coefficient of the i-th |AO| on the jth |MO| in kth kpt
  !
  ! mo_label : Label characterizing the |MOs| (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j, k
  logical                        :: exists
  PROVIDE ezfio_filename

  if (mpi_master) then
    ! Coefs
    call ezfio_has_mo_basis_mo_coef_complex(exists)
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
      stop 'Unable to read mo_coef_kpts with MPI'
    endif
  IRP_ENDIF
  
  if (exists) then
    if (mpi_master) then
      call ezfio_get_mo_basis_mo_coef_kpts(mo_coef_kpts)
      write(*,*) 'Read  mo_coef_kpts'
    endif
    IRP_IF MPI
      call MPI_BCAST( mo_coef_kpts, kpt_num*mo_num_per_kpt*ao_num_per_kpt, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read mo_coef_kpts with MPI'
      endif
    IRP_ENDIF
  else
    ! Orthonormalized AO basis

    do k=1,kpt_num
      do i=1,mo_num_per_kpt
        do j=1,ao_num_per_kpt
          mo_coef_kpts(j,i,k) = ao_ortho_canonical_coef_kpts(j,i,k)
        enddo
      enddo
    enddo
  endif
END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_in_ao_ortho_basis_kpts, (ao_num_per_kpt, mo_num_per_kpt, kpt_num) ]
 implicit none
 BEGIN_DOC
 ! |MO| coefficients in orthogonalized |AO| basis
 !
 ! $C^{-1}.C_{mo}$
 END_DOC
 integer :: k
 do k=1,kpt_num

   call zgemm('N','N',ao_num_per_kpt,mo_num_per_kpt,ao_num_per_kpt,(1.d0,0.d0),                   &
       ao_ortho_cano_coef_inv_kpts(:,:,k), size(ao_ortho_cano_coef_inv_kpts,1),&
       mo_coef_kpts(:,:,k), size(mo_coef_kpts,1), (0.d0,0.d0),                                 &
       mo_coef_in_ao_ortho_basis_kpts(:,:,k), size(mo_coef_in_ao_ortho_basis_kpts,1))
 enddo

END_PROVIDER

 BEGIN_PROVIDER [ complex*16, mo_coef_transp_kpts, (mo_num_per_kpt,ao_num_per_kpt,kpt_num) ]
&BEGIN_PROVIDER [ complex*16, mo_coef_transp_kpts_conjg, (mo_num_per_kpt,ao_num_per_kpt,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! |MO| coefficients on |AO| basis set
  END_DOC
  integer                        :: i, j, k

  do k=1,kpt_num
    do j=1,ao_num_per_kpt
      do i=1,mo_num_per_kpt
        mo_coef_transp_kpts(i,j,k) = mo_coef_kpts(j,i,k)
        mo_coef_transp_kpts_conjg(i,j,k) = dconjg(mo_coef_kpts(j,i,k))
      enddo
    enddo
  enddo

END_PROVIDER

subroutine ao_to_mo_kpts(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  !todo: check this
  BEGIN_DOC
  ! Transform A from the AO basis to the MO basis
  ! where A is complex in the AO basis
  !
  ! C^\dagger.A_ao.C
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)   :: A_ao(LDA_ao,ao_num_per_kpt,kpt_num)
  complex*16, intent(out)  :: A_mo(LDA_mo,mo_num_per_kpt,kpt_num)
  complex*16, allocatable  :: T(:,:)
  
  allocate ( T(ao_num_per_kpt,mo_num_per_kpt) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  integer :: k

  do k=1,kpt_num
    call zgemm('N','N', ao_num_per_kpt, mo_num_per_kpt, ao_num_per_kpt,   &
        (1.d0,0.d0), A_ao(:,:,k),LDA_ao,                                      &
        mo_coef_kpts(:,:,k), size(mo_coef_kpts,1),                                      &
        (0.d0,0.d0), T, size(T,1))
    
    call zgemm('C','N', mo_num_per_kpt, mo_num_per_kpt, ao_num_per_kpt,    &
        (1.d0,0.d0), mo_coef_kpts(:,:,k),size(mo_coef_kpts,1),             &
        T, ao_num_per_kpt,                                                     &   
        (0.d0,0.d0), A_mo(:,:,k), size(A_mo,1))
  enddo
  
  deallocate(T)
end

subroutine ao_to_mo_noconjg_kpts(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the AO basis to the MO basis
  ! where A is complex in the AO basis
  !
  ! C^T.A_ao.C
  ! needed for 4idx tranform in four_idx_novvvv
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)   :: A_ao(LDA_ao,ao_num_per_kpt,kpt_num)
  complex*16, intent(out)  :: A_mo(LDA_mo,mo_num_per_kpt,kpt_num)
  complex*16, allocatable  :: T(:,:)
  
  allocate ( T(ao_num_per_kpt,mo_num_per_kpt) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  integer :: k
  do k=1,kpt_num
    call zgemm('N','N', ao_num_per_kpt, mo_num_per_kpt, ao_num_per_kpt, &
        (1.d0,0.d0), A_ao,LDA_ao,                                      &
        mo_coef_kpts(:,:,k), size(mo_coef_kpts,1),                                      &
        (0.d0,0.d0), T, size(T,1))
    
    call zgemm('T','N', mo_num_per_kpt, mo_num_per_kpt, ao_num_per_kpt, &
        (1.d0,0.d0), mo_coef_kpts(:,:,k),size(mo_coef_kpts,1),            &
        T, ao_num_per_kpt,                                                     &   
        (0.d0,0.d0), A_mo(:,:,k), size(A_mo,1))
  enddo 
  deallocate(T)
end


subroutine ao_ortho_cano_to_ao_kpts(A_ao,LDA_ao,A,LDA)
  implicit none
  !todo: check this; no longer using assumed-size arrays
  BEGIN_DOC
  ! Transform A from the |AO| basis to the orthogonal |AO| basis
  !
  ! $C^{-1}.A_{ao}.C^{\dagger-1}$
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA
  complex*16, intent(in)   :: A_ao(LDA_ao,ao_num_per_kpt,kpt_num)
  complex*16, intent(out)  :: A(LDA,ao_num_per_kpt,kpt_num)
  complex*16, allocatable  :: T(:,:)

  allocate ( T(ao_num_per_kpt,ao_num_per_kpt) )

  integer :: k
  do k=1,kpt_num
    call zgemm('C','N', ao_num_per_kpt, ao_num_per_kpt, ao_num_per_kpt,                        &
        (1.d0,0.d0),                                                     &
        ao_ortho_cano_coef_inv_kpts(:,:,k), size(ao_ortho_cano_coef_inv_kpts,1),&
        A_ao(:,:,k),size(A_ao,1),                                             &
        (0.d0,0.d0), T, size(T,1))

    call zgemm('N','N', ao_num_per_kpt, ao_num_per_kpt, ao_num_per_kpt, (1.d0,0.d0),             &
        T, size(T,1),                                                  &
        ao_ortho_cano_coef_inv_kpts(:,:,k),size(ao_ortho_cano_coef_inv_kpts,1),&
        (0.d0,0.d0), A(:,:,k), size(A,1))
  enddo

  deallocate(T)
end


!============================================!
!                                            !
!               elec kpts                    !
!                                            !
!============================================!

 BEGIN_PROVIDER [ integer, elec_alpha_num_kpts, (kpt_num) ]
&BEGIN_PROVIDER [ integer, elec_beta_num_kpts, (kpt_num) ]
  !todo: reorder? if not integer multiple, use some list of kpts to determine filling order
  implicit none
 
  integer :: i,k,kpt

  PROVIDE elec_alpha_num elec_beta_num

  do k=1,kpt_num
    elec_alpha_num_kpts(k) = 0
    elec_beta_num_kpts(k) = 0
  enddo
  kpt=1
  do i=1,elec_beta_num
    elec_alpha_num_kpts(kpt) += 1
    elec_beta_num_kpts(kpt) += 1
    kpt += 1
    if (kpt > kpt_num) then
      kpt = 1
    endif
  enddo
  do i=elec_beta_num+1,elec_alpha_num
    elec_alpha_num_kpts(kpt) += 1
    kpt += 1
    if (kpt > kpt_num) then
      kpt = 1
    endif
  enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_occ_kpts, (mo_num_per_kpt,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! |MO| occupation numbers
  END_DOC
  PROVIDE ezfio_filename elec_beta_num_kpts elec_alpha_num_kpts
  if (mpi_master) then
    logical :: exists
    call ezfio_has_mo_basis_mo_occ_kpts(exists)
    if (exists) then
      call ezfio_get_mo_basis_mo_occ_kpts(mo_occ_kpts)
    else
      mo_occ_kpts = 0.d0
      integer :: i,k
      do k=1,kpt_num
        do i=1,elec_beta_num_kpts(k)
          mo_occ_kpts(i,k) = 2.d0
        enddo
        do i=elec_beta_num_kpts(k)+1,elec_alpha_num_kpts(k)
          mo_occ_kpts(i,k) = 1.d0
        enddo
      enddo
    endif
    write(*,*) 'Read mo_occ_kpts'
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( mo_occ_kpts, mo_num_per_kpt*kpt_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_occ_kpts with MPI'
    endif
  IRP_ENDIF

END_PROVIDER
