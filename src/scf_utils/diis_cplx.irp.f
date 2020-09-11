
BEGIN_PROVIDER [complex*16, FPS_SPF_Matrix_AO_complex, (AO_num, AO_num)]
  implicit none
  BEGIN_DOC
  !   Commutator FPS - SPF
  END_DOC
  complex*16, allocatable  :: scratch(:,:)
  allocate(                                                          &
      scratch(AO_num, AO_num)                                  &
      )

  ! Compute FP

  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (1.d0,0.d0),                                                          &
      Fock_Matrix_AO_complex,Size(Fock_Matrix_AO_complex,1),                &
      SCF_Density_Matrix_AO_complex,Size(SCF_Density_Matrix_AO_complex,1),  &
      (0.d0,0.d0),                                                          &
      scratch,Size(scratch,1))

  ! Compute FPS

  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (1.d0,0.d0),                                                   &
      scratch,Size(scratch,1),                                       &
      AO_Overlap_complex,Size(AO_Overlap_complex,1),                 &
      (0.d0,0.d0),                                                  &
      FPS_SPF_Matrix_AO_complex,Size(FPS_SPF_Matrix_AO_complex,1))

  ! Compute SP

  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (1.d0,0.d0),                                                          &
      AO_Overlap_complex,Size(AO_Overlap_complex,1),                        &
      SCF_Density_Matrix_AO_complex,Size(SCF_Density_Matrix_AO_complex,1),  &
      (0.d0,0.d0),                                                          &
      scratch,Size(scratch,1))

  ! Compute FPS - SPF

  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (-1.d0,0.d0),                                                  &
      scratch,Size(scratch,1),                                       &
      Fock_Matrix_AO_complex,Size(Fock_Matrix_AO_complex,1),         &
      (1.d0,0.d0),                                                   &
      FPS_SPF_Matrix_AO_complex,Size(FPS_SPF_Matrix_AO_complex,1))

END_PROVIDER

BEGIN_PROVIDER [complex*16, FPS_SPF_Matrix_MO_complex, (mo_num, mo_num)]
  implicit none
  begin_doc
!   Commutator FPS - SPF in MO basis
  end_doc
  call ao_to_mo_complex(FPS_SPF_Matrix_AO_complex, size(FPS_SPF_Matrix_AO_complex,1), &
     FPS_SPF_Matrix_MO_complex, size(FPS_SPF_Matrix_MO_complex,1))
END_PROVIDER


 BEGIN_PROVIDER [ double precision, eigenvalues_Fock_matrix_AO_complex, (AO_num) ]
&BEGIN_PROVIDER [ complex*16, eigenvectors_Fock_matrix_AO_complex, (AO_num,AO_num) ]
   !TODO: finish this provider; write provider for S_half_inv_complex
   BEGIN_DOC
   ! Eigenvalues and eigenvectors of the Fock matrix over the AO basis
   END_DOC

   implicit none

   double precision, allocatable  :: rwork(:)
   integer                        :: lwork,info,lrwork
   complex*16, allocatable        :: scratch(:,:),Xt(:,:),work(:)
   integer                        :: i,j

   allocate(                                                         &
       scratch(AO_num,AO_num),                                 &
       Xt(AO_num,AO_num)                                             &
       )

! Calculate Xt

  do i=1,AO_num
    do j=1,AO_num
      Xt(i,j) = dconjg(S_half_inv_complex(j,i))
    enddo
  enddo

! Calculate Fock matrix in orthogonal basis: F' = Xt.F.X

  call zgemm('N','N',AO_num,AO_num,AO_num,     &
       (1.d0,0.d0),                                   &
       Fock_matrix_AO_complex,size(Fock_matrix_AO_complex,1),  &
       S_half_inv_complex,size(s_half_inv_complex,1),        &
       (0.d0,0.d0),                                   &
       eigenvectors_Fock_matrix_AO_complex, &
       size(eigenvectors_Fock_matrix_AO_complex,1))

  call zgemm('N','N',AO_num,AO_num,AO_num,                              &
       (1.d0,0.d0),                                                            &
       Xt,size(Xt,1),                                                   &
       eigenvectors_Fock_matrix_AO_complex, &
       size(eigenvectors_Fock_matrix_AO_complex,1), &
       (0.d0,0.d0),                                                            &
       scratch,size(scratch,1))

! Diagonalize F' to obtain eigenvectors in orthogonal basis C' and eigenvalues
   lrwork = 3*ao_num - 2
   allocate(rwork(lrwork), work(1))
   lwork = -1

   call zheev('V','U',ao_num, &
        scratch,size(scratch,1), &
        eigenvalues_Fock_matrix_AO_complex, &
        work,lwork,rwork,info)

   lwork = int(work(1))
   deallocate(work)
   allocate(work(lwork))

   call zheev('V','U',ao_num, &
        scratch,size(scratch,1), &
        eigenvalues_Fock_matrix_AO_complex, &
        work,lwork,rwork,info)

   if(info /= 0) then
     print *,  irp_here//' failed : ', info
     stop 1
   endif
   
   deallocate(work,rwork)
! Back-transform eigenvectors: C =X.C'

  call zgemm('N','N',AO_num,AO_num,AO_num,     &
       (1.d0,0.d0),                                   &
       S_half_inv_complex,size(S_half_inv_complex,1),        &
       scratch,size(scratch,1),                &
       (0.d0,0.d0),                                   &
       eigenvectors_Fock_matrix_AO_complex, &
       size(eigenvectors_Fock_matrix_AO_complex,1))

  deallocate(scratch)
END_PROVIDER

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!

BEGIN_PROVIDER [complex*16, FPS_SPF_Matrix_AO_kpts, (AO_num_per_kpt, AO_num_per_kpt,kpt_num)]
  implicit none
  BEGIN_DOC
  !   Commutator FPS - SPF
  END_DOC
  complex*16, allocatable  :: scratch(:,:)
  integer :: k
  allocate(                                                          &
      scratch(ao_num_per_kpt, ao_num_per_kpt)                                  &
      )

  do k=1,kpt_num

    ! Compute FP

    call zgemm('N','N',AO_num_per_kpt,AO_num_per_kpt,AO_num_per_kpt,                           &
        (1.d0,0.d0),                                                          &
        Fock_Matrix_AO_kpts(1,1,k),Size(Fock_Matrix_AO_kpts,1),                &
        scf_density_matrix_ao_kpts(1,1,k),Size(SCF_Density_Matrix_AO_kpts,1),  &
        (0.d0,0.d0),                                                          &
        scratch,Size(scratch,1))

    ! Compute FPS

    call zgemm('N','N',AO_num_per_kpt,AO_num_per_kpt,AO_num_per_kpt,                           &
        (1.d0,0.d0),                                                   &
        scratch,Size(scratch,1),                                       &
        AO_Overlap_kpts(1,1,k),Size(AO_Overlap_kpts,1),                 &
        (0.d0,0.d0),                                                  &
        FPS_SPF_Matrix_AO_kpts(1,1,k),Size(FPS_SPF_Matrix_AO_kpts,1))

    ! Compute SP

    call zgemm('N','N',AO_num_per_kpt,AO_num_per_kpt,AO_num_per_kpt,                           &
        (1.d0,0.d0),                                                          &
        AO_Overlap_kpts(1,1,k),Size(AO_Overlap_kpts,1),                        &
        SCF_Density_Matrix_AO_kpts(1,1,k),Size(SCF_Density_Matrix_AO_kpts,1),  &
        (0.d0,0.d0),                                                          &
        scratch,Size(scratch,1))

    ! Compute FPS - SPF

    call zgemm('N','N',AO_num_per_kpt,AO_num_per_kpt,AO_num_per_kpt,                           &
        (-1.d0,0.d0),                                                  &
        scratch,Size(scratch,1),                                       &
        Fock_Matrix_AO_kpts(1,1,k),Size(Fock_Matrix_AO_kpts,1),         &
        (1.d0,0.d0),                                                   &
        FPS_SPF_Matrix_AO_kpts(1,1,k),Size(FPS_SPF_Matrix_AO_kpts,1))
  enddo
END_PROVIDER

BEGIN_PROVIDER [complex*16, FPS_SPF_Matrix_MO_kpts, (mo_num_per_kpt, mo_num_per_kpt,kpt_num)]
  implicit none
  begin_doc
!   Commutator FPS - SPF in MO basis
  end_doc
  call ao_to_mo_kpts(FPS_SPF_Matrix_AO_kpts, size(FPS_SPF_Matrix_AO_kpts,1), &
     FPS_SPF_Matrix_MO_kpts, size(FPS_SPF_Matrix_MO_kpts,1))
END_PROVIDER


 BEGIN_PROVIDER [ double precision, eigenvalues_fock_matrix_ao_kpts, (ao_num_per_kpt,kpt_num) ]
&BEGIN_PROVIDER [ complex*16, eigenvectors_fock_matrix_ao_kpts, (ao_num_per_kpt,ao_num_per_kpt,kpt_num) ]
   !TODO: finish this provider; write provider for S_half_inv_complex
   BEGIN_DOC
   ! Eigenvalues and eigenvectors of the Fock matrix over the AO basis
   END_DOC

   implicit none

   double precision, allocatable  :: rwork(:)
   integer                        :: lwork,info,lrwork
   complex*16, allocatable        :: scratch(:,:),Xt(:,:),work(:)
   integer                        :: i,j,k


   allocate(                                                         &
       scratch(ao_num_per_kpt,ao_num_per_kpt),                                 &
       Xt(ao_num_per_kpt,ao_num_per_kpt)                                             &
       )

  do k=1,kpt_num
  ! Calculate Xt
  
    do i=1,ao_num_per_kpt
      do j=1,ao_num_per_kpt
!        Xt(i,j) = dconjg(s_half_inv_complex(j,i,k))
        Xt(i,j) = dconjg(S_half_inv_kpts(j,i,k))
      enddo
    enddo
  
  ! Calculate Fock matrix in orthogonal basis: F' = Xt.F.X
  
    call zgemm('N','N',ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt,     &
         (1.d0,0.d0),                                   &
         fock_matrix_ao_kpts(1,1,k),size(fock_matrix_ao_kpts,1),  &
         s_half_inv_kpts(1,1,k),size(s_half_inv_kpts,1),        &
         (0.d0,0.d0),                                   &
         eigenvectors_fock_matrix_ao_kpts(1,1,k), &
         size(eigenvectors_fock_matrix_ao_kpts,1))
  
    call zgemm('N','N',ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt,                              &
         (1.d0,0.d0),                                                            &
         Xt,size(Xt,1),                                                   &
         eigenvectors_fock_matrix_ao_kpts(1,1,k), &
         size(eigenvectors_fock_matrix_ao_kpts,1), &
         (0.d0,0.d0),                                                            &
         scratch,size(scratch,1))
  
  ! Diagonalize F' to obtain eigenvectors in orthogonal basis C' and eigenvalues
     lrwork = 3*ao_num_per_kpt - 2
     allocate(rwork(lrwork), work(1))
     lwork = -1
  
     call zheev('V','U',ao_num_per_kpt, &
          scratch,size(scratch,1), &
          eigenvalues_fock_matrix_ao_kpts(1,k), &
          work,lwork,rwork,info)
  
     lwork = int(work(1))
     deallocate(work)
     allocate(work(lwork))
  
     call zheev('V','U',ao_num_per_kpt, &
          scratch,size(scratch,1), &
          eigenvalues_fock_matrix_ao_kpts(1,k), &
          work,lwork,rwork,info)
  
     if(info /= 0) then
       print *,  irp_here//' failed : ', info
       stop 1
     endif
     
     deallocate(work,rwork)
  ! Back-transform eigenvectors: C =X.C'
  
    call zgemm('N','N',ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt,     &
         (1.d0,0.d0),                                   &
         s_half_inv_kpts(1,1,k),size(s_half_inv_kpts,1),        &
         scratch,size(scratch,1),                &
         (0.d0,0.d0),                                   &
         eigenvectors_fock_matrix_ao_kpts(1,1,k), &
         size(eigenvectors_fock_matrix_ao_kpts,1))
  enddo
  deallocate(scratch)
END_PROVIDER
