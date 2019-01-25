BEGIN_PROVIDER [ double precision, threshold_DIIS_nonzero ]
 implicit none
 BEGIN_DOC
 ! If threshold_DIIS is zero, choose sqrt(thresh_scf)
 END_DOC
 if (threshold_DIIS == 0.d0) then
   threshold_DIIS_nonzero = dsqrt(thresh_scf)
 else
   threshold_DIIS_nonzero = threshold_DIIS
 endif
 ASSERT (threshold_DIIS_nonzero >= 0.d0)

END_PROVIDER

BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_AO, (AO_num, AO_num)]
  implicit none
  BEGIN_DOC
  !   Commutator FPS - SPF
  END_DOC
  double precision, allocatable  :: scratch(:,:)
  allocate(                                                          &
      scratch(AO_num, AO_num)                                  &
      )

  ! Compute FP

  call dgemm('N','N',AO_num,AO_num,AO_num,                           &
      1.d0,                                                          &
      Fock_Matrix_AO,Size(Fock_Matrix_AO,1),                         &
      SCF_Density_Matrix_AO,Size(SCF_Density_Matrix_AO,1),             &
      0.d0,                                                          &
      scratch,Size(scratch,1))

  ! Compute FPS

  call dgemm('N','N',AO_num,AO_num,AO_num,                           &
      1.d0,                                                          &
      scratch,Size(scratch,1),                                       &
      AO_Overlap,Size(AO_Overlap,1),                                 &
      0.d0,                                                          &
      FPS_SPF_Matrix_AO,Size(FPS_SPF_Matrix_AO,1))

  ! Compute SP

  call dgemm('N','N',AO_num,AO_num,AO_num,                           &
      1.d0,                                                          &
      AO_Overlap,Size(AO_Overlap,1),                                 &
      SCF_Density_Matrix_AO,Size(SCF_Density_Matrix_AO,1),             &
      0.d0,                                                          &
      scratch,Size(scratch,1))

  ! Compute FPS - SPF

  call dgemm('N','N',AO_num,AO_num,AO_num,                           &
      -1.d0,                                                         &
      scratch,Size(scratch,1),                                       &
      Fock_Matrix_AO,Size(Fock_Matrix_AO,1),                         &
      1.d0,                                                          &
      FPS_SPF_Matrix_AO,Size(FPS_SPF_Matrix_AO,1))

END_PROVIDER

BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_MO, (mo_num, mo_num)]
  implicit none
  begin_doc
!   Commutator FPS - SPF in MO basis
  end_doc
  call ao_to_mo(FPS_SPF_Matrix_AO, size(FPS_SPF_Matrix_AO,1), &
     FPS_SPF_Matrix_MO, size(FPS_SPF_Matrix_MO,1))
END_PROVIDER


 BEGIN_PROVIDER [ double precision, eigenvalues_Fock_matrix_AO, (AO_num) ]
&BEGIN_PROVIDER [ double precision, eigenvectors_Fock_matrix_AO, (AO_num,AO_num) ]

   BEGIN_DOC
   ! Eigenvalues and eigenvectors of the Fock matrix over the AO basis
   END_DOC

   implicit none

   double precision, allocatable  :: scratch(:,:),work(:),Xt(:,:)
   integer                        :: lwork,info
   integer                        :: i,j

   lwork = 3*AO_num - 1
   allocate(                                                         &
       scratch(AO_num,AO_num),                                 &
       work(lwork),                                                  &
       Xt(AO_num,AO_num)                                             &
       )

! Calculate Xt

  do i=1,AO_num
    do j=1,AO_num
      Xt(i,j) = S_half_inv(j,i)
    enddo
  enddo

! Calculate Fock matrix in orthogonal basis: F' = Xt.F.X

  call dgemm('N','N',AO_num,AO_num,AO_num,     &
       1.d0,                                   &
       Fock_matrix_AO,size(Fock_matrix_AO,1),  &
       S_half_inv,size(S_half_inv,1),        &
       0.d0,                                   &
       eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1))

  call dgemm('N','N',AO_num,AO_num,AO_num,                              &
       1.d0,                                                            &
       Xt,size(Xt,1),                                                   &
       eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1), &
       0.d0,                                                            &
       scratch,size(scratch,1))

! Diagonalize F' to obtain eigenvectors in orthogonal basis C' and eigenvalues

   call dsyev('V','U',AO_num,       &
        scratch,size(scratch,1),    &
        eigenvalues_Fock_matrix_AO, &
        work,lwork,info)

   if(info /= 0) then
     print *,  irp_here//' failed : ', info
     stop 1
   endif

! Back-transform eigenvectors: C =X.C'

  call dgemm('N','N',AO_num,AO_num,AO_num,     &
       1.d0,                                   &
       S_half_inv,size(S_half_inv,1),        &
       scratch,size(scratch,1),                &
       0.d0,                                   &
       eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1))

END_PROVIDER

