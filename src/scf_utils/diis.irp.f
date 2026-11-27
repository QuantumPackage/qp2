use sizet

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

  call dgemm('N','N',int(AO_num,4),int(AO_num,4),int(AO_num,4), &
      1.d0,                                                          &
      Fock_Matrix_AO,int(Size(Fock_Matrix_AO,1),4),                         &
      SCF_Density_Matrix_AO,int(Size(SCF_Density_Matrix_AO,1),4),             &
      0.d0,                                                          &
      scratch,int(Size(scratch,1),4))

  ! Compute FPS

  call dgemm('N','N',int(AO_num,4),int(AO_num,4),int(AO_num,4), &
      1.d0,                                                          &
      scratch,int(Size(scratch,1),4),                                       &
      AO_Overlap,int(Size(AO_Overlap,1),4),                                 &
      0.d0,                                                          &
      FPS_SPF_Matrix_AO,int(Size(FPS_SPF_Matrix_AO,1),4))

  ! Compute SP

  call dgemm('N','N',int(AO_num,4),int(AO_num,4),int(AO_num,4), &
      1.d0,                                                          &
      AO_Overlap,int(Size(AO_Overlap,1),4),                                 &
      SCF_Density_Matrix_AO,int(Size(SCF_Density_Matrix_AO,1),4),             &
      0.d0,                                                          &
      scratch,int(Size(scratch,1),4))

  ! Compute FPS - SPF

  call dgemm('N','N',int(AO_num,4),int(AO_num,4),int(AO_num,4), &
      -1.d0,                                                         &
      scratch,int(Size(scratch,1),4),                                       &
      Fock_Matrix_AO,int(Size(Fock_Matrix_AO,1),4),                         &
      1.d0,                                                          &
      FPS_SPF_Matrix_AO,int(Size(FPS_SPF_Matrix_AO,1),4))

END_PROVIDER

BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_MO, (mo_num, mo_num)]
  implicit none
  begin_doc
!   Commutator FPS - SPF in MO basis
  end_doc
  call ao_to_mo(FPS_SPF_Matrix_AO, int(size(FPS_SPF_Matrix_AO,1),size_t), &
     FPS_SPF_Matrix_MO, int(size(FPS_SPF_Matrix_MO,1),size_t))
END_PROVIDER


 BEGIN_PROVIDER [ double precision, eigenvalues_Fock_matrix_AO, (AO_num) ]
&BEGIN_PROVIDER [ double precision, eigenvectors_Fock_matrix_AO, (AO_num,AO_num) ]

   BEGIN_DOC
   ! Eigenvalues and eigenvectors of the Fock matrix over the AO basis
   END_DOC

   implicit none

   double precision, allocatable  :: scratch(:,:),work(:),Xt(:,:)
   integer                        :: lwork,info
   integer(size_t)                :: i,j

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

  call dgemm('N','N',int(AO_num,4),int(AO_num,4),int(AO_num,4),     &
       1.d0,                                   &
       Fock_matrix_AO,int(size(Fock_matrix_AO,1),4),  &
       S_half_inv,int(size(S_half_inv,1),4),        &
       0.d0,                                   &
       eigenvectors_Fock_matrix_AO,int(size(eigenvectors_Fock_matrix_AO,1),4))

  call dgemm('N','N',int(AO_num,4),int(AO_num,4),int(AO_num,4), &
       1.d0,                                                            &
       Xt,int(size(Xt,1),4),                                                   &
       eigenvectors_Fock_matrix_AO,int(size(eigenvectors_Fock_matrix_AO,1),4), &
       0.d0,                                                            &
       scratch,int(size(scratch,1),4))

! Diagonalize F' to obtain eigenvectors in orthogonal basis C' and eigenvalues

   call dsyev('V','U',int(AO_num,4),       &
        scratch,int(size(scratch,1),4),    &
        eigenvalues_Fock_matrix_AO, &
        work,lwork,info)

   if(info /= 0) then
     print *,  irp_here//' failed : ', info
     stop 1
   endif

! Back-transform eigenvectors: C =X.C'

  call dgemm('N','N',int(AO_num,4),int(AO_num,4),int(AO_num,4),     &
       1.d0,                                   &
       S_half_inv,int(size(S_half_inv,1),4),        &
       scratch,int(size(scratch,1),4),                &
       0.d0,                                   &
       eigenvectors_Fock_matrix_AO,int(size(eigenvectors_Fock_matrix_AO,1),4))

END_PROVIDER

