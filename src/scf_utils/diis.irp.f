! ---

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

! ---

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

! ---

BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_MO, (mo_num, mo_num)]
  implicit none
  begin_doc
!   Commutator FPS - SPF in MO basis
  end_doc
  call ao_to_mo(FPS_SPF_Matrix_AO, size(FPS_SPF_Matrix_AO,1), &
     FPS_SPF_Matrix_MO, size(FPS_SPF_Matrix_MO,1))
END_PROVIDER

! ---

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

! ---

!BEGIN_PROVIDER [double precision, error_diis_Fmo, (ao_num, ao_num)]
!
!  BEGIN_DOC
!  !
!  ! error_diis_Fmo = (S x C) x [F_mo x \eta_occ - \eta_occ x F_mo] x (S x C).T
!  !
!  ! \eta_occ is the matrix of occupation : \eta_occ = \eta_occ(alpha) + \eta_occ(beta)
!  !
!  END_DOC
!
!  implicit none
!  integer                       :: i, j
!  double precision, allocatable :: tmp(:,:)
!
!  provide Fock_matrix_mo
!
!  allocate(tmp(mo_num,mo_num))
!  tmp = 0.d0
!
!  ! F_mo x \eta_occ(alpha) - \eta_occ x F_mo(alpha)
!  do j = 1, elec_alpha_num
!    do i = elec_alpha_num + 1, mo_num
!      tmp(i,j) = Fock_matrix_mo(i,j)
!    enddo
!  enddo
!  do j = elec_alpha_num + 1, mo_num
!    do i = 1, elec_alpha_num
!      tmp(i,j) = -Fock_matrix_mo(i,j)
!    enddo
!  enddo
!
!  ! F_mo x \eta_occ(beta) - \eta_occ x F_mo(beta)
!  do j = 1, elec_beta_num
!    do i = elec_beta_num + 1, mo_num
!      tmp(i,j) += Fock_matrix_mo(i,j)
!    enddo
!  enddo
!  do j = elec_beta_num + 1, mo_num
!    do i = 1, elec_beta_num
!      tmp(i,j) -= Fock_matrix_mo(i,j)
!    enddo
!  enddo
!
!  call mo_to_ao(tmp, size(tmp, 1), error_diis_Fmo, size(error_diis_Fmo, 1))
!  
!  deallocate(tmp)
!
!END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, error_diis_Fmo, (mo_num, mo_num)]

  BEGIN_DOC
  !
  ! error_diis_Fmo = [F_mo x \eta_occ - \eta_occ x F_mo]
  !
  ! \eta_occ is the matrix of occupation : \eta_occ = \eta_occ(alpha) + \eta_occ(beta)
  !
  END_DOC

  implicit none
  integer                       :: i, j
  double precision, allocatable :: tmp(:,:)

  provide Fock_matrix_mo

  error_diis_Fmo = 0.d0

  ! F_mo x \eta_occ(alpha) - \eta_occ x F_mo(alpha)
  do j = 1, elec_alpha_num
    do i = elec_alpha_num + 1, mo_num
      error_diis_Fmo(i,j) += Fock_matrix_mo(i,j)
    enddo
  enddo
  do j = elec_alpha_num + 1, mo_num
    do i = 1, elec_alpha_num
      error_diis_Fmo(i,j) -= Fock_matrix_mo(i,j)
    enddo
  enddo

  ! F_mo x \eta_occ(beta) - \eta_occ x F_mo(beta)
  do j = 1, elec_beta_num
    do i = elec_beta_num + 1, mo_num
      error_diis_Fmo(i,j) += Fock_matrix_mo(i,j)
    enddo
  enddo
  do j = elec_beta_num + 1, mo_num
    do i = 1, elec_beta_num
      error_diis_Fmo(i,j) -= Fock_matrix_mo(i,j)
    enddo
  enddo

  !allocate(tmp(ao_num,ao_num))
  !call mo_to_ao(error_diis_Fmo, size(error_diis_Fmo, 1), tmp, size(tmp, 1))
  !call ao_to_mo(tmp, size(tmp, 1), error_diis_Fmo, size(error_diis_Fmo, 1))
  !deallocate(tmp)

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_AO_a, (AO_num, AO_num)]

  implicit none
  double precision, allocatable  :: scratch(:,:)

  allocate(scratch(AO_num, AO_num))

  call dgemm( 'N', 'N', AO_num, AO_num, AO_num, 1.d0                                                                                 &
            , Fock_Matrix_AO_alpha, size(Fock_Matrix_AO_alpha, 1), SCF_density_matrix_ao_alpha, size(SCF_Density_Matrix_AO_alpha, 1) &
            , 0.d0, scratch, size(scratch, 1) )

  call dgemm( 'N', 'N', AO_num, AO_num, AO_num, 1.d0                     &
            , scratch, size(scratch, 1), AO_Overlap, size(AO_Overlap, 1) &
            , 0.d0, FPS_SPF_Matrix_AO_a, size(FPS_SPF_Matrix_AO_a, 1) )

  call dgemm( 'N', 'N', AO_num, AO_num, AO_num, 1.d0                                                             &
            , AO_Overlap, size(AO_Overlap, 1), SCF_density_matrix_ao_alpha, size(SCF_density_matrix_ao_alpha, 1) & 
            , 0.d0, scratch, size(scratch, 1) )

  call dgemm( 'N', 'N', AO_num, AO_num, AO_num, -1.d0                                        &
            , scratch, size(scratch, 1), Fock_Matrix_AO_alpha, size(Fock_Matrix_AO_alpha, 1) &
            , 1.d0, FPS_SPF_Matrix_AO_a, size(FPS_SPF_Matrix_AO_a, 1) )

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_AO_b, (AO_num, AO_num)]

  implicit none
  double precision, allocatable  :: scratch(:,:)

  allocate(scratch(AO_num, AO_num))

  call dgemm( 'N', 'N', AO_num, AO_num, AO_num, 1.d0                                                                             &
            , Fock_Matrix_AO_beta, size(Fock_Matrix_AO_beta, 1), SCF_density_matrix_ao_beta, size(SCF_Density_Matrix_AO_beta, 1) &
            , 0.d0, scratch, size(scratch, 1) )

  call dgemm( 'N', 'N', AO_num, AO_num, AO_num, 1.d0                     &
            , scratch, size(scratch, 1), AO_Overlap, size(AO_Overlap, 1) &
            , 0.d0, FPS_SPF_Matrix_AO_b, size(FPS_SPF_Matrix_AO_b, 1) )

  call dgemm( 'N', 'N', AO_num, AO_num, AO_num, 1.d0                                                           &
            , AO_Overlap, size(AO_Overlap, 1), SCF_density_matrix_ao_beta, size(SCF_density_matrix_ao_beta, 1) & 
            , 0.d0, scratch, size(scratch, 1) )

  call dgemm( 'N', 'N', AO_num, AO_num, AO_num, -1.d0                                      &
            , scratch, size(scratch, 1), Fock_Matrix_AO_beta, size(Fock_Matrix_AO_beta, 1) &
            , 1.d0, FPS_SPF_Matrix_AO_b, size(FPS_SPF_Matrix_AO_b, 1) )

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_MO_a, (mo_num, mo_num)]
  implicit none
  call ao_to_mo(FPS_SPF_Matrix_AO_a, size(FPS_SPF_Matrix_AO_a, 1), FPS_SPF_Matrix_MO_a, size(FPS_SPF_Matrix_MO_a, 1))
END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_MO_b, (mo_num, mo_num)]
  implicit none
  call ao_to_mo(FPS_SPF_Matrix_AO_b, size(FPS_SPF_Matrix_AO_b, 1), FPS_SPF_Matrix_MO_b, size(FPS_SPF_Matrix_MO_b, 1))
END_PROVIDER

! ---

