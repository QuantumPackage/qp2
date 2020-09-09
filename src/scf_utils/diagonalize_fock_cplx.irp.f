BEGIN_PROVIDER [ complex*16, eigenvectors_Fock_matrix_mo_complex, (ao_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Eigenvectors of the Fock matrix in the |MO| basis obtained with level shift.
   END_DOC

   integer                        :: i,j
   integer                        :: n
   complex*16, allocatable :: F(:,:)
   double precision, allocatable  :: diag(:)


   allocate( F(mo_num,mo_num) )
   allocate (diag(mo_num) )

   do j=1,mo_num
     do i=1,mo_num
       F(i,j) = fock_matrix_mo_complex(i,j)
     enddo
   enddo

   if(frozen_orb_scf)then
     integer                        :: iorb,jorb
     do i = 1, n_core_orb
      iorb = list_core(i)
      do j = 1, n_act_orb
       jorb = list_act(j)
       F(iorb,jorb) = (0.d0,0.d0)
       F(jorb,iorb) = (0.d0,0.d0)
      enddo
     enddo
   endif

   ! Insert level shift here
   do i = elec_beta_num+1, elec_alpha_num
     F(i,i) += 0.5d0*level_shift
   enddo

   do i = elec_alpha_num+1, mo_num
     F(i,i) += level_shift
   enddo

   n = mo_num
   call lapack_diagd_diag_in_place_complex(diag,F,n,n)

   call zgemm('N','N',ao_num,mo_num,mo_num, (1.d0,0.d0),            &
       mo_coef_complex, size(mo_coef_complex,1), F, size(F,1),                       &
       (0.d0,0.d0), eigenvectors_Fock_matrix_mo_complex, size(eigenvectors_Fock_matrix_mo_complex,1))
   deallocate(F, diag)


END_PROVIDER

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!
BEGIN_PROVIDER [ complex*16, eigenvectors_Fock_matrix_mo_kpts, (ao_num_per_kpt,mo_num_per_kpt,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Eigenvectors of the Fock matrix in the |MO| basis obtained with level shift.
  END_DOC

  integer                        :: i,j,k
  integer                        :: n
  complex*16, allocatable :: F(:,:)
  double precision, allocatable  :: diag(:)


  allocate( F(mo_num_per_kpt,mo_num_per_kpt) )
  allocate (diag(mo_num_per_kpt) )

  do k=1,kpt_num 
    do j=1,mo_num_per_kpt
      do i=1,mo_num_per_kpt
        !F(i,j) = fock_matrix_mo_complex(i,j)
        F(i,j) = fock_matrix_mo_kpts(i,j,k)
      enddo
    enddo
  
    if(frozen_orb_scf)then
      integer                        :: iorb,jorb
      !todo: core/act per kpt
      do i = 1, n_core_orb
       iorb = list_core(i)
       do j = 1, n_act_orb
        jorb = list_act(j)
        F(iorb,jorb) = (0.d0,0.d0)
        F(jorb,iorb) = (0.d0,0.d0)
       enddo
      enddo
    endif
  
    ! Insert level shift here
    !todo: elec per kpt
    do i = elec_beta_num_kpts(k)+1, elec_alpha_num_kpts(k)
      F(i,i) += 0.5d0*level_shift
    enddo
  
    do i = elec_alpha_num_kpts(k)+1, mo_num_per_kpt
      F(i,i) += level_shift
    enddo
  
    n = mo_num_per_kpt
    call lapack_diagd_diag_in_place_complex(diag,F,n,n)
  
    call zgemm('N','N',ao_num_per_kpt,mo_num_per_kpt,mo_num_per_kpt, (1.d0,0.d0),            &
        mo_coef_kpts(:,:,k), size(mo_coef_kpts,1), F, size(F,1),                       &
        (0.d0,0.d0), eigenvectors_Fock_matrix_mo_kpts(:,:,k), size(eigenvectors_Fock_matrix_mo_kpts,1))
  enddo
  deallocate(F, diag)


END_PROVIDER
BEGIN_PROVIDER [ complex*16, eigenvectors_Fock_matrix_mo_kpts_real, (ao_num_per_kpt,mo_num_per_kpt,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Eigenvectors of the Fock matrix in the |MO| basis obtained with level shift.
  END_DOC

  integer                        :: i,j,k
  integer                        :: n
  !complex*16, allocatable :: F(:,:)
  double precision, allocatable :: F(:,:)
  double precision, allocatable  :: diag(:), mo_coef_tmp(:,:), eigvecs_tmp(:,:)

  allocate( F(mo_num_per_kpt,mo_num_per_kpt) )
  allocate (diag(mo_num_per_kpt) )
  allocate (mo_coef_tmp(ao_num_per_kpt,mo_num_per_kpt) )
  allocate (eigvecs_tmp(ao_num_per_kpt,mo_num_per_kpt) )

  do k=1,kpt_num 
    do j=1,mo_num_per_kpt
      do i=1,mo_num_per_kpt
        !F(i,j) = fock_matrix_mo_complex(i,j)
        F(i,j) = dble(fock_matrix_mo_kpts(i,j,k))
      enddo
    enddo
  
    if(frozen_orb_scf)then
      integer                        :: iorb,jorb
      !todo: core/act per kpt
      do i = 1, n_core_orb
       iorb = list_core(i)
       do j = 1, n_act_orb
        jorb = list_act(j)
        F(iorb,jorb) = 0.d0
        F(jorb,iorb) = 0.d0
       enddo
      enddo
    endif
  
    ! Insert level shift here
    !todo: elec per kpt
    do i = elec_beta_num_kpts(k)+1, elec_alpha_num_kpts(k)
      F(i,i) += 0.5d0*level_shift
    enddo
  
    do i = elec_alpha_num_kpts(k)+1, mo_num_per_kpt
      F(i,i) += level_shift
    enddo
  
    n = mo_num_per_kpt
    call lapack_diagd_diag_in_place(diag,F,n,n)
    
    mo_coef_tmp = dble(mo_coef_kpts(:,:,k))
    call dgemm('N','N',ao_num_per_kpt,mo_num_per_kpt,mo_num_per_kpt, 1.d0,            &
        mo_coef_tmp, size(mo_coef_tmp,1), F, size(F,1),                       &
        0.d0, eigvecs_tmp, size(eigvecs_tmp,1))

    call zlacp2('X',ao_num_per_kpt,mo_num_per_kpt,eigvecs_tmp,size(eigvecs_tmp,1), &
        eigenvectors_fock_matrix_mo_kpts_real(:,:,k), size(eigenvectors_Fock_matrix_mo_kpts_real,1))

!    call zgemm('N','N',ao_num_per_kpt,mo_num_per_kpt,mo_num_per_kpt, (1.d0,0.d0),            &
!        mo_coef_kpts(:,:,k), size(mo_coef_kpts,1), F, size(F,1),                       &
!        (0.d0,0.d0), eigenvectors_Fock_matrix_mo_kpts(:,:,k), size(eigenvectors_Fock_matrix_mo_kpts,1))
  enddo
  deallocate(F, diag,mo_coef_tmp,eigvecs_tmp)


END_PROVIDER
