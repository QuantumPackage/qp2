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
       F(i,j) = Fock_matrix_mo_complex(i,j)
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

