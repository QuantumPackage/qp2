BEGIN_PROVIDER [ complex*16, scf_density_matrix_ao_alpha_complex, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\alpha$ MOs
   END_DOC
   
   complex*16, allocatable :: mo_coef_alpha_tmp(:,:)
   integer :: occ(N_int*bit_kind_size)
   integer :: na, i

   call bitstring_to_list(hf_bitmask(1,1), occ, na, n_int)
   allocate(mo_coef_alpha_tmp(ao_num,na))
   do i=1,na
     mo_coef_alpha_tmp(:,i) = mo_coef_complex(:,occ(i))
   enddo


   call zgemm('N','C',ao_num,ao_num,elec_alpha_num,(1.d0,0.d0), &
        mo_coef_alpha_tmp, size(mo_coef_alpha_tmp,1), &
        mo_coef_alpha_tmp, size(mo_coef_alpha_tmp,1), (0.d0,0.d0), &
        scf_density_matrix_ao_alpha_complex, size(scf_density_matrix_ao_alpha_complex,1))

   deallocate(mo_coef_alpha_tmp)
   !call zgemm('N','C',ao_num,ao_num,elec_alpha_num,(1.d0,0.d0), &
   !     mo_coef_complex, size(mo_coef_complex,1), &
   !     mo_coef_complex, size(mo_coef_complex,1), (0.d0,0.d0), &
   !     scf_density_matrix_ao_alpha_complex, size(scf_density_matrix_ao_alpha_complex,1))

END_PROVIDER

BEGIN_PROVIDER [ complex*16, scf_density_matrix_ao_beta_complex,  (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\beta$ MOs
   END_DOC

   complex*16, allocatable :: mo_coef_beta_tmp(:,:)
   integer :: occ(N_int*bit_kind_size)
   integer :: nb, i

   call bitstring_to_list(hf_bitmask(1,2), occ, nb, n_int)
   allocate(mo_coef_beta_tmp(ao_num,nb))
   do i=1,nb
     mo_coef_beta_tmp(:,i) = mo_coef_complex(:,occ(i))
   enddo


   call zgemm('N','C',ao_num,ao_num,elec_beta_num,(1.d0,0.d0), &
        mo_coef_beta_tmp, size(mo_coef_beta_tmp,1), &
        mo_coef_beta_tmp, size(mo_coef_beta_tmp,1), (0.d0,0.d0), &
        scf_density_matrix_ao_beta_complex, size(scf_density_matrix_ao_beta_complex,1))

   deallocate(mo_coef_beta_tmp)
   !call zgemm('N','C',ao_num,ao_num,elec_beta_num,(1.d0,0.d0), &
   !     mo_coef_complex, size(mo_coef_complex,1), &
   !     mo_coef_complex, size(mo_coef_complex,1), (0.d0,0.d0), &
   !     scf_density_matrix_ao_beta_complex, size(scf_density_matrix_ao_beta_complex,1))

END_PROVIDER

BEGIN_PROVIDER [ complex*16, scf_density_matrix_ao_complex, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! Sum of $\alpha$ and $\beta$ density matrices
   END_DOC
   ASSERT (size(scf_density_matrix_ao_complex,1) == size(scf_density_matrix_ao_alpha_complex,1))
   if (elec_alpha_num== elec_beta_num) then
     scf_density_matrix_ao_complex = scf_density_matrix_ao_alpha_complex + scf_density_matrix_ao_alpha_complex
   else
     ASSERT (size(scf_density_matrix_ao_complex,1) == size(scf_density_matrix_ao_beta_complex ,1))
     scf_density_matrix_ao_complex = scf_density_matrix_ao_alpha_complex + scf_density_matrix_ao_beta_complex
   endif

END_PROVIDER

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!

BEGIN_PROVIDER [ complex*16, scf_density_matrix_ao_alpha_kpts, (ao_num_per_kpt,ao_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\alpha$ MOs
   END_DOC
   
   integer :: k
   do k=1,kpt_num
     call zgemm('N','C',ao_num_per_kpt,ao_num_per_kpt,elec_alpha_num_kpts(k),(1.d0,0.d0), &
          mo_coef_kpts(1,1,k), size(mo_coef_kpts,1), &
          mo_coef_kpts(1,1,k), size(mo_coef_kpts,1), (0.d0,0.d0), &
          scf_density_matrix_ao_alpha_kpts(1,1,k), size(scf_density_matrix_ao_alpha_kpts,1))
   enddo
END_PROVIDER

BEGIN_PROVIDER [ complex*16, scf_density_matrix_ao_beta_kpts, (ao_num_per_kpt,ao_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\beta$ MOs
   END_DOC
   
   integer :: k
   do k=1,kpt_num
     call zgemm('N','C',ao_num_per_kpt,ao_num_per_kpt,elec_beta_num_kpts(k),(1.d0,0.d0), &
          mo_coef_kpts(1,1,k), size(mo_coef_kpts,1), &
          mo_coef_kpts(1,1,k), size(mo_coef_kpts,1), (0.d0,0.d0), &
          scf_density_matrix_ao_beta_kpts(1,1,k), size(scf_density_matrix_ao_beta_kpts,1))
   enddo
END_PROVIDER

BEGIN_PROVIDER [ complex*16, scf_density_matrix_ao_kpts, (ao_num_per_kpt,ao_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! Sum of $\alpha$ and $\beta$ density matrices
   END_DOC
   ASSERT (size(scf_density_matrix_ao_kpts,1) == size(scf_density_matrix_ao_alpha_kpts,1))
   if (elec_alpha_num== elec_beta_num) then
     scf_density_matrix_ao_kpts = scf_density_matrix_ao_alpha_kpts + scf_density_matrix_ao_alpha_kpts
   else
     ASSERT (size(scf_density_matrix_ao_kpts,1) == size(scf_density_matrix_ao_beta_kpts ,1))
     scf_density_matrix_ao_kpts = scf_density_matrix_ao_alpha_kpts + scf_density_matrix_ao_beta_kpts
   endif

END_PROVIDER

