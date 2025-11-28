BEGIN_PROVIDER [ logical, all_shells_closed ]
 implicit none
 BEGIN_DOC
 ! 
 END_DOC
 all_shells_closed = (elec_alpha_num == elec_beta_num)
END_PROVIDER

BEGIN_PROVIDER [double precision, SCF_density_matrix_ao_alpha, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\alpha$ MOs
   END_DOC

   if (use_optimal_damping) then
     double precision, allocatable :: mo_coef_tmp(:,:)
     allocate(mo_coef_tmp(ao_num,mo_num))
     integer :: j
     do j=1,mo_num
       mo_coef_tmp(:,j) = mo_coef(:,j) * partial_occupation(j) * 0.5d0
     enddo
     call dgemm('N','T',ao_num,ao_num,mo_num,1.d0, &
        mo_coef_tmp, size(mo_coef_tmp,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        SCF_density_matrix_ao_alpha, size(SCF_density_matrix_ao_alpha,1))
   else
     call dgemm('N','T',ao_num,ao_num,elec_alpha_num,1.d0, &
        mo_coef, size(mo_coef,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        SCF_density_matrix_ao_alpha, size(SCF_density_matrix_ao_alpha,1))
   endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, SCF_density_matrix_ao_beta,  (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\beta$ MOs
   END_DOC

   if (use_optimal_damping) then
     double precision, allocatable :: mo_coef_tmp(:,:)
     allocate(mo_coef_tmp(ao_num,mo_num))
     integer :: j
     do j=1,mo_num
       mo_coef_tmp(:,j) = mo_coef(:,j) * partial_occupation(j) * 0.5d0
     enddo
     call dgemm('N','T',ao_num,ao_num,mo_num,1.d0, &
        mo_coef_tmp, size(mo_coef_tmp,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        SCF_density_matrix_ao_beta, size(SCF_density_matrix_ao_beta,1))
   else
     call dgemm('N','T',ao_num,ao_num,elec_beta_num,1.d0, &
        mo_coef, size(mo_coef,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        SCF_density_matrix_ao_beta, size(SCF_density_matrix_ao_beta,1))
    endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, SCF_density_matrix_ao, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! Sum of $\alpha$ and $\beta$ density matrices
   END_DOC
   ASSERT (size(SCF_density_matrix_ao,1) == size(SCF_density_matrix_ao_alpha,1))
   if (all_shells_closed) then
     SCF_density_matrix_ao = SCF_density_matrix_ao_alpha + SCF_density_matrix_ao_alpha
   else
     ASSERT (size(SCF_density_matrix_ao,1) == size(SCF_density_matrix_ao_beta ,1))
     SCF_density_matrix_ao = SCF_density_matrix_ao_alpha + SCF_density_matrix_ao_beta
   endif

END_PROVIDER

