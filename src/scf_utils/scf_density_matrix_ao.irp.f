BEGIN_PROVIDER [double precision, SCF_density_matrix_ao_alpha, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\alpha$ MOs
   END_DOC
   if(elec_alpha_num > 0)then
     call dgemm('N','T',ao_num,ao_num,elec_alpha_num,1.d0, &
          mo_coef, size(mo_coef,1), &
          mo_coef, size(mo_coef,1), 0.d0, &
          SCF_density_matrix_ao_alpha, size(SCF_density_matrix_ao_alpha,1))
   else
     SCF_density_matrix_ao_alpha = 0.d0
   endif


END_PROVIDER

BEGIN_PROVIDER [ double precision, SCF_density_matrix_ao_beta,  (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\beta$ MOs
   END_DOC
   if(elec_beta_num > 0)then
     call dgemm('N','T',ao_num,ao_num,elec_beta_num,1.d0, &
          mo_coef, size(mo_coef,1), &
          mo_coef, size(mo_coef,1), 0.d0, &
          SCF_density_matrix_ao_beta, size(SCF_density_matrix_ao_beta,1))
   else
     SCF_density_matrix_ao_beta = 0.d0
   endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, SCF_density_matrix_ao, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! Sum of $\alpha$ and $\beta$ density matrices
   END_DOC
   ASSERT (size(SCF_density_matrix_ao,1) == size(SCF_density_matrix_ao_alpha,1))
   if (elec_alpha_num== elec_beta_num) then
     SCF_density_matrix_ao = SCF_density_matrix_ao_alpha + SCF_density_matrix_ao_alpha
   else
     ASSERT (size(SCF_density_matrix_ao,1) == size(SCF_density_matrix_ao_beta ,1))
     SCF_density_matrix_ao = SCF_density_matrix_ao_alpha + SCF_density_matrix_ao_beta
   endif

END_PROVIDER

