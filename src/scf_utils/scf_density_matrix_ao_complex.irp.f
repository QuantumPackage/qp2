BEGIN_PROVIDER [complex*16, SCF_density_matrix_ao_alpha_complex, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\alpha$ MOs
   END_DOC

   call zgemm('N','C',ao_num,ao_num,elec_alpha_num,(1.d0,0.d0), &
        mo_coef_complex, size(mo_coef_complex,1), &
        mo_coef_complex, size(mo_coef_complex,1), (0.d0,0.d0), &
        SCF_density_matrix_ao_alpha_complex, size(SCF_density_matrix_ao_alpha_complex,1))

END_PROVIDER

BEGIN_PROVIDER [ complex*16, SCF_density_matrix_ao_beta_complex,  (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\beta$ MOs
   END_DOC

   call zgemm('N','C',ao_num,ao_num,elec_beta_num,(1.d0,0.d0), &
        mo_coef_complex, size(mo_coef_complex,1), &
        mo_coef_complex, size(mo_coef_complex,1), (0.d0,0.d0), &
        SCF_density_matrix_ao_beta_complex, size(SCF_density_matrix_ao_beta_complex,1))

END_PROVIDER

BEGIN_PROVIDER [ complex*16, SCF_density_matrix_ao_complex, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! Sum of $\alpha$ and $\beta$ density matrices
   END_DOC
   ASSERT (size(SCF_density_matrix_ao_complex,1) == size(SCF_density_matrix_ao_alpha_complex,1))
   if (elec_alpha_num== elec_beta_num) then
     SCF_density_matrix_ao_complex = SCF_density_matrix_ao_alpha_complex + SCF_density_matrix_ao_alpha_complex
   else
     ASSERT (size(SCF_density_matrix_ao_complex,1) == size(SCF_density_matrix_ao_beta_complex ,1))
     SCF_density_matrix_ao_complex = SCF_density_matrix_ao_alpha_complex + SCF_density_matrix_ao_beta_complex
   endif

END_PROVIDER

