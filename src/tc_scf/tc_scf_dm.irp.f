! ---

BEGIN_PROVIDER [ double precision, TCSCF_density_matrix_ao_beta, (ao_num, ao_num) ]

  implicit none

  if(bi_ortho) then
    PROVIDE mo_l_coef mo_r_coef
    TCSCF_density_matrix_ao_beta = TCSCF_bi_ort_dm_ao_beta
  else
    TCSCF_density_matrix_ao_beta = SCF_density_matrix_ao_beta
  endif
END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, TCSCF_density_matrix_ao_alpha, (ao_num, ao_num) ]

  implicit none

  if(bi_ortho) then
    PROVIDE mo_l_coef mo_r_coef
    TCSCF_density_matrix_ao_alpha = TCSCF_bi_ort_dm_ao_alpha
  else
    TCSCF_density_matrix_ao_alpha = SCF_density_matrix_ao_alpha
  endif
END_PROVIDER 


! ---

BEGIN_PROVIDER [ double precision, TCSCF_density_matrix_ao_tot, (ao_num, ao_num) ]
  implicit none
  TCSCF_density_matrix_ao_tot = TCSCF_density_matrix_ao_beta + TCSCF_density_matrix_ao_alpha
END_PROVIDER 


