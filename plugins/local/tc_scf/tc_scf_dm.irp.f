! ---

BEGIN_PROVIDER [double precision, TCSCF_density_matrix_ao_beta, (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! TC-SCF transition density matrix on the AO basis for BETA electrons 
  !
  END_DOC

  implicit none

  PROVIDE mo_l_coef mo_r_coef
  TCSCF_density_matrix_ao_beta = TCSCF_bi_ort_dm_ao_beta

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, TCSCF_density_matrix_ao_alpha, (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! TC-SCF transition density matrix on the AO basis for ALPHA electrons 
  !
  END_DOC

  implicit none

  PROVIDE mo_l_coef mo_r_coef
  TCSCF_density_matrix_ao_alpha = TCSCF_bi_ort_dm_ao_alpha

END_PROVIDER 


! ---

BEGIN_PROVIDER [double precision, TCSCF_density_matrix_ao_tot, (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! TC-SCF transition density matrix on the AO basis for ALPHA+BETA electrons 
  !
  END_DOC

  implicit none

  TCSCF_density_matrix_ao_tot = TCSCF_density_matrix_ao_beta + TCSCF_density_matrix_ao_alpha

END_PROVIDER 


