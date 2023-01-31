
! ---

BEGIN_PROVIDER [double precision, TCSCF_bi_ort_dm_ao_alpha, (ao_num, ao_num) ]

  BEGIN_DOC
  ! TCSCF_bi_ort_dm_ao_alpha(i,j) = <Chi_0| a^dagger_i,alpha a_j,alpha |Phi_0> where i,j are AO basis. 
  !
  ! This is the equivalent of the alpha density of the HF Slater determinant, but with a couple of bi-orthonormal Slater determinant |Chi_0> and |Phi_0>
  END_DOC

  implicit none

  PROVIDE mo_l_coef mo_r_coef

  call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0               &
            , mo_l_coef, size(mo_l_coef, 1), mo_r_coef, size(mo_r_coef, 1) &
            !, mo_r_coef, size(mo_r_coef, 1), mo_l_coef, size(mo_l_coef, 1) &
            , 0.d0, TCSCF_bi_ort_dm_ao_alpha, size(TCSCF_bi_ort_dm_ao_alpha, 1) )

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, TCSCF_bi_ort_dm_ao_beta, (ao_num, ao_num) ]

  BEGIN_DOC
  ! TCSCF_bi_ort_dm_ao_beta(i,j) = <Chi_0| a^dagger_i,beta a_j,beta |Phi_0> where i,j are AO basis. 
  !
  ! This is the equivalent of the beta density of the HF Slater determinant, but with a couple of bi-orthonormal Slater determinant |Chi_0> and |Phi_0>
  END_DOC

  implicit none

  PROVIDE mo_l_coef mo_r_coef

  call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0                &
            , mo_l_coef, size(mo_l_coef, 1), mo_r_coef, size(mo_r_coef, 1) &
            !, mo_r_coef, size(mo_r_coef, 1), mo_l_coef, size(mo_l_coef, 1) &
            , 0.d0, TCSCF_bi_ort_dm_ao_beta, size(TCSCF_bi_ort_dm_ao_beta, 1) )

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, TCSCF_bi_ort_dm_ao, (ao_num, ao_num) ]

  BEGIN_DOC
  ! TCSCF_bi_ort_dm_ao(i,j) = <Chi_0| a^dagger_i,beta+alpha a_j,beta+alpha |Phi_0> where i,j are AO basis. 
  !
  ! This is the equivalent of the total electronic density of the HF Slater determinant, but with a couple of bi-orthonormal Slater determinant |Chi_0> and |Phi_0>
  END_DOC

  implicit none

  PROVIDE mo_l_coef mo_r_coef

  ASSERT(size(TCSCF_bi_ort_dm_ao, 1) == size(TCSCF_bi_ort_dm_ao_alpha, 1))

  if(elec_alpha_num==elec_beta_num) then
    TCSCF_bi_ort_dm_ao = TCSCF_bi_ort_dm_ao_alpha + TCSCF_bi_ort_dm_ao_alpha
  else
    ASSERT(size(TCSCF_bi_ort_dm_ao, 1) == size(TCSCF_bi_ort_dm_ao_beta, 1))
    TCSCF_bi_ort_dm_ao = TCSCF_bi_ort_dm_ao_alpha + TCSCF_bi_ort_dm_ao_beta
  endif

END_PROVIDER

! ---

