
! ---

BEGIN_PROVIDER [double precision, TCSCF_bi_ort_dm_ao_alpha, (ao_num, ao_num) ]

  BEGIN_DOC
  ! TCSCF_bi_ort_dm_ao_alpha(i,j) = <Chi_0| a^dagger_i,alpha a_j,alpha |Phi_0> where i,j are AO basis. 
  !
  ! This is the equivalent of the alpha density of the HF Slater determinant, but with a couple of bi-orthonormal Slater determinant |Chi_0> and |Phi_0>
  ! 
  ! WARNING : if core_tc_op then only the valence orbitals are considered
  END_DOC

  implicit none

  PROVIDE mo_l_coef mo_r_coef
  if(.not.core_tc_op)then
   call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0               &
             , mo_l_coef, size(mo_l_coef, 1), mo_r_coef, size(mo_r_coef, 1) &
             , 0.d0, TCSCF_bi_ort_dm_ao_alpha, size(TCSCF_bi_ort_dm_ao_alpha, 1) )
  else ! take only the valence orbitals, which are not the core orbitals 
   TCSCF_bi_ort_dm_ao_alpha =  TCSCF_bi_ort_valence_dm_ao_alpha
  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, TCSCF_bi_ort_dm_ao_beta, (ao_num, ao_num) ]

  BEGIN_DOC
  ! TCSCF_bi_ort_dm_ao_beta(i,j) = <Chi_0| a^dagger_i,beta a_j,beta |Phi_0> where i,j are AO basis. 
  !
  ! This is the equivalent of the beta density of the HF Slater determinant, but with a couple of bi-orthonormal Slater determinant |Chi_0> and |Phi_0>
  ! 
  ! WARNING : if core_tc_op then only the valence orbitals are considered
  END_DOC

  implicit none

  PROVIDE mo_l_coef mo_r_coef

  if(.not.core_tc_op)then
   call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0                &
             , mo_l_coef, size(mo_l_coef, 1), mo_r_coef, size(mo_r_coef, 1) &
             , 0.d0, TCSCF_bi_ort_dm_ao_beta, size(TCSCF_bi_ort_dm_ao_beta, 1) )
  else ! take only the valence orbitals, which are not the core orbitals 
   TCSCF_bi_ort_dm_ao_beta =  TCSCF_bi_ort_valence_dm_ao_beta
  endif

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


BEGIN_PROVIDER [ double precision, TCSCF_bi_ort_core_dm_ao, (ao_num, ao_num)]
 BEGIN_DOC
  ! 
  ! core density matrix = density matrix built with only the doubly occupied core orbitals
  ! 
  ! contains both alpha and beta electrons 
  ! 
 END_DOC
 PROVIDE mo_r_core_coef mo_l_core_coef
 call dgemm( 'N', 'T', ao_num, ao_num, n_core_orb, 1.d0                &
           , mo_l_core_coef, size(mo_l_core_coef, 1), mo_r_core_coef, size(mo_r_core_coef, 1) &
           , 0.d0, TCSCF_bi_ort_core_dm_ao, size(TCSCF_bi_ort_core_dm_ao, 1) )
 ! factor two because it is the total DM 
 TCSCF_bi_ort_core_dm_ao = 2.d0 * TCSCF_bi_ort_core_dm_ao
END_PROVIDER 

BEGIN_PROVIDER [ double precision, TCSCF_bi_ort_valence_dm_ao_beta, (ao_num, ao_num)]
 BEGIN_DOC
  ! 
  ! VALENCE ONLY beta density matrix = density matrix built with only the NOT CORE BETA occupied orbitals
  ! 
 END_DOC
 PROVIDE mo_r_valence_coef mo_l_valence_coef
   call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0                &
             , mo_l_valence_coef, size(mo_l_valence_coef, 1), mo_r_valence_coef, size(mo_r_valence_coef, 1) &
             , 0.d0, TCSCF_bi_ort_valence_dm_ao_beta, size(TCSCF_bi_ort_valence_dm_ao_beta, 1) )
END_PROVIDER 

BEGIN_PROVIDER [ double precision, TCSCF_bi_ort_valence_dm_ao_alpha, (ao_num, ao_num)]
 BEGIN_DOC
  ! 
  ! VALENCE ONLY alpha density matrix = density matrix built with only the NOT CORE ALPHA occupied orbitals
  ! 
 END_DOC
   PROVIDE mo_r_valence_coef mo_l_valence_coef
   call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0               &
             , mo_l_valence_coef, size(mo_l_valence_coef, 1), mo_r_valence_coef, size(mo_r_valence_coef, 1) &
             , 0.d0, TCSCF_bi_ort_valence_dm_ao_alpha, size(TCSCF_bi_ort_valence_dm_ao_alpha, 1) )
END_PROVIDER 

BEGIN_PROVIDER [ double precision, TCSCF_bi_ort_valence_dm_ao_full, (ao_num, ao_num)]
 BEGIN_DOC
  ! 
  ! VALENCE ONLY alpha+beta density matrix = density matrix built with only the NOT CORE ALPHA+BETA occupied orbitals
  ! 
 END_DOC
 TCSCF_bi_ort_valence_dm_ao_full = TCSCF_bi_ort_valence_dm_ao_alpha + TCSCF_bi_ort_valence_dm_ao_beta

END_PROVIDER 
