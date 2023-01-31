! ---

BEGIN_PROVIDER [ double precision, threshold_DIIS_nonzero_TCSCF ]
 
  implicit none
 
  if(threshold_DIIS_TCSCF == 0.d0) then
    threshold_DIIS_nonzero_TCSCF = dsqrt(thresh_tcscf)
  else
    threshold_DIIS_nonzero_TCSCF = threshold_DIIS_TCSCF
  endif
  ASSERT(threshold_DIIS_nonzero_TCSCF >= 0.d0)

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, Q_alpha, (ao_num, ao_num) ]

  BEGIN_DOC
  !
  ! Q_alpha = mo_r_coef x eta_occ_alpha x mo_l_coef.T
  !
  ! [Q_alpha]_ij = \sum_{k=1}^{elec_alpha_num} [mo_r_coef]_ik [mo_l_coef]_jk
  !
  END_DOC

  implicit none

  Q_alpha = 0.d0
  call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0               &
            , mo_r_coef, size(mo_r_coef, 1), mo_l_coef, size(mo_l_coef, 1) &
            , 0.d0, Q_alpha, size(Q_alpha, 1) )

END_PROVIDER
  
! ---
    
BEGIN_PROVIDER [ double precision, Q_beta, (ao_num, ao_num) ]

  BEGIN_DOC
  !
  ! Q_beta = mo_r_coef x eta_occ_beta x mo_l_coef.T
  !
  ! [Q_beta]_ij = \sum_{k=1}^{elec_beta_num} [mo_r_coef]_ik [mo_l_coef]_jk
  !
  END_DOC

  implicit none

  Q_beta = 0.d0
  call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0                &
            , mo_r_coef, size(mo_r_coef, 1), mo_l_coef, size(mo_l_coef, 1) &
            , 0.d0, Q_beta, size(Q_beta, 1) )

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Q_matrix, (ao_num, ao_num) ]

  BEGIN_DOC
  !
  ! Q_matrix = 2 mo_r_coef x eta_occ x mo_l_coef.T
  ! 
  ! with: 
  !                        | 1   if i = j = 1, ..., nb of occ orbitals
  !        [eta_occ]_ij =  |     
  !                        | 0   otherwise
  !
  ! the diis error is defines as:
  !                         e = F_ao x Q x ao_overlap - ao_overlap x Q x F_ao
  ! with: 
  !       mo_l_coef.T x ao_overlap x mo_r_coef = I
  !       F_mo = mo_l_coef.T x F_ao x mo_r_coef
  !       F_ao = (ao_overlap x mo_r_coef) x F_mo x (ao_overlap x mo_l_coef).T
  !
  ! ==> e = 2 ao_overlap x mo_r_coef x [ F_mo x eta_occ - eta_occ x F_mo ] x (ao_overlap x mo_l_coef).T
  !
  !      at convergence:
  !                                      F_mo x eta_occ - eta_occ x F_mo = 0
  !                                  ==> [F_mo]_ij ([eta_occ]_ii - [eta_occ]_jj) = 0  
  !                                  ==> [F_mo]_ia = [F_mo]_ai = 0 where: i = occ and a = vir
  !                                  ==> Brillouin conditions
  !
  END_DOC

  implicit none

  if(elec_alpha_num == elec_beta_num) then
    Q_matrix = Q_alpha + Q_alpha
  else
    Q_matrix = Q_alpha + Q_beta
  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, FQS_SQF_ao, (ao_num, ao_num)]

  implicit none
  double precision, allocatable :: tmp(:,:)

  allocate(tmp(ao_num,ao_num))

  ! F x Q
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, 1.d0                                             &
            , Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot, 1), Q_matrix, size(Q_matrix, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  ! F x Q x S
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, 1.d0             &
            , tmp, size(tmp, 1), ao_overlap, size(ao_overlap, 1) &
            , 0.d0, FQS_SQF_ao, size(FQS_SQF_ao, 1) )

  ! S x Q
  tmp = 0.d0
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, 1.d0                       &
            , ao_overlap, size(ao_overlap, 1), Q_matrix, size(Q_matrix, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  ! F x Q x S - S x Q x F
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, -1.d0                                  &
            , tmp, size(tmp, 1), Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot, 1) &
            , 1.d0, FQS_SQF_ao, size(FQS_SQF_ao, 1) )

  deallocate(tmp)

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, FQS_SQF_mo, (mo_num, mo_num)]

  implicit none

  call ao_to_mo_bi_ortho( FQS_SQF_ao, size(FQS_SQF_ao, 1) &
                        , FQS_SQF_mo, size(FQS_SQF_mo, 1) )

END_PROVIDER

! ---

! BEGIN_PROVIDER [ double precision, eigenval_Fock_tc_ao, (ao_num) ]
!&BEGIN_PROVIDER [ double precision, eigenvec_Fock_tc_ao, (ao_num,ao_num) ]
!
!  BEGIN_DOC
!  !
!  ! Eigenvalues and eigenvectors of the Fock matrix over the ao basis
!  !
!  ! F' = X.T x F x X   where X = ao_overlap^(-1/2)
!  ! 
!  ! F'   x Cr' = Cr' x E ==> F   Cr = Cr x E with Cr = X x Cr'
!  ! F'.T x Cl' = Cl' x E ==> F.T Cl = Cl x E with Cl = X x Cl'
!  !
!  END_DOC
!
!  implicit none
!  double precision, allocatable :: tmp1(:,:), tmp2(:,:)
!
!  ! ---
!  ! Fock matrix in orthogonal basis: F' = X.T x F x X
!
!  allocate(tmp1(ao_num,ao_num))
!  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, 1.d0                                                 &
!            , Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot, 1), S_half_inv, size(S_half_inv, 1) &
!            , 0.d0, tmp1, size(tmp1, 1) )
!
!  allocate(tmp2(ao_num,ao_num))
!  call dgemm( 'T', 'N', ao_num, ao_num, ao_num, 1.d0               &
!            , S_half_inv, size(S_half_inv, 1), tmp1, size(tmp1, 1) &
!            , 0.d0, tmp2, size(tmp2, 1) )
!
!  ! ---
!
!  ! Diagonalize F' to obtain eigenvectors in orthogonal basis C' and eigenvalues
!  ! TODO
!
!  ! Back-transform eigenvectors: C =X.C'
!
!END_PROVIDER

! ---

~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
