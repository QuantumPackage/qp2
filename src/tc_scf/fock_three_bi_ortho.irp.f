
! ---

BEGIN_PROVIDER [ double precision, fock_3_mat_a_op_sh_bi_orth_old, (mo_num, mo_num)]

  BEGIN_DOC
  ! Fock matrix for opposite spin contribution for bi ortho 
  END_DOC
 
  implicit none
  integer          :: j, m, i, a
  double precision :: direct_int, exch_int 
 
  fock_3_mat_a_op_sh_bi_orth_old = 0.d0 

  do i = 1, mo_num ! alpha  single excitation 
    do a = 1, mo_num ! alpha  single excitation 

      ! ---

      do j = 1, elec_beta_num 
        do m = 1, elec_beta_num
          call  give_integrals_3_body_bi_ort(a, m, j, i, m, j, direct_int)
          fock_3_mat_a_op_sh_bi_orth_old(a,i) += 1.d0 * direct_int
          call  give_integrals_3_body_bi_ort(a, m, j, j, m, i, exch_int)
          fock_3_mat_a_op_sh_bi_orth_old(a,i) += -1.d0 * exch_int
        enddo
      enddo

      ! ---

      do j = 1, elec_beta_num ! beta 
        do m = j+1, elec_beta_num ! beta 
          call  give_integrals_3_body_bi_ort(a, m, j, i, m, j, direct_int)
          fock_3_mat_a_op_sh_bi_orth_old(a,i) += 1.d0 * direct_int
          call  give_integrals_3_body_bi_ort(a, m, j, i, j, m, exch_int)
          fock_3_mat_a_op_sh_bi_orth_old(a,i) += -1.d0 * exch_int
        enddo
      enddo

      ! ---

   enddo
  enddo
 
  fock_3_mat_a_op_sh_bi_orth_old = - fock_3_mat_a_op_sh_bi_orth_old

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, fock_3_mat_a_op_sh_bi_orth, (mo_num, mo_num)]

  BEGIN_DOC
  ! Fock matrix for opposite spin contribution for bi ortho 
  END_DOC

  implicit none
  integer          :: i, a
  double precision :: integral1, integral2, integral3

  fock_3_mat_a_op_sh_bi_orth = 0.d0 

 !$OMP PARALLEL                                        &
 !$OMP DEFAULT (NONE)                                  &
 !$OMP PRIVATE (i, a, integral1, integral2, integral3) & 
 !$OMP SHARED (mo_num, fock_3_mat_a_op_sh_bi_orth)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num ! alpha  single excitation 
    do a = 1, mo_num ! alpha  single excitation 

      call direct_term_imj_bi_ortho(a, i, integral1)
      call exch_term_jmi_bi_ortho  (a, i, integral2)
      call exch_term_ijm_bi_ortho  (a, i, integral3)

      fock_3_mat_a_op_sh_bi_orth(a,i) += 1.5d0 * integral1 - integral2 - 0.5d0 * integral3
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  fock_3_mat_a_op_sh_bi_orth = - fock_3_mat_a_op_sh_bi_orth

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, fock_3_mat_a_sa_sh_bi_orth_old, (mo_num, mo_num)]

  BEGIN_DOC
  ! Fock matrix for same spin contribution for bi ortho 
  END_DOC

  implicit none
  integer          :: j, m, i, a
  double precision :: direct_int, cyclic_1, cyclic_2, non_cyclic_1, non_cyclic_2, non_cyclic_3

  fock_3_mat_a_sa_sh_bi_orth_old = 0.d0 

  do i = 1, mo_num
    do a = 1, mo_num
      do j = 1, elec_beta_num
        do m = j+1, elec_beta_num

          call give_integrals_3_body_bi_ort(a, m, j, i, m, j, direct_int)
          call give_integrals_3_body_bi_ort(a, m, j, j, i, m, cyclic_1)
          call give_integrals_3_body_bi_ort(a, m, j, m, j, i, cyclic_2)
          fock_3_mat_a_sa_sh_bi_orth_old(a,i) += direct_int + cyclic_1 + cyclic_2

          call give_integrals_3_body_bi_ort(a, m, j, j, m, i, non_cyclic_1)
          call give_integrals_3_body_bi_ort(a, m, j, i, j, m, non_cyclic_2)
          call give_integrals_3_body_bi_ort(a, m, j, m, i, j, non_cyclic_3)
          fock_3_mat_a_sa_sh_bi_orth_old(a,i) += -1.d0 * (non_cyclic_1 + non_cyclic_2 + non_cyclic_3)

        enddo
      enddo
    enddo
  enddo

  fock_3_mat_a_sa_sh_bi_orth_old = -fock_3_mat_a_sa_sh_bi_orth_old

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, fock_3_mat_a_sa_sh_bi_orth, (mo_num, mo_num)]

  BEGIN_DOC
  ! Fock matrix for same spin contribution for bi ortho 
  END_DOC

  implicit none
  integer          :: j, m, i, a
  double precision :: integral1, integral2, integral3, integral4

  fock_3_mat_a_sa_sh_bi_orth = 0.d0 

 !$OMP PARALLEL                                                   &
 !$OMP DEFAULT (NONE)                                             &
 !$OMP PRIVATE (i, a, integral1, integral2, integral3, integral4) & 
 !$OMP SHARED (mo_num, fock_3_mat_a_sa_sh_bi_orth)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do a = 1, mo_num
      call direct_term_imj_bi_ortho(a, i, integral1)
      call cyclic_term_jim_bi_ortho(a, i, integral2)
      call exch_term_jmi_bi_ortho  (a, i, integral3)
      call exch_term_ijm_bi_ortho  (a, i, integral4)
      fock_3_mat_a_sa_sh_bi_orth(a,i) += 0.5d0 * (integral1 - integral4) + integral2 - integral3
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  fock_3_mat_a_sa_sh_bi_orth = -fock_3_mat_a_sa_sh_bi_orth

END_PROVIDER 

! ---

