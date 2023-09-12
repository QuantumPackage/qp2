
! ---

BEGIN_PROVIDER [double precision, noL_0e]

  implicit none
  integer                       :: i, j, k
  double precision              :: I_ijk_ijk, I_ijk_kij, I_ijk_jik, I_ijk_jki, I_ijk_ikj, I_ijk_kji
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:)

  call wall_time(t0)
  print*, " Providing noL_0e ..."

  if(elec_alpha_num .eq. elec_beta_num) then

    allocate(tmp(elec_beta_num))

    !$OMP PARALLEL                                 &
    !$OMP DEFAULT (NONE)                           &
    !$OMP PRIVATE (i, j, k,                        &
    !$OMP         I_ijk_ijk, I_ijk_kij, I_ijk_jik) &
    !$OMP SHARED (elec_beta_num, tmp)

    !$OMP DO
    do i = 1, elec_beta_num

      tmp(i) = 0.d0
      do j = 1, elec_beta_num
        do k = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, k, i, j, I_ijk_kij)
          call give_integrals_3_body_bi_ort(i, j, k, j, i, k, I_ijk_jik)

          tmp(i) = tmp(i) + 4.d0 * (2.d0 * I_ijk_ijk + I_ijk_kij - 3.d0 * I_ijk_jik)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    noL_0e = -1.d0 * (-sum(tmp)) / 6.d0

    deallocate(tmp)

  else

    allocate(tmp(elec_alpha_num))

    !$OMP PARALLEL                                 &
    !$OMP DEFAULT (NONE)                           &
    !$OMP PRIVATE (i, j, k,                        &
    !$OMP         I_ijk_ijk, I_ijk_kij, I_ijk_jik, &
    !$OMP         I_ijk_jki, I_ijk_ikj, I_ijk_kji) &
    !$OMP SHARED (elec_beta_num, elec_alpha_num, tmp)

    !$OMP DO
    do i = 1, elec_beta_num

      tmp(i) = 0.d0
      do j = 1, elec_beta_num
        do k = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, k, i, j, I_ijk_kij)
          call give_integrals_3_body_bi_ort(i, j, k, j, i, k, I_ijk_jik)

          tmp(i) = tmp(i) + 4.d0 * (2.d0 * I_ijk_ijk + I_ijk_kij - 3.d0 * I_ijk_jik)
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP DO
    do i = elec_beta_num+1, elec_alpha_num

      tmp(i) = 0.d0
      do j = elec_beta_num+1, elec_alpha_num
        do k = elec_beta_num+1, elec_alpha_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, k, i, j, I_ijk_kij)
          call give_integrals_3_body_bi_ort(i, j, k, j, i, k, I_ijk_jik)

          tmp(i) = tmp(i) + I_ijk_ijk + 2.d0 * I_ijk_kij - 3.d0 * I_ijk_jik
        enddo ! k
      enddo ! j

      do j = 1, elec_beta_num
        do k = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, j, k, i, I_ijk_jki)
          call give_integrals_3_body_bi_ort(i, j, k, i, k, j, I_ijk_ikj)
          call give_integrals_3_body_bi_ort(i, j, k, j, i, k, I_ijk_jik)
          call give_integrals_3_body_bi_ort(i, j, k, k, j, i, I_ijk_kji)

          tmp(i) = tmp(i) + 6.d0 * (2.d0 * I_ijk_ijk + I_ijk_jki - I_ijk_ikj - I_ijk_jik - I_ijk_kji)
        enddo ! k

        do k = elec_beta_num+1, elec_alpha_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, j, k, i, I_ijk_jki)
          call give_integrals_3_body_bi_ort(i, j, k, i, k, j, I_ijk_ikj)
          call give_integrals_3_body_bi_ort(i, j, k, j, i, k, I_ijk_jik)
          call give_integrals_3_body_bi_ort(i, j, k, k, j, i, I_ijk_kji)

          tmp(i) = tmp(i) + 3.d0 * (2.d0 * I_ijk_ijk + 2.d0 * I_ijk_jki - I_ijk_ikj - I_ijk_jik - 2.d0 * I_ijk_kji)
        enddo ! k
      enddo ! j
    enddo ! i
    !$OMP END DO
    !$OMP END PARALLEL

    noL_0e = -1.d0 * (-sum(tmp)) / 6.d0

    deallocate(tmp)

  endif

  call wall_time(t1)
  print*, " Wall time for noL_0e (min) = ", (t1 - t0)/60.d0

  print*, " noL_0e = ", noL_0e

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, noL_1e, (mo_num, mo_num)]

  BEGIN_DOC
  !
  ! x (-1) because integrals are over -L
  !
  END_DOC

  implicit none
  integer          :: p, s, i, j
  double precision :: I_pij_sij, I_pij_isj, I_pij_ijs, I_pij_sji, I_pij_jsi, I_pij_jis
  double precision :: t0, t1

  call wall_time(t0)
  print*, " Providing noL_1e ..."

  if(elec_alpha_num .eq. elec_beta_num) then

    !$OMP PARALLEL                                 &
    !$OMP DEFAULT (NONE)                           &
    !$OMP PRIVATE (p, s, i, j,                     &
    !$OMP         I_pij_sij, I_pij_isj, I_pij_ijs, &
    !$OMP         I_pij_sji)                       &
    !$OMP SHARED (mo_num, elec_beta_num, noL_1e)

    !$OMP DO COLLAPSE(2)
    do s = 1, mo_num
      do p = 1, mo_num
  
        noL_1e(p,s) = 0.d0
        do i = 1, elec_beta_num
          do j = 1, elec_beta_num

            call give_integrals_3_body_bi_ort(p, i, j, s, i, j, I_pij_sij)
            call give_integrals_3_body_bi_ort(p, i, j, i, s, j, I_pij_isj)
            call give_integrals_3_body_bi_ort(p, i, j, i, j, s, I_pij_ijs)
            call give_integrals_3_body_bi_ort(p, i, j, s, j, i, I_pij_sji)
          
            noL_1e(p,s) = noL_1e(p,s) - (2.d0*I_pij_sij - 2.d0*I_pij_isj + I_pij_ijs - I_pij_sji)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  else

    !$OMP PARALLEL                                 &
    !$OMP DEFAULT (NONE)                           &
    !$OMP PRIVATE (p, s, i, j,                     &
    !$OMP         I_pij_sij, I_pij_isj, I_pij_ijs, &
    !$OMP         I_pij_sji, I_pij_jsi, I_pij_jis) &
    !$OMP SHARED (mo_num, elec_beta_num, elec_alpha_num, noL_1e)

    !$OMP DO COLLAPSE(2)
    do s = 1, mo_num
      do p = 1, mo_num
  
        noL_1e(p,s) = 0.d0
        do i = 1, elec_beta_num
          do j = 1, elec_beta_num

            call give_integrals_3_body_bi_ort(p, i, j, s, i, j, I_pij_sij)
            call give_integrals_3_body_bi_ort(p, i, j, i, s, j, I_pij_isj)
            call give_integrals_3_body_bi_ort(p, i, j, i, j, s, I_pij_ijs)
            call give_integrals_3_body_bi_ort(p, i, j, s, j, i, I_pij_sji)
          
            noL_1e(p,s) = noL_1e(p,s) - (2.d0*I_pij_sij - 2.d0*I_pij_isj + I_pij_ijs - I_pij_sji)
          enddo ! j
        enddo ! i

        do i = elec_beta_num+1, elec_alpha_num
          do j = 1, elec_beta_num

            call give_integrals_3_body_bi_ort(p, i, j, s, j, i, I_pij_sji)
            call give_integrals_3_body_bi_ort(p, i, j, j, s, i, I_pij_jsi)
            call give_integrals_3_body_bi_ort(p, i, j, j, i, s, I_pij_jis)
            call give_integrals_3_body_bi_ort(p, i, j, s, i, j, I_pij_sij)
            call give_integrals_3_body_bi_ort(p, i, j, i, s, j, I_pij_isj)
            call give_integrals_3_body_bi_ort(p, i, j, i, j, s, I_pij_ijs)
          
            noL_1e(p,s) = noL_1e(p,s) + 0.5d0 * (2.d0*I_pij_sji - I_pij_jsi + 2.d0*I_pij_jis - 4.d0*I_pij_sij + 2.d0*I_pij_isj - I_pij_ijs)
          enddo ! j

          do j = elec_beta_num+1, elec_alpha_num

            call give_integrals_3_body_bi_ort(p, i, j, s, i, j, I_pij_sij)
            call give_integrals_3_body_bi_ort(p, i, j, i, s, j, I_pij_isj)
            call give_integrals_3_body_bi_ort(p, i, j, i, j, s, I_pij_ijs)
            call give_integrals_3_body_bi_ort(p, i, j, s, j, i, I_pij_sji)
          
            noL_1e(p,s) = noL_1e(p,s) - 0.5d0 * (I_pij_sij - I_pij_isj + I_pij_ijs - I_pij_sji)
          enddo ! j
        enddo ! i

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

  endif

  call wall_time(t1)
  print*, " Wall time for noL_1e (min) = ", (t1 - t0)/60.d0

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, noL_2e, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! x (-1) because integrals are over -L
  !
  END_DOC

  implicit none
  integer          :: p, q, s, t, i
  double precision :: I_ipq_sit, I_ipq_tsi, I_ipq_ist
  double precision :: t0, t1

  call wall_time(t0)
  print*, " Providing noL_2e ..."

  if(elec_alpha_num .eq. elec_beta_num) then

    !$OMP PARALLEL                                  &
    !$OMP DEFAULT (NONE)                            &
    !$OMP PRIVATE (p, q, s, t, i,                   &
    !$OMP          I_ipq_sit, I_ipq_tsi, I_ipq_ist) &
    !$OMP SHARED (mo_num, elec_beta_num, noL_2e)

    !$OMP DO COLLAPSE(4)
    do t = 1, mo_num
      do s = 1, mo_num
        do q = 1, mo_num
          do p = 1, mo_num
  
            noL_2e(p,q,s,t) = 0.d0
            do i = 1, elec_beta_num

              call give_integrals_3_body_bi_ort(i, p, q, s, i, t, I_ipq_sit)
              call give_integrals_3_body_bi_ort(i, p, q, t, s, i, I_ipq_tsi)
              call give_integrals_3_body_bi_ort(i, p, q, i, s, t, I_ipq_ist)
          
              noL_2e(p,q,s,t) = noL_2e(p,q,s,t) - 0.5d0 * (I_ipq_sit + I_ipq_tsi - 2.d0*I_ipq_ist)
            enddo
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  else

    !$OMP PARALLEL                                  &
    !$OMP DEFAULT (NONE)                            &
    !$OMP PRIVATE (p, q, s, t, i,                   &
    !$OMP          I_ipq_sit, I_ipq_tsi, I_ipq_ist) &
    !$OMP SHARED (mo_num, elec_beta_num, elec_alpha_num, noL_2e)

    !$OMP DO COLLAPSE(4)
    do t = 1, mo_num
      do s = 1, mo_num
        do q = 1, mo_num
          do p = 1, mo_num
  
            noL_2e(p,q,s,t) = 0.d0
            do i = 1, elec_beta_num

              call give_integrals_3_body_bi_ort(i, p, q, s, i, t, I_ipq_sit)
              call give_integrals_3_body_bi_ort(i, p, q, t, s, i, I_ipq_tsi)
              call give_integrals_3_body_bi_ort(i, p, q, i, s, t, I_ipq_ist)
            
              noL_2e(p,q,s,t) = noL_2e(p,q,s,t) - 0.5d0 * (I_ipq_sit + I_ipq_tsi - 2.d0*I_ipq_ist)
            enddo ! i

            do i = elec_beta_num+1, elec_alpha_num

              call give_integrals_3_body_bi_ort(i, p, q, s, i, t, I_ipq_sit)
              call give_integrals_3_body_bi_ort(i, p, q, t, s, i, I_ipq_tsi)
              call give_integrals_3_body_bi_ort(i, p, q, i, s, t, I_ipq_ist)
            
              noL_2e(p,q,s,t) = noL_2e(p,q,s,t) - 0.25d0 * (I_ipq_sit + I_ipq_tsi - 2.d0*I_ipq_ist)
            enddo ! i

          enddo ! p
        enddo ! q
      enddo ! s
    enddo ! t
    !$OMP END DO
    !$OMP END PARALLEL

  endif

  call wall_time(t1)
  print*, " Wall time for noL_2e (min) = ", (t1 - t0)/60.d0

END_PROVIDER

! ---

