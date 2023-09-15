
! ---

BEGIN_PROVIDER [double precision, noL_0e_naive]

  implicit none
  integer                       :: ii, jj, kk
  integer                       :: i, j, k
  double precision              :: sigma_i, sigma_j, sigma_k
  double precision              :: I_ijk_ijk, I_ijk_kij, I_ijk_jki, I_ijk_jik, I_ijk_kji, I_ijk_ikj
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:)

  print*, " Providing noL_0e_naive ..."
  call wall_time(t0)

  allocate(tmp(elec_num))

  !$OMP PARALLEL                                                 &
  !$OMP DEFAULT (NONE)                                           &
  !$OMP PRIVATE (ii, i, sigma_i, jj, j, sigma_j, kk, k, sigma_k, &
  !$OMP          I_ijk_ijk, I_ijk_kij, I_ijk_jki, I_ijk_jik,     &
  !$OMP          I_ijk_kji, I_ijk_ikj)                           &
  !$OMP SHARED (elec_beta_num, elec_num, tmp)
  !$OMP DO

  do ii = 1, elec_num

    if(ii .le. elec_beta_num) then
      i       = ii 
      sigma_i = -1.d0
    else
      i       = ii - elec_beta_num
      sigma_i = +1.d0
    endif

    tmp(ii) = 0.d0

    do jj = 1, elec_num

      if(jj .le. elec_beta_num) then
        j       = jj
        sigma_j = -1.d0
      else
        j       = jj - elec_beta_num
        sigma_j = +1.d0
      endif

      do kk = 1, elec_num

        if(kk .le. elec_beta_num) then
          k       = kk
          sigma_k = -1.d0
        else
          k       = kk - elec_beta_num
          sigma_k = +1.d0
        endif

        call give_integrals_3_body_bi_ort_spin( i, sigma_i, j, sigma_j, k, sigma_k &
                                              , i, sigma_i, j, sigma_j, k, sigma_k &
                                              , I_ijk_ijk)

        call give_integrals_3_body_bi_ort_spin( i, sigma_i, j, sigma_j, k, sigma_k &
                                              , k, sigma_k, i, sigma_i, j, sigma_j &
                                              , I_ijk_kij)

        call give_integrals_3_body_bi_ort_spin( i, sigma_i, j, sigma_j, k, sigma_k &
                                              , j, sigma_j, k, sigma_k, i, sigma_i &
                                              , I_ijk_jki)

        call give_integrals_3_body_bi_ort_spin( i, sigma_i, j, sigma_j, k, sigma_k &
                                              , j, sigma_j, i, sigma_i, k, sigma_k &
                                              , I_ijk_jik)

        call give_integrals_3_body_bi_ort_spin( i, sigma_i, j, sigma_j, k, sigma_k &
                                              , k, sigma_k, j, sigma_j, i, sigma_i &
                                              , I_ijk_kji)

        call give_integrals_3_body_bi_ort_spin( i, sigma_i, j, sigma_j, k, sigma_k &
                                              , i, sigma_i, k, sigma_k, j, sigma_j &
                                              , I_ijk_ikj)


        tmp(ii) = tmp(ii) + I_ijk_ijk + I_ijk_kij + I_ijk_jki - I_ijk_jik - I_ijk_kji - I_ijk_ikj
        !       = tmp(ii) + I_ijk_ijk + 2.d0 * I_ijk_kij - 3.d0 * I_ijk_jik
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  noL_0e_naive = -1.d0 * (sum(tmp)) / 6.d0

  deallocate(tmp)

  call wall_time(t1)
  print*, " Wall time for noL_0e_naive (min) = ", (t1 - t0)/60.d0

  print*, " noL_0e_naive = ", noL_0e_naive
  
END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, noL_1e_naive, (mo_num, mo_num)]

  BEGIN_DOC
  ! 
  ! < p | H(1) | s > is dressed with noL_1e_naive(p,s)
  !
  END_DOC

  implicit none
  integer          :: ii, jj
  integer          :: i, j, p, s
  double precision :: sigma_i, sigma_j, sigma_p, sigma_s
  double precision :: I_pij_sji, I_pij_sij, I_pij_jis, I_pij_ijs, I_pij_isj, I_pij_jsi
  double precision :: t0, t1

  print*, " Providing noL_1e_naive ..."
  call wall_time(t0)

  ! ----
  ! up-up part

  sigma_p = +1.d0
  sigma_s = +1.d0

  !$OMP PARALLEL                                   &
  !$OMP DEFAULT (NONE)                             &
  !$OMP PRIVATE (ii, i, sigma_i, jj, j, sigma_j,   &
  !$OMP          I_pij_sji, I_pij_sij, I_pij_jis,  & 
  !$OMP          I_pij_ijs, I_pij_isj, I_pij_jsi ) &
  !$OMP SHARED (mo_num, elec_beta_num, elec_num,   &
  !$OMP         sigma_p, sigma_s, noL_1e_naive)

  !$OMP DO COLLAPSE (2)

  do s = 1, mo_num
    do p = 1, mo_num

      noL_1e_naive(p,s) = 0.d0
      do ii = 1, elec_num
        if(ii .le. elec_beta_num) then
          i       = ii 
          sigma_i = -1.d0
        else
          i       = ii - elec_beta_num
          sigma_i = +1.d0
        endif

        do jj = 1, elec_num
          if(jj .le. elec_beta_num) then
            j       = jj
            sigma_j = -1.d0
          else
            j       = jj - elec_beta_num
            sigma_j = +1d0
          endif

          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , s, sigma_s, j, sigma_j, i, sigma_i &
                                                , I_pij_sji)
  
          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , s, sigma_s, i, sigma_i, j, sigma_j &
                                                , I_pij_sij)
  
          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , j, sigma_j, i, sigma_i, s, sigma_s &
                                                , I_pij_jis)
  
          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , i, sigma_i, j, sigma_j, s, sigma_s &
                                                , I_pij_ijs)
  
          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , i, sigma_i, s, sigma_s, j, sigma_j &
                                                , I_pij_isj)
  
          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , j, sigma_j, s, sigma_s, i, sigma_i &
                                                , I_pij_jsi)

          ! x 0.5  because we consider 0.5 (up + down)
          noL_1e_naive(p,s) = noL_1e_naive(p,s) - 0.25d0 * (I_pij_sji - I_pij_sij + I_pij_jis - I_pij_ijs + I_pij_isj - I_pij_jsi)
        enddo ! j
      enddo ! i
    enddo ! s
  enddo ! p
  !$OMP END DO
  !$OMP END PARALLEL


  ! ----
  ! down-down part

  sigma_p = -1.d0
  sigma_s = -1.d0

  !$OMP PARALLEL                                   &
  !$OMP DEFAULT (NONE)                             &
  !$OMP PRIVATE (ii, i, sigma_i, jj, j, sigma_j,   &
  !$OMP          I_pij_sji, I_pij_sij, I_pij_jis,  & 
  !$OMP          I_pij_ijs, I_pij_isj, I_pij_jsi ) &
  !$OMP SHARED (mo_num, elec_beta_num, elec_num,   &
  !$OMP         sigma_p, sigma_s, noL_1e_naive)

  !$OMP DO COLLAPSE (2)

  do s = 1, mo_num
    do p = 1, mo_num

      do ii = 1, elec_num
        if(ii .le. elec_beta_num) then
          i       = ii 
          sigma_i = -1.d0
        else
          i       = ii - elec_beta_num
          sigma_i = +1.d0
        endif

        do jj = 1, elec_num
          if(jj .le. elec_beta_num) then
            j       = jj
            sigma_j = -1.d0
          else
            j       = jj - elec_beta_num
            sigma_j = +1d0
          endif

          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , s, sigma_s, j, sigma_j, i, sigma_i &
                                                , I_pij_sji)
  
          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , s, sigma_s, i, sigma_i, j, sigma_j &
                                                , I_pij_sij)
  
          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , j, sigma_j, i, sigma_i, s, sigma_s &
                                                , I_pij_jis)
  
          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , i, sigma_i, j, sigma_j, s, sigma_s &
                                                , I_pij_ijs)
  
          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , i, sigma_i, s, sigma_s, j, sigma_j &
                                                , I_pij_isj)
  
          call give_integrals_3_body_bi_ort_spin( p, sigma_p, i, sigma_i, j, sigma_j &
                                                , j, sigma_j, s, sigma_s, i, sigma_i &
                                                , I_pij_jsi)

          ! x 0.5  because we consider 0.5 (up + down)
          noL_1e_naive(p,s) = noL_1e_naive(p,s) - 0.25d0 * (I_pij_sji - I_pij_sij + I_pij_jis - I_pij_ijs + I_pij_isj - I_pij_jsi)
        enddo ! j
      enddo ! i
    enddo ! s
  enddo ! p
  !$OMP END DO
  !$OMP END PARALLEL

  ! ---

  call wall_time(t1)
  print*, " Wall time for noL_1e_naive (min) = ", (t1 - t0)/60.d0

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, noL_2e_naive, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  ! 
  ! < p q | H(2) | s t > is dressed with noL_2e_naive(p,q,s,t)
  !
  END_DOC

  implicit none
  integer          :: ii
  integer          :: i, p, q, s, t
  double precision :: sigma_i, sigma_p, sigma_q, sigma_s, sigma_t
  double precision :: I_ipq_ist, I_ipq_sit, I_ipq_tsi
  double precision :: t0, t1

  print*, " Providing noL_2e_naive ..."
  call wall_time(t0)

  ! ----
  ! up-up & up-up part

  sigma_p = +1.d0
  sigma_s = +1.d0
  sigma_q = +1.d0
  sigma_t = +1.d0

  !$OMP PARALLEL                                    &
  !$OMP DEFAULT (NONE)                              &
  !$OMP PRIVATE (ii, i, sigma_i, p, q, s, t,        &
  !$OMP          I_ipq_ist, I_ipq_sit, I_ipq_tsi)   &
  !$OMP SHARED (mo_num, elec_beta_num, elec_num,    &
  !$OMP         sigma_p, sigma_q, sigma_s, sigma_t, &
  !$OMP         noL_2e_naive)

  !$OMP DO COLLAPSE (4)
  do t = 1, mo_num
    do s = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num

          noL_2e_naive(p,q,s,t) = 0.d0
          do ii = 1, elec_num
            if(ii .le. elec_beta_num) then
              i       = ii 
              sigma_i = -1.d0
            else
              i       = ii - elec_beta_num
              sigma_i = +1.d0
            endif

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , i, sigma_i, s, sigma_s, t, sigma_t &
                                                  , I_ipq_ist)

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , s, sigma_s, i, sigma_i, t, sigma_t &
                                                  , I_ipq_sit)

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , t, sigma_t, s, sigma_s, i, sigma_i &
                                                  , I_ipq_tsi)

            ! x 0.25  because we consider 0.25 (up-up + up-down + down-up + down-down)
            noL_2e_naive(p,q,s,t) = noL_2e_naive(p,q,s,t) - 0.125d0 * (I_ipq_ist - I_ipq_sit - I_ipq_tsi)
          enddo ! i
        enddo ! p
      enddo ! q
    enddo ! s
  enddo ! t
  !$OMP END DO
  !$OMP END PARALLEL

  ! ----
  ! up-up & down-down part

  sigma_p = +1.d0
  sigma_s = +1.d0
  sigma_q = -1.d0
  sigma_t = -1.d0

  !$OMP PARALLEL                                    &
  !$OMP DEFAULT (NONE)                              &
  !$OMP PRIVATE (ii, i, sigma_i, p, q, s, t,        &
  !$OMP          I_ipq_ist, I_ipq_sit, I_ipq_tsi)   &
  !$OMP SHARED (mo_num, elec_beta_num, elec_num,    &
  !$OMP         sigma_p, sigma_q, sigma_s, sigma_t, &
  !$OMP         noL_2e_naive)

  !$OMP DO COLLAPSE (4)
  do t = 1, mo_num
    do s = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num

          do ii = 1, elec_num
            if(ii .le. elec_beta_num) then
              i       = ii 
              sigma_i = -1.d0
            else
              i       = ii - elec_beta_num
              sigma_i = +1.d0
            endif

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , i, sigma_i, s, sigma_s, t, sigma_t &
                                                  , I_ipq_ist)

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , s, sigma_s, i, sigma_i, t, sigma_t &
                                                  , I_ipq_sit)

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , t, sigma_t, s, sigma_s, i, sigma_i &
                                                  , I_ipq_tsi)

            ! x 0.25  because we consider 0.25 (up-up + up-down + down-up + down-down)
            noL_2e_naive(p,q,s,t) = noL_2e_naive(p,q,s,t) - 0.125d0 * (I_ipq_ist - I_ipq_sit - I_ipq_tsi)
          enddo ! i
        enddo ! p
      enddo ! q
    enddo ! s
  enddo ! t
  !$OMP END DO
  !$OMP END PARALLEL

  ! ----
  ! down-down & up-up part

  sigma_p = -1.d0
  sigma_s = -1.d0
  sigma_q = +1.d0
  sigma_t = +1.d0

  !$OMP PARALLEL                                    &
  !$OMP DEFAULT (NONE)                              &
  !$OMP PRIVATE (ii, i, sigma_i, p, q, s, t,        &
  !$OMP          I_ipq_ist, I_ipq_sit, I_ipq_tsi)   &
  !$OMP SHARED (mo_num, elec_beta_num, elec_num,    &
  !$OMP         sigma_p, sigma_q, sigma_s, sigma_t, &
  !$OMP         noL_2e_naive)

  !$OMP DO COLLAPSE (4)
  do t = 1, mo_num
    do s = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num

          do ii = 1, elec_num
            if(ii .le. elec_beta_num) then
              i       = ii 
              sigma_i = -1.d0
            else
              i       = ii - elec_beta_num
              sigma_i = +1.d0
            endif

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , i, sigma_i, s, sigma_s, t, sigma_t &
                                                  , I_ipq_ist)

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , s, sigma_s, i, sigma_i, t, sigma_t &
                                                  , I_ipq_sit)

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , t, sigma_t, s, sigma_s, i, sigma_i &
                                                  , I_ipq_tsi)

            ! x 0.25  because we consider 0.25 (up-up + up-down + down-up + down-down)
            noL_2e_naive(p,q,s,t) = noL_2e_naive(p,q,s,t) - 0.125d0 * (I_ipq_ist - I_ipq_sit - I_ipq_tsi)
          enddo ! i
        enddo ! p
      enddo ! q
    enddo ! s
  enddo ! t
  !$OMP END DO
  !$OMP END PARALLEL

  ! ----
  ! down-down & down-down part

  sigma_p = -1.d0
  sigma_s = -1.d0
  sigma_q = -1.d0
  sigma_t = -1.d0

  !$OMP PARALLEL                                    &
  !$OMP DEFAULT (NONE)                              &
  !$OMP PRIVATE (ii, i, sigma_i, p, q, s, t,        &
  !$OMP          I_ipq_ist, I_ipq_sit, I_ipq_tsi)   &
  !$OMP SHARED (mo_num, elec_beta_num, elec_num,    &
  !$OMP         sigma_p, sigma_q, sigma_s, sigma_t, &
  !$OMP         noL_2e_naive)

  !$OMP DO COLLAPSE (4)
  do t = 1, mo_num
    do s = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num

          do ii = 1, elec_num
            if(ii .le. elec_beta_num) then
              i       = ii 
              sigma_i = -1.d0
            else
              i       = ii - elec_beta_num
              sigma_i = +1.d0
            endif

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , i, sigma_i, s, sigma_s, t, sigma_t &
                                                  , I_ipq_ist)

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , s, sigma_s, i, sigma_i, t, sigma_t &
                                                  , I_ipq_sit)

            call give_integrals_3_body_bi_ort_spin( i, sigma_i, p, sigma_p, q, sigma_q &
                                                  , t, sigma_t, s, sigma_s, i, sigma_i &
                                                  , I_ipq_tsi)

            ! x 0.25  because we consider 0.25 (up-up + up-down + down-up + down-down)
            noL_2e_naive(p,q,s,t) = noL_2e_naive(p,q,s,t) - 0.125d0 * (I_ipq_ist - I_ipq_sit - I_ipq_tsi)
          enddo ! i
        enddo ! p
      enddo ! q
    enddo ! s
  enddo ! t
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(t1)
  print*, " Wall time for noL_2e_naive (min) = ", (t1 - t0)/60.d0

END_PROVIDER

! ---


