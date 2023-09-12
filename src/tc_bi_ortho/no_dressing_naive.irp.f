
! ---

BEGIN_PROVIDER [double precision, no_0_naive]

  implicit none
  integer           :: ii, jj, kk
  integer           :: i, j, k
  double precision  :: sigma_i, sigma_j, sigma_k
  double precision  :: tmp
  double precision  :: I_ijk_ijk, I_ijk_kij, I_ijk_jki, I_ijk_jik, I_ijk_kji, I_ijk_ikj
  double precision  :: t0, t1
  logical, external :: is_same_spin

  print*, " Providing no_0_naive ..."
  call wall_time(t0)


  tmp = 0.d0
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


        tmp = tmp + I_ijk_ijk + I_ijk_kij + I_ijk_jki - I_ijk_jik - I_ijk_kji - I_ijk_ikj
        !tmp = tmp + I_ijk_ijk + 2.d0 * I_ijk_kij - 3.d0 * I_ijk_jik
      enddo
    enddo
  enddo

  no_0_naive = -1.d0 * (-tmp) / 6.d0

  call wall_time(t1)
  print*, " Wall time for no_0_naive (sec) = ", (t1 - t0)/60.d0

  print*, " no_0_naive = ", no_0_naive
  
END_PROVIDER

! ---


