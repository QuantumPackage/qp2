
! ---

BEGIN_PROVIDER [double precision, noL_0e_v0]

  implicit none
  integer                       :: i, j, k
  double precision              :: I_ijk_ijk, I_ijk_kij, I_ijk_jik, I_ijk_jki, I_ijk_ikj, I_ijk_kji
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:)

  call wall_time(t0)
  print*, " Providing noL_0e_v0 ..."

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

    noL_0e_v0 = -1.d0 * (sum(tmp)) / 6.d0

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

          tmp(i) = tmp(i) + 6.d0 * (2.d0 * I_ijk_ijk + I_ijk_jki - I_ijk_ikj - 2.d0 * I_ijk_jik)
        enddo ! k

        do k = elec_beta_num+1, elec_alpha_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, j, k, i, I_ijk_jki)
          call give_integrals_3_body_bi_ort(i, j, k, i, k, j, I_ijk_ikj)
          call give_integrals_3_body_bi_ort(i, j, k, k, j, i, I_ijk_kji)

          tmp(i) = tmp(i) + 6.d0 * (I_ijk_ijk + I_ijk_jki - I_ijk_ikj - I_ijk_kji)
        enddo ! k
      enddo ! j
    enddo ! i
    !$OMP END DO
    !$OMP END PARALLEL

    noL_0e_v0 = -1.d0 * (sum(tmp)) / 6.d0

    deallocate(tmp)

  endif

  call wall_time(t1)
  print*, " Wall time for noL_0e_v0 (min) = ", (t1 - t0)/60.d0

  print*, " noL_0e_v0 = ", noL_0e_v0

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, noL_1e_v0, (mo_num, mo_num)]

  implicit none
  integer          :: p, s, i, j
  double precision :: I_pij_sij, I_pij_isj, I_pij_ijs, I_pij_sji, I_pij_jsi, I_pij_jis
  double precision :: t0, t1

  call wall_time(t0)
  print*, " Providing noL_1e_v0 ..."

  if(elec_alpha_num .eq. elec_beta_num) then

    !$OMP PARALLEL                                 &
    !$OMP DEFAULT (NONE)                           &
    !$OMP PRIVATE (p, s, i, j,                     &
    !$OMP         I_pij_sij, I_pij_isj, I_pij_ijs, &
    !$OMP         I_pij_sji)                       &
    !$OMP SHARED (mo_num, elec_beta_num, noL_1e_v0)

    !$OMP DO COLLAPSE(2)
    do s = 1, mo_num
      do p = 1, mo_num
  
        noL_1e_v0(p,s) = 0.d0
        do i = 1, elec_beta_num
          do j = 1, elec_beta_num

            call give_integrals_3_body_bi_ort(p, i, j, s, i, j, I_pij_sij)
            call give_integrals_3_body_bi_ort(p, i, j, i, s, j, I_pij_isj)
            call give_integrals_3_body_bi_ort(p, i, j, i, j, s, I_pij_ijs)
            call give_integrals_3_body_bi_ort(p, i, j, s, j, i, I_pij_sji)
          
            noL_1e_v0(p,s) = noL_1e_v0(p,s) + (2.d0*I_pij_sij - 2.d0*I_pij_isj + I_pij_ijs - I_pij_sji)
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
    !$OMP SHARED (mo_num, elec_beta_num, elec_alpha_num, noL_1e_v0)

    !$OMP DO COLLAPSE(2)
    do s = 1, mo_num
      do p = 1, mo_num
  
        noL_1e_v0(p,s) = 0.d0
        do i = 1, elec_beta_num
          do j = 1, elec_beta_num

            call give_integrals_3_body_bi_ort(p, i, j, s, i, j, I_pij_sij)
            call give_integrals_3_body_bi_ort(p, i, j, i, s, j, I_pij_isj)
            call give_integrals_3_body_bi_ort(p, i, j, i, j, s, I_pij_ijs)
            call give_integrals_3_body_bi_ort(p, i, j, s, j, i, I_pij_sji)
          
            noL_1e_v0(p,s) = noL_1e_v0(p,s) + (2.d0*I_pij_sij - 2.d0*I_pij_isj + I_pij_ijs - I_pij_sji)
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
          
            noL_1e_v0(p,s) = noL_1e_v0(p,s) - 0.5d0 * (2.d0*I_pij_sji - I_pij_jsi + 2.d0*I_pij_jis - 4.d0*I_pij_sij + 2.d0*I_pij_isj - I_pij_ijs)
          enddo ! j

          do j = elec_beta_num+1, elec_alpha_num

            call give_integrals_3_body_bi_ort(p, i, j, s, i, j, I_pij_sij)
            call give_integrals_3_body_bi_ort(p, i, j, i, s, j, I_pij_isj)
            call give_integrals_3_body_bi_ort(p, i, j, i, j, s, I_pij_ijs)
            call give_integrals_3_body_bi_ort(p, i, j, s, j, i, I_pij_sji)
          
            noL_1e_v0(p,s) = noL_1e_v0(p,s) + 0.5d0 * (I_pij_sij - I_pij_isj + I_pij_ijs - I_pij_sji)
          enddo ! j
        enddo ! i

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

  endif

  call wall_time(t1)
  print*, " Wall time for noL_1e_v0 (min) = ", (t1 - t0)/60.d0

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, noL_2e_v0, (mo_num, mo_num, mo_num, mo_num)]

  implicit none
  integer          :: p, q, s, t, i
  double precision :: I_ipq_sit, I_ipq_tsi, I_ipq_ist
  double precision :: t0, t1

  call wall_time(t0)
  print*, " Providing noL_2e_v0 ..."

  if(elec_alpha_num .eq. elec_beta_num) then

    !$OMP PARALLEL                                  &
    !$OMP DEFAULT (NONE)                            &
    !$OMP PRIVATE (p, q, s, t, i,                   &
    !$OMP          I_ipq_sit, I_ipq_tsi, I_ipq_ist) &
    !$OMP SHARED (mo_num, elec_beta_num, noL_2e_v0)

    !$OMP DO COLLAPSE(4)
    do t = 1, mo_num
      do s = 1, mo_num
        do q = 1, mo_num
          do p = 1, mo_num
  
            noL_2e_v0(p,q,s,t) = 0.d0
            do i = 1, elec_beta_num

              call give_integrals_3_body_bi_ort(i, p, q, s, i, t, I_ipq_sit)
              call give_integrals_3_body_bi_ort(i, p, q, t, s, i, I_ipq_tsi)
              call give_integrals_3_body_bi_ort(i, p, q, i, s, t, I_ipq_ist)
          
              noL_2e_v0(p,q,s,t) = noL_2e_v0(p,q,s,t) + 0.5d0 * (I_ipq_sit + I_ipq_tsi - 2.d0*I_ipq_ist)
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
    !$OMP SHARED (mo_num, elec_beta_num, elec_alpha_num, noL_2e_v0)

    !$OMP DO COLLAPSE(4)
    do t = 1, mo_num
      do s = 1, mo_num
        do q = 1, mo_num
          do p = 1, mo_num
  
            noL_2e_v0(p,q,s,t) = 0.d0
            do i = 1, elec_beta_num

              call give_integrals_3_body_bi_ort(i, p, q, s, i, t, I_ipq_sit)
              call give_integrals_3_body_bi_ort(i, p, q, t, s, i, I_ipq_tsi)
              call give_integrals_3_body_bi_ort(i, p, q, i, s, t, I_ipq_ist)
            
              noL_2e_v0(p,q,s,t) = noL_2e_v0(p,q,s,t) + 0.5d0 * (I_ipq_sit + I_ipq_tsi - 2.d0*I_ipq_ist)
            enddo ! i

            do i = elec_beta_num+1, elec_alpha_num

              call give_integrals_3_body_bi_ort(i, p, q, s, i, t, I_ipq_sit)
              call give_integrals_3_body_bi_ort(i, p, q, t, s, i, I_ipq_tsi)
              call give_integrals_3_body_bi_ort(i, p, q, i, s, t, I_ipq_ist)
            
              noL_2e_v0(p,q,s,t) = noL_2e_v0(p,q,s,t) + 0.25d0 * (I_ipq_sit + I_ipq_tsi - 2.d0*I_ipq_ist)
            enddo ! i

          enddo ! p
        enddo ! q
      enddo ! s
    enddo ! t
    !$OMP END DO
    !$OMP END PARALLEL

  endif

  call wall_time(t1)
  print*, " Wall time for noL_2e_v0 (min) = ", (t1 - t0)/60.d0

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, noL_0e]

  implicit none
  integer                       :: i, j, k, ipoint
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:)
  double precision, allocatable :: tmp_L(:,:), tmp_R(:,:)
  double precision, allocatable :: tmp_M(:,:), tmp_S(:), tmp_O(:), tmp_J(:,:)
  double precision, allocatable :: tmp_M_priv(:,:), tmp_S_priv(:), tmp_O_priv(:), tmp_J_priv(:,:)


  call wall_time(t0)
  print*, " Providing noL_0e ..."

  if(elec_alpha_num .eq. elec_beta_num) then

    allocate(tmp(elec_beta_num))
    allocate(tmp_L(n_points_final_grid,3), tmp_R(n_points_final_grid,3))

    !$OMP PARALLEL                                                 &
    !$OMP DEFAULT(NONE)                                            &
    !$OMP PRIVATE(j, i, ipoint, tmp_L, tmp_R)                      &
    !$OMP SHARED(elec_beta_num, n_points_final_grid,               & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
    !$OMP        int2_grad1_u12_bimo_t, tmp, final_weight_at_r_vector)

    !$OMP DO
    do j = 1, elec_beta_num

      tmp_L = 0.d0
      tmp_R = 0.d0
      do i = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          tmp_L(ipoint,1) = tmp_L(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,2) = tmp_L(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,3) = tmp_L(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i)

          tmp_R(ipoint,1) = tmp_R(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,2) = tmp_R(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,3) = tmp_R(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_r_in_r_array_transp(ipoint,i)
        enddo
      enddo

      tmp(j) = 0.d0
      do ipoint = 1, n_points_final_grid
        tmp(j) = tmp(j) + final_weight_at_r_vector(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
      enddo
    enddo ! j
    !$OMP END DO
    !$OMP END PARALLEL

    noL_0e = -2.d0 * sum(tmp)

    deallocate(tmp)
    deallocate(tmp_L, tmp_R)

    ! ---

    allocate(tmp_O(n_points_final_grid), tmp_J(n_points_final_grid,3))
    tmp_O = 0.d0
    tmp_J = 0.d0

    !$OMP PARALLEL                                                  &
    !$OMP DEFAULT(NONE)                                             &
    !$OMP PRIVATE(i, ipoint, tmp_O_priv, tmp_J_priv)                &
    !$OMP SHARED(elec_beta_num, n_points_final_grid,                & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,  &
    !$OMP        int2_grad1_u12_bimo_t, tmp_O, tmp_J)

    allocate(tmp_O_priv(n_points_final_grid), tmp_J_priv(n_points_final_grid,3))
    tmp_O_priv = 0.d0
    tmp_J_priv = 0.d0
  
    !$OMP DO 
    do i = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid
        tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,i)
        tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,i)
        tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_O = tmp_O + tmp_O_priv
    tmp_J = tmp_J + tmp_J_priv
    !$OMP END CRITICAL

    deallocate(tmp_O_priv, tmp_J_priv)
    !$OMP END PARALLEL

    allocate(tmp_M(n_points_final_grid,3), tmp_S(n_points_final_grid))
    tmp_M = 0.d0
    tmp_S = 0.d0

    !$OMP PARALLEL                                                 &
    !$OMP DEFAULT(NONE)                                            &
    !$OMP PRIVATE(i, j, ipoint, tmp_M_priv, tmp_S_priv)            &
    !$OMP SHARED(elec_beta_num, n_points_final_grid,               & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
    !$OMP        int2_grad1_u12_bimo_t, tmp_M, tmp_S)

    allocate(tmp_M_priv(n_points_final_grid,3), tmp_S_priv(n_points_final_grid))
    tmp_M_priv = 0.d0
    tmp_S_priv = 0.d0
  
    !$OMP DO COLLAPSE(2)
    do i = 1, elec_beta_num
      do j = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_M = tmp_M + tmp_M_priv
    tmp_S = tmp_S + tmp_S_priv
    !$OMP END CRITICAL

    deallocate(tmp_M_priv, tmp_S_priv)
    !$OMP END PARALLEL

    allocate(tmp(n_points_final_grid))

    do ipoint = 1, n_points_final_grid

      tmp_S(ipoint) = 2.d0 * (tmp_J(ipoint,1)*tmp_J(ipoint,1) + tmp_J(ipoint,2)*tmp_J(ipoint,2) + tmp_J(ipoint,3)*tmp_J(ipoint,3)) - tmp_S(ipoint)

      tmp(ipoint) = final_weight_at_r_vector(ipoint) * ( tmp_O(ipoint) * tmp_S(ipoint)              &
                                                       - 2.d0 * ( tmp_J(ipoint,1) * tmp_M(ipoint,1) &
                                                                + tmp_J(ipoint,2) * tmp_M(ipoint,2) &
                                                                + tmp_J(ipoint,3) * tmp_M(ipoint,3)))
    enddo

    noL_0e = noL_0e -2.d0 * (sum(tmp))

    deallocate(tmp)

  else

    allocate(tmp(elec_alpha_num))
    allocate(tmp_L(n_points_final_grid,3), tmp_R(n_points_final_grid,3))

    !$OMP PARALLEL                                                   &
    !$OMP DEFAULT(NONE)                                              &
    !$OMP PRIVATE(j, i, ipoint, tmp_L, tmp_R)                        &
    !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
    !$OMP        int2_grad1_u12_bimo_t, tmp, final_weight_at_r_vector)

    !$OMP DO
    do j = 1, elec_beta_num

      tmp_L = 0.d0
      tmp_R = 0.d0
      do i = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid

          tmp_L(ipoint,1) = tmp_L(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,2) = tmp_L(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,3) = tmp_L(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i)

          tmp_R(ipoint,1) = tmp_R(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,2) = tmp_R(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,3) = tmp_R(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_r_in_r_array_transp(ipoint,i)
        enddo
      enddo

      tmp(j) = 0.d0
      do ipoint = 1, n_points_final_grid
        tmp(j) = tmp(j) + final_weight_at_r_vector(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
      enddo

      do i = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          tmp_L(ipoint,1) = tmp_L(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,2) = tmp_L(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,3) = tmp_L(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i)

          tmp_R(ipoint,1) = tmp_R(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,2) = tmp_R(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,3) = tmp_R(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_r_in_r_array_transp(ipoint,i)
        enddo
      enddo

      do ipoint = 1, n_points_final_grid
        tmp(j) = tmp(j) + final_weight_at_r_vector(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
      enddo
    enddo ! j
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    !$OMP PARALLEL                                                   &
    !$OMP DEFAULT(NONE)                                              &
    !$OMP PRIVATE(j, i, ipoint, tmp_L, tmp_R)                        &
    !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
    !$OMP        int2_grad1_u12_bimo_t, tmp, final_weight_at_r_vector)

    !$OMP DO
    do j = elec_beta_num+1, elec_alpha_num

      tmp_L = 0.d0
      tmp_R = 0.d0
      do i = 1, elec_alpha_num
        do ipoint = 1, n_points_final_grid
          tmp_L(ipoint,1) = tmp_L(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,2) = tmp_L(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,3) = tmp_L(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i)

          tmp_R(ipoint,1) = tmp_R(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,2) = tmp_R(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,3) = tmp_R(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_r_in_r_array_transp(ipoint,i)
        enddo
      enddo

      tmp(j) = 0.d0
      do ipoint = 1, n_points_final_grid
        tmp(j) = tmp(j) + 0.5d0 * final_weight_at_r_vector(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
      enddo
    enddo ! j
    !$OMP END DO
    !$OMP END PARALLEL

    noL_0e = -2.d0 * sum(tmp)

    deallocate(tmp)
    deallocate(tmp_L, tmp_R)

    ! ---

    allocate(tmp_O(n_points_final_grid), tmp_J(n_points_final_grid,3))
    tmp_O = 0.d0
    tmp_J = 0.d0

    !$OMP PARALLEL                                                   &
    !$OMP DEFAULT(NONE)                                              &
    !$OMP PRIVATE(i, ipoint, tmp_O_priv, tmp_J_priv)                 &
    !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
    !$OMP        int2_grad1_u12_bimo_t, tmp_O, tmp_J)

    allocate(tmp_O_priv(n_points_final_grid), tmp_J_priv(n_points_final_grid,3))
    tmp_O_priv = 0.d0
    tmp_J_priv = 0.d0
  
    !$OMP DO 
    do i = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid
        tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,i)
        tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,i)
        tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO 
    do i = elec_beta_num+1, elec_alpha_num
      do ipoint = 1, n_points_final_grid
        tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + 0.5d0 * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,i)
        tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,i)
        tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_O = tmp_O + tmp_O_priv
    tmp_J = tmp_J + tmp_J_priv
    !$OMP END CRITICAL

    deallocate(tmp_O_priv, tmp_J_priv)
    !$OMP END PARALLEL

    ! ---

    allocate(tmp_M(n_points_final_grid,3), tmp_S(n_points_final_grid))
    tmp_M = 0.d0
    tmp_S = 0.d0

    !$OMP PARALLEL                                                   &
    !$OMP DEFAULT(NONE)                                              &
    !$OMP PRIVATE(i, j, ipoint, tmp_M_priv, tmp_S_priv)              &
    !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
    !$OMP        int2_grad1_u12_bimo_t, tmp_M, tmp_S)

    allocate(tmp_M_priv(n_points_final_grid,3), tmp_S_priv(n_points_final_grid))
    tmp_M_priv = 0.d0
    tmp_S_priv = 0.d0
  
    !$OMP DO COLLAPSE(2)
    do i = 1, elec_beta_num
      do j = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO COLLAPSE(2)
    do i = elec_beta_num+1, elec_alpha_num
      do j = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO COLLAPSE(2)
    do i = elec_beta_num+1, elec_alpha_num
      do j = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                      + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                      + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_M = tmp_M + tmp_M_priv
    tmp_S = tmp_S + tmp_S_priv
    !$OMP END CRITICAL

    deallocate(tmp_M_priv, tmp_S_priv)
    !$OMP END PARALLEL

    allocate(tmp(n_points_final_grid))

    do ipoint = 1, n_points_final_grid

      tmp_S(ipoint) = 2.d0 * (tmp_J(ipoint,1)*tmp_J(ipoint,1) + tmp_J(ipoint,2)*tmp_J(ipoint,2) + tmp_J(ipoint,3)*tmp_J(ipoint,3)) - tmp_S(ipoint)

      tmp(ipoint) = final_weight_at_r_vector(ipoint) * ( tmp_O(ipoint) * tmp_S(ipoint)              &
                                                       - 2.d0 * ( tmp_J(ipoint,1) * tmp_M(ipoint,1) &
                                                                + tmp_J(ipoint,2) * tmp_M(ipoint,2) &
                                                                + tmp_J(ipoint,3) * tmp_M(ipoint,3)))
    enddo

    noL_0e = noL_0e -2.d0 * (sum(tmp))

    deallocate(tmp)

  endif

  call wall_time(t1)
  print*, " Wall time for noL_0e (min) = ", (t1 - t0)/60.d0

  print*, " noL_0e = ", noL_0e

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, noL_1e, (mo_num, mo_num)]

  implicit none
  integer                       :: p, s, i, j, ipoint
  double precision              :: t0, t1
  double precision, allocatable :: tmp1(:,:,:,:), tmp2(:,:), tmp3(:,:,:), tmp4(:,:,:)
  double precision, allocatable :: tmp_L(:,:,:), tmp_R(:,:,:), tmp_M(:,:), tmp_S(:), tmp_O(:), tmp_J(:,:)
  double precision, allocatable :: tmp_L0(:,:,:), tmp_R0(:,:,:)
  double precision, allocatable :: tmp_M_priv(:,:), tmp_S_priv(:), tmp_O_priv(:), tmp_J_priv(:,:)


  PROVIDE int2_grad1_u12_bimo_t
  PROVIDE mos_l_in_r_array_transp mos_r_in_r_array_transp

  call wall_time(t0)
  print*, " Providing noL_1e ..."

  if(elec_alpha_num .eq. elec_beta_num) then

    allocate(tmp_O(n_points_final_grid), tmp_J(n_points_final_grid,3))
    tmp_O = 0.d0
    tmp_J = 0.d0

    !$OMP PARALLEL                                                  &
    !$OMP DEFAULT(NONE)                                             &
    !$OMP PRIVATE(i, ipoint, tmp_O_priv, tmp_J_priv)                &
    !$OMP SHARED(elec_beta_num, n_points_final_grid,                & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,  &
    !$OMP        int2_grad1_u12_bimo_t, tmp_O, tmp_J)

    allocate(tmp_O_priv(n_points_final_grid), tmp_J_priv(n_points_final_grid,3))
    tmp_O_priv = 0.d0
    tmp_J_priv = 0.d0
  
    !$OMP DO 
    do i = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid
        tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,i)
        tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,i)
        tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_O = tmp_O + tmp_O_priv
    tmp_J = tmp_J + tmp_J_priv
    !$OMP END CRITICAL

    deallocate(tmp_O_priv, tmp_J_priv)
    !$OMP END PARALLEL

    ! ---

    allocate(tmp_M(n_points_final_grid,3), tmp_S(n_points_final_grid))
    tmp_M = 0.d0
    tmp_S = 0.d0

    !$OMP PARALLEL                                                 &
    !$OMP DEFAULT(NONE)                                            &
    !$OMP PRIVATE(i, j, ipoint, tmp_M_priv, tmp_S_priv)            &
    !$OMP SHARED(elec_beta_num, n_points_final_grid,               & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
    !$OMP        int2_grad1_u12_bimo_t, tmp_M, tmp_S)

    allocate(tmp_M_priv(n_points_final_grid,3), tmp_S_priv(n_points_final_grid))
    tmp_M_priv = 0.d0
    tmp_S_priv = 0.d0
  
    !$OMP DO COLLAPSE(2)
    do i = 1, elec_beta_num
      do j = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_M = tmp_M + tmp_M_priv
    tmp_S = tmp_S + tmp_S_priv
    !$OMP END CRITICAL

    deallocate(tmp_M_priv, tmp_S_priv)
    !$OMP END PARALLEL

    ! ---

    allocate(tmp2(n_points_final_grid,4))
    allocate(tmp1(n_points_final_grid,4,mo_num,mo_num))

    do ipoint = 1, n_points_final_grid

      tmp2(ipoint,1) = final_weight_at_r_vector(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,1) - tmp_M(ipoint,1))
      tmp2(ipoint,2) = final_weight_at_r_vector(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,2) - tmp_M(ipoint,2))
      tmp2(ipoint,3) = final_weight_at_r_vector(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,3) - tmp_M(ipoint,3))
      tmp2(ipoint,4) = -final_weight_at_r_vector(ipoint) * tmp_O(ipoint)

      tmp_S(ipoint) = 2.d0 * (tmp_J(ipoint,1) * tmp_J(ipoint,1) + tmp_J(ipoint,2) * tmp_J(ipoint,2) + tmp_J(ipoint,3) * tmp_J(ipoint,3)) - tmp_S(ipoint)
    enddo

    deallocate(tmp_O, tmp_M)

    !$OMP PARALLEL                                           &
    !$OMP DEFAULT(NONE)                                      &
    !$OMP PRIVATE(p, s, i, ipoint)                           &
    !$OMP SHARED(mo_num, elec_beta_num, n_points_final_grid, & 
    !$OMP        int2_grad1_u12_bimo_t, tmp1)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, mo_num
      do p = 1, mo_num

        do ipoint = 1, n_points_final_grid
          tmp1(ipoint,1,p,s) = int2_grad1_u12_bimo_t(ipoint,1,p,s)
          tmp1(ipoint,2,p,s) = int2_grad1_u12_bimo_t(ipoint,2,p,s)
          tmp1(ipoint,3,p,s) = int2_grad1_u12_bimo_t(ipoint,3,p,s)
        enddo

        tmp1(:,4,p,s) = 0.d0
        do i = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,4,p,s) = tmp1(ipoint,4,p,s) + int2_grad1_u12_bimo_t(ipoint,1,p,i) * int2_grad1_u12_bimo_t(ipoint,1,i,s) &
                                                    + int2_grad1_u12_bimo_t(ipoint,2,p,i) * int2_grad1_u12_bimo_t(ipoint,2,i,s) &
                                                    + int2_grad1_u12_bimo_t(ipoint,3,p,i) * int2_grad1_u12_bimo_t(ipoint,3,i,s)
          enddo
        enddo

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemv( 'T', 4*n_points_final_grid, mo_num*mo_num, 2.d0 &
              , tmp1(1,1,1,1), size(tmp1, 1) * size(tmp1, 2)    &
              , tmp2(1,1), 1                                    &
              , 0.d0, noL_1e(1,1), 1)

    deallocate(tmp1, tmp2)

    ! ---

    allocate(tmp_L(n_points_final_grid,3,mo_num))
    allocate(tmp_R(n_points_final_grid,3,mo_num))

    !$OMP PARALLEL                                                 &
    !$OMP DEFAULT(NONE)                                            &
    !$OMP PRIVATE(p, i, ipoint)                                    &
    !$OMP SHARED(elec_beta_num, n_points_final_grid, mo_num,       & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
    !$OMP        int2_grad1_u12_bimo_t, tmp_L, tmp_R)

    !$OMP DO
    do p = 1, mo_num

      tmp_L(:,1:3,p) = 0.d0
      tmp_R(:,1:3,p) = 0.d0

      do i = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          tmp_L(ipoint,1,p) = tmp_L(ipoint,1,p) + int2_grad1_u12_bimo_t(ipoint,1,p,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,2,p) = tmp_L(ipoint,2,p) + int2_grad1_u12_bimo_t(ipoint,2,p,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,3,p) = tmp_L(ipoint,3,p) + int2_grad1_u12_bimo_t(ipoint,3,p,i) * mos_l_in_r_array_transp(ipoint,i)

          tmp_R(ipoint,1,p) = tmp_R(ipoint,1,p) + int2_grad1_u12_bimo_t(ipoint,1,i,p) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,2,p) = tmp_R(ipoint,2,p) + int2_grad1_u12_bimo_t(ipoint,2,i,p) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,3,p) = tmp_R(ipoint,3,p) + int2_grad1_u12_bimo_t(ipoint,3,i,p) * mos_r_in_r_array_transp(ipoint,i)
        enddo
      enddo
    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    allocate(tmp3(n_points_final_grid,5,mo_num))
    allocate(tmp4(n_points_final_grid,5,mo_num))

    !$OMP PARALLEL                                                 &
    !$OMP DEFAULT(NONE)                                            &
    !$OMP PRIVATE(p, i, j, ipoint)                                 &
    !$OMP SHARED(elec_beta_num, n_points_final_grid, mo_num,       & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
    !$OMP        int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
    !$OMP        tmp_L, tmp_R, tmp_J, tmp_S, tmp3, tmp4)

    !$OMP DO
    do p = 1, mo_num

      do ipoint = 1, n_points_final_grid

        tmp3(ipoint,1,p) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,p)
        tmp3(ipoint,2,p) = -2.d0 * (tmp_L(ipoint,1,p) * tmp_J(ipoint,1) + tmp_L(ipoint,2,p) * tmp_J(ipoint,2) + tmp_L(ipoint,3,p) * tmp_J(ipoint,3))
        tmp3(ipoint,3,p) = final_weight_at_r_vector(ipoint) * tmp_L(ipoint,1,p)
        tmp3(ipoint,4,p) = final_weight_at_r_vector(ipoint) * tmp_L(ipoint,2,p)
        tmp3(ipoint,5,p) = final_weight_at_r_vector(ipoint) * tmp_L(ipoint,3,p)

        tmp4(ipoint,1,p) = -2.d0 * (tmp_R(ipoint,1,p) * tmp_J(ipoint,1) + tmp_R(ipoint,2,p) * tmp_J(ipoint,2) + tmp_R(ipoint,3,p) * tmp_J(ipoint,3)) &
                         + mos_r_in_r_array_transp(ipoint,p) * tmp_S(ipoint)
        tmp4(ipoint,2,p) = final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,p)
        tmp4(ipoint,3,p) = tmp_R(ipoint,1,p)
        tmp4(ipoint,4,p) = tmp_R(ipoint,2,p)
        tmp4(ipoint,5,p) = tmp_R(ipoint,3,p)
      enddo

      do i = 1, elec_beta_num
        do j = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid

            tmp3(ipoint,2,p) = tmp3(ipoint,2,p) + mos_l_in_r_array_transp(ipoint,j) * ( int2_grad1_u12_bimo_t(ipoint,1,p,i) * int2_grad1_u12_bimo_t(ipoint,1,i,j) & 
                                                                                      + int2_grad1_u12_bimo_t(ipoint,2,p,i) * int2_grad1_u12_bimo_t(ipoint,2,i,j) &                                                                                     
                                                                                      + int2_grad1_u12_bimo_t(ipoint,3,p,i) * int2_grad1_u12_bimo_t(ipoint,3,i,j) )

            tmp4(ipoint,1,p) = tmp4(ipoint,1,p) + mos_r_in_r_array_transp(ipoint,i) * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,p) & 
                                                                                      + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,p) &                                                                                     
                                                                                      + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmp_L, tmp_R, tmp_J, tmp_S)

    call dgemm( 'T', 'N', mo_num, mo_num, 5*n_points_final_grid, 1.d0                  &
              , tmp3(1,1,1), 5*n_points_final_grid, tmp4(1,1,1), 5*n_points_final_grid &
              , 1.d0, noL_1e(1,1), mo_num)
   
    deallocate(tmp3, tmp4)

    ! ---

  else

    allocate(tmp_O(n_points_final_grid), tmp_J(n_points_final_grid,3))
    tmp_O = 0.d0
    tmp_J = 0.d0

    !$OMP PARALLEL                                                   &
    !$OMP DEFAULT(NONE)                                              &
    !$OMP PRIVATE(i, ipoint, tmp_O_priv, tmp_J_priv)                 &
    !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
    !$OMP        int2_grad1_u12_bimo_t, tmp_O, tmp_J)

    allocate(tmp_O_priv(n_points_final_grid), tmp_J_priv(n_points_final_grid,3))
    tmp_O_priv = 0.d0
    tmp_J_priv = 0.d0
  
    !$OMP DO 
    do i = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid
        tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,i)
        tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,i)
        tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO 
    do i = elec_beta_num+1, elec_alpha_num
      do ipoint = 1, n_points_final_grid
        tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + 0.5d0 * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,i)
        tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,i)
        tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_O = tmp_O + tmp_O_priv
    tmp_J = tmp_J + tmp_J_priv
    !$OMP END CRITICAL

    deallocate(tmp_O_priv, tmp_J_priv)
    !$OMP END PARALLEL

    ! ---

    allocate(tmp_M(n_points_final_grid,3), tmp_S(n_points_final_grid))
    tmp_M = 0.d0
    tmp_S = 0.d0

    !$OMP PARALLEL                                                   &
    !$OMP DEFAULT(NONE)                                              &
    !$OMP PRIVATE(i, j, ipoint, tmp_M_priv, tmp_S_priv)              &
    !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
    !$OMP        int2_grad1_u12_bimo_t, tmp_M, tmp_S)

    allocate(tmp_M_priv(n_points_final_grid,3), tmp_S_priv(n_points_final_grid))
    tmp_M_priv = 0.d0
    tmp_S_priv = 0.d0
  
    !$OMP DO COLLAPSE(2)
    do i = 1, elec_beta_num
      do j = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO COLLAPSE(2)
    do i = elec_beta_num+1, elec_alpha_num
      do j = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                      + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO COLLAPSE(2)
    do i = elec_beta_num+1, elec_alpha_num
      do j = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                      + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                      + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_M = tmp_M + tmp_M_priv
    tmp_S = tmp_S + tmp_S_priv
    !$OMP END CRITICAL

    deallocate(tmp_M_priv, tmp_S_priv)
    !$OMP END PARALLEL

    ! ---

    allocate(tmp2(n_points_final_grid,4))
    allocate(tmp1(n_points_final_grid,4,mo_num,mo_num))

    do ipoint = 1, n_points_final_grid

      tmp2(ipoint,1) = final_weight_at_r_vector(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,1) - tmp_M(ipoint,1))
      tmp2(ipoint,2) = final_weight_at_r_vector(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,2) - tmp_M(ipoint,2))
      tmp2(ipoint,3) = final_weight_at_r_vector(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,3) - tmp_M(ipoint,3))
      tmp2(ipoint,4) = -final_weight_at_r_vector(ipoint) * tmp_O(ipoint)

      tmp_S(ipoint) = 2.d0 * (tmp_J(ipoint,1) * tmp_J(ipoint,1) + tmp_J(ipoint,2) * tmp_J(ipoint,2) + tmp_J(ipoint,3) * tmp_J(ipoint,3)) - tmp_S(ipoint)
    enddo

    deallocate(tmp_O, tmp_M)

    !$OMP PARALLEL                                           &
    !$OMP DEFAULT(NONE)                                      &
    !$OMP PRIVATE(p, s, i, ipoint)                           &
    !$OMP SHARED(mo_num, elec_beta_num, n_points_final_grid, & 
    !$OMP        elec_alpha_num, int2_grad1_u12_bimo_t, tmp1)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, mo_num
      do p = 1, mo_num

        do ipoint = 1, n_points_final_grid
          tmp1(ipoint,1,p,s) = int2_grad1_u12_bimo_t(ipoint,1,p,s)
          tmp1(ipoint,2,p,s) = int2_grad1_u12_bimo_t(ipoint,2,p,s)
          tmp1(ipoint,3,p,s) = int2_grad1_u12_bimo_t(ipoint,3,p,s)
        enddo

        tmp1(:,4,p,s) = 0.d0
        do i = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,4,p,s) = tmp1(ipoint,4,p,s) + int2_grad1_u12_bimo_t(ipoint,1,p,i) * int2_grad1_u12_bimo_t(ipoint,1,i,s) &
                                                    + int2_grad1_u12_bimo_t(ipoint,2,p,i) * int2_grad1_u12_bimo_t(ipoint,2,i,s) &
                                                    + int2_grad1_u12_bimo_t(ipoint,3,p,i) * int2_grad1_u12_bimo_t(ipoint,3,i,s)
          enddo
        enddo
        do i = elec_beta_num+1, elec_alpha_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,4,p,s) = tmp1(ipoint,4,p,s) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,p,i) * int2_grad1_u12_bimo_t(ipoint,1,i,s) &
                                                    + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,p,i) * int2_grad1_u12_bimo_t(ipoint,2,i,s) &
                                                    + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,p,i) * int2_grad1_u12_bimo_t(ipoint,3,i,s)
          enddo
        enddo

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemv( 'T', 4*n_points_final_grid, mo_num*mo_num, 2.d0 &
              , tmp1(1,1,1,1), size(tmp1, 1) * size(tmp1, 2)    &
              , tmp2(1,1), 1                                    &
              , 0.d0, noL_1e(1,1), 1)

    deallocate(tmp1, tmp2)

    ! ---

    allocate(tmp_L(n_points_final_grid,3,mo_num), tmp_L0(n_points_final_grid,3,mo_num))
    allocate(tmp_R(n_points_final_grid,3,mo_num), tmp_R0(n_points_final_grid,3,mo_num))

    !$OMP PARALLEL                                                           &
    !$OMP DEFAULT(NONE)                                                      &
    !$OMP PRIVATE(p, i, ipoint)                                              &
    !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, mo_num, & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,           &
    !$OMP        int2_grad1_u12_bimo_t, tmp_L0, tmp_R0, tmp_L, tmp_R)

    !$OMP DO
    do p = 1, mo_num

      tmp_L0(:,1:3,p) = 0.d0
      tmp_R0(:,1:3,p) = 0.d0
      do i = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid

          tmp_L0(ipoint,1,p) = tmp_L0(ipoint,1,p) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,p,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L0(ipoint,2,p) = tmp_L0(ipoint,2,p) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,p,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L0(ipoint,3,p) = tmp_L0(ipoint,3,p) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,p,i) * mos_l_in_r_array_transp(ipoint,i)
                                     
          tmp_R0(ipoint,1,p) = tmp_R0(ipoint,1,p) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,p) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R0(ipoint,2,p) = tmp_R0(ipoint,2,p) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,p) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R0(ipoint,3,p) = tmp_R0(ipoint,3,p) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,p) * mos_r_in_r_array_transp(ipoint,i)
        enddo
      enddo

      tmp_L(:,1:3,p) = tmp_L0(:,1:3,p)
      tmp_R(:,1:3,p) = tmp_R0(:,1:3,p)
      do i = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          tmp_L(ipoint,1,p) = tmp_L(ipoint,1,p) + int2_grad1_u12_bimo_t(ipoint,1,p,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,2,p) = tmp_L(ipoint,2,p) + int2_grad1_u12_bimo_t(ipoint,2,p,i) * mos_l_in_r_array_transp(ipoint,i)
          tmp_L(ipoint,3,p) = tmp_L(ipoint,3,p) + int2_grad1_u12_bimo_t(ipoint,3,p,i) * mos_l_in_r_array_transp(ipoint,i)
                                   
          tmp_R(ipoint,1,p) = tmp_R(ipoint,1,p) + int2_grad1_u12_bimo_t(ipoint,1,i,p) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,2,p) = tmp_R(ipoint,2,p) + int2_grad1_u12_bimo_t(ipoint,2,i,p) * mos_r_in_r_array_transp(ipoint,i)
          tmp_R(ipoint,3,p) = tmp_R(ipoint,3,p) + int2_grad1_u12_bimo_t(ipoint,3,i,p) * mos_r_in_r_array_transp(ipoint,i)
        enddo
      enddo

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    allocate(tmp3(n_points_final_grid,8,mo_num))
    allocate(tmp4(n_points_final_grid,8,mo_num))

    !$OMP PARALLEL                                                           &
    !$OMP DEFAULT(NONE)                                                      &
    !$OMP PRIVATE(p, i, j, ipoint)                                           &
    !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, mo_num, & 
    !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,           &
    !$OMP        int2_grad1_u12_bimo_t, final_weight_at_r_vector,            &
    !$OMP        tmp_L, tmp_L0, tmp_R, tmp_R0, tmp_J, tmp_S, tmp3, tmp4)

    !$OMP DO
    do p = 1, mo_num

      do ipoint = 1, n_points_final_grid

        tmp3(ipoint,1,p) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,p)
        tmp3(ipoint,2,p) = -2.d0 * (tmp_L(ipoint,1,p) * tmp_J(ipoint,1) + tmp_L(ipoint,2,p) * tmp_J(ipoint,2) + tmp_L(ipoint,3,p) * tmp_J(ipoint,3))
        tmp3(ipoint,3,p) = final_weight_at_r_vector(ipoint) * tmp_L(ipoint,1,p)
        tmp3(ipoint,4,p) = final_weight_at_r_vector(ipoint) * tmp_L(ipoint,2,p)
        tmp3(ipoint,5,p) = final_weight_at_r_vector(ipoint) * tmp_L(ipoint,3,p)
        tmp3(ipoint,6,p) = final_weight_at_r_vector(ipoint) * tmp_L0(ipoint,1,p)
        tmp3(ipoint,7,p) = final_weight_at_r_vector(ipoint) * tmp_L0(ipoint,2,p)
        tmp3(ipoint,8,p) = final_weight_at_r_vector(ipoint) * tmp_L0(ipoint,3,p)

        tmp4(ipoint,1,p) = -2.d0 * (tmp_R(ipoint,1,p) * tmp_J(ipoint,1) + tmp_R(ipoint,2,p) * tmp_J(ipoint,2) + tmp_R(ipoint,3,p) * tmp_J(ipoint,3)) &
                         + mos_r_in_r_array_transp(ipoint,p) * tmp_S(ipoint)
        tmp4(ipoint,2,p) = final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,p)
        tmp4(ipoint,3,p) = tmp_R(ipoint,1,p)
        tmp4(ipoint,4,p) = tmp_R(ipoint,2,p)
        tmp4(ipoint,5,p) = tmp_R(ipoint,3,p)
        tmp4(ipoint,6,p) = tmp_R0(ipoint,1,p)
        tmp4(ipoint,7,p) = tmp_R0(ipoint,2,p)
        tmp4(ipoint,8,p) = tmp_R0(ipoint,3,p)
      enddo

      do i = 1, elec_beta_num
        do j = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid

            tmp3(ipoint,2,p) = tmp3(ipoint,2,p) + mos_l_in_r_array_transp(ipoint,j) * ( int2_grad1_u12_bimo_t(ipoint,1,p,i) * int2_grad1_u12_bimo_t(ipoint,1,i,j) & 
                                                                                      + int2_grad1_u12_bimo_t(ipoint,2,p,i) * int2_grad1_u12_bimo_t(ipoint,2,i,j) &                                                                                     
                                                                                      + int2_grad1_u12_bimo_t(ipoint,3,p,i) * int2_grad1_u12_bimo_t(ipoint,3,i,j) )

            tmp4(ipoint,1,p) = tmp4(ipoint,1,p) + mos_r_in_r_array_transp(ipoint,i) * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,p) & 
                                                                                      + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,p) &                                                                                     
                                                                                      + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

      do i = elec_beta_num+1, elec_alpha_num
        do j = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid

            tmp3(ipoint,2,p) = tmp3(ipoint,2,p) + 0.5d0 * mos_l_in_r_array_transp(ipoint,j) * ( int2_grad1_u12_bimo_t(ipoint,1,p,i) * int2_grad1_u12_bimo_t(ipoint,1,i,j) & 
                                                                                              + int2_grad1_u12_bimo_t(ipoint,2,p,i) * int2_grad1_u12_bimo_t(ipoint,2,i,j) &                                                                                     
                                                                                              + int2_grad1_u12_bimo_t(ipoint,3,p,i) * int2_grad1_u12_bimo_t(ipoint,3,i,j) )
            tmp3(ipoint,2,p) = tmp3(ipoint,2,p) + 0.5d0 * mos_l_in_r_array_transp(ipoint,i) * ( int2_grad1_u12_bimo_t(ipoint,1,p,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) & 
                                                                                              + int2_grad1_u12_bimo_t(ipoint,2,p,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &                                                                                     
                                                                                              + int2_grad1_u12_bimo_t(ipoint,3,p,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i) )

            tmp4(ipoint,1,p) = tmp4(ipoint,1,p) + 0.5d0 * mos_r_in_r_array_transp(ipoint,i) * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,p) & 
                                                                                              + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,p) &                                                                                     
                                                                                              + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,p) )
            tmp4(ipoint,1,p) = tmp4(ipoint,1,p) + 0.5d0 * mos_r_in_r_array_transp(ipoint,j) * ( int2_grad1_u12_bimo_t(ipoint,1,j,i) * int2_grad1_u12_bimo_t(ipoint,1,i,p) & 
                                                                                              + int2_grad1_u12_bimo_t(ipoint,2,j,i) * int2_grad1_u12_bimo_t(ipoint,2,i,p) &                                                                                     
                                                                                              + int2_grad1_u12_bimo_t(ipoint,3,j,i) * int2_grad1_u12_bimo_t(ipoint,3,i,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

      do i = elec_beta_num+1, elec_alpha_num
        do j = elec_beta_num+1, elec_alpha_num
          do ipoint = 1, n_points_final_grid

            tmp3(ipoint,2,p) = tmp3(ipoint,2,p) + 0.5d0 * mos_l_in_r_array_transp(ipoint,j) * ( int2_grad1_u12_bimo_t(ipoint,1,p,i) * int2_grad1_u12_bimo_t(ipoint,1,i,j) & 
                                                                                              + int2_grad1_u12_bimo_t(ipoint,2,p,i) * int2_grad1_u12_bimo_t(ipoint,2,i,j) &                                                                                     
                                                                                              + int2_grad1_u12_bimo_t(ipoint,3,p,i) * int2_grad1_u12_bimo_t(ipoint,3,i,j) )

            tmp4(ipoint,1,p) = tmp4(ipoint,1,p) + 0.5d0 * mos_r_in_r_array_transp(ipoint,i) * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,p) & 
                                                                                              + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,p) &                                                                                     
                                                                                              + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmp_L0, tmp_L, tmp_R0, tmp_R, tmp_J, tmp_S)

    call dgemm( 'T', 'N', mo_num, mo_num, 8*n_points_final_grid, 1.d0                  &
              , tmp3(1,1,1), 8*n_points_final_grid, tmp4(1,1,1), 8*n_points_final_grid &
              , 1.d0, noL_1e(1,1), mo_num)
   
    deallocate(tmp3, tmp4)

  endif

  call wall_time(t1)
  print*, " Wall time for noL_1e (min) = ", (t1 - t0)/60.d0

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, noL_2e, (mo_num, mo_num, mo_num, mo_num)]

  implicit none
  integer                       :: p, q, s, t, i, ipoint
  double precision              :: t0, t1
  double precision, allocatable :: tmp_O(:), tmp_J(:,:)
  double precision, allocatable :: tmp_A(:,:,:), tmp_B(:,:,:)
  double precision, allocatable :: tmp1(:,:,:,:), tmp2(:,:,:,:)
  double precision, allocatable :: tmp(:,:,:,:)

  PROVIDE int2_grad1_u12_bimo_t
  PROVIDE mos_l_in_r_array_transp mos_r_in_r_array_transp

  call wall_time(t0)
  print*, " Providing noL_2e ..."

  if(elec_alpha_num .eq. elec_beta_num) then

    allocate(tmp_O(n_points_final_grid), tmp_J(n_points_final_grid,3))
    allocate(tmp_A(n_points_final_grid,3,mo_num), tmp_B(n_points_final_grid,3,mo_num))
    allocate(tmp1(n_points_final_grid,4,mo_num,mo_num), tmp2(n_points_final_grid,4,mo_num,mo_num))
    allocate(tmp(mo_num,mo_num,mo_num,mo_num))

    tmp_O = 0.d0
    tmp_J = 0.d0
    do i = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid
        tmp_O(ipoint)   = tmp_O(ipoint)   + final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp_J(ipoint,1) = tmp_J(ipoint,1) + final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i)
        tmp_J(ipoint,2) = tmp_J(ipoint,2) + final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i)
        tmp_J(ipoint,3) = tmp_J(ipoint,3) + final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i)
      enddo
    enddo

    !$OMP PARALLEL                                                  &
    !$OMP DEFAULT(NONE)                                             &
    !$OMP PRIVATE(p, i, ipoint)                                     &
    !$OMP SHARED(mo_num, elec_beta_num, n_points_final_grid,        & 
    !$OMP        final_weight_at_r_vector, mos_l_in_r_array_transp, &
    !$OMP        mos_r_in_r_array_transp, int2_grad1_u12_bimo_t,    &
    !$OMP        tmp_A, tmp_B)
  
    !$OMP DO
    do p = 1, mo_num

      tmp_A(:,:,p) = 0.d0
      tmp_B(:,:,p) = 0.d0
      do i = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid
          tmp_A(ipoint,1,p) = tmp_A(ipoint,1,p) + final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,1,p,i)
          tmp_A(ipoint,2,p) = tmp_A(ipoint,2,p) + final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,2,p,i)
          tmp_A(ipoint,3,p) = tmp_A(ipoint,3,p) + final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,3,p,i)
          tmp_B(ipoint,1,p) = tmp_B(ipoint,1,p) + final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,1,i,p)
          tmp_B(ipoint,2,p) = tmp_B(ipoint,2,p) + final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,2,i,p)
          tmp_B(ipoint,3,p) = tmp_B(ipoint,3,p) + final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,3,i,p)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


    !$OMP PARALLEL                                                  &
    !$OMP DEFAULT(NONE)                                             &
    !$OMP PRIVATE(p, s, i, ipoint)                                  &
    !$OMP SHARED(mo_num, elec_beta_num, n_points_final_grid,        & 
    !$OMP        final_weight_at_r_vector, mos_l_in_r_array_transp, &
    !$OMP        mos_r_in_r_array_transp, int2_grad1_u12_bimo_t,    &
    !$OMP        tmp_A, tmp_B, tmp_O, tmp_J, tmp1, tmp2)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, mo_num
      do p = 1, mo_num

        do ipoint = 1, n_points_final_grid

          tmp1(ipoint,1,p,s) = mos_r_in_r_array_transp(ipoint,s) * tmp_A(ipoint,1,p) &
                             + mos_l_in_r_array_transp(ipoint,p) * tmp_B(ipoint,1,s) &
                             - tmp_O(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p,s)   &
                             - 2.d0 * mos_l_in_r_array_transp(ipoint,p) * mos_r_in_r_array_transp(ipoint,s) * tmp_J(ipoint,1)
          tmp1(ipoint,2,p,s) = mos_r_in_r_array_transp(ipoint,s) * tmp_A(ipoint,2,p) &
                             + mos_l_in_r_array_transp(ipoint,p) * tmp_B(ipoint,2,s) &
                             - tmp_O(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p,s)   &
                             - 2.d0 * mos_l_in_r_array_transp(ipoint,p) * mos_r_in_r_array_transp(ipoint,s) * tmp_J(ipoint,2)
          tmp1(ipoint,3,p,s) = mos_r_in_r_array_transp(ipoint,s) * tmp_A(ipoint,3,p) &
                             + mos_l_in_r_array_transp(ipoint,p) * tmp_B(ipoint,3,s) &
                             - tmp_O(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p,s)   &
                             - 2.d0 * mos_l_in_r_array_transp(ipoint,p) * mos_r_in_r_array_transp(ipoint,s) * tmp_J(ipoint,3)

          tmp2(ipoint,1,p,s) = int2_grad1_u12_bimo_t(ipoint,1,p,s)
          tmp2(ipoint,2,p,s) = int2_grad1_u12_bimo_t(ipoint,2,p,s)
          tmp2(ipoint,3,p,s) = int2_grad1_u12_bimo_t(ipoint,3,p,s)
          tmp2(ipoint,4,p,s) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,p) * mos_r_in_r_array_transp(ipoint,s)

        enddo ! ipoint

        tmp1(:,4,p,s) = 0.d0
        do i = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,4,p,s) = tmp1(ipoint,4,p,s) + int2_grad1_u12_bimo_t(ipoint,1,p,i) * int2_grad1_u12_bimo_t(ipoint,1,i,s) &
                                                    + int2_grad1_u12_bimo_t(ipoint,2,p,i) * int2_grad1_u12_bimo_t(ipoint,2,i,s) &
                                                    + int2_grad1_u12_bimo_t(ipoint,3,p,i) * int2_grad1_u12_bimo_t(ipoint,3,i,s)
          enddo ! ipoint
        enddo ! i

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmp_O, tmp_J, tmp_A, tmp_B)


    call dgemm( 'T', 'N', mo_num*mo_num, mo_num*mo_num, 4*n_points_final_grid, 0.5d0       &
              , tmp1(1,1,1,1), 4*n_points_final_grid, tmp2(1,1,1,1), 4*n_points_final_grid &
              , 0.d0, tmp(1,1,1,1), mo_num*mo_num)

    deallocate(tmp1, tmp2)

    call sum_a_at(tmp, mo_num*mo_num)

    !$OMP PARALLEL            &
    !$OMP DEFAULT(NONE)       &
    !$OMP PRIVATE(t, s, q, p) &
    !$OMP SHARED(mo_num, tmp, noL_2e)
  
    !$OMP DO COLLAPSE(3)
    do t = 1, mo_num
      do s = 1, mo_num
        do q = 1, mo_num
          do p = 1, mo_num
            noL_2e(p,q,s,t) = tmp(p,s,q,t)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    deallocate(tmp)

  else

    allocate(tmp_O(n_points_final_grid), tmp_J(n_points_final_grid,3))
    allocate(tmp_A(n_points_final_grid,3,mo_num), tmp_B(n_points_final_grid,3,mo_num))
    allocate(tmp1(n_points_final_grid,4,mo_num,mo_num), tmp2(n_points_final_grid,4,mo_num,mo_num))
    allocate(tmp(mo_num,mo_num,mo_num,mo_num))

    tmp_O = 0.d0
    tmp_J = 0.d0
    do i = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid
        tmp_O(ipoint)   = tmp_O(ipoint)   + final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp_J(ipoint,1) = tmp_J(ipoint,1) + final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i)
        tmp_J(ipoint,2) = tmp_J(ipoint,2) + final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i)
        tmp_J(ipoint,3) = tmp_J(ipoint,3) + final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i)
      enddo
    enddo
    do i = elec_beta_num+1, elec_alpha_num
      do ipoint = 1, n_points_final_grid
        tmp_O(ipoint)   = tmp_O(ipoint)   + 0.5d0 * final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp_J(ipoint,1) = tmp_J(ipoint,1) + 0.5d0 * final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i)
        tmp_J(ipoint,2) = tmp_J(ipoint,2) + 0.5d0 * final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i)
        tmp_J(ipoint,3) = tmp_J(ipoint,3) + 0.5d0 * final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i)
      enddo
    enddo

    !$OMP PARALLEL                                                           &
    !$OMP DEFAULT(NONE)                                                      &
    !$OMP PRIVATE(p, i, ipoint)                                              &
    !$OMP SHARED(mo_num, elec_alpha_num, elec_beta_num, n_points_final_grid, & 
    !$OMP        final_weight_at_r_vector, mos_l_in_r_array_transp,          &
    !$OMP        mos_r_in_r_array_transp, int2_grad1_u12_bimo_t,             &
    !$OMP        tmp_A, tmp_B)
  
    !$OMP DO
    do p = 1, mo_num

      tmp_A(:,:,p) = 0.d0
      tmp_B(:,:,p) = 0.d0
      do i = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid
          tmp_A(ipoint,1,p) = tmp_A(ipoint,1,p) + final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,1,p,i)
          tmp_A(ipoint,2,p) = tmp_A(ipoint,2,p) + final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,2,p,i)
          tmp_A(ipoint,3,p) = tmp_A(ipoint,3,p) + final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,3,p,i)
          tmp_B(ipoint,1,p) = tmp_B(ipoint,1,p) + final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,1,i,p)
          tmp_B(ipoint,2,p) = tmp_B(ipoint,2,p) + final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,2,i,p)
          tmp_B(ipoint,3,p) = tmp_B(ipoint,3,p) + final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,3,i,p)
        enddo
      enddo
      do i = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid
          tmp_A(ipoint,1,p) = tmp_A(ipoint,1,p) + 0.5d0 * final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,1,p,i)
          tmp_A(ipoint,2,p) = tmp_A(ipoint,2,p) + 0.5d0 * final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,2,p,i)
          tmp_A(ipoint,3,p) = tmp_A(ipoint,3,p) + 0.5d0 * final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,3,p,i)
          tmp_B(ipoint,1,p) = tmp_B(ipoint,1,p) + 0.5d0 * final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,1,i,p)
          tmp_B(ipoint,2,p) = tmp_B(ipoint,2,p) + 0.5d0 * final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,2,i,p)
          tmp_B(ipoint,3,p) = tmp_B(ipoint,3,p) + 0.5d0 * final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,i) * int2_grad1_u12_bimo_t(ipoint,3,i,p)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


    !$OMP PARALLEL                                                           &
    !$OMP DEFAULT(NONE)                                                      &
    !$OMP PRIVATE(p, s, i, ipoint)                                           &
    !$OMP SHARED(mo_num, elec_alpha_num, elec_beta_num, n_points_final_grid, & 
    !$OMP        final_weight_at_r_vector, mos_l_in_r_array_transp,          &
    !$OMP        mos_r_in_r_array_transp, int2_grad1_u12_bimo_t,             &
    !$OMP        tmp_A, tmp_B, tmp_O, tmp_J, tmp1, tmp2)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, mo_num
      do p = 1, mo_num

        do ipoint = 1, n_points_final_grid

          tmp1(ipoint,1,p,s) = mos_r_in_r_array_transp(ipoint,s) * tmp_A(ipoint,1,p) &
                             + mos_l_in_r_array_transp(ipoint,p) * tmp_B(ipoint,1,s) &
                             - tmp_O(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p,s)   &
                             - 2.d0 * mos_l_in_r_array_transp(ipoint,p) * mos_r_in_r_array_transp(ipoint,s) * tmp_J(ipoint,1)
          tmp1(ipoint,2,p,s) = mos_r_in_r_array_transp(ipoint,s) * tmp_A(ipoint,2,p) &
                             + mos_l_in_r_array_transp(ipoint,p) * tmp_B(ipoint,2,s) &
                             - tmp_O(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p,s)   &
                             - 2.d0 * mos_l_in_r_array_transp(ipoint,p) * mos_r_in_r_array_transp(ipoint,s) * tmp_J(ipoint,2)
          tmp1(ipoint,3,p,s) = mos_r_in_r_array_transp(ipoint,s) * tmp_A(ipoint,3,p) &
                             + mos_l_in_r_array_transp(ipoint,p) * tmp_B(ipoint,3,s) &
                             - tmp_O(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p,s)   &
                             - 2.d0 * mos_l_in_r_array_transp(ipoint,p) * mos_r_in_r_array_transp(ipoint,s) * tmp_J(ipoint,3)

          tmp2(ipoint,1,p,s) = int2_grad1_u12_bimo_t(ipoint,1,p,s)
          tmp2(ipoint,2,p,s) = int2_grad1_u12_bimo_t(ipoint,2,p,s)
          tmp2(ipoint,3,p,s) = int2_grad1_u12_bimo_t(ipoint,3,p,s)
          tmp2(ipoint,4,p,s) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,p) * mos_r_in_r_array_transp(ipoint,s)

        enddo ! ipoint

        tmp1(:,4,p,s) = 0.d0
        do i = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,4,p,s) = tmp1(ipoint,4,p,s) + int2_grad1_u12_bimo_t(ipoint,1,p,i) * int2_grad1_u12_bimo_t(ipoint,1,i,s) &
                                                    + int2_grad1_u12_bimo_t(ipoint,2,p,i) * int2_grad1_u12_bimo_t(ipoint,2,i,s) &
                                                    + int2_grad1_u12_bimo_t(ipoint,3,p,i) * int2_grad1_u12_bimo_t(ipoint,3,i,s)
          enddo ! ipoint
        enddo ! i
        do i = elec_beta_num+1, elec_alpha_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,4,p,s) = tmp1(ipoint,4,p,s) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,p,i) * int2_grad1_u12_bimo_t(ipoint,1,i,s) &
                                                    + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,p,i) * int2_grad1_u12_bimo_t(ipoint,2,i,s) &
                                                    + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,p,i) * int2_grad1_u12_bimo_t(ipoint,3,i,s)
          enddo ! ipoint
        enddo ! i

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmp_O, tmp_J, tmp_A, tmp_B)


    call dgemm( 'T', 'N', mo_num*mo_num, mo_num*mo_num, 4*n_points_final_grid, 0.5d0       &
              , tmp1(1,1,1,1), 4*n_points_final_grid, tmp2(1,1,1,1), 4*n_points_final_grid &
              , 0.d0, tmp(1,1,1,1), mo_num*mo_num)

    deallocate(tmp1, tmp2)

    call sum_a_at(tmp, mo_num*mo_num)

    !$OMP PARALLEL            &
    !$OMP DEFAULT(NONE)       &
    !$OMP PRIVATE(t, s, q, p) &
    !$OMP SHARED(mo_num, tmp, noL_2e)
  
    !$OMP DO COLLAPSE(3)
    do t = 1, mo_num
      do s = 1, mo_num
        do q = 1, mo_num
          do p = 1, mo_num
            noL_2e(p,q,s,t) = tmp(p,s,q,t)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    deallocate(tmp)

  endif

  call wall_time(t1)
  print*, " Wall time for noL_2e (min) = ", (t1 - t0)/60.d0

END_PROVIDER

! ---

