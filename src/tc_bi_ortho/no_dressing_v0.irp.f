
! ---

BEGIN_PROVIDER [double precision, no_0_v0]

  implicit none
  integer          :: i, j, k
  double precision :: tmp
  double precision :: I_ijk_ijk, I_ijk_kij, I_ijk_jik, I_ijk_jki, I_ijk_ikj, I_ijk_kji
  double precision :: t0, t1

  call wall_time(t0)
  print*, " Providing no_0_v0 ..."

  if(elec_alpha_num .eq. elec_beta_num) then

    tmp = 0.d0

    do i = 1, elec_beta_num
      do j = 1, elec_beta_num
        do k = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, k, i, j, I_ijk_kij)
          call give_integrals_3_body_bi_ort(i, j, k, j, i, k, I_ijk_jik)

          tmp = tmp + 4.d0 * (2.d0 * I_ijk_ijk + I_ijk_kij - 3.d0 * I_ijk_jik)
        enddo
      enddo
    enddo

    no_0_v0 = -1.d0 * (-tmp) / 6.d0

  else

    tmp = 0.d0

    do i = 1, elec_beta_num
      do j = 1, elec_beta_num
        do k = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, k, i, j, I_ijk_kij)
          call give_integrals_3_body_bi_ort(i, j, k, j, i, k, I_ijk_jik)

          tmp = tmp + 4.d0 * (2.d0 * I_ijk_ijk + I_ijk_kij - 3.d0 * I_ijk_jik)
        enddo
      enddo
    enddo

    do i = elec_beta_num+1, elec_alpha_num
      do j = elec_beta_num+1, elec_alpha_num
        do k = elec_beta_num+1, elec_alpha_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, k, i, j, I_ijk_kij)
          call give_integrals_3_body_bi_ort(i, j, k, j, i, k, I_ijk_jik)

          tmp = tmp + I_ijk_ijk + 2.d0 * I_ijk_kij - 3.d0 * I_ijk_jik
        enddo
      enddo
    enddo

    do i = elec_beta_num+1, elec_alpha_num
      do j = 1, elec_beta_num

        do k = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, j, k, i, I_ijk_jki)
          call give_integrals_3_body_bi_ort(i, j, k, i, k, j, I_ijk_ikj)
          call give_integrals_3_body_bi_ort(i, j, k, j, i, k, I_ijk_jik)
          call give_integrals_3_body_bi_ort(i, j, k, k, j, i, I_ijk_kji)

          tmp = tmp + 6.d0 * (2.d0 * I_ijk_ijk + I_ijk_jki - I_ijk_ikj - I_ijk_jik - I_ijk_kji)
        enddo

        do k = elec_beta_num+1, elec_alpha_num

          call give_integrals_3_body_bi_ort(i, j, k, i, j, k, I_ijk_ijk)
          call give_integrals_3_body_bi_ort(i, j, k, j, k, i, I_ijk_jki)
          call give_integrals_3_body_bi_ort(i, j, k, i, k, j, I_ijk_ikj)
          call give_integrals_3_body_bi_ort(i, j, k, j, i, k, I_ijk_jik)
          call give_integrals_3_body_bi_ort(i, j, k, k, j, i, I_ijk_kji)

          tmp = tmp + 3.d0 * (2.d0 * I_ijk_ijk + 2.d0 * I_ijk_jki - I_ijk_ikj - I_ijk_jik - 2.d0 * I_ijk_kji)
        enddo

      enddo
    enddo

    no_0_v0 = -1.d0 * (-tmp) / 6.d0

  endif

  call wall_time(t1)
  print*, " Wall time for no_0_v0 (sec) = ", (t1 - t0)/60.d0

  print*, " no_0_v0 = ", no_0_v0

END_PROVIDER

! ---


