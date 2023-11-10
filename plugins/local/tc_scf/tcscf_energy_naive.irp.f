
! ---

BEGIN_PROVIDER [double precision, tcscf_energy_3e_naive]

  implicit none
  integer          :: i, j, k
  integer          :: neu, ned, D(elec_num)
  integer          :: ii, jj, kk
  integer          :: si, sj, sk
  double precision :: I_ijk, I_jki, I_kij, I_jik, I_ikj, I_kji
  double precision :: I_tot

  PROVIDE mo_l_coef mo_r_coef

  neu = elec_alpha_num
  ned = elec_beta_num
  if (neu > 0) D(1:neu) = [(2*i-1, i = 1, neu)]
  if (ned > 0) D(neu+1:neu+ned) = [(2*i, i = 1, ned)]

  !print*, "D = "
  !do i = 1, elec_num
  !  ii = (D(i) - 1) / 2 + 1
  !  si = mod(D(i), 2)
  !  print*, i, D(i), ii, si
  !enddo

  tcscf_energy_3e_naive = 0.d0

  do i = 1, elec_num - 2
    ii = (D(i) - 1) / 2 + 1
    si = mod(D(i), 2)

    do j = i + 1, elec_num - 1
      jj = (D(j) - 1) / 2 + 1
      sj = mod(D(j), 2)

      do k = j + 1, elec_num
        kk = (D(k) - 1) / 2 + 1
        sk = mod(D(k), 2)

        call give_integrals_3_body_bi_ort(ii, jj, kk, ii, jj, kk, I_ijk)
        I_tot = I_ijk

        if(sj==si .and. sk==sj) then
          call give_integrals_3_body_bi_ort(ii, jj, kk, jj, kk, ii, I_jki)
          I_tot += I_jki
        endif

        if(sk==si .and. si==sj) then
          call give_integrals_3_body_bi_ort(ii, jj, kk, kk, ii, jj, I_kij)
          I_tot += I_kij
        endif

        if(sj==si) then
          call give_integrals_3_body_bi_ort(ii, jj, kk, jj, ii, kk, I_jik)
          I_tot -= I_jik
        endif

        if(sk==sj) then
          call give_integrals_3_body_bi_ort(ii, jj, kk, ii, kk, jj, I_ikj)
          I_tot -= I_ikj
        endif

        if(sk==si) then
          call give_integrals_3_body_bi_ort(ii, jj, kk, kk, jj, ii, I_kji)
          I_tot -= I_kji
        endif

        tcscf_energy_3e_naive += I_tot
      enddo
    enddo
  enddo

  tcscf_energy_3e_naive = -tcscf_energy_3e_naive

END_PROVIDER

! ---

