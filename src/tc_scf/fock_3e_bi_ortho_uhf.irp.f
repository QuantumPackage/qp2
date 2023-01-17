
! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_cs, (mo_num, mo_num)]

  implicit none
  integer          :: a, b, i, j
  double precision :: I_bij_aij, I_bij_ija, I_bij_jai, I_bij_aji, I_bij_iaj, I_bij_jia
  double precision :: ti, tf

  PROVIDE mo_l_coef mo_r_coef

  !print *, ' PROVIDING fock_3e_uhf_mo_cs ...'
  call wall_time(ti)

  fock_3e_uhf_mo_cs = 0.d0

  do a = 1, mo_num
    do b = 1, mo_num
   
      do j = 1, elec_beta_num
        do i = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          fock_3e_uhf_mo_cs(b,a) -= 0.5d0 * ( 4.d0 * I_bij_aij &
                                            +        I_bij_ija &
                                            +        I_bij_jai &
                                            - 2.d0 * I_bij_aji &
                                            - 2.d0 * I_bij_iaj &
                                            - 2.d0 * I_bij_jia )

        enddo
      enddo
    enddo
  enddo

  call wall_time(tf)
  !print *, ' total Wall time for fock_3e_uhf_mo_cs =', tf - ti

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_a, (mo_num, mo_num)]

  implicit none
  integer          :: a, b, i, j, o
  double precision :: I_bij_aij, I_bij_ija, I_bij_jai, I_bij_aji, I_bij_iaj, I_bij_jia
  double precision :: ti, tf

  PROVIDE mo_l_coef mo_r_coef

  !print *, ' PROVIDING fock_3e_uhf_mo_a ...'
  call wall_time(ti)

  o = elec_beta_num + 1

  fock_3e_uhf_mo_a = fock_3e_uhf_mo_cs

  do a = 1, mo_num
    do b = 1, mo_num

      ! ---

      do j = o, elec_alpha_num
        do i = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          fock_3e_uhf_mo_a(b,a) -= 0.5d0 * ( 2.d0 * I_bij_aij &
                                           +        I_bij_ija &
                                           +        I_bij_jai &
                                           -        I_bij_aji &
                                           -        I_bij_iaj &
                                           - 2.d0 * I_bij_jia )

        enddo
      enddo

      ! ---

      do j = 1, elec_beta_num
        do i = o, elec_alpha_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          fock_3e_uhf_mo_a(b,a) -= 0.5d0 * ( 2.d0 * I_bij_aij &
                                           +        I_bij_ija &
                                           +        I_bij_jai &
                                           -        I_bij_aji &
                                           - 2.d0 * I_bij_iaj &
                                           -        I_bij_jia )

        enddo
      enddo

      ! ---

      do j = o, elec_alpha_num
        do i = o, elec_alpha_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          fock_3e_uhf_mo_a(b,a) -= 0.5d0 * ( I_bij_aij &
                                           + I_bij_ija &
                                           + I_bij_jai &
                                           - I_bij_aji &
                                           - I_bij_iaj &
                                           - I_bij_jia )

        enddo
      enddo

      ! ---

    enddo
  enddo

  call wall_time(tf)
  !print *, ' total Wall time for fock_3e_uhf_mo_a =', tf - ti

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_b, (mo_num, mo_num)]

  implicit none
  integer          :: a, b, i, j, o
  double precision :: I_bij_aij, I_bij_ija, I_bij_jai, I_bij_aji, I_bij_iaj, I_bij_jia
  double precision :: ti, tf

  PROVIDE mo_l_coef mo_r_coef

  !print *, ' PROVIDING fock_3e_uhf_mo_b ...'
  call wall_time(ti)

  o = elec_beta_num + 1

  fock_3e_uhf_mo_b = fock_3e_uhf_mo_cs

  do a = 1, mo_num
    do b = 1, mo_num

      ! ---

      do j = o, elec_alpha_num
        do i = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          fock_3e_uhf_mo_b(b,a) -= 0.5d0 * ( 2.d0 * I_bij_aij &
                                           -        I_bij_aji &
                                           -        I_bij_iaj )

        enddo
      enddo

      ! ---

      do j = 1, elec_beta_num
        do i = o, elec_alpha_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          fock_3e_uhf_mo_b(b,a) -= 0.5d0 * ( 2.d0 * I_bij_aij &
                                           -        I_bij_aji &
                                           -        I_bij_jia )

        enddo
      enddo

      ! ---

      do j = o, elec_alpha_num
        do i = o, elec_alpha_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          fock_3e_uhf_mo_b(b,a) -= 0.5d0 * ( I_bij_aij &
                                           - I_bij_aji )

        enddo
      enddo

      ! ---

    enddo
  enddo

  call wall_time(tf)
  !print *, ' total Wall time for fock_3e_uhf_mo_b =', tf - ti

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_ao_a, (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! Equations (B6) and (B7)
  !
  ! g <--> gamma
  ! d <--> delta
  ! e <--> eta
  ! k <--> kappa
  !
  END_DOC

  implicit none
  integer                       :: g, d, e, k, mu, nu
  double precision              :: dm_ge_a, dm_ge_b, dm_ge
  double precision              :: dm_dk_a, dm_dk_b, dm_dk
  double precision              :: i_mugd_nuek, i_mugd_eknu, i_mugd_knue, i_mugd_nuke, i_mugd_enuk, i_mugd_kenu
  double precision              :: ti, tf
  double precision, allocatable :: f_tmp(:,:)

  print *, ' PROVIDING fock_3e_uhf_ao_a ...'
  call wall_time(ti)

  fock_3e_uhf_ao_a = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                                &
 !$OMP PRIVATE (g, e, d, k, mu, nu, dm_ge_a, dm_ge_b, dm_ge, dm_dk_a, dm_dk_b, dm_dk, f_tmp,  &
 !$OMP          i_mugd_nuek, i_mugd_eknu, i_mugd_knue, i_mugd_nuke, i_mugd_enuk, i_mugd_kenu) &
 !$OMP SHARED  (ao_num, TCSCF_bi_ort_dm_ao_alpha, TCSCF_bi_ort_dm_ao_beta, fock_3e_uhf_ao_a)

  allocate(f_tmp(ao_num,ao_num))
  f_tmp = 0.d0

 !$OMP DO
  do g = 1, ao_num
    do e = 1, ao_num
      dm_ge_a = TCSCF_bi_ort_dm_ao_alpha(g,e)
      dm_ge_b = TCSCF_bi_ort_dm_ao_beta (g,e)
      dm_ge   = dm_ge_a + dm_ge_b
      do d = 1, ao_num
        do k = 1, ao_num
          dm_dk_a = TCSCF_bi_ort_dm_ao_alpha(d,k)
          dm_dk_b = TCSCF_bi_ort_dm_ao_beta (d,k)
          dm_dk   = dm_dk_a + dm_dk_b
          do mu = 1, ao_num
            do nu = 1, ao_num
              call give_integrals_3_body_bi_ort_ao(mu, g, d, nu, e, k, i_mugd_nuek)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, e, k, nu, i_mugd_eknu)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, k, nu, e, i_mugd_knue)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, nu, k, e, i_mugd_nuke)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, e, nu, k, i_mugd_enuk)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, k, e, nu, i_mugd_kenu)
              f_tmp(mu,nu) -= 0.5d0 * ( dm_ge   * dm_dk   * i_mugd_nuek &
                                      + dm_ge_a * dm_dk_a * i_mugd_eknu &
                                      + dm_ge_a * dm_dk_a * i_mugd_knue &
                                      - dm_ge_a * dm_dk   * i_mugd_enuk &
                                      - dm_ge   * dm_dk_a * i_mugd_kenu &
                                      - dm_ge_a * dm_dk_a * i_mugd_nuke &
                                      - dm_ge_b * dm_dk_b * i_mugd_nuke )
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO NOWAIT

 !$OMP CRITICAL
  do mu = 1, ao_num
    do nu = 1, ao_num
      fock_3e_uhf_ao_a(mu,nu) += f_tmp(mu,nu)
    enddo
  enddo
 !$OMP END CRITICAL

  deallocate(f_tmp)
 !$OMP END PARALLEL

  call wall_time(tf)
  print *, ' total Wall time for fock_3e_uhf_ao_a =', tf - ti

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_ao_b, (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! Equations (B6) and (B7)
  !
  ! g <--> gamma
  ! d <--> delta
  ! e <--> eta
  ! k <--> kappa
  !
  END_DOC

  implicit none
  integer                       :: g, d, e, k, mu, nu
  double precision              :: dm_ge_a, dm_ge_b, dm_ge
  double precision              :: dm_dk_a, dm_dk_b, dm_dk
  double precision              :: i_mugd_nuek, i_mugd_eknu, i_mugd_knue, i_mugd_nuke, i_mugd_enuk, i_mugd_kenu
  double precision              :: ti, tf
  double precision, allocatable :: f_tmp(:,:)

  print *, ' PROVIDING fock_3e_uhf_ao_b ...'
  call wall_time(ti)

  fock_3e_uhf_ao_b = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                                &
 !$OMP PRIVATE (g, e, d, k, mu, nu, dm_ge_a, dm_ge_b, dm_ge, dm_dk_a, dm_dk_b, dm_dk, f_tmp,  &
 !$OMP          i_mugd_nuek, i_mugd_eknu, i_mugd_knue, i_mugd_nuke, i_mugd_enuk, i_mugd_kenu) &
 !$OMP SHARED  (ao_num, TCSCF_bi_ort_dm_ao_alpha, TCSCF_bi_ort_dm_ao_beta, fock_3e_uhf_ao_b)

  allocate(f_tmp(ao_num,ao_num))
  f_tmp = 0.d0

 !$OMP DO
  do g = 1, ao_num
    do e = 1, ao_num
      dm_ge_a = TCSCF_bi_ort_dm_ao_alpha(g,e)
      dm_ge_b = TCSCF_bi_ort_dm_ao_beta (g,e)
      dm_ge   = dm_ge_a + dm_ge_b
      do d = 1, ao_num
        do k = 1, ao_num
          dm_dk_a = TCSCF_bi_ort_dm_ao_alpha(d,k)
          dm_dk_b = TCSCF_bi_ort_dm_ao_beta (d,k)
          dm_dk   = dm_dk_a + dm_dk_b
          do mu = 1, ao_num
            do nu = 1, ao_num
              call give_integrals_3_body_bi_ort_ao(mu, g, d, nu, e, k, i_mugd_nuek)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, e, k, nu, i_mugd_eknu)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, k, nu, e, i_mugd_knue)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, nu, k, e, i_mugd_nuke)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, e, nu, k, i_mugd_enuk)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, k, e, nu, i_mugd_kenu)
              f_tmp(mu,nu) -= 0.5d0 * ( dm_ge   * dm_dk   * i_mugd_nuek &
                                      + dm_ge_b * dm_dk_b * i_mugd_eknu &
                                      + dm_ge_b * dm_dk_b * i_mugd_knue &
                                      - dm_ge_b * dm_dk   * i_mugd_enuk &
                                      - dm_ge   * dm_dk_b * i_mugd_kenu &
                                      - dm_ge_b * dm_dk_b * i_mugd_nuke &
                                      - dm_ge_a * dm_dk_a * i_mugd_nuke )
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO NOWAIT

 !$OMP CRITICAL
  do mu = 1, ao_num
    do nu = 1, ao_num
      fock_3e_uhf_ao_b(mu,nu) += f_tmp(mu,nu)
    enddo
  enddo
 !$OMP END CRITICAL

  deallocate(f_tmp)
 !$OMP END PARALLEL

  call wall_time(tf)
  print *, ' total Wall time for fock_3e_uhf_ao_b =', tf - ti

END_PROVIDER 

! ---

