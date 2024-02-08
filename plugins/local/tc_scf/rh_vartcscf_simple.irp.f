! ---

subroutine rh_vartcscf_simple()

  implicit none
  integer          :: i, j, it, dim_DIIS
  double precision :: t0, t1
  double precision :: e_save, e_delta, rho_delta
  double precision :: etc_tot, etc_1e, etc_2e, etc_3e
  double precision :: er_DIIS


  it       = 0
  e_save   = 0.d0
  dim_DIIS = 0

  ! ---

  PROVIDE level_shift_tcscf
  PROVIDE mo_r_coef

  write(6, '(A4,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A4, 1X, A8)')                      &
    '====', '================', '================', '================', '================', '================' &
          , '================', '================', '====', '========'
  write(6, '(A4,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A4, 1X, A8)')                      &
    ' it ', '  SCF TC Energy ', '      E(1e)     ', '      E(2e)     ', '      E(3e)     ', '   energy diff  ' &
          , '    DIIS error  ', '  level shift   ', 'DIIS', '  WT (m)'
  write(6, '(A4,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A4, 1X, A8)')                      &
    '====', '================', '================', '================', '================', '================' &
          , '================', '================', '====', '========'


  ! first iteration (HF orbitals)
  call wall_time(t0)

  etc_tot = VARTC_HF_energy
  etc_1e  = VARTC_HF_one_e_energy
  etc_2e  = VARTC_HF_two_e_energy
  etc_3e  = 0.d0
  if(three_body_h_tc) then
    etc_3e = diag_three_elem_hf
  endif
  er_DIIS = maxval(abs(FQS_SQF_mo))
  e_delta = dabs(etc_tot - e_save)
  e_save  = etc_tot

  call wall_time(t1)
  write(6, '(I4,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, I4,1X, F8.2)')  &
    it, etc_tot, etc_1e, etc_2e, etc_3e, e_delta, er_DIIS, level_shift_tcscf, dim_DIIS, (t1-t0)/60.d0

  do while(er_DIIS .gt. dsqrt(thresh_tcscf))
    call wall_time(t0)

    it += 1
    if(it > n_it_tcscf_max) then
      print *, ' max of TCSCF iterations is reached ', n_it_TCSCF_max
      stop
    endif

    mo_r_coef = fock_vartc_eigvec_ao
    mo_l_coef = mo_r_coef
    call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
    call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
    TOUCH mo_l_coef mo_r_coef

    etc_tot = VARTC_HF_energy
    etc_1e  = VARTC_HF_one_e_energy
    etc_2e  = VARTC_HF_two_e_energy
    etc_3e  = 0.d0
    if(three_body_h_tc) then
      etc_3e = diag_three_elem_hf
    endif
    er_DIIS = maxval(abs(FQS_SQF_mo))
    e_delta = dabs(etc_tot - e_save)
    e_save  = etc_tot

    call ezfio_set_tc_scf_bitc_energy(etc_tot)

    call wall_time(t1)
    write(6, '(I4,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, I4,1X, F8.2)')  &
      it, etc_tot, etc_1e, etc_2e, etc_3e, e_delta, er_DIIS, level_shift_tcscf, dim_DIIS, (t1-t0)/60.d0
  enddo

  print *, ' VAR-TCSCF Simple converged !'

end

! ---

