! ---

subroutine rh_tcscf_simple()

  implicit none
  integer                       :: i, j, it, dim_DIIS
  double precision              :: t0, t1
  double precision              :: e_save, e_delta, rho_delta
  double precision              :: etc_tot, etc_1e, etc_2e, etc_3e, tc_grad
  double precision              :: er_DIIS
  double precision, allocatable :: rho_old(:,:), rho_new(:,:)

  allocate(rho_old(ao_num,ao_num), rho_new(ao_num,ao_num))

  it       = 0
  e_save   = 0.d0
  dim_DIIS = 0

  ! ---

  if(.not. bi_ortho) then
   print *, ' grad_hermit = ', grad_hermit
   call save_good_hermit_tc_eigvectors
   TOUCH mo_coef 
   call save_mos
  endif

  ! ---

  if(bi_ortho) then

    PROVIDE level_shift_tcscf
    PROVIDE mo_l_coef mo_r_coef

    write(6, '(A4,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A4, 1X, A8)')              &
      '====', '================', '================', '================', '================', '================' &
            , '================', '================', '================', '====', '========'

    write(6, '(A4,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A4, 1X, A8)')              &
      ' it ', '  SCF TC Energy ', '      E(1e)     ', '      E(2e)     ', '      E(3e)     ', '   energy diff  ' &
            , '    gradient    ', '    DIIS error  ', '  level shift   ', 'DIIS', '  WT (m)'

    write(6, '(A4,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A4, 1X, A8)')              &
      '====', '================', '================', '================', '================', '================' &
            , '================', '================', '================', '====', '========'


    ! first iteration (HF orbitals)
    call wall_time(t0)

    etc_tot = TC_HF_energy
    etc_1e  = TC_HF_one_e_energy
    etc_2e  = TC_HF_two_e_energy
    etc_3e  = 0.d0
    if(three_body_h_tc) then
      etc_3e = diag_three_elem_hf
    endif
    tc_grad = grad_non_hermit
    er_DIIS = maxval(abs(FQS_SQF_mo))
    e_delta = dabs(etc_tot - e_save)
    e_save  = etc_tot

    call wall_time(t1)
    write(6, '(I4,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, I4,1X, F8.2)')  &
      it, etc_tot, etc_1e, etc_2e, etc_3e, e_delta, tc_grad, er_DIIS, level_shift_tcscf, dim_DIIS, (t1-t0)/60.d0

    do while(tc_grad .gt. dsqrt(thresh_tcscf))
      call wall_time(t0)

      it += 1
      if(it > n_it_tcscf_max) then
        print *, ' max of TCSCF iterations is reached ', n_it_TCSCF_max
        stop
      endif

      mo_l_coef = fock_tc_leigvec_ao
      mo_r_coef = fock_tc_reigvec_ao
      call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
      call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
      TOUCH mo_l_coef mo_r_coef

      etc_tot = TC_HF_energy
      etc_1e  = TC_HF_one_e_energy
      etc_2e  = TC_HF_two_e_energy
      etc_3e  = 0.d0
      if(three_body_h_tc) then
        etc_3e = diag_three_elem_hf
      endif
      tc_grad = grad_non_hermit
      er_DIIS = maxval(abs(FQS_SQF_mo))
      e_delta = dabs(etc_tot - e_save)
      e_save  = etc_tot

      call ezfio_set_tc_scf_bitc_energy(etc_tot)

      call wall_time(t1)
      write(6, '(I4,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, I4,1X, F8.2)')  &
        it, etc_tot, etc_1e, etc_2e, etc_3e, e_delta, tc_grad, er_DIIS, level_shift_tcscf, dim_DIIS, (t1-t0)/60.d0
    enddo

  else

   do while( (grad_hermit.gt.dsqrt(thresh_tcscf)) .and. (it.lt.n_it_tcscf_max) )
      print*,'grad_hermit = ',grad_hermit
      it += 1
      print *, 'iteration = ', it
      print *, '***'
      print *, 'TC HF total energy = ', TC_HF_energy
      print *, 'TC HF 1 e   energy = ', TC_HF_one_e_energy
      print *, 'TC HF 2 e   energy = ', TC_HF_two_e_energy
      print *, 'TC HF 3 body       = ', diag_three_elem_hf
      print *, '***'
      print *, ''
      call save_good_hermit_tc_eigvectors
      TOUCH mo_coef 
      call save_mos
    enddo

  endif

  print *, ' TCSCF Simple converged !'
  call print_energy_and_mos()

  deallocate(rho_old, rho_new)

end

! ---

