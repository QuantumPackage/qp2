! ---

subroutine rh_tcscf_diis()

  implicit none

  integer                       :: i, j, it
  integer                       :: dim_DIIS, index_dim_DIIS
  double precision              :: etc_tot, etc_1e, etc_2e, etc_3e, e_save, e_delta
  double precision              :: tc_grad, g_save, g_delta, g_delta_th
  double precision              :: level_shift_save, rate_th
  double precision              :: t0, t1
  double precision              :: er_DIIS, er_delta, er_save, er_delta_th
  double precision, allocatable :: F_DIIS(:,:,:), E_DIIS(:,:,:)
  double precision, allocatable :: mo_r_coef_save(:,:), mo_l_coef_save(:,:)

  logical, external             :: qp_stop

  it          = 0
  e_save      = 0.d0
  dim_DIIS    = 0
  g_delta_th  = 1d0
  er_delta_th = 1d0
  rate_th     = 100.d0 !0.01d0 !0.2d0

  allocate(mo_r_coef_save(ao_num,mo_num), mo_l_coef_save(ao_num,mo_num))
  mo_l_coef_save = 0.d0
  mo_r_coef_save = 0.d0

  allocate(F_DIIS(ao_num,ao_num,max_dim_DIIS_TCSCF), E_DIIS(ao_num,ao_num,max_dim_DIIS_TCSCF))
  F_DIIS = 0.d0
  E_DIIS = 0.d0

  call write_time(6)

  ! ---

  PROVIDE level_shift_TCSCF
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
  g_save  = tc_grad
  er_save = er_DIIS

  call wall_time(t1)
  write(6, '(I4,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, I4,1X, F8.2)')  &
    it, etc_tot, etc_1e, etc_2e, etc_3e, e_delta, tc_grad, er_DIIS, level_shift_tcscf, dim_DIIS, (t1-t0)/60.d0

  ! ---

  PROVIDE FQS_SQF_ao Fock_matrix_tc_ao_tot

  do while((tc_grad .gt. dsqrt(thresh_tcscf)) .and. (er_DIIS .gt. threshold_DIIS_nonzero_TCSCF))

    call wall_time(t0)

    it += 1
    if(it > n_it_TCSCF_max) then
      print *, ' max of TCSCF iterations is reached ', n_it_TCSCF_max
      stop
    endif

    dim_DIIS = min(dim_DIIS+1, max_dim_DIIS_TCSCF)

    ! ---

    if(dabs(e_delta) > 1.d-12) then

      index_dim_DIIS = mod(dim_DIIS-1, max_dim_DIIS_TCSCF) + 1
      do j = 1, ao_num
        do i = 1, ao_num
          F_DIIS(i,j,index_dim_DIIS) = Fock_matrix_tc_ao_tot(i,j)
          E_DIIS(i,j,index_dim_DIIS) = FQS_SQF_ao           (i,j)
        enddo
      enddo

      call extrapolate_TC_Fock_matrix(E_DIIS, F_DIIS, Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot, 1), it, dim_DIIS)

      call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot, 1) &
                            , Fock_matrix_tc_mo_tot, size(Fock_matrix_tc_mo_tot, 1) )
      TOUCH Fock_matrix_tc_mo_tot fock_matrix_tc_diag_mo_tot
    endif

    ! ---

    mo_l_coef(1:ao_num,1:mo_num) = fock_tc_leigvec_ao(1:ao_num,1:mo_num)
    mo_r_coef(1:ao_num,1:mo_num) = fock_tc_reigvec_ao(1:ao_num,1:mo_num)
    !call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
    !call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
    TOUCH mo_l_coef mo_r_coef

    ! ---

    g_delta  = grad_non_hermit         -  g_save
    er_delta = maxval(abs(FQS_SQF_mo)) - er_save

    !if((g_delta > rate_th * g_delta_th) .and. (er_delta > rate_th * er_delta_th) .and. (it > 1)) then
    if((g_delta > rate_th * g_delta_th) .and. (it > 1)) then
    !if((g_delta > 0.d0) .and. (it > 1)) then

      Fock_matrix_tc_ao_tot(1:ao_num,1:ao_num) = F_DIIS(1:ao_num,1:ao_num,index_dim_DIIS)
      call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot, 1) &
                            , Fock_matrix_tc_mo_tot, size(Fock_matrix_tc_mo_tot, 1) )
      TOUCH Fock_matrix_tc_mo_tot fock_matrix_tc_diag_mo_tot

      mo_l_coef(1:ao_num,1:mo_num) = fock_tc_leigvec_ao(1:ao_num,1:mo_num)
      mo_r_coef(1:ao_num,1:mo_num) = fock_tc_reigvec_ao(1:ao_num,1:mo_num)
      !call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
      !call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
      TOUCH mo_l_coef mo_r_coef

    endif

    ! ---

    g_delta  = grad_non_hermit         -  g_save
    er_delta = maxval(abs(FQS_SQF_mo)) - er_save

    mo_l_coef_save(1:ao_num,1:mo_num) = mo_l_coef(1:ao_num,1:mo_num)
    mo_r_coef_save(1:ao_num,1:mo_num) = mo_r_coef(1:ao_num,1:mo_num)

    !do while((g_delta > rate_th * g_delta_th) .and. (er_delta > rate_th * er_delta_th) .and. (it > 1))
    do while((g_delta > rate_th * g_delta_th) .and. (it > 1))
      print *, ' big or bad step : ', g_delta, rate_th * g_delta_th

      mo_l_coef(1:ao_num,1:mo_num) = mo_l_coef_save(1:ao_num,1:mo_num) 
      mo_r_coef(1:ao_num,1:mo_num) = mo_r_coef_save(1:ao_num,1:mo_num) 
      if(level_shift_TCSCF <= .1d0) then
        level_shift_TCSCF = 1.d0
      else
        level_shift_TCSCF = level_shift_TCSCF * 3.0d0
      endif
      TOUCH mo_l_coef mo_r_coef level_shift_TCSCF

      mo_l_coef(1:ao_num,1:mo_num) = fock_tc_leigvec_ao(1:ao_num,1:mo_num)
      mo_r_coef(1:ao_num,1:mo_num) = fock_tc_reigvec_ao(1:ao_num,1:mo_num)
      !call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
      !call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
      TOUCH mo_l_coef mo_r_coef

      g_delta  = grad_non_hermit         -  g_save
      er_delta = maxval(abs(FQS_SQF_mo)) - er_save

      if(level_shift_TCSCF - level_shift_save > 40.d0) then
        level_shift_TCSCF = level_shift_save * 4.d0
        SOFT_TOUCH level_shift_TCSCF
        exit
      endif

      dim_DIIS = 0
    enddo

    ! ---

    level_shift_TCSCF = level_shift_TCSCF * 0.5d0
    SOFT_TOUCH level_shift_TCSCF

    etc_tot = TC_HF_energy
    etc_1e  = TC_HF_one_e_energy
    etc_2e  = TC_HF_two_e_energy
    etc_3e  = 0.d0
    if(three_body_h_tc) then
      etc_3e = diag_three_elem_hf
    endif
    tc_grad  = grad_non_hermit
    er_DIIS  = maxval(abs(FQS_SQF_mo))
    e_delta  = dabs(etc_tot - e_save)
    g_delta  = tc_grad - g_save
    er_delta = er_DIIS - er_save
    
    e_save           = etc_tot
    g_save           = tc_grad
    level_shift_save = level_shift_TCSCF
    er_save          = er_DIIS

    g_delta_th  = dabs(tc_grad) ! g_delta)
    er_delta_th = dabs(er_DIIS) !er_delta)

    call wall_time(t1)
    write(6, '(I4,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, I4,1X, F8.2)')  &
      it, etc_tot, etc_1e, etc_2e, etc_3e, e_delta, tc_grad, er_DIIS, level_shift_tcscf, dim_DIIS, (t1-t0)/60.d0

    if(g_delta .lt. 0.d0) then
      call ezfio_set_tc_scf_bitc_energy(etc_tot)
      call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
      call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
    endif

    if(qp_stop()) exit
  enddo

  ! ---

  print *, ' TCSCF DIIS converged !'
  call print_energy_and_mos()

  call write_time(6)

  deallocate(mo_r_coef_save, mo_l_coef_save, F_DIIS, E_DIIS)

  call ezfio_set_tc_scf_bitc_energy(TC_HF_energy)
  call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
  call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)

end

! ---

subroutine extrapolate_TC_Fock_matrix(E_DIIS, F_DIIS, F_ao, size_F_ao, it, dim_DIIS)

  BEGIN_DOC
  !
  ! Compute the extrapolated Fock matrix using the DIIS procedure
  !
  ! e = \sum_i c_i e_i and \sum_i c_i = 1 
  ! ==> lagrange multiplier with L = |e|^2 - \lambda (\sum_i c_i = 1)
  !
  END_DOC

  implicit none

  integer,          intent(in)    :: it, size_F_ao
  integer,          intent(inout) :: dim_DIIS
  double precision, intent(in)    :: F_DIIS(ao_num,ao_num,dim_DIIS)
  double precision, intent(in)    :: E_DIIS(ao_num,ao_num,dim_DIIS)
  double precision, intent(inout) :: F_ao(size_F_ao,ao_num)

  double precision, allocatable   :: B_matrix_DIIS(:,:), X_vector_DIIS(:), C_vector_DIIS(:)

  integer                         :: i, j, k, l, i_DIIS, j_DIIS
  integer                         :: lwork
  double precision                :: rcond, ferr, berr
  integer,          allocatable   :: iwork(:)
  double precision, allocatable   :: scratch(:,:)

  if(dim_DIIS < 1) then
    return
  endif

  allocate( B_matrix_DIIS(dim_DIIS+1,dim_DIIS+1), X_vector_DIIS(dim_DIIS+1) &
          , C_vector_DIIS(dim_DIIS+1), scratch(ao_num,ao_num) )

  ! Compute the matrices B and X
  B_matrix_DIIS(:,:) = 0.d0
  do j = 1, dim_DIIS
    j_DIIS = min(dim_DIIS, mod(it-j, max_dim_DIIS_TCSCF)+1)

    do i = 1, dim_DIIS
      i_DIIS = min(dim_DIIS, mod(it-i, max_dim_DIIS_TCSCF)+1)

      ! Compute product of two errors vectors
      do l = 1, ao_num
        do k = 1, ao_num
          B_matrix_DIIS(i,j) = B_matrix_DIIS(i,j) + E_DIIS(k,l,i_DIIS) * E_DIIS(k,l,j_DIIS)
        enddo
      enddo

    enddo
  enddo

  ! Pad B matrix and build the X matrix

  C_vector_DIIS(:) = 0.d0
  do i = 1, dim_DIIS
    B_matrix_DIIS(i,dim_DIIS+1) = -1.d0
    B_matrix_DIIS(dim_DIIS+1,i) = -1.d0
  enddo
  C_vector_DIIS(dim_DIIS+1) = -1.d0

  deallocate(scratch)

  ! Estimate condition number of B
  integer                       :: info
  double precision              :: anorm
  integer,          allocatable :: ipiv(:)
  double precision, allocatable :: AF(:,:)
  double precision, external :: dlange

  lwork = max((dim_DIIS+1)**2, (dim_DIIS+1)*5)
  allocate(AF(dim_DIIS+1,dim_DIIS+1))
  allocate(ipiv(2*(dim_DIIS+1)), iwork(2*(dim_DIIS+1)) )
  allocate(scratch(lwork,1))
  scratch(:,1) = 0.d0

  anorm = dlange('1', dim_DIIS+1, dim_DIIS+1, B_matrix_DIIS, size(B_matrix_DIIS, 1), scratch(1,1))

  AF(:,:) = B_matrix_DIIS(:,:)
  call dgetrf(dim_DIIS+1, dim_DIIS+1, AF, size(AF, 1), ipiv, info)
  if(info /= 0) then
    dim_DIIS = 0
    return
  endif

  call dgecon('1', dim_DIIS+1, AF, size(AF, 1), anorm, rcond, scratch, iwork, info)
  if(info /= 0) then
    dim_DIIS = 0
    return
  endif

  if(rcond < 1.d-14) then
    dim_DIIS = 0
    return
  endif

  ! solve the linear system C = B x X

  X_vector_DIIS = C_vector_DIIS
  call dgesv(dim_DIIS+1, 1, B_matrix_DIIS, size(B_matrix_DIIS, 1), ipiv , X_vector_DIIS, size(X_vector_DIIS, 1), info)

  deallocate(scratch, AF, iwork)
  if(info < 0) then
    stop ' bug in TC-DIIS'
  endif

  ! Compute extrapolated Fock matrix

  !$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(SHARED) if (ao_num > 200)
  do j = 1, ao_num
    do i = 1, ao_num
      F_ao(i,j) = 0.d0
    enddo
    do k = 1, dim_DIIS
      if(dabs(X_vector_DIIS(k)) < 1.d-10) cycle
      do i = 1,ao_num
        ! FPE here
        F_ao(i,j) = F_ao(i,j) + X_vector_DIIS(k) * F_DIIS(i,j,dim_DIIS-k+1)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

end

! ---

