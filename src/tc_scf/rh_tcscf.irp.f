! ---

subroutine rh_tcscf()

  BEGIN_DOC
  !
  ! Roothaan-Hall algorithm for TC-SCF calculation
  !
  END_DOC

  implicit none

  integer                       :: i, j
  integer                       :: iteration_TCSCF, dim_DIIS, index_dim_DIIS
  double precision              :: energy_TCSCF, energy_TCSCF_1e, energy_TCSCF_2e, energy_TCSCF_3e, gradie_TCSCF
  double precision              :: energy_TCSCF_previous, delta_energy_TCSCF
  double precision              :: gradie_TCSCF_previous, delta_gradie_TCSCF
  double precision              :: max_error_DIIS_TCSCF
  double precision              :: level_shift_save
  double precision              :: delta_energy_tmp, delta_gradie_tmp
  double precision, allocatable :: F_DIIS(:,:,:), e_DIIS(:,:,:)
  double precision, allocatable :: mo_r_coef_save(:,:), mo_l_coef_save(:,:)

  logical, external             :: qp_stop


  !PROVIDE ao_md5 mo_occ
  PROVIDE level_shift_TCSCF

  allocate( mo_r_coef_save(ao_num,mo_num), mo_l_coef_save(ao_num,mo_num) &
          , F_DIIS(ao_num,ao_num,max_dim_DIIS_TCSCF), e_DIIS(ao_num,ao_num,max_dim_DIIS_TCSCF) )

  F_DIIS         = 0.d0
  e_DIIS         = 0.d0
  mo_l_coef_save = 0.d0
  mo_r_coef_save = 0.d0

  call write_time(6)

  ! ---
  ! Initialize energies and density matrices

  energy_TCSCF_previous = TC_HF_energy
  energy_TCSCF_1e       = TC_HF_one_e_energy
  energy_TCSCF_2e       = TC_HF_two_e_energy
  energy_TCSCF_3e       = 0.d0
  if(three_body_h_tc) then
    energy_TCSCF_3e     = diag_three_elem_hf
  endif
  gradie_TCSCF_previous = grad_non_hermit
  delta_energy_TCSCF    = 1.d0
  delta_gradie_TCSCF    = 1.d0
  iteration_TCSCF       = 0
  dim_DIIS              = 0
  max_error_DIIS_TCSCF  = 1.d0

  ! ---

  ! Start of main SCF loop

  PROVIDE FQS_SQF_ao Fock_matrix_tc_ao_tot

  do while( (max_error_DIIS_TCSCF > threshold_DIIS_nonzero_TCSCF) .or. &
            !(dabs(delta_energy_TCSCF) > thresh_TCSCF)             .or. &
            (dabs(gradie_TCSCF_previous) > dsqrt(thresh_TCSCF))        )

    iteration_TCSCF += 1
    if(iteration_TCSCF > n_it_TCSCF_max) then
      print *, ' max of TCSCF iterations is reached ', n_it_TCSCF_max
      stop
    endif

    dim_DIIS = min(dim_DIIS+1, max_dim_DIIS_TCSCF)

    ! ---

    if((tcscf_algorithm == 'DIIS') .and. (dabs(delta_energy_TCSCF) > 1.d-6))  then

      ! store Fock and error matrices at each iteration
      index_dim_DIIS = mod(dim_DIIS-1, max_dim_DIIS_TCSCF) + 1
      do j = 1, ao_num
        do i = 1, ao_num
          F_DIIS(i,j,index_dim_DIIS) = Fock_matrix_tc_ao_tot(i,j)
          e_DIIS(i,j,index_dim_DIIS) = FQS_SQF_ao(i,j)
        enddo
      enddo

      call extrapolate_TC_Fock_matrix(e_DIIS, F_DIIS, Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot, 1), iteration_TCSCF, dim_DIIS)

      Fock_matrix_tc_ao_alpha = 0.5d0 * Fock_matrix_tc_ao_tot
      Fock_matrix_tc_ao_beta  = 0.5d0 * Fock_matrix_tc_ao_tot
      !TOUCH Fock_matrix_tc_ao_alpha Fock_matrix_tc_ao_beta

      call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_alpha, size(Fock_matrix_tc_ao_alpha, 1) &
                            , Fock_matrix_tc_mo_alpha, size(Fock_matrix_tc_mo_alpha, 1) )
      call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_beta , size(Fock_matrix_tc_ao_beta , 1) &
                            , Fock_matrix_tc_mo_beta , size(Fock_matrix_tc_mo_beta , 1) )
      TOUCH Fock_matrix_tc_mo_alpha Fock_matrix_tc_mo_beta
    endif

    ! ---

    mo_l_coef(1:ao_num,1:mo_num) = fock_tc_leigvec_ao(1:ao_num,1:mo_num)
    mo_r_coef(1:ao_num,1:mo_num) = fock_tc_reigvec_ao(1:ao_num,1:mo_num)
    TOUCH mo_l_coef mo_r_coef

    ! ---

    ! calculate error vectors
    max_error_DIIS_TCSCF = maxval(abs(FQS_SQF_mo))

    ! ---

    delta_energy_tmp = TC_HF_energy    - energy_TCSCF_previous
    delta_gradie_tmp = grad_non_hermit - gradie_TCSCF_previous

    ! ---

    do while((delta_gradie_tmp > 1.d-7) .and. (iteration_TCSCF > 1))
    !do while((dabs(delta_energy_tmp) > 0.5d0) .and. (iteration_TCSCF > 1))
      print *, ' very big or bad step  : ', delta_energy_tmp, delta_gradie_tmp
      print *, ' TC level shift = ', level_shift_TCSCF

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
      TOUCH mo_l_coef mo_r_coef

      delta_energy_tmp = TC_HF_energy    - energy_TCSCF_previous
      delta_gradie_tmp = grad_non_hermit - gradie_TCSCF_previous

      if(level_shift_TCSCF - level_shift_save > 40.d0) then
        level_shift_TCSCF = level_shift_save * 4.d0
        SOFT_TOUCH level_shift_TCSCF
        exit
      endif

      dim_DIIS = 0
    enddo
!    print *, ' very big step  : ', delta_energy_tmp
!    print *, ' TC level shift = ', level_shift_TCSCF

    ! ---

    level_shift_TCSCF = 0.d0
    !level_shift_TCSCF = level_shift_TCSCF * 0.5d0
    SOFT_TOUCH level_shift_TCSCF

    gradie_TCSCF       = grad_non_hermit
    energy_TCSCF       = TC_HF_energy
    energy_TCSCF_1e    = TC_HF_one_e_energy
    energy_TCSCF_2e    = TC_HF_two_e_energy
    energy_TCSCF_3e    = 0.d0
    if(three_body_h_tc) then
      energy_TCSCF_3e  = diag_three_elem_hf
    endif
    delta_energy_TCSCF = energy_TCSCF - energy_TCSCF_previous
    delta_gradie_TCSCF = gradie_TCSCF - gradie_TCSCF_previous

    energy_TCSCF_previous = energy_TCSCF
    gradie_TCSCF_previous = gradie_TCSCF


    level_shift_save = level_shift_TCSCF
    mo_l_coef_save(1:ao_num,1:mo_num) = mo_l_coef(1:ao_num,1:mo_num)
    mo_r_coef_save(1:ao_num,1:mo_num) = mo_r_coef(1:ao_num,1:mo_num)


    print *, ' iteration         = ', iteration_TCSCF
    print *, ' total TC energy   = ', energy_TCSCF 
    print *, ' 1-e   TC energy   = ', energy_TCSCF_1e
    print *, ' 2-e   TC energy   = ', energy_TCSCF_2e
    print *, ' 3-e   TC energy   = ', energy_TCSCF_3e
    print *, ' |delta TC energy| = ', dabs(delta_energy_TCSCF)
    print *, ' TC gradient       = ', gradie_TCSCF
    print *, ' delta TC gradient = ', delta_gradie_TCSCF
    print *, ' max TC DIIS error = ', max_error_DIIS_TCSCF 
    print *, ' TC DIIS dim       = ', dim_DIIS
    print *, ' TC level shift    = ', level_shift_TCSCF
    print *, ' '

    call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
    call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)

    if(qp_stop()) exit
  enddo

  ! ---

  print *, ' TCSCF DIIS converged !'
  call print_energy_and_mos()

  call write_time(6)

  deallocate(mo_r_coef_save, mo_l_coef_save, F_DIIS, e_DIIS)

end

! ---

subroutine extrapolate_TC_Fock_matrix(e_DIIS, F_DIIS, F_ao, size_F_ao, iteration_TCSCF, dim_DIIS)

  BEGIN_DOC
  !
  ! Compute the extrapolated Fock matrix using the DIIS procedure
  !
  ! e = \sum_i c_i e_i and \sum_i c_i = 1 
  ! ==> lagrange multiplier with L = |e|^2 - \lambda (\sum_i c_i = 1)
  !
  END_DOC

  implicit none

  integer,          intent(in)    :: iteration_TCSCF, size_F_ao
  integer,          intent(inout) :: dim_DIIS
  double precision, intent(in)    :: F_DIIS(ao_num,ao_num,dim_DIIS)
  double precision, intent(in)    :: e_DIIS(ao_num,ao_num,dim_DIIS)
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
    j_DIIS = min(dim_DIIS, mod(iteration_TCSCF-j, max_dim_DIIS_TCSCF)+1)

    do i = 1, dim_DIIS
      i_DIIS = min(dim_DIIS, mod(iteration_TCSCF-i, max_dim_DIIS_TCSCF)+1)

      ! Compute product of two errors vectors
      do l = 1, ao_num
        do k = 1, ao_num
          B_matrix_DIIS(i,j) = B_matrix_DIIS(i,j) + e_DIIS(k,l,i_DIIS) * e_DIIS(k,l,j_DIIS)
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

