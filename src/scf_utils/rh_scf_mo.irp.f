! ---

subroutine Roothaan_Hall_SCF_MO()

  BEGIN_DOC
  !
  ! Roothaan-Hall algorithm for SCF Hartree-Fock calculation
  !
  END_DOC

  implicit none

  double precision              :: energy_SCF, energy_SCF_previous, Delta_energy_SCF
  double precision              :: max_error_DIIS
  double precision, allocatable :: Fock_matrix_DIIS(:,:,:), error_matrix_DIIS(:,:,:)

  integer                       :: iteration_SCF, dim_DIIS, index_dim_DIIS

  integer                       :: i, j
  double precision              :: level_shift_save
  double precision, allocatable :: mo_coef_save(:,:)

  logical, external             :: qp_stop

  PROVIDE ao_md5 mo_occ level_shift

  allocate( mo_coef_save(ao_num,mo_num)                   &
          , Fock_matrix_DIIS (mo_num,mo_num,max_dim_DIIS) &
          , error_matrix_DIIS(mo_num,mo_num,max_dim_DIIS) )

  Fock_matrix_DIIS  = 0.d0
  error_matrix_DIIS = 0.d0
  mo_coef_save      = 0.d0

  call write_time(6)

  print*,'energy of the guess = ',SCF_energy
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '  N ', 'energy  ', 'energy diff  ',  'DIIS error  ', 'Level shift   '
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'

! Initialize energies and density matrices
  energy_SCF_previous = SCF_energy
  Delta_energy_SCF    = 1.d0
  iteration_SCF       = 0
  dim_DIIS            = 0
  max_error_DIIS      = 1.d0


!
! Start of main SCF loop
!
  PROVIDE Fock_matrix_mo error_diis_Fmo

  do while ( &
    ( (max_error_DIIS > threshold_DIIS_nonzero) .or. &
      (dabs(Delta_energy_SCF) > thresh_SCF) &
    ) .and. (iteration_SCF < n_it_SCF_max) )

    iteration_SCF += 1
    if(frozen_orb_scf) then
     call initialize_mo_coef_begin_iteration
    endif

    dim_DIIS = min(dim_DIIS+1, max_dim_DIIS)

    if( (scf_algorithm == 'DIIS_MO').and.(dabs(Delta_energy_SCF) > 1.d-6))  then
    !if(scf_algorithm == 'DIIS_MO') then

      index_dim_DIIS = mod(dim_DIIS-1, max_dim_DIIS) + 1
      do j = 1, mo_num
        do i = 1, mo_num
          Fock_matrix_DIIS (i,j,index_dim_DIIS) = Fock_matrix_mo(i,j)
          error_matrix_DIIS(i,j,index_dim_DIIS) = error_diis_Fmo(i,j)
        enddo
      enddo

      call extrapolate_Fock_matrix_mo(error_matrix_DIIS, Fock_matrix_DIIS, Fock_matrix_mo, size(Fock_matrix_mo, 1), iteration_SCF, dim_DIIS)
      do i = 1, mo_num
        Fock_matrix_diag_mo(i) = Fock_matrix_mo(i,i)
      enddo
      TOUCH Fock_matrix_mo fock_matrix_diag_mo
    endif

    mo_coef = eigenvectors_Fock_matrix_mo
    if(frozen_orb_scf) then
      call reorder_core_orb
      call initialize_mo_coef_begin_iteration
    endif

    TOUCH mo_coef

    max_error_DIIS = maxval(Abs(error_diis_Fmo))

    energy_SCF = SCF_energy
    Delta_energy_SCF = energy_SCF - energy_SCF_previous

    if( (SCF_algorithm == 'DIIS_MO') .and. (Delta_energy_SCF > 0.d0) ) then
      Fock_matrix_MO(1:mo_num,1:mo_num) = Fock_matrix_DIIS(1:mo_num,1:mo_num,index_dim_DIIS)
      do i = 1, mo_num
        Fock_matrix_diag_mo(i) = Fock_matrix_mo(i,i)
      enddo
      TOUCH Fock_matrix_mo fock_matrix_diag_mo
      mo_coef = eigenvectors_Fock_matrix_mo
      max_error_DIIS = maxval(Abs(error_diis_Fmo))
      energy_SCF = SCF_energy
      Delta_energy_SCF = energy_SCF - energy_SCF_previous
    endif

    level_shift_save = level_shift
    mo_coef_save(1:ao_num,1:mo_num) = mo_coef(1:ao_num,1:mo_num)
    do while(Delta_energy_SCF > 0.d0)
      mo_coef(1:ao_num,1:mo_num) = mo_coef_save(1:ao_num,1:mo_num)
      if(level_shift <= .1d0) then
        level_shift = 1.d0
      else
        level_shift = level_shift * 3.0d0
      endif
      TOUCH mo_coef level_shift
      mo_coef(1:ao_num,1:mo_num) = eigenvectors_Fock_matrix_mo(1:ao_num,1:mo_num)
      if(frozen_orb_scf) then
        call reorder_core_orb
        call initialize_mo_coef_begin_iteration
      endif
      TOUCH mo_coef
      Delta_energy_SCF = SCF_energy - energy_SCF_previous
      energy_SCF = SCF_energy
      if(level_shift-level_shift_save > 40.d0) then
        level_shift = level_shift_save * 4.d0
        SOFT_TOUCH level_shift
        exit
      endif

      dim_DIIS=0
    enddo

    level_shift = level_shift * 0.5d0
    SOFT_TOUCH level_shift
    energy_SCF_previous = energy_SCF

!   Print results at the end of each iteration

    write(6,'(I4, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, I3)')  &
      iteration_SCF, energy_SCF, Delta_energy_SCF, max_error_DIIS, level_shift, dim_DIIS

    if(Delta_energy_SCF < 0.d0) then
      call save_mos
    endif

    if(qp_stop()) exit
  enddo

!
! End of Main SCF loop
!

  if(iteration_SCF < n_it_SCF_max) then
    mo_label = 'Canonical'
  endif

  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,*)

  if(.not.frozen_orb_scf)then
   call mo_as_eigvectors_of_mo_matrix(Fock_matrix_mo, size(Fock_matrix_mo, 1), size(Fock_matrix_mo, 2), mo_label, 1, .true.)
   call restore_symmetry(ao_num, mo_num, mo_coef, size(mo_coef, 1), 1.d-10)
   call orthonormalize_mos
   call save_mos
  endif

  call write_double(6, energy_SCF, 'SCF energy')

  call write_time(6)

end

! ---

subroutine extrapolate_Fock_matrix_mo(error_matrix_DIIS, Fock_matrix_DIIS, Fock_matrix_MO_, size_Fock_matrix_MO, iteration_SCF, dim_DIIS)

  BEGIN_DOC
  ! Compute the extrapolated Fock matrix using the DIIS procedure
  END_DOC

  implicit none

  integer,intent(inout)         :: dim_DIIS
  double precision,intent(in)   :: Fock_matrix_DIIS(mo_num,mo_num,dim_DIIS), error_matrix_DIIS(mo_num,mo_num,dim_DIIS)
  integer,intent(in)            :: iteration_SCF, size_Fock_matrix_MO
  double precision,intent(inout):: Fock_matrix_MO_(size_Fock_matrix_MO,mo_num)

  double precision,allocatable  :: B_matrix_DIIS(:,:),X_vector_DIIS(:)
  double precision,allocatable  :: C_vector_DIIS(:)

  double precision,allocatable  :: scratch(:,:)
  integer                       :: i,j,k,l,i_DIIS,j_DIIS
  double precision :: rcond, ferr, berr
  integer, allocatable :: iwork(:)
  integer :: lwork

  if(dim_DIIS < 1) then
    return
  endif

  allocate(                               &
    B_matrix_DIIS(dim_DIIS+1,dim_DIIS+1), &
    X_vector_DIIS(dim_DIIS+1),            &
    C_vector_DIIS(dim_DIIS+1),            &
    scratch(mo_num,mo_num)                &
  )

  ! Compute the matrices B and X
  B_matrix_DIIS(:,:) = 0.d0
  do j = 1, dim_DIIS
    j_DIIS = min(dim_DIIS, mod(iteration_SCF-j, max_dim_DIIS) + 1)

    do i = 1, dim_DIIS
      i_DIIS = min(dim_DIIS, mod(iteration_SCF-i, max_dim_DIIS) + 1)

      ! Compute product of two errors vectors
      do l = 1, mo_num
        do k = 1, mo_num
          B_matrix_DIIS(i,j) = B_matrix_DIIS(i,j) + error_matrix_DIIS(k,l,i_DIIS) * error_matrix_DIIS(k,l,j_DIIS)
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
  double precision :: anorm
  integer              :: info
  integer,allocatable  :: ipiv(:)
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

  call dgecon( '1', dim_DIIS+1, AF, size(AF, 1), anorm, rcond, scratch, iwork, info)
  if(info /= 0) then
    dim_DIIS = 0
    return
  endif

  if(rcond < 1.d-14) then
    dim_DIIS = 0
    return
  endif

  ! solve the linear system C = B.X

  X_vector_DIIS = C_vector_DIIS
  call dgesv(dim_DIIS+1 , 1, B_matrix_DIIS, size(B_matrix_DIIS, 1), ipiv, X_vector_DIIS, size(X_vector_DIIS, 1), info)

  deallocate(scratch, AF, iwork)

  if(info < 0) then
    stop 'bug in DIIS_MO'
  endif

  ! Compute extrapolated Fock matrix


  !$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(SHARED) if (mo_num > 200)
  do j = 1, mo_num
    do i = 1, mo_num
      Fock_matrix_MO_(i,j) = 0.d0
    enddo
    do k = 1, dim_DIIS
      if(dabs(X_vector_DIIS(k)) < 1.d-10) cycle
      do i = 1, mo_num
        ! FPE here
        Fock_matrix_MO_(i,j) = Fock_matrix_MO_(i,j) + X_vector_DIIS(k) * Fock_matrix_DIIS(i,j,dim_DIIS-k+1)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

end

