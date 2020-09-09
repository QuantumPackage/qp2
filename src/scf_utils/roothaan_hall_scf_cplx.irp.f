subroutine Roothaan_Hall_SCF_complex

BEGIN_DOC
! Roothaan-Hall algorithm for SCF Hartree-Fock calculation
END_DOC

  implicit none

  double precision               :: energy_SCF,energy_SCF_previous,Delta_energy_SCF
  double precision               :: max_error_DIIS,max_error_DIIS_alpha,max_error_DIIS_beta
  complex*16, allocatable  :: Fock_matrix_DIIS(:,:,:),error_matrix_DIIS(:,:,:)

  integer                        :: iteration_SCF,dim_DIIS,index_dim_DIIS

  integer                        :: i,j
  logical, external              :: qp_stop
  complex*16, allocatable :: mo_coef_save(:,:)

  PROVIDE ao_md5 mo_occ level_shift

  allocate(mo_coef_save(ao_num,mo_num),                          &
      Fock_matrix_DIIS (ao_num,ao_num,max_dim_DIIS),                 &
      error_matrix_DIIS(ao_num,ao_num,max_dim_DIIS)                  &
      )

  call write_time(6)

  print*,'Energy of the guess = ',SCF_energy
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '  N ', 'Energy  ', 'Energy diff  ',  'DIIS error  ', 'Level shift   '
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
  PROVIDE FPS_SPF_matrix_AO_complex Fock_matrix_AO_complex 

  do while ( &
    ( (max_error_DIIS > threshold_DIIS_nonzero) .or. &
      (dabs(Delta_energy_SCF) > thresh_SCF) &
    ) .and. (iteration_SCF < n_it_SCF_max) )

! Increment cycle number

    iteration_SCF += 1
    if(frozen_orb_scf)then
     call initialize_mo_coef_begin_iteration
    endif

! Current size of the DIIS space

    dim_DIIS = min(dim_DIIS+1,max_dim_DIIS)

    if (scf_algorithm == 'DIIS') then

      ! Store Fock and error matrices at each iteration
      do j=1,ao_num
        do i=1,ao_num
          index_dim_DIIS = mod(dim_DIIS-1,max_dim_DIIS)+1
          Fock_matrix_DIIS (i,j,index_dim_DIIS) = Fock_matrix_AO_complex(i,j)
          error_matrix_DIIS(i,j,index_dim_DIIS) = FPS_SPF_matrix_AO_complex(i,j)
        enddo
      enddo

      ! Compute the extrapolated Fock matrix

      call extrapolate_Fock_matrix_complex(                            &
          error_matrix_DIIS,Fock_matrix_DIIS,                          &
          Fock_matrix_AO_complex,size(Fock_matrix_AO_complex,1),       &
          iteration_SCF,dim_DIIS                                       &
          )

      Fock_matrix_AO_alpha_complex = Fock_matrix_AO_complex*0.5d0
      Fock_matrix_AO_beta_complex  = Fock_matrix_AO_complex*0.5d0
      TOUCH Fock_matrix_AO_alpha_complex Fock_matrix_AO_beta_complex

    endif

    mo_coef_complex = eigenvectors_fock_matrix_mo_complex
    if(frozen_orb_scf)then
     call reorder_core_orb
     call initialize_mo_coef_begin_iteration
    endif

    TOUCH mo_coef_complex

!   Calculate error vectors

    max_error_DIIS = maxval(cdabs(FPS_SPF_Matrix_MO_complex))

!   SCF energy
!    call print_debug_scf_complex
    energy_SCF = scf_energy
    Delta_Energy_SCF = energy_SCF - energy_SCF_previous
    if ( (SCF_algorithm == 'DIIS').and.(Delta_Energy_SCF > 0.d0) ) then
      Fock_matrix_AO_complex(1:ao_num,1:ao_num) = Fock_matrix_DIIS (1:ao_num,1:ao_num,index_dim_DIIS)
      Fock_matrix_AO_alpha_complex = Fock_matrix_AO_complex*0.5d0
      Fock_matrix_AO_beta_complex  = Fock_matrix_AO_complex*0.5d0
      TOUCH Fock_matrix_AO_alpha_complex Fock_matrix_AO_beta_complex
    endif

    double precision :: level_shift_save
    level_shift_save = level_shift
    mo_coef_save(1:ao_num,1:mo_num) = mo_coef_complex(1:ao_num,1:mo_num)
    do while (Delta_energy_SCF > 0.d0)
      mo_coef_complex(1:ao_num,1:mo_num) = mo_coef_save
      if (level_shift <= .1d0) then
        level_shift = 1.d0
      else
        level_shift = level_shift * 3.0d0
      endif
      TOUCH mo_coef_complex level_shift
      mo_coef_complex(1:ao_num,1:mo_num) = eigenvectors_fock_matrix_mo_complex(1:ao_num,1:mo_num)
      if(frozen_orb_scf)then
        call reorder_core_orb
        call initialize_mo_coef_begin_iteration
      endif
      TOUCH mo_coef_complex
      Delta_Energy_SCF = SCF_energy - energy_SCF_previous
      energy_SCF = SCF_energy
      if (level_shift-level_shift_save > 40.d0) then
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
      iteration_SCF, energy_scf, Delta_energy_SCF, max_error_DIIS, level_shift, dim_DIIS

    if (Delta_energy_SCF < 0.d0) then
      call save_mos
    endif
    if (qp_stop()) exit

  enddo

 if (iteration_SCF < n_it_SCF_max) then
   mo_label = "Canonical"
 endif
!
! End of Main SCF loop
!

  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,*)

  if(.not.frozen_orb_scf)then
   call mo_as_eigvectors_of_mo_matrix_complex(Fock_matrix_mo_complex,size(Fock_matrix_mo_complex,1),size(Fock_matrix_mo_complex,2),mo_label,1,.true.)
   call save_mos
  endif

  call write_double(6, Energy_SCF, 'SCF energy')

  call write_time(6)

end

subroutine extrapolate_Fock_matrix_complex(       &
           error_matrix_DIIS,Fock_matrix_DIIS,    &
           Fock_matrix_AO_,size_Fock_matrix_AO,   &
           iteration_SCF,dim_DIIS                 &
           )

BEGIN_DOC
! Compute the extrapolated Fock matrix using the DIIS procedure
END_DOC

  implicit none

  complex*16,intent(in)   :: Fock_matrix_DIIS(ao_num,ao_num,*),error_matrix_DIIS(ao_num,ao_num,*)
  integer,intent(in)            :: iteration_SCF, size_Fock_matrix_AO
  complex*16,intent(inout):: Fock_matrix_AO_(size_Fock_matrix_AO,ao_num)
  integer,intent(inout)         :: dim_DIIS

  double precision,allocatable  :: B_matrix_DIIS(:,:),X_vector_DIIS(:)
  double precision,allocatable  :: C_vector_DIIS(:)
  double precision              :: accum_im, thr_im
  complex*16,allocatable  :: scratch(:,:)
  integer                       :: i,j,k,i_DIIS,j_DIIS
  thr_im = 1.0d-10
  allocate(                               &
    B_matrix_DIIS(dim_DIIS+1,dim_DIIS+1), &
    X_vector_DIIS(dim_DIIS+1),            &
    C_vector_DIIS(dim_DIIS+1),            &
    scratch(ao_num,ao_num)                &
  )

! Compute the matrices B and X
  do j=1,dim_DIIS
    do i=1,dim_DIIS

      j_DIIS = mod(iteration_SCF-j,max_dim_DIIS)+1
      i_DIIS = mod(iteration_SCF-i,max_dim_DIIS)+1

! Compute product of two errors vectors

      call zgemm('N','N',ao_num,ao_num,ao_num,                      &
           (1.d0,0.d0),                                             &
           error_matrix_DIIS(1,1,i_DIIS),size(error_matrix_DIIS,1), &
           error_matrix_DIIS(1,1,j_DIIS),size(error_matrix_DIIS,1), &
           (0.d0,0.d0),                                             &
           scratch,size(scratch,1))

! Compute Trace

      B_matrix_DIIS(i,j) = 0.d0
      accum_im = 0.d0
      do k=1,ao_num
        B_matrix_DIIS(i,j) = B_matrix_DIIS(i,j) + dble(scratch(k,k))
        accum_im = accum_im + dimag(scratch(k,k))
      enddo
      if (dabs(accum_im) .gt. thr_im) then
        !stop 'problem with imaginary parts in DIIS B_matrix?'
        print*, 'problem with imaginary parts in DIIS B_matrix?',accum_im
      endif
    enddo
  enddo
  deallocate(scratch)
! Pad B matrix and build the X matrix

  do i=1,dim_DIIS
    B_matrix_DIIS(i,dim_DIIS+1) = -1.d0
    B_matrix_DIIS(dim_DIIS+1,i) = -1.d0
    C_vector_DIIS(i) = 0.d0
  enddo
  B_matrix_DIIS(dim_DIIS+1,dim_DIIS+1) = 0.d0
  C_vector_DIIS(dim_DIIS+1) = -1.d0

! Solve the linear system C = B.X

  integer              :: info
  integer,allocatable  :: ipiv(:)

  allocate(          &
    ipiv(dim_DIIS+1) &
  )

  double precision, allocatable :: AF(:,:),scratch_d1(:)
  allocate (AF(dim_DIIS+1,dim_DIIS+1),scratch_d1(1))
  double precision :: rcond, ferr, berr
  integer :: iwork(dim_DIIS+1), lwork

  call dsysvx('N','U',dim_DIIS+1,1,      &
    B_matrix_DIIS,size(B_matrix_DIIS,1), &
    AF, size(AF,1),                      &
    ipiv,                                &
    C_vector_DIIS,size(C_vector_DIIS,1), &
    X_vector_DIIS,size(X_vector_DIIS,1), &
    rcond,                               &
    ferr,                                &
    berr,                                &
    scratch_d1,-1,                       &
    iwork,                               &
    info                                 &
  )
  lwork = int(scratch_d1(1))
  deallocate(scratch_d1)
  allocate(scratch_d1(lwork))

  call dsysvx('N','U',dim_DIIS+1,1,      &
    B_matrix_DIIS,size(B_matrix_DIIS,1), &
    AF, size(AF,1),                      &
    ipiv,                                &
    C_vector_DIIS,size(C_vector_DIIS,1), &
    X_vector_DIIS,size(X_vector_DIIS,1), &
    rcond,                               &
    ferr,                                &
    berr,                                &
    scratch_d1,size(scratch_d1),         &
    iwork,                               &
    info                                 &
  )
  deallocate(scratch_d1,ipiv)

 if(info < 0) then
   stop 'bug in DIIS'
 endif

 if (rcond > 1.d-12) then

  ! Compute extrapolated Fock matrix


      !$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(SHARED) if (ao_num > 200)
      do j=1,ao_num
        do i=1,ao_num
          Fock_matrix_AO_(i,j) = (0.d0,0.d0)
        enddo
        do k=1,dim_DIIS
          do i=1,ao_num
            Fock_matrix_AO_(i,j) = Fock_matrix_AO_(i,j) +            &
                X_vector_DIIS(k)*Fock_matrix_DIIS(i,j,dim_DIIS-k+1)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

  else
    dim_DIIS = 0
  endif

end

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!

subroutine Roothaan_Hall_SCF_kpts

BEGIN_DOC
! Roothaan-Hall algorithm for SCF Hartree-Fock calculation
END_DOC

  implicit none

  double precision               :: energy_SCF,energy_SCF_previous,Delta_energy_SCF
  double precision               :: max_error_DIIS,max_error_DIIS_alpha,max_error_DIIS_beta
  complex*16, allocatable  :: Fock_matrix_DIIS(:,:,:,:),error_matrix_DIIS(:,:,:,:)

  integer                        :: iteration_SCF,dim_DIIS,index_dim_DIIS

  integer                        :: i,j,k,kk
  logical, external              :: qp_stop
  complex*16, allocatable :: mo_coef_save(:,:,:)

  PROVIDE ao_md5 mo_occ_kpts level_shift

  allocate(mo_coef_save(ao_num_per_kpt,mo_num_per_kpt,kpt_num),                          &
      Fock_matrix_DIIS (ao_num_per_kpt,ao_num_per_kpt,max_dim_DIIS,kpt_num),                 &
      error_matrix_DIIS(ao_num_per_kpt,ao_num_per_kpt,max_dim_DIIS,kpt_num)                  &
      )
      !todo: add kpt_num dim to diis mats? (3 or 4)
  call write_time(6)

  print*,'Energy of the guess = ',scf_energy
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '  N ', 'Energy  ', 'Energy diff  ',  'DIIS error  ', 'Level shift   '
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
  !PROVIDE fps_spf_matrix_ao_complex fock_matrix_ao_complex 
  PROVIDE fps_spf_matrix_ao_kpts fock_matrix_ao_kpts 

  do while ( &
    ( (max_error_DIIS > threshold_DIIS_nonzero) .or. &
      (dabs(Delta_energy_SCF) > thresh_SCF) &
    ) .and. (iteration_SCF < n_it_SCF_max) )

! Increment cycle number

    iteration_SCF += 1
    if(frozen_orb_scf)then
     call initialize_mo_coef_begin_iteration
    endif

! Current size of the DIIS space

    dim_DIIS = min(dim_DIIS+1,max_dim_DIIS)

    if (scf_algorithm == 'DIIS') then

      do kk=1,kpt_num
        ! Store Fock and error matrices at each iteration
        do j=1,ao_num_per_kpt
          do i=1,ao_num_per_kpt
            index_dim_DIIS = mod(dim_DIIS-1,max_dim_DIIS)+1
            Fock_matrix_DIIS (i,j,index_dim_DIIS,kk) = fock_matrix_ao_kpts(i,j,kk)
            error_matrix_DIIS(i,j,index_dim_DIIS,kk) = fps_spf_matrix_ao_kpts(i,j,kk)
          enddo
        enddo

        ! Compute the extrapolated Fock matrix

        call extrapolate_fock_matrix_kpts(                            &
            error_matrix_DIIS(1,1,1,kk),Fock_matrix_DIIS(1,1,1,kk),        &
            Fock_matrix_AO_kpts(1,1,kk),size(Fock_matrix_AO_kpts,1),       &
            iteration_SCF,dim_DIIS                                       &
            )
      enddo
      Fock_matrix_AO_alpha_kpts = Fock_matrix_AO_kpts*0.5d0
      Fock_matrix_AO_beta_kpts  = Fock_matrix_AO_kpts*0.5d0
      TOUCH Fock_matrix_AO_alpha_kpts Fock_matrix_AO_beta_kpts

    endif

    mo_coef_kpts = eigenvectors_fock_matrix_mo_kpts
    if(frozen_orb_scf)then
     call reorder_core_orb
     call initialize_mo_coef_begin_iteration
    endif

    TOUCH mo_coef_kpts

!   Calculate error vectors

    max_error_DIIS = maxval(cdabs(FPS_SPF_Matrix_MO_kpts))

!   SCF energy
!    call print_debug_scf_complex
    energy_SCF = scf_energy
    Delta_Energy_SCF = energy_SCF - energy_SCF_previous
    if ( (SCF_algorithm == 'DIIS').and.(Delta_Energy_SCF > 0.d0) ) then
      do kk=1,kpt_num
        Fock_matrix_AO_kpts(1:ao_num_per_kpt,1:ao_num_per_kpt,kk) = &
            Fock_matrix_DIIS (1:ao_num_per_kpt,1:ao_num_per_kpt,index_dim_DIIS,kk)
      enddo
      Fock_matrix_AO_alpha_kpts = Fock_matrix_AO_kpts*0.5d0
      Fock_matrix_AO_beta_kpts  = Fock_matrix_AO_kpts*0.5d0
      TOUCH Fock_matrix_AO_alpha_kpts Fock_matrix_AO_beta_kpts
    endif

    double precision :: level_shift_save
    level_shift_save = level_shift
    mo_coef_save(1:ao_num_per_kpt,1:mo_num_per_kpt,1:kpt_num) = mo_coef_kpts(1:ao_num_per_kpt,1:mo_num_per_kpt,1:kpt_num)
    do while (Delta_energy_SCF > 0.d0)
      mo_coef_kpts(1:ao_num_per_kpt,1:mo_num_per_kpt,1:kpt_num) = mo_coef_save
      if (level_shift <= .1d0) then
        level_shift = 1.d0
      else
        level_shift = level_shift * 3.0d0
      endif
      TOUCH mo_coef_kpts level_shift
      mo_coef_kpts(1:ao_num_per_kpt,1:mo_num_per_kpt,1:kpt_num) = &
          eigenvectors_fock_matrix_mo_kpts(1:ao_num_per_kpt,1:mo_num_per_kpt,1:kpt_num)
      if(frozen_orb_scf)then
        call reorder_core_orb
        call initialize_mo_coef_begin_iteration
      endif
      TOUCH mo_coef_kpts
      Delta_Energy_SCF = SCF_energy - energy_SCF_previous
      energy_SCF = SCF_energy
      if (level_shift-level_shift_save > 40.d0) then
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
      iteration_SCF, energy_scf, Delta_energy_SCF, max_error_DIIS, level_shift, dim_DIIS

    if (Delta_energy_SCF < 0.d0) then
      call save_mos
    endif
    if (qp_stop()) exit

  enddo

 if (iteration_SCF < n_it_SCF_max) then
   mo_label = "Canonical"
 endif
!
! End of Main SCF loop
!

  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,*)

  if(.not.frozen_orb_scf)then
   call mo_as_eigvectors_of_mo_matrix_kpts(Fock_matrix_mo_kpts,size(Fock_matrix_mo_kpts,1),size(Fock_matrix_mo_kpts,2),size(Fock_matrix_mo_kpts,3),mo_label,1,.true.)
   call save_mos
  endif

  call write_double(6, Energy_SCF, 'SCF energy')

  call write_time(6)

end

subroutine extrapolate_Fock_matrix_kpts(       &
           error_matrix_DIIS,Fock_matrix_DIIS,    &
           Fock_matrix_AO_,size_Fock_matrix_AO,   &
           iteration_SCF,dim_DIIS                 &
           )

BEGIN_DOC
! Compute the extrapolated Fock matrix using the DIIS procedure
END_DOC

  implicit none

  complex*16,intent(in)   :: Fock_matrix_DIIS(ao_num_per_kpt,ao_num_per_kpt,*),error_matrix_DIIS(ao_num_per_kpt,ao_num_per_kpt,*)
  integer,intent(in)            :: iteration_SCF, size_Fock_matrix_AO
  complex*16,intent(inout):: Fock_matrix_AO_(size_Fock_matrix_AO,ao_num_per_kpt)
  integer,intent(inout)         :: dim_DIIS

  double precision,allocatable  :: B_matrix_DIIS(:,:),X_vector_DIIS(:)
  double precision,allocatable  :: C_vector_DIIS(:)
  double precision              :: accum_im, thr_im
  complex*16,allocatable  :: scratch(:,:)
  integer                       :: i,j,k,i_DIIS,j_DIIS
  thr_im = 1.0d-10
  allocate(                               &
    B_matrix_DIIS(dim_DIIS+1,dim_DIIS+1), &
    X_vector_DIIS(dim_DIIS+1),            &
    C_vector_DIIS(dim_DIIS+1),            &
    scratch(ao_num,ao_num)                &
  )

! Compute the matrices B and X
  do j=1,dim_DIIS
    do i=1,dim_DIIS

      j_DIIS = mod(iteration_SCF-j,max_dim_DIIS)+1
      i_DIIS = mod(iteration_SCF-i,max_dim_DIIS)+1

! Compute product of two errors vectors

      call zgemm('N','N',ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt,                      &
           (1.d0,0.d0),                                             &
           error_matrix_DIIS(1,1,i_DIIS),size(error_matrix_DIIS,1), &
           error_matrix_DIIS(1,1,j_DIIS),size(error_matrix_DIIS,1), &
           (0.d0,0.d0),                                             &
           scratch,size(scratch,1))

! Compute Trace

      B_matrix_DIIS(i,j) = 0.d0
      accum_im = 0.d0
      do k=1,ao_num_per_kpt
        B_matrix_DIIS(i,j) = B_matrix_DIIS(i,j) + dble(scratch(k,k))
        accum_im = accum_im + dimag(scratch(k,k))
      enddo
      if (dabs(accum_im) .gt. thr_im) then
        !stop 'problem with imaginary parts in DIIS B_matrix?'
        print*, 'problem with imaginary parts in DIIS B_matrix?',accum_im
      endif
    enddo
  enddo
  deallocate(scratch)
! Pad B matrix and build the X matrix

  do i=1,dim_DIIS
    B_matrix_DIIS(i,dim_DIIS+1) = -1.d0
    B_matrix_DIIS(dim_DIIS+1,i) = -1.d0
    C_vector_DIIS(i) = 0.d0
  enddo
  B_matrix_DIIS(dim_DIIS+1,dim_DIIS+1) = 0.d0
  C_vector_DIIS(dim_DIIS+1) = -1.d0

! Solve the linear system C = B.X

  integer              :: info
  integer,allocatable  :: ipiv(:)

  allocate(          &
    ipiv(dim_DIIS+1) &
  )

  double precision, allocatable :: AF(:,:),scratch_d1(:)
  allocate (AF(dim_DIIS+1,dim_DIIS+1),scratch_d1(1))
  double precision :: rcond, ferr, berr
  integer :: iwork(dim_DIIS+1), lwork

  call dsysvx('N','U',dim_DIIS+1,1,      &
    B_matrix_DIIS,size(B_matrix_DIIS,1), &
    AF, size(AF,1),                      &
    ipiv,                                &
    C_vector_DIIS,size(C_vector_DIIS,1), &
    X_vector_DIIS,size(X_vector_DIIS,1), &
    rcond,                               &
    ferr,                                &
    berr,                                &
    scratch_d1,-1,                       &
    iwork,                               &
    info                                 &
  )
  lwork = int(scratch_d1(1))
  deallocate(scratch_d1)
  allocate(scratch_d1(lwork))

  call dsysvx('N','U',dim_DIIS+1,1,      &
    B_matrix_DIIS,size(B_matrix_DIIS,1), &
    AF, size(AF,1),                      &
    ipiv,                                &
    C_vector_DIIS,size(C_vector_DIIS,1), &
    X_vector_DIIS,size(X_vector_DIIS,1), &
    rcond,                               &
    ferr,                                &
    berr,                                &
    scratch_d1,size(scratch_d1),         &
    iwork,                               &
    info                                 &
  )
  deallocate(scratch_d1,ipiv)

 if(info < 0) then
   stop 'bug in DIIS'
 endif

 if (rcond > 1.d-12) then

  ! Compute extrapolated Fock matrix


      !$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(SHARED) if (ao_num_per_kpt > 200)
      do j=1,ao_num_per_kpt
        do i=1,ao_num_per_kpt
          Fock_matrix_AO_(i,j) = (0.d0,0.d0)
        enddo
        do k=1,dim_DIIS
          do i=1,ao_num_per_kpt
            Fock_matrix_AO_(i,j) = Fock_matrix_AO_(i,j) +            &
                X_vector_DIIS(k)*Fock_matrix_DIIS(i,j,dim_DIIS-k+1)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

  else
    dim_DIIS = 0
  endif

end

!============================================!
!                                            !
!                    kpts_real               !
!                                            !
!============================================!

subroutine Roothaan_Hall_SCF_kpts_real

BEGIN_DOC
! Roothaan-Hall algorithm for SCF Hartree-Fock calculation
END_DOC

  implicit none

  double precision               :: energy_SCF,energy_SCF_previous,Delta_energy_SCF
  double precision               :: max_error_DIIS,max_error_DIIS_alpha,max_error_DIIS_beta
  complex*16, allocatable  :: Fock_matrix_DIIS(:,:,:,:),error_matrix_DIIS(:,:,:,:)

  integer                        :: iteration_SCF,dim_DIIS,index_dim_DIIS

  integer                        :: i,j,k,kk
  logical, external              :: qp_stop
  complex*16, allocatable :: mo_coef_save(:,:,:)

  PROVIDE ao_md5 mo_occ_kpts level_shift

  allocate(mo_coef_save(ao_num_per_kpt,mo_num_per_kpt,kpt_num),                          &
      Fock_matrix_DIIS (ao_num_per_kpt,ao_num_per_kpt,max_dim_DIIS,kpt_num),                 &
      error_matrix_DIIS(ao_num_per_kpt,ao_num_per_kpt,max_dim_DIIS,kpt_num)                  &
      )
      !todo: add kpt_num dim to diis mats? (3 or 4)
  call write_time(6)

  print*,'Energy of the guess = ',scf_energy
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '  N ', 'Energy  ', 'Energy diff  ',  'DIIS error  ', 'Level shift   '
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
  !PROVIDE fps_spf_matrix_ao_complex fock_matrix_ao_complex 
  PROVIDE fps_spf_matrix_ao_kpts fock_matrix_ao_kpts 

  do while ( &
    ( (max_error_DIIS > threshold_DIIS_nonzero) .or. &
      (dabs(Delta_energy_SCF) > thresh_SCF) &
    ) .and. (iteration_SCF < n_it_SCF_max) )

! Increment cycle number

    iteration_SCF += 1
    if(frozen_orb_scf)then
     call initialize_mo_coef_begin_iteration
    endif

! Current size of the DIIS space

    dim_DIIS = min(dim_DIIS+1,max_dim_DIIS)

    if (scf_algorithm == 'DIIS') then

      do kk=1,kpt_num
        ! Store Fock and error matrices at each iteration
        do j=1,ao_num_per_kpt
          do i=1,ao_num_per_kpt
            index_dim_DIIS = mod(dim_DIIS-1,max_dim_DIIS)+1
            Fock_matrix_DIIS (i,j,index_dim_DIIS,kk) = fock_matrix_ao_kpts(i,j,kk)
            error_matrix_DIIS(i,j,index_dim_DIIS,kk) = fps_spf_matrix_ao_kpts(i,j,kk)
          enddo
        enddo

        ! Compute the extrapolated Fock matrix

        call extrapolate_fock_matrix_kpts(                            &
            error_matrix_DIIS(1,1,1,kk),Fock_matrix_DIIS(1,1,1,kk),        &
            Fock_matrix_AO_kpts(1,1,kk),size(Fock_matrix_AO_kpts,1),       &
            iteration_SCF,dim_DIIS                                       &
            )
      enddo
      Fock_matrix_AO_alpha_kpts = Fock_matrix_AO_kpts*0.5d0
      Fock_matrix_AO_beta_kpts  = Fock_matrix_AO_kpts*0.5d0
      TOUCH Fock_matrix_AO_alpha_kpts Fock_matrix_AO_beta_kpts

    endif

    mo_coef_kpts = eigenvectors_fock_matrix_mo_kpts
    if(frozen_orb_scf)then
     call reorder_core_orb
     call initialize_mo_coef_begin_iteration
    endif

    TOUCH mo_coef_kpts

!   Calculate error vectors

    max_error_DIIS = maxval(cdabs(FPS_SPF_Matrix_MO_kpts))

!   SCF energy
!    call print_debug_scf_complex
    energy_SCF = scf_energy
    Delta_Energy_SCF = energy_SCF - energy_SCF_previous
    if ( (SCF_algorithm == 'DIIS').and.(Delta_Energy_SCF > 0.d0) ) then
      do kk=1,kpt_num
        Fock_matrix_AO_kpts(1:ao_num_per_kpt,1:ao_num_per_kpt,kk) = &
            Fock_matrix_DIIS (1:ao_num_per_kpt,1:ao_num_per_kpt,index_dim_DIIS,kk)
      enddo
      Fock_matrix_AO_alpha_kpts = Fock_matrix_AO_kpts*0.5d0
      Fock_matrix_AO_beta_kpts  = Fock_matrix_AO_kpts*0.5d0
      TOUCH fock_matrix_ao_alpha_kpts Fock_matrix_AO_beta_kpts
    endif

    double precision :: level_shift_save
    level_shift_save = level_shift
    mo_coef_save(1:ao_num_per_kpt,1:mo_num_per_kpt,1:kpt_num) = mo_coef_kpts(1:ao_num_per_kpt,1:mo_num_per_kpt,1:kpt_num)
    do while (Delta_energy_SCF > 0.d0)
      mo_coef_kpts(1:ao_num_per_kpt,1:mo_num_per_kpt,1:kpt_num) = mo_coef_save
      if (level_shift <= .1d0) then
        level_shift = 1.d0
      else
        level_shift = level_shift * 3.0d0
      endif
      TOUCH mo_coef_kpts level_shift
      mo_coef_kpts(1:ao_num_per_kpt,1:mo_num_per_kpt,1:kpt_num) = &
          eigenvectors_fock_matrix_mo_kpts_real(1:ao_num_per_kpt,1:mo_num_per_kpt,1:kpt_num)
      if(frozen_orb_scf)then
        call reorder_core_orb
        call initialize_mo_coef_begin_iteration
      endif
      TOUCH mo_coef_kpts
      Delta_Energy_SCF = SCF_energy - energy_SCF_previous
      energy_SCF = SCF_energy
      if (level_shift-level_shift_save > 40.d0) then
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
      iteration_SCF, energy_scf, Delta_energy_SCF, max_error_DIIS, level_shift, dim_DIIS

    if (Delta_energy_SCF < 0.d0) then
      call save_mos
    endif
    if (qp_stop()) exit

  enddo

 if (iteration_SCF < n_it_SCF_max) then
   mo_label = "Canonical"
 endif
!
! End of Main SCF loop
!

  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,*)

  if(.not.frozen_orb_scf)then
   call mo_as_eigvectors_of_mo_matrix_kpts_real(fock_matrix_mo_kpts_real,size(Fock_matrix_mo_kpts_real,1),size(Fock_matrix_mo_kpts_real,2),size(Fock_matrix_mo_kpts_real,3),mo_label,1,.true.)
   call save_mos
  endif

  call write_double(6, Energy_SCF, 'SCF energy')

  call write_time(6)

end

