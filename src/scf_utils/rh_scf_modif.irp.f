subroutine Roothaan_Hall_SCF_MODIF

BEGIN_DOC
! Roothaan-Hall algorithm for SCF Hartree-Fock calculation
END_DOC

  implicit none

  double precision               :: energy_SCF,energy_SCF_previous,Delta_energy_SCF
  double precision               :: max_error_DIIS,max_error_DIIS_alpha,max_error_DIIS_beta
  double precision, allocatable  :: Fock_matrix_DIIS(:,:,:),error_matrix_DIIS(:,:,:)

  integer                        :: iteration_SCF,dim_DIIS,index_dim_DIIS

  integer                        :: i,j
  logical, external              :: qp_stop
  double precision, allocatable :: mo_coef_save(:,:)

  PROVIDE ao_md5 mo_occ level_shift

  allocate(mo_coef_save(ao_num,mo_num),                          &
      Fock_matrix_DIIS (ao_num,ao_num,max_dim_DIIS),                 &
      error_matrix_DIIS(ao_num,ao_num,max_dim_DIIS)                  &
      )

  Fock_matrix_DIIS = 0.d0
  error_matrix_DIIS = 0.d0
  mo_coef_save = 0.d0

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
  PROVIDE FPS_SPF_matrix_AO Fock_matrix_AO 

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

    if( (scf_algorithm == 'DIIS_MODIF') .and. (dabs(Delta_energy_SCF) > 1.d-6) )  then
    !if(scf_algorithm == 'DIIS_MODIF') then

      ! Store Fock and error matrices at each iteration
      index_dim_DIIS = mod(dim_DIIS-1,max_dim_DIIS)+1
      do j=1,ao_num
        do i=1,ao_num
          Fock_matrix_DIIS (i,j,index_dim_DIIS) = Fock_matrix_AO(i,j)
          error_matrix_DIIS(i,j,index_dim_DIIS) = FPS_SPF_matrix_AO(i,j)
        enddo
      enddo

      ! Compute the extrapolated Fock matrix

      call extrapolate_Fock_matrix(                                    &
          error_matrix_DIIS,Fock_matrix_DIIS,                          &
          Fock_matrix_AO,size(Fock_matrix_AO,1),                       &
          iteration_SCF,dim_DIIS                                       &
          )
      call ao_to_mo(Fock_matrix_AO, size(Fock_matrix_AO, 1), Fock_matrix_MO, size(Fock_matrix_MO, 1))
      do i = 1, mo_num
        Fock_matrix_diag_MO(i) = Fock_matrix_MO(i,i)
      enddo
      TOUCH Fock_matrix_MO Fock_matrix_diag_MO

      !Fock_matrix_AO_alpha = Fock_matrix_AO*0.5d0
      !Fock_matrix_AO_beta  = Fock_matrix_AO*0.5d0
      !TOUCH Fock_matrix_AO_alpha Fock_matrix_AO_beta
    endif

    MO_coef = eigenvectors_Fock_matrix_MO
    if(frozen_orb_scf)then
     call reorder_core_orb
     call initialize_mo_coef_begin_iteration
    endif

    TOUCH MO_coef

!   Calculate error vectors

    max_error_DIIS = maxval(Abs(FPS_SPF_Matrix_MO))

!   SCF energy

    energy_SCF = SCF_energy
    Delta_energy_SCF = energy_SCF - energy_SCF_previous
    if( (SCF_algorithm == 'DIIS_MODIF') .and. (Delta_energy_SCF > 0.d0) ) then
      Fock_matrix_AO(1:ao_num,1:ao_num) = Fock_matrix_DIIS(1:ao_num,1:ao_num,index_dim_DIIS)
      call ao_to_mo(Fock_matrix_AO, size(Fock_matrix_AO, 1), Fock_matrix_MO, size(Fock_matrix_MO, 1))
      do i = 1, mo_num
        Fock_matrix_diag_MO(i) = Fock_matrix_MO(i,i)
      enddo
      TOUCH Fock_matrix_MO Fock_matrix_diag_MO

      !Fock_matrix_AO_alpha = Fock_matrix_AO*0.5d0
      !Fock_matrix_AO_beta  = Fock_matrix_AO*0.5d0
      !TOUCH Fock_matrix_AO_alpha Fock_matrix_AO_beta
    endif

    double precision :: level_shift_save
    level_shift_save = level_shift
    mo_coef_save(1:ao_num,1:mo_num) = mo_coef(1:ao_num,1:mo_num)
    do while (Delta_energy_SCF > 0.d0)
      mo_coef(1:ao_num,1:mo_num) = mo_coef_save
      if (level_shift <= .1d0) then
        level_shift = 1.d0
      else
        level_shift = level_shift * 3.0d0
      endif
      TOUCH mo_coef level_shift
      mo_coef(1:ao_num,1:mo_num) = eigenvectors_Fock_matrix_MO(1:ao_num,1:mo_num)
      if(frozen_orb_scf)then
        call reorder_core_orb
        call initialize_mo_coef_begin_iteration
      endif
      TOUCH mo_coef
      Delta_energy_SCF = SCF_energy - energy_SCF_previous
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
      iteration_SCF, energy_SCF, Delta_energy_SCF, max_error_DIIS, level_shift, dim_DIIS

    if (Delta_energy_SCF < 0.d0) then
      call save_mos
    endif
    if (qp_stop()) exit

  enddo

 if (iteration_SCF < n_it_SCF_max) then
   mo_label = 'Canonical'
 endif
!
! End of Main SCF loop
!

  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,*)

  if(.not.frozen_orb_scf)then
   call mo_as_eigvectors_of_mo_matrix(Fock_matrix_mo,size(Fock_matrix_mo,1), &
      size(Fock_matrix_mo,2),mo_label,1,.true.)
   call restore_symmetry(ao_num, mo_num, mo_coef, size(mo_coef,1), 1.d-10)
   call orthonormalize_mos
   call save_mos
  endif

  call write_double(6, energy_SCF, 'SCF energy')

  call write_time(6)

end

