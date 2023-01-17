subroutine Roothaan_Hall_SCF_Simple

BEGIN_DOC
! Roothaan-Hall algorithm for SCF Hartree-Fock calculation
END_DOC

  implicit none

  integer                        :: iteration_SCF, dim_DIIS
  double precision               :: energy_SCF,energy_SCF_previous,Delta_energy_SCF
  double precision               :: max_error_DIIS

  integer                        :: i,j
  logical, external              :: qp_stop
  double precision, allocatable :: mo_coef_save(:,:)

  PROVIDE ao_md5 mo_occ level_shift

  allocate(mo_coef_save(ao_num,mo_num))


  dim_DIIS     = 0
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
  max_error_DIIS      = 1.d0

  do while ( &
    ( (max_error_DIIS > threshold_DIIS_nonzero) .or. &
      (dabs(Delta_energy_SCF) > thresh_SCF) &
    ) .and. (iteration_SCF < n_it_SCF_max) )

    iteration_SCF += 1
    if(frozen_orb_scf)then
     call initialize_mo_coef_begin_iteration
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

    !double precision :: level_shift_save
    !level_shift_save = level_shift
    !mo_coef_save(1:ao_num,1:mo_num) = mo_coef(1:ao_num,1:mo_num)
    !do while (Delta_energy_SCF > 0.d0)
    !  mo_coef(1:ao_num,1:mo_num) = mo_coef_save
    !  if (level_shift <= .1d0) then
    !    level_shift = 1.d0
    !  else
    !    level_shift = level_shift * 3.0d0
    !  endif
    !  TOUCH mo_coef level_shift
    !  mo_coef(1:ao_num,1:mo_num) = eigenvectors_Fock_matrix_MO(1:ao_num,1:mo_num)
    !  if(frozen_orb_scf)then
    !    call reorder_core_orb
    !    call initialize_mo_coef_begin_iteration
    !  endif
    !  TOUCH mo_coef
    !  Delta_energy_SCF = SCF_energy - energy_SCF_previous
    !  energy_SCF = SCF_energy
    !  if (level_shift-level_shift_save > 40.d0) then
    !    level_shift = level_shift_save * 4.d0
    !    SOFT_TOUCH level_shift
    !    exit
    !  endif
    !enddo
    !level_shift = level_shift * 0.5d0
    !SOFT_TOUCH level_shift

    energy_SCF_previous = energy_SCF

!   Print results at the end of each iteration

    write(6,'(I4, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, I3)')  &
      iteration_SCF, energy_SCF, Delta_energy_SCF, max_error_DIIS, level_shift, dim_DIIS

    if(Delta_energy_SCF < 0.d0) then
      call save_mos()
    endif
    if(qp_stop()) exit

  enddo

  if (iteration_SCF < n_it_SCF_max) then
    mo_label = 'Canonical'
  endif

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

