! ---

program tc_scf

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'starting ...'

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
!  my_n_pt_r_grid = 10 ! small grid for quick debug
!  my_n_pt_a_grid = 26 ! small grid for quick debug
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  PROVIDE mu_erf 
  print *, ' mu = ', mu_erf
  PROVIDE j1b_type
  print *, ' j1b_type = ', j1b_type
  print *, j1b_pen

  !call create_guess()
  !call orthonormalize_mos()

  PROVIDE tcscf_algorithm
  if(tcscf_algorithm == 'DIIS') then
    call rh_tcscf_diis()
  elseif(tcscf_algorithm == 'Simple') then
    call rh_tcscf_simple()
  else
    print *, ' not implemented yet', tcscf_algorithm
    stop
  endif

  call minimize_tc_orb_angles()
  call print_energy_and_mos()

end

! ---

subroutine create_guess()

  implicit none
  logical :: exists

  PROVIDE ezfio_filename
  !call ezfio_has_mo_basis_mo_coef(exists)
  exists = .false.

  if(.not.exists) then
    mo_label = 'Guess'
    if(mo_guess_type == "HCore") then
      mo_coef = ao_ortho_lowdin_coef
      call restore_symmetry(ao_num, mo_num, mo_coef, size(mo_coef, 1), 1.d-10)
      TOUCH mo_coef
      call mo_as_eigvectors_of_mo_matrix(mo_one_e_integrals, size(mo_one_e_integrals, 1), size(mo_one_e_integrals, 2), mo_label, 1, .false.)
      call restore_symmetry(ao_num, mo_num, mo_coef, size(mo_coef, 1), 1.d-10)
      SOFT_TOUCH mo_coef
    elseif (mo_guess_type == "Huckel") then
      call huckel_guess
    else
      print *,  'Unrecognized MO guess type : '//mo_guess_type
      stop 1
    endif
    SOFT_TOUCH mo_label
  endif

end subroutine create_guess

! ---
