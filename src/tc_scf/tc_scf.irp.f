! ---

program tc_scf

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  write(json_unit,json_array_open_fmt) 'tc-scf'

  print *, ' starting ...'

  my_grid_becke  = .True.

  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a

  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  PROVIDE mu_erf 
  print *, ' mu = ', mu_erf
  PROVIDE j1b_type
  print *, ' j1b_type = ', j1b_type
  print *, j1b_pen

  !call create_guess()
  !call orthonormalize_mos()

  PROVIDE tcscf_algorithm
  PROVIDE var_tc

  if(var_tc) then

    print *, ' VAR-TC'

    if(tcscf_algorithm == 'DIIS') then
      print*, ' NOT implemented yet'
    elseif(tcscf_algorithm == 'Simple') then
      call rh_vartcscf_simple()
    else
      print *, ' not implemented yet', tcscf_algorithm
      stop
    endif

  else

    if(tcscf_algorithm == 'DIIS') then
      call rh_tcscf_diis()
    elseif(tcscf_algorithm == 'Simple') then
      call rh_tcscf_simple()
    else
      print *, ' not implemented yet', tcscf_algorithm
      stop
    endif

    call minimize_tc_orb_angles()

  endif

  write(json_unit,json_array_close_fmtx)
  call json_close

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

