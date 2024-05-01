! ---

program tc_scf

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none
  integer :: i
  logical :: good_angles

  PROVIDE j1e_type
  PROVIDE j2e_type
  PROVIDE tcscf_algorithm

  print *, ' TC-SCF with:'
  print *, ' j1e_type = ', j1e_type
  print *, ' j2e_type = ', j2e_type

  write(json_unit,json_array_open_fmt) 'tc-scf'

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call write_int(6, my_n_pt_r_grid, 'radial  external grid over')
  call write_int(6, my_n_pt_a_grid, 'angular external grid over')


  if(tc_integ_type .eq. "numeric") then
    my_extra_grid_becke  = .True.
    PROVIDE tc_grid2_a tc_grid2_r
    my_n_pt_r_extra_grid = tc_grid2_r
    my_n_pt_a_extra_grid = tc_grid2_a
    touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid

    call write_int(6, my_n_pt_r_extra_grid, 'radial  internal grid over')
    call write_int(6, my_n_pt_a_extra_grid, 'angular internal grid over')
  endif

  !call create_guess()
  !call orthonormalize_mos()

  if(tcscf_algorithm == 'DIIS') then
    call rh_tcscf_diis()
  elseif(tcscf_algorithm == 'Simple') then
    call rh_tcscf_simple()
  else
    print *, ' not implemented yet', tcscf_algorithm
    stop
  endif

  PROVIDE Fock_matrix_tc_diag_mo_tot
  print*, ' Eigenvalues:' 
  do i = 1, mo_num
    print*, i, Fock_matrix_tc_diag_mo_tot(i)
  enddo

  ! TODO 
  ! rotate angles in separate code only if necessary
  if(minimize_lr_angles)then
   call minimize_tc_orb_angles()
  endif
  call print_energy_and_mos(good_angles)


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

