
! ---

program fci_tc_bi

  BEGIN_DOC
  ! Selected Full Configuration Interaction with stochastic selection
  ! and PT2.
  !
  ! This program performs a |CIPSI|-like selected |CI| using a
  ! stochastic scheme for both the selection of the important Slater
  ! determinants and the computation of the |PT2| correction. This
  ! |CIPSI|-like algorithm will be performed for the lowest states of
  ! the variational space (see :option:`determinants n_states`). The
  ! |FCI| program will stop when reaching at least one the two following
  ! conditions:
  !
  ! * number of Slater determinants > :option:`determinants n_det_max`
  ! * abs(|PT2|) less than :option:`perturbation pt2_max`
  !
  ! The following other options can be of interest:
  !
  ! :option:`determinants read_wf`
  !   When set to |false|, the program starts with a ROHF-like Slater
  !   determinant as a guess wave function. When set to |true|, the
  !   program starts with the wave function(s) stored in the |EZFIO|
  !   directory as guess wave function(s).
  !
  ! :option:`determinants s2_eig`
  !   When set to |true|, the selection will systematically add all the
  !   necessary Slater determinants in order to have a pure spin wave
  !   function with an |S^2| value corresponding to
  !   :option:`determinants expected_s2`.
  !
  ! For excited states calculations, it is recommended to start with
  ! :ref:`cis` or :ref:`cisd` guess wave functions, eventually in
  ! a restricted set of |MOs|, and to set :option:`determinants s2_eig`
  ! to |true|.
  !
  END_DOC

  implicit none

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 

  pruning = -1.d0
  touch pruning

!  pt2_relative_error = 0.01d0
!  touch pt2_relative_error

  call run_cipsi_tc()

end

! ---

subroutine run_cipsi_tc()

  implicit none

  if (.not. is_zmq_slave) then

    PROVIDE psi_det psi_coef mo_bi_ortho_tc_two_e mo_bi_ortho_tc_one_e

    if((elec_alpha_num+elec_beta_num) .ge. 3) then
      if(three_body_h_tc) then
        call provide_all_three_ints_bi_ortho()
      endif
    endif

    FREE int2_grad1_u12_ao int2_grad1_u12_ao_t int2_grad1_u12_ao_transp
    FREE int2_grad1_u12_bimo_transp 

    write(json_unit,json_array_open_fmt) 'fci_tc'

    if(do_pt2) then
      call run_stochastic_cipsi()
    else
      call run_cipsi()
    endif

    write(json_unit,json_dict_uopen_fmt)
    write(json_unit,json_dict_close_fmtx)
    write(json_unit,json_array_close_fmtx)
    call json_close

  else

    PROVIDE mo_bi_ortho_tc_one_e mo_bi_ortho_tc_two_e pt2_min_parallel_tasks

    if((elec_alpha_num+elec_beta_num) .ge. 3) then
      if(three_body_h_tc) then
        call provide_all_three_ints_bi_ortho()
      endif
    endif

    FREE int2_grad1_u12_ao int2_grad1_u12_ao_t int2_grad1_u12_ao_transp
    FREE int2_grad1_u12_bimo_transp 

    call run_slave_cipsi

  endif

end
