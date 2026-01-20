.. _module_ezfio_files: 
 
.. program:: ezfio_files 
 
.. default-role:: option 
 
===========
ezfio_files
===========

This modules essentially contains the name of the |EZFIO| directory in the
:c:data:`ezfio_filename` variable. This is read as the first argument of the
command-line, or as the :envvar:`QP_INPUT` environment variable.

 
 
 
Providers 
--------- 
 
.. c:var:: ezfio_filename


    File : :file:`ezfio_files/ezfio.irp.f`

    .. code:: fortran

        character*(1024)	:: ezfio_filename	


    Name of EZFIO file. It is obtained from the QPACKAGE_INPUT environment
    variable if it is set, or as the 1st argument of the command line.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`file_lock`
       * :c:data:`mpi_initialized`
       * :c:data:`output_wall_time_0`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`absolute_eig`
       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`act_mos_opt`
       * :c:data:`adaptive_pt2_max`
       * :c:data:`ao_cartesian`
       * :c:data:`ao_cholesky_threshold`
       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_expo_im`
       * :c:data:`ao_expo_phase`
       * :c:data:`ao_expo_pw`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_md5`
       * :c:data:`ao_normalized`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_threshold`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`avoid_saddle`
       * :c:data:`basis`
       * :c:data:`basis_nucleus_index`
       * :c:data:`calc_dipole_moment`
       * :c:data:`calc_energy_components`
       * :c:data:`calc_osc_str`
       * :c:data:`calc_tr_dipole_moment`
       * :c:data:`correlation_energy_ratio_max`
       * :c:data:`correlation_functional`
       * :c:data:`criterion_casscf`
       * :c:data:`csf_based`
       * :c:data:`damping_for_rs_dft`
       * :c:data:`data_energy_proj`
       * :c:data:`data_energy_var`
       * :c:data:`data_one_e_dm_alpha_ao`
       * :c:data:`data_one_e_dm_alpha_mo`
       * :c:data:`data_one_e_dm_beta_ao`
       * :c:data:`data_one_e_dm_beta_mo`
       * :c:data:`davidson_sze_max`
       * :c:data:`density_for_dft`
       * :c:data:`diag_hess_cas`
       * :c:data:`disk_based_davidson`
       * :c:data:`distributed_davidson`
       * :c:data:`do_ao_cholesky`
       * :c:data:`do_mom`
       * :c:data:`do_ormas`
       * :c:data:`do_pseudo`
       * :c:data:`do_pt2`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`exchange_functional`
       * :c:data:`excitation_alpha_max`
       * :c:data:`excitation_beta_max`
       * :c:data:`excitation_max`
       * :c:data:`excitation_ref`
       * :c:data:`ezfio_work_dir`
       * :c:data:`fast_2rdm`
       * :c:data:`frozen_orb_scf`
       * :c:data:`generators_bitmask`
       * :c:data:`h0_type`
       * :c:data:`hess_cv_cv`
       * :c:data:`hf_exchange`
       * :c:data:`io_ao_cholesky`
       * :c:data:`io_ao_integrals_kinetic`
       * :c:data:`io_ao_integrals_n_e`
       * :c:data:`io_ao_integrals_overlap`
       * :c:data:`io_ao_integrals_pseudo`
       * :c:data:`io_ao_one_e_integrals`
       * :c:data:`io_ao_two_e_integrals`
       * :c:data:`io_ao_two_e_integrals_erf`
       * :c:data:`io_mo_cholesky`
       * :c:data:`io_mo_integrals_kinetic`
       * :c:data:`io_mo_integrals_n_e`
       * :c:data:`io_mo_integrals_pseudo`
       * :c:data:`io_mo_one_e_integrals`
       * :c:data:`io_mo_two_e_integrals`
       * :c:data:`io_mo_two_e_integrals_erf`
       * :c:data:`io_nuclear_repulsion`
       * :c:data:`io_two_body_rdm_aa`
       * :c:data:`io_two_body_rdm_ab`
       * :c:data:`io_two_body_rdm_bb`
       * :c:data:`io_two_body_rdm_spin_trace`
       * :c:data:`is_periodic`
       * :c:data:`json_filename`
       * :c:data:`level_shift`
       * :c:data:`level_shift_casscf`
       * :c:data:`lin_dep_cutoff`
       * :c:data:`max_dim_diis`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_aux`
       * :c:data:`mo_coef_imag`
       * :c:data:`mo_guess_type`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_label`
       * :c:data:`mo_num`
       * :c:data:`mo_occ`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mu_dft_type`
       * :c:data:`mu_erf`
       * :c:data:`n_big_act_orb`
       * :c:data:`n_det`
       * :c:data:`n_det_max`
       * :c:data:`n_det_max_full`
       * :c:data:`n_det_max_opt`
       * :c:data:`n_det_print_wf`
       * :c:data:`n_it_scf_max`
       * :c:data:`n_pts_charge`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nb_it_max_lambda`
       * :c:data:`nb_it_max_pre_search`
       * :c:data:`no_core_density`
       * :c:data:`no_oa_or_av_opt`
       * :c:data:`normalize_dm`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_charge_remove`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`
       * :c:data:`nucleus_shell_num`
       * :c:data:`only_expected_s2`
       * :c:data:`optimization_max_nb_iter`
       * :c:data:`optimization_method`
       * :c:data:`ormas_max_e`
       * :c:data:`ormas_min_e`
       * :c:data:`ormas_mstart`
       * :c:data:`ormas_n_space`
       * :c:data:`point_charges`
       * :c:data:`prim_coef`
       * :c:data:`prim_expo`
       * :c:data:`prim_normalization_factor`
       * :c:data:`prim_num`
       * :c:data:`primitives_normalized`
       * :c:data:`print_all_transitions`
       * :c:data:`pruning`
       * :c:data:`pseudo_dz_k`
       * :c:data:`pseudo_dz_kl`
       * :c:data:`pseudo_grid_rmax`
       * :c:data:`pseudo_grid_size`
       * :c:data:`pseudo_klocmax`
       * :c:data:`pseudo_kmax`
       * :c:data:`pseudo_lmax`
       * :c:data:`pseudo_n_k`
       * :c:data:`pseudo_n_kl`
       * :c:data:`pseudo_sym`
       * :c:data:`pseudo_v_k`
       * :c:data:`pseudo_v_kl`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`pt2_max`
       * :c:data:`pt2_min_casscf`
       * :c:data:`pt2_min_parallel_tasks`
       * :c:data:`pt2_relative_error`
       * :c:data:`pts_charge_coord`
       * :c:data:`pts_charge_z`
       * :c:data:`qp_stop_filename`
       * :c:data:`read_wf`
       * :c:data:`restore_symm`
       * :c:data:`s2_eig`
       * :c:data:`save_threshold`
       * :c:data:`save_wf_after_selection`
       * :c:data:`scf_algorithm`
       * :c:data:`selection_factor`
       * :c:data:`seniority_max`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_index`
       * :c:data:`shell_normalization_factor`
       * :c:data:`shell_num`
       * :c:data:`shell_prim_num`
       * :c:data:`small_active_space`
       * :c:data:`state_following`
       * :c:data:`state_following_casscf`
       * :c:data:`target_energy`
       * :c:data:`thresh_casscf`
       * :c:data:`thresh_cc`
       * :c:data:`thresh_delta`
       * :c:data:`thresh_eig`
       * :c:data:`thresh_model`
       * :c:data:`thresh_model_2`
       * :c:data:`thresh_opt_max_elem_grad`
       * :c:data:`thresh_rho`
       * :c:data:`thresh_rho_2`
       * :c:data:`thresh_scf`
       * :c:data:`thresh_sym`
       * :c:data:`thresh_wtg`
       * :c:data:`thresh_wtg2`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_davidson_from_pt2`
       * :c:data:`threshold_diis`
       * :c:data:`threshold_generators`
       * :c:data:`threshold_nonsym_davidson`
       * :c:data:`twice_hierarchy_max`
       * :c:data:`typ`
       * :c:data:`use_cgtos`
       * :c:data:`use_only_lr`
       * :c:data:`variance_max`
       * :c:data:`version_avoid_saddle`
       * :c:data:`version_lambda_search`
       * :c:data:`weight_one_e_dm`
       * :c:data:`weight_selection`
       * :c:data:`without_diagonal`

 
.. c:var:: ezfio_work_dir


    File : :file:`ezfio_files/ezfio.irp.f`

    .. code:: fortran

        character*(1024)	:: ezfio_work_dir	


    EZFIO/work/

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`ezfio_work_dir_pid`

 
.. c:var:: ezfio_work_dir_pid


    File : :file:`ezfio_files/ezfio.irp.f`

    .. code:: fortran

        character*(1024)	:: ezfio_work_dir_pid	


    EZFIO/work/pid_

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_work_dir`


 
.. c:var:: file_lock


    File : :file:`ezfio_files/lock.irp.f`

    .. code:: fortran

        integer(omp_lock_kind)	:: file_lock	


    OpenMP Lock for I/O

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`json_filename`
       * :c:data:`json_unit`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`
       * :c:data:`qp_max_mem`

 
.. c:var:: nthreads_pt2


    File : :file:`ezfio_files/environment.irp.f`

    .. code:: fortran

        integer	:: nthreads_pt2	


    Number of threads for Davidson

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`file_lock`
       * :c:data:`mpi_master`
       * :c:data:`nproc`


 
.. c:var:: output_cpu_time_0


    File : :file:`ezfio_files/output.irp.f`

    .. code:: fortran

        double precision	:: output_wall_time_0	
        double precision	:: output_cpu_time_0	


    Initial CPU and wall times when printing in the output files

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy`
       * :c:data:`ezfio_filename`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_pts_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nuclear_repulsion`
       * :c:data:`prim_normalization_factor`
       * :c:data:`shell_normalization_factor`

 
.. c:var:: output_wall_time_0


    File : :file:`ezfio_files/output.irp.f`

    .. code:: fortran

        double precision	:: output_wall_time_0	
        double precision	:: output_cpu_time_0	


    Initial CPU and wall times when printing in the output files

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy`
       * :c:data:`ezfio_filename`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_pts_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nuclear_repulsion`
       * :c:data:`prim_normalization_factor`
       * :c:data:`shell_normalization_factor`

 
.. c:var:: qp_kill_filename


    File : :file:`ezfio_files/qp_stop.irp.f`

    .. code:: fortran

        character*(256)	:: qp_stop_filename	
        character*(256)	:: qp_kill_filename	
        integer	:: qp_stop_variable	


    Name of the file to check for qp stop

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`


 
.. c:var:: qp_stop_filename


    File : :file:`ezfio_files/qp_stop.irp.f`

    .. code:: fortran

        character*(256)	:: qp_stop_filename	
        character*(256)	:: qp_kill_filename	
        integer	:: qp_stop_variable	


    Name of the file to check for qp stop

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`


 
.. c:var:: qp_stop_variable


    File : :file:`ezfio_files/qp_stop.irp.f`

    .. code:: fortran

        character*(256)	:: qp_stop_filename	
        character*(256)	:: qp_kill_filename	
        integer	:: qp_stop_variable	


    Name of the file to check for qp stop

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: getunitandopen:


    File : :file:`ezfio_files/get_unit_and_open.irp.f`

    .. code:: fortran

        integer function getUnitAndOpen(f,mode)


    :f:
       file name
    
    :mode:
       'R' : READ, UNFORMATTED
       'W' : WRITE, UNFORMATTED
       'r' : READ, FORMATTED
       'w' : WRITE, FORMATTED
       'a' : APPEND, FORMATTED
       'x' : READ/WRITE, FORMATTED
    

 
.. c:function:: lock_io:


    File : :file:`ezfio_files/lock.irp.f`

    .. code:: fortran

        subroutine lock_io()


    Needs to be called because before doing I/O because internal read and write
    are not thread safe.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`file_lock`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_slave_work`
       * :c:func:`format_w_error`
       * :c:func:`json_close`
       * :c:data:`json_filename`
       * :c:data:`json_unit`
       * :c:func:`load_mo_integrals`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`
       * :c:data:`qp_max_mem`
       * :c:func:`read_array_two_rdm`
       * :c:func:`read_array_two_trans_rdm`
       * :c:func:`resident_memory`
       * :c:func:`roothaan_hall_scf`
       * :c:func:`total_memory`
       * :c:func:`write_array_two_rdm`
       * :c:func:`write_array_two_trans_rdm`
       * :c:func:`write_cipsi_json`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`

 
.. c:function:: qp_stop:


    File : :file:`ezfio_files/qp_stop.irp.f`

    .. code:: fortran

        logical function qp_stop()


    Checks if the qp_stop command was invoked for the clean termination of the program

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_stop_filename`

 
.. c:function:: unlock_io:


    File : :file:`ezfio_files/lock.irp.f`

    .. code:: fortran

        subroutine unlock_io()


    Needs to be called because afterdoing I/O because internal read and write
    are not thread safe.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`file_lock`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_slave_work`
       * :c:func:`format_w_error`
       * :c:func:`json_close`
       * :c:data:`json_filename`
       * :c:data:`json_unit`
       * :c:func:`load_mo_integrals`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`
       * :c:data:`qp_max_mem`
       * :c:func:`read_array_two_rdm`
       * :c:func:`read_array_two_trans_rdm`
       * :c:func:`resident_memory`
       * :c:func:`roothaan_hall_scf`
       * :c:func:`total_memory`
       * :c:func:`write_array_two_rdm`
       * :c:func:`write_array_two_trans_rdm`
       * :c:func:`write_cipsi_json`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_unset_lock`

 
.. c:function:: write_bool:


    File : :file:`ezfio_files/output.irp.f`

    .. code:: fortran

        subroutine write_bool(iunit,value,label)


    Write an logical value in output

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

 
.. c:function:: write_double:


    File : :file:`ezfio_files/output.irp.f`

    .. code:: fortran

        subroutine write_double(iunit,value,label)


    Write a double precision value in output

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ci_energy`
       * :c:func:`damping_scf`
       * :c:func:`dav_double_dressed`
       * :c:func:`davidson_diag_csf_hjj`
       * :c:func:`davidson_diag_hjj`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`davidson_general`
       * :c:func:`davidson_general_diag_dressed_ext_rout_nonsym_b1space`
       * :c:func:`davidson_general_ext_rout`
       * :c:func:`davidson_general_ext_rout_diag_dressed`
       * :c:func:`davidson_general_ext_rout_dressed`
       * :c:func:`davidson_general_ext_rout_nonsym_b1space`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_coef_max`
       * :c:data:`pt2_e0_denominator`
       * :c:func:`roothaan_hall_scf`
       * :c:func:`run_cipsi`
       * :c:func:`run_slave_main`
       * :c:func:`run_stochastic_cipsi`
       * :c:func:`zmq_pt2`
       * :c:func:`zmq_selection`

 
.. c:function:: write_int:


    File : :file:`ezfio_files/output.irp.f`

    .. code:: fortran

        subroutine write_int(iunit,value,label)


    Write an integer value in output

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`dav_double_dressed`
       * :c:func:`davidson_diag_csf_hjj`
       * :c:func:`davidson_diag_hjj`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`davidson_general`
       * :c:func:`davidson_general_diag_dressed_ext_rout_nonsym_b1space`
       * :c:func:`davidson_general_ext_rout`
       * :c:func:`davidson_general_ext_rout_diag_dressed`
       * :c:func:`davidson_general_ext_rout_dressed`
       * :c:func:`davidson_general_ext_rout_nonsym_b1space`
       * :c:func:`make_s2_eigenfunction`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_del_orb`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_inact_orb`
       * :c:data:`n_int`
       * :c:data:`n_virt_orb`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`
       * :c:data:`psi_cas`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_size`
       * :c:data:`pt2_f`
       * :c:data:`pt2_n_teeth`
       * :c:data:`qp_max_mem`
       * :c:func:`remove_small_contributions`
       * :c:func:`save_wavefunction_general`
       * :c:func:`save_wavefunction_general_unormalized`
       * :c:func:`save_wavefunction_specified`
       * :c:func:`zmq_pt2`

 
.. c:function:: write_time:


    File : :file:`ezfio_files/output.irp.f`

    .. code:: fortran

        subroutine write_time(iunit)


    Write a time stamp in the output for chronological reconstruction

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`output_wall_time_0`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy`
       * :c:func:`damping_scf`
       * :c:func:`dav_double_dressed`
       * :c:func:`davidson_diag_csf_hjj`
       * :c:func:`davidson_diag_hjj`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`davidson_general`
       * :c:func:`davidson_general_diag_dressed_ext_rout_nonsym_b1space`
       * :c:func:`davidson_general_ext_rout`
       * :c:func:`davidson_general_ext_rout_diag_dressed`
       * :c:func:`davidson_general_ext_rout_dressed`
       * :c:func:`davidson_general_ext_rout_nonsym_b1space`
       * :c:func:`make_s2_eigenfunction`
       * :c:func:`mo_as_eigvectors_of_mo_matrix`
       * :c:func:`mo_as_svd_vectors_of_mo_matrix`
       * :c:func:`mo_as_svd_vectors_of_mo_matrix_eig`
       * :c:func:`mo_coef_new_as_svd_vectors_of_mo_matrix_eig`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_pts_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nuclear_repulsion`
       * :c:data:`prim_normalization_factor`
       * :c:func:`roothaan_hall_scf`
       * :c:data:`shell_normalization_factor`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`cpu_time`
       * :c:func:`print_memory_usage`
       * :c:func:`wall_time`

