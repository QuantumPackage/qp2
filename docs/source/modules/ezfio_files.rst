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

        character*(128)	:: ezfio_filename	


    Name of EZFIO file. It is obtained from the QPACKAGE_INPUT environment
    variable if it is set, or as the 1st argument of the command line.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_initialized`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cartesian`
       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_md5`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`cas_bitmask`
       * :c:data:`correlation_energy_ratio_max`
       * :c:data:`data_energy_proj`
       * :c:data:`data_energy_var`
       * :c:data:`data_one_e_dm_alpha_ao`
       * :c:data:`data_one_e_dm_alpha_mo`
       * :c:data:`data_one_e_dm_beta_ao`
       * :c:data:`data_one_e_dm_beta_mo`
       * :c:data:`davidson_sze_max`
       * :c:data:`disk_access_nuclear_repulsion`
       * :c:data:`disk_based_davidson`
       * :c:data:`distributed_davidson`
       * :c:data:`do_direct_integrals`
       * :c:data:`do_pseudo`
       * :c:data:`do_pt2`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`energy_iterations`
       * :c:data:`ezfio_work_dir`
       * :c:data:`frozen_orb_scf`
       * :c:data:`generators_bitmask`
       * :c:data:`generators_bitmask_restart`
       * :c:data:`h0_type`
       * :c:data:`io_ao_integrals_e_n`
       * :c:data:`io_ao_integrals_kinetic`
       * :c:data:`io_ao_integrals_overlap`
       * :c:data:`io_ao_integrals_pseudo`
       * :c:data:`io_ao_one_e_integrals`
       * :c:data:`io_ao_two_e_integrals`
       * :c:data:`io_ao_two_e_integrals_erf`
       * :c:data:`io_mo_integrals_e_n`
       * :c:data:`io_mo_integrals_kinetic`
       * :c:data:`io_mo_integrals_pseudo`
       * :c:data:`io_mo_one_e_integrals`
       * :c:data:`io_mo_two_e_integrals`
       * :c:data:`io_mo_two_e_integrals_erf`
       * :c:data:`level_shift`
       * :c:data:`max_dim_diis`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_guess_type`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_label`
       * :c:data:`mo_num`
       * :c:data:`mo_occ`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mu_erf`
       * :c:data:`n_cas_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_det_iterations`
       * :c:data:`n_det_max`
       * :c:data:`n_det_max_full`
       * :c:data:`n_det_print_wf`
       * :c:data:`n_generators_bitmask`
       * :c:data:`n_generators_bitmask_restart`
       * :c:data:`n_it_scf_max`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`no_ivvv_integrals`
       * :c:data:`no_vvv_integrals`
       * :c:data:`no_vvvv_integrals`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_charge_remove`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`
       * :c:data:`only_expected_s2`
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
       * :c:data:`pt2_iterations`
       * :c:data:`pt2_max`
       * :c:data:`pt2_relative_error`
       * :c:data:`qp_stop_filename`
       * :c:data:`read_wf`
       * :c:data:`s2_eig`
       * :c:data:`scf_algorithm`
       * :c:data:`state_following`
       * :c:data:`target_energy`
       * :c:data:`thresh_scf`
       * :c:data:`thresh_sym`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_diis`
       * :c:data:`threshold_generators`
       * :c:data:`used_weight`
       * :c:data:`variance_max`

 
.. c:var:: ezfio_work_dir


    File : :file:`ezfio_files/ezfio.irp.f`

    .. code:: fortran

        character*(128)	:: ezfio_work_dir	


    EZFIO/work/

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`


 
.. c:var:: file_lock


    File : :file:`ezfio_files/lock.irp.f`

    .. code:: fortran

        integer(omp_lock_kind)	:: file_lock	


    OpenMP Lock for I/O


 
.. c:var:: output_cpu_time_0


    File : :file:`ezfio_files/output.irp.f`

    .. code:: fortran

        double precision	:: output_wall_time_0	
        double precision	:: output_cpu_time_0	


    Initial CPU and wall times when printing in the output files

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cartesian`
       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_md5`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ci_energy`
       * :c:data:`correlation_energy_ratio_max`
       * :c:data:`data_energy_proj`
       * :c:data:`data_energy_var`
       * :c:data:`data_one_e_dm_alpha_ao`
       * :c:data:`data_one_e_dm_alpha_mo`
       * :c:data:`data_one_e_dm_beta_ao`
       * :c:data:`data_one_e_dm_beta_mo`
       * :c:data:`davidson_sze_max`
       * :c:data:`disk_access_nuclear_repulsion`
       * :c:data:`disk_based_davidson`
       * :c:data:`distributed_davidson`
       * :c:data:`do_direct_integrals`
       * :c:data:`do_pseudo`
       * :c:data:`do_pt2`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`energy_iterations`
       * :c:data:`frozen_orb_scf`
       * :c:data:`h0_type`
       * :c:data:`io_ao_integrals_e_n`
       * :c:data:`io_ao_integrals_kinetic`
       * :c:data:`io_ao_integrals_overlap`
       * :c:data:`io_ao_integrals_pseudo`
       * :c:data:`io_ao_one_e_integrals`
       * :c:data:`io_ao_two_e_integrals`
       * :c:data:`io_ao_two_e_integrals_erf`
       * :c:data:`io_mo_integrals_e_n`
       * :c:data:`io_mo_integrals_kinetic`
       * :c:data:`io_mo_integrals_pseudo`
       * :c:data:`io_mo_one_e_integrals`
       * :c:data:`io_mo_two_e_integrals`
       * :c:data:`io_mo_two_e_integrals_erf`
       * :c:data:`level_shift`
       * :c:data:`max_dim_diis`
       * :c:data:`mo_class`
       * :c:data:`mo_guess_type`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mu_erf`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_iterations`
       * :c:data:`n_det_max`
       * :c:data:`n_det_max_full`
       * :c:data:`n_det_print_wf`
       * :c:data:`n_det_selectors`
       * :c:data:`n_it_scf_max`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`no_ivvv_integrals`
       * :c:data:`no_vvv_integrals`
       * :c:data:`no_vvvv_integrals`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_charge_remove`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`
       * :c:data:`nuclear_repulsion`
       * :c:data:`only_expected_s2`
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
       * :c:data:`pt2_iterations`
       * :c:data:`pt2_max`
       * :c:data:`pt2_relative_error`
       * :c:data:`read_wf`
       * :c:data:`s2_eig`
       * :c:data:`scf_algorithm`
       * :c:data:`state_following`
       * :c:data:`target_energy`
       * :c:data:`thresh_scf`
       * :c:data:`thresh_sym`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_diis`
       * :c:data:`threshold_generators`
       * :c:data:`used_weight`
       * :c:data:`variance_max`

 
.. c:var:: output_wall_time_0


    File : :file:`ezfio_files/output.irp.f`

    .. code:: fortran

        double precision	:: output_wall_time_0	
        double precision	:: output_cpu_time_0	


    Initial CPU and wall times when printing in the output files

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cartesian`
       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_md5`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ci_energy`
       * :c:data:`correlation_energy_ratio_max`
       * :c:data:`data_energy_proj`
       * :c:data:`data_energy_var`
       * :c:data:`data_one_e_dm_alpha_ao`
       * :c:data:`data_one_e_dm_alpha_mo`
       * :c:data:`data_one_e_dm_beta_ao`
       * :c:data:`data_one_e_dm_beta_mo`
       * :c:data:`davidson_sze_max`
       * :c:data:`disk_access_nuclear_repulsion`
       * :c:data:`disk_based_davidson`
       * :c:data:`distributed_davidson`
       * :c:data:`do_direct_integrals`
       * :c:data:`do_pseudo`
       * :c:data:`do_pt2`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`energy_iterations`
       * :c:data:`frozen_orb_scf`
       * :c:data:`h0_type`
       * :c:data:`io_ao_integrals_e_n`
       * :c:data:`io_ao_integrals_kinetic`
       * :c:data:`io_ao_integrals_overlap`
       * :c:data:`io_ao_integrals_pseudo`
       * :c:data:`io_ao_one_e_integrals`
       * :c:data:`io_ao_two_e_integrals`
       * :c:data:`io_ao_two_e_integrals_erf`
       * :c:data:`io_mo_integrals_e_n`
       * :c:data:`io_mo_integrals_kinetic`
       * :c:data:`io_mo_integrals_pseudo`
       * :c:data:`io_mo_one_e_integrals`
       * :c:data:`io_mo_two_e_integrals`
       * :c:data:`io_mo_two_e_integrals_erf`
       * :c:data:`level_shift`
       * :c:data:`max_dim_diis`
       * :c:data:`mo_class`
       * :c:data:`mo_guess_type`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mu_erf`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_iterations`
       * :c:data:`n_det_max`
       * :c:data:`n_det_max_full`
       * :c:data:`n_det_print_wf`
       * :c:data:`n_det_selectors`
       * :c:data:`n_it_scf_max`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`no_ivvv_integrals`
       * :c:data:`no_vvv_integrals`
       * :c:data:`no_vvvv_integrals`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_charge_remove`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`
       * :c:data:`nuclear_repulsion`
       * :c:data:`only_expected_s2`
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
       * :c:data:`pt2_iterations`
       * :c:data:`pt2_max`
       * :c:data:`pt2_relative_error`
       * :c:data:`read_wf`
       * :c:data:`s2_eig`
       * :c:data:`scf_algorithm`
       * :c:data:`state_following`
       * :c:data:`target_energy`
       * :c:data:`thresh_scf`
       * :c:data:`thresh_sym`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_diis`
       * :c:data:`threshold_generators`
       * :c:data:`used_weight`
       * :c:data:`variance_max`

 
.. c:var:: qp_kill_filename


    File : :file:`ezfio_files/qp_stop.irp.f`

    .. code:: fortran

        character*(128)	:: qp_stop_filename	
        character*(128)	:: qp_kill_filename	
        integer	:: qp_stop_variable	


    Name of the file to check for qp stop

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`


 
.. c:var:: qp_stop_filename


    File : :file:`ezfio_files/qp_stop.irp.f`

    .. code:: fortran

        character*(128)	:: qp_stop_filename	
        character*(128)	:: qp_kill_filename	
        integer	:: qp_stop_variable	


    Name of the file to check for qp stop

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`


 
.. c:var:: qp_stop_variable


    File : :file:`ezfio_files/qp_stop.irp.f`

    .. code:: fortran

        character*(128)	:: qp_stop_filename	
        character*(128)	:: qp_kill_filename	
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
    

 
.. c:function:: qp_stop:


    File : :file:`ezfio_files/qp_stop.irp.f`

    .. code:: fortran

        logical function qp_stop()


    Checks if the qp_stop command was invoked for the clean termination of the program

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_stop_filename`

 
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

       * :c:data:`ci_energy`
       * :c:func:`damping_scf`
       * :c:func:`davidson_diag_hjj_sjj`
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

       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`make_s2_eigenfunction`
       * :c:data:`mo_num`
       * :c:data:`n_cas_bitmask`
       * :c:data:`n_core_orb`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_generators_bitmask`
       * :c:data:`n_generators_bitmask_restart`
       * :c:data:`n_int`
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

       * :c:data:`output_wall_time_0`
       * :c:data:`mpi_master`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cartesian`
       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_md5`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ci_energy`
       * :c:data:`correlation_energy_ratio_max`
       * :c:func:`damping_scf`
       * :c:data:`data_energy_proj`
       * :c:data:`data_energy_var`
       * :c:data:`data_one_e_dm_alpha_ao`
       * :c:data:`data_one_e_dm_alpha_mo`
       * :c:data:`data_one_e_dm_beta_ao`
       * :c:data:`data_one_e_dm_beta_mo`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:data:`davidson_sze_max`
       * :c:data:`disk_access_nuclear_repulsion`
       * :c:data:`disk_based_davidson`
       * :c:data:`distributed_davidson`
       * :c:data:`do_direct_integrals`
       * :c:data:`do_pseudo`
       * :c:data:`do_pt2`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`energy_iterations`
       * :c:data:`frozen_orb_scf`
       * :c:data:`h0_type`
       * :c:data:`io_ao_integrals_e_n`
       * :c:data:`io_ao_integrals_kinetic`
       * :c:data:`io_ao_integrals_overlap`
       * :c:data:`io_ao_integrals_pseudo`
       * :c:data:`io_ao_one_e_integrals`
       * :c:data:`io_ao_two_e_integrals`
       * :c:data:`io_ao_two_e_integrals_erf`
       * :c:data:`io_mo_integrals_e_n`
       * :c:data:`io_mo_integrals_kinetic`
       * :c:data:`io_mo_integrals_pseudo`
       * :c:data:`io_mo_one_e_integrals`
       * :c:data:`io_mo_two_e_integrals`
       * :c:data:`io_mo_two_e_integrals_erf`
       * :c:data:`level_shift`
       * :c:func:`make_s2_eigenfunction`
       * :c:data:`max_dim_diis`
       * :c:func:`mo_as_eigvectors_of_mo_matrix`
       * :c:func:`mo_as_svd_vectors_of_mo_matrix`
       * :c:func:`mo_as_svd_vectors_of_mo_matrix_eig`
       * :c:data:`mo_class`
       * :c:data:`mo_guess_type`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mu_erf`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_iterations`
       * :c:data:`n_det_max`
       * :c:data:`n_det_max_full`
       * :c:data:`n_det_print_wf`
       * :c:data:`n_det_selectors`
       * :c:data:`n_it_scf_max`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`no_ivvv_integrals`
       * :c:data:`no_vvv_integrals`
       * :c:data:`no_vvvv_integrals`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_charge_remove`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`
       * :c:data:`nuclear_repulsion`
       * :c:data:`only_expected_s2`
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
       * :c:data:`pt2_iterations`
       * :c:data:`pt2_max`
       * :c:data:`pt2_relative_error`
       * :c:data:`read_wf`
       * :c:func:`roothaan_hall_scf`
       * :c:data:`s2_eig`
       * :c:data:`scf_algorithm`
       * :c:data:`state_following`
       * :c:data:`target_energy`
       * :c:data:`thresh_scf`
       * :c:data:`thresh_sym`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_diis`
       * :c:data:`threshold_generators`
       * :c:data:`used_weight`
       * :c:data:`variance_max`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`cpu_time`
       * :c:func:`print_memory_usage`
       * :c:func:`wall_time`

