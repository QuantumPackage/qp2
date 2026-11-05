.. _module_ezfio_files: 
 
.. program:: ezfio_files 
 
.. default-role:: option 
 
===========
ezfio_files
===========

Access to the |EZFIO| database and to the basic execution
information.

The provider :c:data:`ezfio_filename` is the name of the |EZFIO|
directory: it is read from the :envvar:`QPACKAGE_INPUT` environment
variable, or from the first argument of the command line.

The other providers and subroutines of this module are:

* the process and thread counts: :c:data:`nproc`,
  :c:data:`nproc_max`, :c:data:`nthreads_pt2`,
  :c:data:`nthreads_davidson`
* the file locking: :c:func:`lock_io`, :c:func:`unlock_io`,
  :c:data:`file_lock`, and :c:func:`getunitandopen` to open a file
  with a free unit number
* the output helpers: :c:func:`write_time`, :c:func:`wall_time`,
  :c:data:`output_wall_time_0`, :c:data:`output_cpu_time_0`, and
  :c:func:`write_bool`, :c:func:`write_int`, :c:func:`write_double`
* the clean stopping of the calculation: :c:data:`qp_stop`,
  :c:data:`qp_stop_variable`, :c:data:`qp_stop_filename` and
  :c:data:`qp_kill_filename`
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: ezfio_convention
 
    Version of the EZFIO conventions
 
    Default: 20210101
 
 
Providers 
--------- 
 
.. c:var:: char ezfio_filename


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
       * :c:data:`nproc`
       * :c:data:`output_wall_time_0`

    Needed by:

    .. hlist::
       :columns: 3

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
       * :c:data:`auto_symmetry`
       * :c:data:`basis`
       * :c:data:`basis_nucleus_index`
       * :c:data:`basis_type`
       * :c:data:`calc_dipole_moment`
       * :c:data:`calc_energy_components`
       * :c:data:`calc_osc_str`
       * :c:data:`calc_tr_dipole_moment`
       * :c:data:`correlation_functional`
       * :c:data:`csf_based`
       * :c:data:`damping_for_rs_dft`
       * :c:data:`data_energy_proj`
       * :c:data:`data_energy_var`
       * :c:data:`data_one_e_dm_alpha_ao`
       * :c:data:`data_one_e_dm_alpha_mo`
       * :c:data:`data_one_e_dm_beta_ao`
       * :c:data:`data_one_e_dm_beta_mo`
       * :c:data:`data_one_e_dm_tot_ao`
       * :c:data:`davidson_sze_max`
       * :c:data:`density_for_dft`
       * :c:data:`dft_grid_mo_dipole`
       * :c:data:`disk_based_davidson`
       * :c:data:`distributed_davidson`
       * :c:data:`do_ao_cholesky`
       * :c:data:`do_mo_cholesky`
       * :c:data:`do_ormas`
       * :c:data:`do_pseudo`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`exchange_functional`
       * :c:data:`extra_grid_type_sgn`
       * :c:data:`extra_rad_grid_type`
       * :c:data:`ezfio_convention`
       * :c:data:`ezfio_work_dir`
       * :c:data:`generators_bitmask`
       * :c:data:`grid_type_sgn`
       * :c:data:`hf_exchange`
       * :c:data:`io_ao_cholesky`
       * :c:data:`io_ao_integrals_kinetic`
       * :c:data:`io_ao_integrals_n_e`
       * :c:data:`io_ao_integrals_overlap`
       * :c:data:`io_ao_integrals_pseudo`
       * :c:data:`io_ao_one_e_integrals`
       * :c:data:`io_ao_two_e_integrals`
       * :c:data:`io_ao_two_e_integrals_erf`
       * :c:data:`io_extra_grid`
       * :c:data:`io_mo_cholesky`
       * :c:data:`io_mo_integrals_kinetic`
       * :c:data:`io_mo_integrals_n_e`
       * :c:data:`io_mo_integrals_pseudo`
       * :c:data:`io_mo_one_e_integrals`
       * :c:data:`io_mo_two_e_integrals`
       * :c:data:`io_mo_two_e_integrals_erf`
       * :c:data:`io_nuclear_repulsion`
       * :c:data:`irrep_states`
       * :c:data:`is_periodic`
       * :c:data:`lin_dep_cutoff`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_aux`
       * :c:data:`mo_coef_imag`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_label`
       * :c:data:`mo_num`
       * :c:data:`mo_occ`
       * :c:data:`mo_symmetry`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mu_dft_type`
       * :c:data:`mu_erf`
       * :c:data:`my_extra_grid_becke`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_extra_grid`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_extra_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`n_det`
       * :c:data:`n_det_max`
       * :c:data:`n_det_max_full`
       * :c:data:`n_det_print_wf`
       * :c:data:`n_pts_charge`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`no_core_density`
       * :c:data:`normalize_dm`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_charge_remove`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`
       * :c:data:`nucleus_shell_num`
       * :c:data:`only_expected_s2`
       * :c:data:`ormas_max_e`
       * :c:data:`ormas_min_e`
       * :c:data:`ormas_mstart`
       * :c:data:`ormas_n_space`
       * :c:data:`pg_thresh`
       * :c:data:`pg_tol`
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
       * :c:data:`pts_charge_coord`
       * :c:data:`pts_charge_z`
       * :c:data:`qp_stop_filename`
       * :c:data:`rad_grid_type`
       * :c:data:`read_wf`
       * :c:data:`restore_symm`
       * :c:data:`s2_eig`
       * :c:data:`save_threshold`
       * :c:data:`selection_factor`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_index`
       * :c:data:`shell_num`
       * :c:data:`shell_prim_num`
       * :c:data:`slater_r_power`
       * :c:data:`state_following`
       * :c:data:`target_energy`
       * :c:data:`thresh_extra_grid`
       * :c:data:`thresh_grid`
       * :c:data:`thresh_sym`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_davidson_from_pt2`
       * :c:data:`threshold_generators`
       * :c:data:`threshold_nonsym_davidson`
       * :c:data:`use_cgtos`
       * :c:data:`use_only_lr`
       * :c:data:`weight_one_e_dm`
       * :c:data:`weight_selection`

 
.. c:var:: char ezfio_work_dir


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

 
.. c:var:: char ezfio_work_dir_pid


    File : :file:`ezfio_files/ezfio.irp.f`

    .. code:: fortran

        character*(1024)	:: ezfio_work_dir_pid	


    EZFIO/work/pid\_

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_work_dir`


 
.. c:var:: int file_lock


    File : :file:`ezfio_files/lock.irp.f`

    .. code:: fortran

        integer(omp_lock_kind)	:: file_lock	


    OpenMP Lock for I/O

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`
       * :c:data:`qp_max_mem`

 
.. c:var:: int nproc


    File : :file:`ezfio_files/environment.irp.f`

    .. code:: fortran

        integer	:: nproc	


    Number of default OpenMP threads

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_ao_num`
       * :c:data:`ezfio_filename`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`n_det`
       * :c:data:`nproc_max`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`

 
.. c:var:: int nproc_max


    File : :file:`ezfio_files/environment.irp.f`

    .. code:: fortran

        integer	:: nproc_max	


    Max number of threads

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_integral_buffer1`

 
.. c:var:: int nthreads_davidson


    File : :file:`ezfio_files/environment.irp.f`

    .. code:: fortran

        integer	:: nthreads_davidson	


    Number of threads for Davidson

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`file_lock`
       * :c:data:`mpi_master`
       * :c:data:`nproc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`nproc_max`

 
.. c:var:: int nthreads_pt2


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_busy`
       * :c:data:`gpu_num`
       * :c:data:`nproc_max`

 
.. c:var:: double output_cpu_time_0


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
       * :c:data:`n_pts_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nuclear_repulsion`
       * :c:data:`prim_normalization_factor`

 
.. c:var:: double output_wall_time_0


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
       * :c:data:`n_pts_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nuclear_repulsion`
       * :c:data:`prim_normalization_factor`

 
.. c:var:: char qp_kill_filename


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


 
.. c:var:: char qp_stop_filename


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


 
.. c:var:: int qp_stop_variable


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
 
.. c:function:: void getunitandopen(void)


    File : :file:`ezfio_files/get_unit_and_open.irp.f`

    .. code:: fortran

        integer function getUnitAndOpen(f,mode)


    Returns a free Fortran unit number, searched from 20 upward,
    on which file `f` is opened in the mode given by `mode`.
    In modes 'r' and 'R', a missing file is created as empty
    before the file is opened for reading.
    
    :f:
       file name
    
    :mode:
       'R' : READ, UNFORMATTED
       'W' : WRITE, UNFORMATTED
       'A' : APPEND, UNFORMATTED
       'r' : READ, FORMATTED
       'w' : WRITE, FORMATTED
       'a' : APPEND, FORMATTED
       'x' : READ/WRITE, FORMATTED
    

 
.. c:function:: void lock_io(void)


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
       * :c:func:`load_mo_integrals`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`
       * :c:data:`qp_max_mem`
       * :c:func:`resident_memory`
       * :c:func:`total_memory`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`

 
.. c:function:: void qp_stop(void)


    File : :file:`ezfio_files/qp_stop.irp.f`

    .. code:: fortran

        logical function qp_stop()


    Checks if the qp_stop command was invoked for the clean termination of the program

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_stop_filename`

 
.. c:function:: void unlock_io(void)


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
       * :c:func:`load_mo_integrals`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`
       * :c:data:`qp_max_mem`
       * :c:func:`resident_memory`
       * :c:func:`total_memory`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_unset_lock`

 
.. c:function:: void write_bool(void)


    File : :file:`ezfio_files/output.irp.f`

    .. code:: fortran

        subroutine write_bool(iunit,value,label)


    Write an logical value in output

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

 
.. c:function:: void write_double(void)


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
       * :c:func:`davidson_diag_csf_hjj`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_coef_max`

 
.. c:function:: void write_int(void)


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

       * :c:func:`davidson_diag_csf_hjj`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`make_s2_eigenfunction`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_del_orb`
       * :c:data:`n_det`
       * :c:data:`n_inact_orb`
       * :c:data:`n_int`
       * :c:data:`n_virt_orb`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`
       * :c:data:`psi_cas`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_size`
       * :c:data:`qp_max_mem`
       * :c:func:`save_wavefunction_general`
       * :c:func:`save_wavefunction_general_unormalized`
       * :c:func:`save_wavefunction_specified`

 
.. c:function:: void write_time(void)


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
       * :c:func:`davidson_diag_csf_hjj`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`make_s2_eigenfunction`
       * :c:func:`mo_as_eigvectors_of_mo_matrix`
       * :c:func:`mo_as_svd_vectors_of_mo_matrix`
       * :c:func:`mo_as_svd_vectors_of_mo_matrix_eig`
       * :c:func:`mo_coef_new_as_svd_vectors_of_mo_matrix_eig`
       * :c:data:`n_pts_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nuclear_repulsion`
       * :c:data:`prim_normalization_factor`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`cpu_time`
       * :c:func:`print_memory_usage`
       * :c:func:`wall_time`

