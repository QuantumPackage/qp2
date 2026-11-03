.. _module_mpi: 
 
.. program:: mpi 
 
.. default-role:: option 
 
===
mpi
===

Initialization and basic routines for the parallelization with |MPI|.

The provider :c:data:`mpi_initialized` starts |MPI|. The providers
:c:data:`mpi_rank`, :c:data:`mpi_size` and :c:data:`mpi_master` give
the rank of the process, the number of processes, and the boolean that
is |true| only on the master rank.

:c:func:`mpi_print` prints a string only on the master rank. The
subroutines :c:func:`broadcast_chunks_double`,
:c:func:`broadcast_chunks_integer` and
:c:func:`broadcast_chunks_integer8` broadcast an array from the master
rank in chunks of about 2 GB.
 
 
 
Providers 
--------- 
 
.. c:var:: bool mpi_initialized


    File : :file:`mpi/mpi.irp.f`

    .. code:: fortran

        logical	:: mpi_initialized	


    Always true. Initialized MPI

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`

 
.. c:var:: bool mpi_master


    File : :file:`mpi/mpi.irp.f`

    .. code:: fortran

        logical	:: mpi_master	


    If true, rank is zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_rank`

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
       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`auto_symmetry`
       * :c:data:`basis`
       * :c:data:`basis_nucleus_index`
       * :c:data:`basis_type`
       * :c:data:`calc_dipole_moment`
       * :c:data:`calc_energy_components`
       * :c:data:`calc_osc_str`
       * :c:data:`calc_tr_dipole_moment`
       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
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
       * :c:data:`element_name`
       * :c:data:`exchange_functional`
       * :c:data:`extra_grid_type_sgn`
       * :c:data:`extra_rad_grid_type`
       * :c:data:`ezfio_convention`
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
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_del_orb`
       * :c:data:`n_det`
       * :c:data:`n_det_max`
       * :c:data:`n_det_max_full`
       * :c:data:`n_det_print_wf`
       * :c:data:`n_inact_orb`
       * :c:data:`n_int`
       * :c:data:`n_pts_charge`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`n_virt_orb`
       * :c:data:`no_core_density`
       * :c:data:`normalize_dm`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_charge_remove`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`
       * :c:data:`nuclear_repulsion`
       * :c:data:`nucleus_shell_num`
       * :c:data:`only_expected_s2`
       * :c:data:`ormas_max_e`
       * :c:data:`ormas_min_e`
       * :c:data:`ormas_mstart`
       * :c:data:`ormas_n_space`
       * :c:data:`pg_ao_ok`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_irrep_filter`
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
       * :c:data:`psi_cas`
       * :c:data:`psi_coef`
       * :c:data:`psi_coef_max`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_size`
       * :c:data:`pts_charge_coord`
       * :c:data:`pts_charge_z`
       * :c:data:`qp_max_mem`
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

 
.. c:var:: int mpi_rank


    File : :file:`mpi/mpi.irp.f`

    .. code:: fortran

        integer	:: mpi_rank	
        integer	:: mpi_size	


    Rank of MPI process and number of MPI processes

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

 
.. c:var:: int mpi_size


    File : :file:`mpi/mpi.irp.f`

    .. code:: fortran

        integer	:: mpi_rank	
        integer	:: mpi_size	


    Rank of MPI process and number of MPI processes

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void broadcast_chunks_double(void)


    File : :file:`mpi/mpi.irp.f_template_97`

    .. code:: fortran

        subroutine broadcast_chunks_double(A, LDA)


    Broadcast with chunks of ~2GB

 
.. c:function:: void broadcast_chunks_integer(void)


    File : :file:`mpi/mpi.irp.f_template_97`

    .. code:: fortran

        subroutine broadcast_chunks_integer(A, LDA)


    Broadcast with chunks of ~2GB

 
.. c:function:: void broadcast_chunks_integer8(void)


    File : :file:`mpi/mpi.irp.f_template_97`

    .. code:: fortran

        subroutine broadcast_chunks_integer8(A, LDA)


    Broadcast with chunks of ~2GB

 
.. c:function:: void mpi_print(void)


    File : :file:`mpi/mpi.irp.f`

    .. code:: fortran

        subroutine mpi_print(string)


    Print string to stdout if the MPI rank is zero.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

