.. _module_mpi: 
 
.. program:: mpi 
 
.. default-role:: option 
 
===
mpi
===

Contains all the functions and providers for parallelization with |MPI|.
 
 
 
Providers 
--------- 
 
.. c:var:: mpi_initialized


    File : :file:`mpi/mpi.irp.f`

    .. code:: fortran

        logical	:: mpi_initialized	


    Always true. Initialized MPI

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`

 
.. c:var:: mpi_master


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
       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_md5`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`cas_bitmask`
       * :c:data:`ci_energy`
       * :c:data:`core_bitmask`
       * :c:data:`correlation_energy_ratio_max`
       * :c:data:`data_energy_proj`
       * :c:data:`data_energy_var`
       * :c:data:`data_one_e_dm_alpha_mo`
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
       * :c:data:`element_name`
       * :c:data:`energy_iterations`
       * :c:data:`frozen_orb_scf`
       * :c:data:`generators_bitmask`
       * :c:data:`generators_bitmask_restart`
       * :c:data:`inact_bitmask`
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
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mu_erf`
       * :c:data:`n_act_orb`
       * :c:data:`n_cas_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_iterations`
       * :c:data:`n_det_max`
       * :c:data:`n_det_max_full`
       * :c:data:`n_det_print_wf`
       * :c:data:`n_det_selectors`
       * :c:data:`n_generators_bitmask`
       * :c:data:`n_generators_bitmask_restart`
       * :c:data:`n_int`
       * :c:data:`n_it_scf_max`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`no_ivvv_integrals`
       * :c:data:`no_vvv_integrals`
       * :c:data:`no_vvvv_integrals`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`
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
       * :c:data:`pseudo_v_k`
       * :c:data:`pseudo_v_kl`
       * :c:data:`psi_cas`
       * :c:data:`psi_coef`
       * :c:data:`psi_coef_max`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_size`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_iterations`
       * :c:data:`pt2_max`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_relative_error`
       * :c:data:`qp_max_mem`
       * :c:data:`read_wf`
       * :c:data:`s2_eig`
       * :c:data:`scf_algorithm`
       * :c:data:`state_following`
       * :c:data:`target_energy`
       * :c:data:`thresh_scf`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_diis`
       * :c:data:`threshold_generators`
       * :c:data:`used_weight`

 
.. c:var:: mpi_rank


    File : :file:`mpi/mpi.irp.f`

    .. code:: fortran

        integer	:: mpi_rank	
        integer	:: mpi_size	


    Rank of MPI process and number of MPI processes

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

 
.. c:var:: mpi_size


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
 
.. c:function:: broadcast_chunks_double:


    File : :file:`mpi/mpi.irp.f_template_97`

    .. code:: fortran

        subroutine broadcast_chunks_double(A, LDA)


    Broadcast with chunks of ~2GB

 
.. c:function:: broadcast_chunks_integer:


    File : :file:`mpi/mpi.irp.f_template_97`

    .. code:: fortran

        subroutine broadcast_chunks_integer(A, LDA)


    Broadcast with chunks of ~2GB

 
.. c:function:: broadcast_chunks_integer8:


    File : :file:`mpi/mpi.irp.f_template_97`

    .. code:: fortran

        subroutine broadcast_chunks_integer8(A, LDA)


    Broadcast with chunks of ~2GB

 
.. c:function:: mpi_print:


    File : :file:`mpi/mpi.irp.f`

    .. code:: fortran

        subroutine mpi_print(string)


    Print string to stdout if the MPI rank is zero.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_slave_main`

