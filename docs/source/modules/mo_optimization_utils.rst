.. _module_mo_optimization_utils: 
 
.. program:: mo_optimization_utils 
 
.. default-role:: option 
 
mo_optimization_utils
=====================

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
EZFIO parameters 
---------------- 
 
.. option:: optimization_method
 
    Define the kind of hessian for the orbital optimization full : full hessian, diag : diagonal hessian, none : no hessian
 
    Default: diag
 
.. option:: n_det_max_opt
 
    Maximal number of the determinants in the wf for the orbital optimization (to stop the optimization if n_det > n_det_max_opt)
 
    Default: 200000
 
.. option:: optimization_max_nb_iter
 
    Maximal number of iterations for the orbital optimization
 
    Default: 10
 
.. option:: thresh_opt_max_elem_grad
 
    Threshold for the convergence, the optimization exits when the biggest element in the gradient is smaller than thresh_optimization_max_elem_grad
 
    Default: 1.e-5
 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void diag_hessian_list_opt(void)


    File : :file:`mo_optimization_utils/diagonal_hessian_list_opt.irp.f`

    .. code:: fortran

        subroutine diag_hessian_list_opt(n, m, list, H)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`omp_set_max_active_levels`
       * :c:func:`vec_to_mat_index`
       * :c:func:`wall_time`

 
.. c:function:: void diag_hessian_opt(void)


    File : :file:`mo_optimization_utils/diagonal_hessian_opt.irp.f`

    .. code:: fortran

        subroutine diag_hessian_opt(n,H)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`omp_set_max_active_levels`
       * :c:func:`vec_to_mat_index`
       * :c:func:`wall_time`

 
.. c:function:: void diagonalization_hessian(void)


    File : :file:`mo_optimization_utils/diagonalization_hessian.irp.f`

    .. code:: fortran

        subroutine diagonalization_hessian(n,H,e_val,w)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`optimization_method`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`abort`
       * :c:func:`dsort`
       * :c:func:`dsyev`
       * :c:func:`wall_time`

 
.. c:function:: void first_diag_hessian_list_opt(void)


    File : :file:`mo_optimization_utils/first_diagonal_hessian_list_opt.irp.f`

    .. code:: fortran

        subroutine first_diag_hessian_list_opt(tmp_n,m,list,H)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`vec_to_mat_index`
       * :c:func:`wall_time`

 
.. c:function:: void first_diag_hessian_opt(void)


    File : :file:`mo_optimization_utils/first_diagonal_hessian_opt.irp.f`

    .. code:: fortran

        subroutine first_diag_hessian_opt(n,H, h_tmpr)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`vec_to_mat_index`
       * :c:func:`wall_time`

 
.. c:function:: void first_gradient_list_opt(void)


    File : :file:`mo_optimization_utils/first_gradient_list_opt.irp.f`

    .. code:: fortran

        subroutine first_gradient_list_opt(tmp_n,m,list,v_grad)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`vec_to_mat_index`

 
.. c:function:: void first_gradient_opt(void)


    File : :file:`mo_optimization_utils/first_gradient_opt.irp.f`

    .. code:: fortran

        subroutine first_gradient_opt(n,v_grad)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`vec_to_mat_index`

 
.. c:function:: void first_hessian_list_opt(void)


    File : :file:`mo_optimization_utils/first_hessian_list_opt.irp.f`

    .. code:: fortran

        subroutine first_hessian_list_opt(tmp_n,m,list,H,h_tmpr)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`vec_to_mat_index`
       * :c:func:`wall_time`

 
.. c:function:: void first_hessian_opt(void)


    File : :file:`mo_optimization_utils/first_hessian_opt.irp.f`

    .. code:: fortran

        subroutine first_hessian_opt(n,H,h_tmpr)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`vec_to_mat_index`
       * :c:func:`wall_time`

 
.. c:function:: void gradient_list_opt(void)


    File : :file:`mo_optimization_utils/gradient_list_opt.irp.f`

    .. code:: fortran

        subroutine gradient_list_opt(n,m,list,v_grad,max_elem,norm)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`omp_set_max_active_levels`
       * :c:func:`vec_to_mat_index`
       * :c:func:`wall_time`

 
.. c:function:: void gradient_opt(void)


    File : :file:`mo_optimization_utils/gradient_opt.irp.f`

    .. code:: fortran

        subroutine gradient_opt(n,v_grad,max_elem)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`omp_set_max_active_levels`
       * :c:func:`vec_to_mat_index`
       * :c:func:`wall_time`

 
.. c:function:: void hessian_list_opt(void)


    File : :file:`mo_optimization_utils/hessian_list_opt.irp.f`

    .. code:: fortran

        subroutine hessian_list_opt(n,m,list,H,h_tmpr)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`omp_set_max_active_levels`
       * :c:func:`vec_to_mat_index`
       * :c:func:`wall_time`

 
.. c:function:: void hessian_opt(void)


    File : :file:`mo_optimization_utils/hessian_opt.irp.f`

    .. code:: fortran

        subroutine hessian_opt(n,H,h_tmpr)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`omp_set_max_active_levels`
       * :c:func:`vec_to_mat_index`
       * :c:func:`wall_time`

 
.. c:function:: void run_optimization_mos_cipsi(void)


    File : :file:`mo_optimization_utils/routine_opt_mos.irp.f`

    .. code:: fortran

        subroutine run_optimization_mos_CIPSI



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`ao_pseudo_integrals`
       * :c:data:`ezfio_filename`
       * :c:data:`n_det`
       * :c:data:`n_det_max`
       * :c:data:`n_det_max_opt`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`
       * :c:func:`run_stochastic_cipsi`
       * :c:func:`state_average_energy`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`psi_configuration`
       * :c:data:`n_det`
       * :c:data:`n_det_max`
       * :c:data:`c0_weight`
       * :c:data:`ci_energy`
       * :c:data:`mo_coef`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy`
       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`threshold_generators`
       * :c:data:`two_e_dm_mo`
       * :c:data:`variance_match_weight`

 
.. c:function:: void run_orb_opt_trust_v2(void)


    File : :file:`mo_optimization_utils/run_orb_opt_trust_v2.irp.f`

    .. code:: fortran

        subroutine run_orb_opt_trust_v2


    Orbital optimization

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`absolute_eig`
       * :c:data:`all_mo_integrals`
       * :c:data:`ao_num`
       * :c:data:`ci_energy`
       * :c:data:`dim_list_act_orb`
       * :c:data:`list_act`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`optimization_max_nb_iter`
       * :c:data:`optimization_method`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`thresh_opt_max_elem_grad`
       * :c:data:`two_e_dm_mo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_optimization_mos_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`abort`
       * :c:func:`apply_mo_rotation`
       * :c:func:`clear_mo_map`
       * :c:func:`diag_hessian_list_opt`
       * :c:func:`diagonalization_hessian`
       * :c:func:`diagonalize_ci`
       * :c:func:`dsort`
       * :c:func:`gradient_list_opt`
       * :c:func:`hessian_list_opt`
       * :c:func:`rotation_matrix`
       * :c:func:`save_mos`
       * :c:func:`save_wavefunction_unsorted`
       * :c:func:`state_average_energy`
       * :c:func:`state_weight_normalization`
       * :c:func:`sub_to_full_rotation_matrix`
       * :c:func:`trust_region_is_step_cancelled`
       * :c:func:`trust_region_step_w_expected_e`
       * :c:func:`vec_to_mat_v2`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
       * :c:data:`mo_coef`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy`
       * :c:data:`state_average_weight`
       * :c:data:`two_e_dm_mo`

 
.. c:function:: void state_average_energy(void)


    File : :file:`mo_optimization_utils/state_average_energy.irp.f`

    .. code:: fortran

        subroutine state_average_energy(energy)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`nuclear_repulsion`
       * :c:data:`one_e_dm_mo`
       * :c:data:`two_e_dm_mo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_optimization_mos_cipsi`
       * :c:func:`run_orb_opt_trust_v2`

 
.. c:function:: void state_weight_normalization(void)


    File : :file:`mo_optimization_utils/state_weight_normalization.irp.f`

    .. code:: fortran

        subroutine state_weight_normalization


    Renormalization of the state weights or enforcing state average
    weights for orbital optimization

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`state_average_weight`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`state_average_weight`

 
.. c:function:: void update_parameters(void)


    File : :file:`mo_optimization_utils/update_parameters.irp.f`

    .. code:: fortran

        subroutine update_parameters()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_coef`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`clear_mo_map`
       * :c:func:`diagonalize_ci`
       * :c:func:`save_wavefunction_unsorted`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`mo_coef`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy`

 
.. c:function:: void update_st_av_ci_energy(void)


    File : :file:`mo_optimization_utils/update_st_av_ci_energy.irp.f`

    .. code:: fortran

        subroutine update_st_av_ci_energy(energy)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy`
       * :c:data:`n_states`
       * :c:data:`state_average_weight`

