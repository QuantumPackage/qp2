.. _module_cipsi_tc_bi_ortho: 
 
.. program:: cipsi_tc_bi_ortho 
 
.. default-role:: option 
 
cipsi_tc_bi_ortho
=================

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
EZFIO parameters 
---------------- 
 
.. option:: save_wf_after_selection
 
    If true, saves the wave function after the selection, before the diagonalization
 
    Default: False
 
.. option:: seniority_max
 
    Maximum number of allowed open shells. Using -1 selects all determinants
 
    Default: -1
 
.. option:: excitation_ref
 
    1: Hartree-Fock determinant, 2:All determinants of the dominant configuration
 
    Default: 1
 
.. option:: excitation_max
 
    Maximum number of excitation with respect to the Hartree-Fock determinant. Using -1 selects all determinants
 
    Default: -1
 
.. option:: excitation_alpha_max
 
    Maximum number of excitation for alpha determinants with respect to the Hartree-Fock determinant. Using -1 selects all determinants
 
    Default: -1
 
.. option:: excitation_beta_max
 
    Maximum number of excitation for beta determinants with respect to the Hartree-Fock determinant. Using -1 selects all determinants
 
    Default: -1
 
 
Providers 
--------- 
 
.. c:var:: bool initialize_pt2_e0_denominator


    File : :file:`cipsi_tc_bi_ortho/energy.irp.f`

    .. code:: fortran

        logical	:: initialize_pt2_e0_denominator	


    If true, initialize pt2_E0_denominator


 
.. c:var:: double pt2_e0_denominator


    File : :file:`cipsi_tc_bi_ortho/energy.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pt2_e0_denominator	(N_states)


    E0 in the denominator of the PT2

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`n_states`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void bitstring_to_list_in_selection(void)


    File : :file:`cipsi_tc_bi_ortho/selection.irp.f`

    .. code:: fortran

        subroutine bitstring_to_list_in_selection( string, list, n_elements, Nint)


    Gives the inidices(+1) of the bits set to 1 in the bit string

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`splash_pq`
       * :c:func:`spot_isinwf`

 
.. c:function:: void build_fock_tmp_tc(void)


    File : :file:`cipsi_tc_bi_ortho/fock_diag.irp.f`

    .. code:: fortran

        subroutine build_fock_tmp_tc(fock_diag_tmp,det_ref,Nint)


    Build the diagonal of the Fock matrix corresponding to a generator
    determinant. $F_{00}$ is $\langle i|H|i \rangle = E_0$.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`select_connected`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`debug_det`

 
.. c:function:: void fill_buffer_double(void)


    File : :file:`cipsi_tc_bi_ortho/selection.irp.f`

    .. code:: fortran

        subroutine fill_buffer_double(i_generator, sp, h1, h2, bannedOrb, banned, fock_diag_tmp, E0, pt2_data, mat, buf, mat_l, mat_r)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`c0_weight`
       * :c:data:`det_to_configuration`
       * :c:data:`do_only_1h1p`
       * :c:data:`do_ormas`
       * :c:data:`dominant_dets_of_cfgs`
       * :c:data:`elec_alpha_num`
       * :c:data:`excitation_alpha_max`
       * :c:data:`excitation_beta_max`
       * :c:data:`excitation_max`
       * :c:data:`excitation_ref`
       * :c:data:`h0_type`
       * :c:data:`hf_bitmask`
       * :c:data:`mo_num`
       * :c:data:`n_dominant_dets_of_cfgs`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_configuration_hii`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_det_hii`
       * :c:data:`selection_tc`
       * :c:data:`selection_weight`
       * :c:data:`seniority_max`
       * :c:data:`weight_selection`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`select_singles_and_doubles`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`add_to_selection_buffer`
       * :c:func:`apply_holes`
       * :c:func:`apply_particles`
       * :c:func:`configuration_to_dets_size`
       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_excitation_degree_spin`

 
.. c:function:: void get_d0(void)


    File : :file:`cipsi_tc_bi_ortho/get_d.irp.f`

    .. code:: fortran

        subroutine get_d0(gen, phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_particles`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_no_3e`

 
.. c:function:: void get_d0_new(void)


    File : :file:`cipsi_tc_bi_ortho/get_d0_good.irp.f`

    .. code:: fortran

        subroutine get_d0_new(gen, phasemask, bannedOrb, banned, mat_l, mat_r, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`splash_pq`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_particles`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_no_3e`

 
.. c:function:: void get_d0_reference(void)


    File : :file:`cipsi_tc_bi_ortho/get_d.irp.f`

    .. code:: fortran

        subroutine get_d0_reference(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_particles`
       * :c:func:`i_h_j`

 
.. c:function:: void get_d0_transp(void)


    File : :file:`cipsi_tc_bi_ortho/get_d0_transp.irp.f`

    .. code:: fortran

        subroutine get_d0_transp(gen, phasemask, bannedOrb, banned, mat_l, mat_r, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_bi_ortho_tc_two_e_transp`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`splash_pq`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_particles`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_no_3e_both`

 
.. c:function:: void get_d1(void)


    File : :file:`cipsi_tc_bi_ortho/get_d.irp.f`

    .. code:: fortran

        subroutine get_d1(gen, phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_particles`
       * :c:func:`debug_det`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_no_3e`

 
.. c:function:: void get_d1_new(void)


    File : :file:`cipsi_tc_bi_ortho/get_d1_good.irp.f`

    .. code:: fortran

        subroutine get_d1_new(gen, phasemask, bannedOrb, banned, mat_l, mat_r, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`splash_pq`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_particles`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_no_3e`

 
.. c:function:: void get_d1_reference(void)


    File : :file:`cipsi_tc_bi_ortho/get_d.irp.f`

    .. code:: fortran

        subroutine get_d1_reference(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_particles`
       * :c:func:`i_h_j`

 
.. c:function:: void get_d1_transp(void)


    File : :file:`cipsi_tc_bi_ortho/get_d1_transp.irp.f`

    .. code:: fortran

        subroutine get_d1_transp(gen, phasemask, bannedOrb, banned, mat_l, mat_r, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_bi_ortho_tc_two_e_transp`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`splash_pq`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_particles`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_no_3e_both`

 
.. c:function:: void get_d2(void)


    File : :file:`cipsi_tc_bi_ortho/get_d.irp.f`

    .. code:: fortran

        subroutine get_d2(gen, phasemask, bannedOrb, banned, mat_l, mat_r, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

 
.. c:function:: void get_d2_new(void)


    File : :file:`cipsi_tc_bi_ortho/get_d2_good.irp.f`

    .. code:: fortran

        subroutine get_d2_new(gen, phasemask, bannedOrb, banned, mat_l, mat_r, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`splash_pq`

 
.. c:function:: void get_d2_new_transp(void)


    File : :file:`cipsi_tc_bi_ortho/get_d2_transp.irp.f`

    .. code:: fortran

        subroutine get_d2_new_transp(gen, phasemask, bannedOrb, banned, mat_l, mat_r, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_bi_ortho_tc_two_e_transp`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`splash_pq`

 
.. c:function:: void get_d2_reference(void)


    File : :file:`cipsi_tc_bi_ortho/get_d.irp.f`

    .. code:: fortran

        subroutine get_d2_reference(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

 
.. c:function:: void get_d3_h(void)


    File : :file:`cipsi_tc_bi_ortho/get_d.irp.f`

    .. code:: fortran

        subroutine get_d3_h(gen, bannedOrb, banned, mat, mask, p, sp, coefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_particles`
       * :c:func:`i_h_j`

 
.. c:function:: void get_d3_htc(void)


    File : :file:`cipsi_tc_bi_ortho/get_d.irp.f`

    .. code:: fortran

        subroutine get_d3_htc(gen, bannedOrb, banned, mat_m, mat_p, mask, p, sp, rcoefs, lcoefs)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_particles`
       * :c:func:`htilde_mu_mat_bi_ortho_tot`

 
.. c:function:: void get_mask_phase(void)


    File : :file:`cipsi_tc_bi_ortho/selection.irp.f`

    .. code:: fortran

        subroutine get_mask_phase(det1, pm, Nint)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`splash_pq`

 
.. c:function:: void get_phase_bi(void)


    File : :file:`cipsi_tc_bi_ortho/selection.irp.f`

    .. code:: fortran

        double precision function get_phase_bi(phasemask, s1, s2, h1, p1, h2, p2, Nint)



 
.. c:function:: void past_d1(void)


    File : :file:`cipsi_tc_bi_ortho/selection.irp.f`

    .. code:: fortran

        subroutine past_d1(bannedOrb, p)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`splash_pq`

 
.. c:function:: void past_d2(void)


    File : :file:`cipsi_tc_bi_ortho/selection.irp.f`

    .. code:: fortran

        subroutine past_d2(banned, p, sp)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`splash_pq`

 
.. c:function:: void provide_for_tc_pt2(void)


    File : :file:`cipsi_tc_bi_ortho/pt2_stoch_routines.irp.f`

    .. code:: fortran

        subroutine provide_for_tc_pt2



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dominant_dets_of_cfgs`
       * :c:data:`n_dominant_dets_of_cfgs`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_selectors_coef_transp`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`
       * :c:func:`tc_pt2`

 
.. c:function:: void run_cipsi(void)


    File : :file:`cipsi_tc_bi_ortho/cipsi.irp.f`

    .. code:: fortran

        subroutine run_cipsi


    Selected Full Configuration Interaction with deterministic selection and
    stochastic PT2.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`distributed_davidson`
       * :c:data:`do_pt2`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_det_max`
       * :c:data:`n_int`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_det_sorted_gen`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`pt2_max`
       * :c:data:`pt2_relative_error`
       * :c:data:`ref_bitmask_energy`
       * :c:data:`s2_eig`
       * :c:data:`save_wf_after_selection`
       * :c:data:`selection_factor`
       * :c:data:`threshold_generators`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi_tc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`copy_h_apply_buffer_to_wf`
       * :c:func:`diagonalize_ci_tc_bi_ortho`
       * :c:func:`ezfio_get_hartree_fock_energy`
       * :c:func:`ezfio_has_hartree_fock_energy`
       * :c:func:`make_s2_eigenfunction`
       * :c:func:`omp_selection`
       * :c:func:`provide_for_tc_pt2`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`pt2_stochastic`
       * :c:func:`save_wavefunction`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`psi_configuration`
       * :c:data:`n_det`
       * :c:data:`c0_weight`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_energy`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`psi_energy`
       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`threshold_generators`
       * :c:data:`variance_match_weight`

 
.. c:function:: void run_pouet(void)


    File : :file:`cipsi_tc_bi_ortho/stochastic_cipsi.irp.f`

    .. code:: fortran

        subroutine run_pouet


    Selected Full Configuration Interaction with Stochastic selection and PT2.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`distributed_davidson`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

 
.. c:function:: void run_stochastic_cipsi(void)


    File : :file:`cipsi_tc_bi_ortho/stochastic_cipsi.irp.f`

    .. code:: fortran

        subroutine run_stochastic_cipsi


    Selected Full Configuration Interaction with Stochastic selection and PT2.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`angle_left_right`
       * :c:data:`distributed_davidson`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`mo_bi_ortho_tc_two_e_transp`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_det`
       * :c:data:`n_det_max`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy_with_nucl_rep`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`pt2_max`
       * :c:data:`pt2_relative_error`
       * :c:data:`selection_factor`
       * :c:data:`tc_2e_3idx_coulomb_integrals_transp`
       * :c:data:`threshold_generators`
       * :c:data:`transpose_two_e_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi_tc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`copy_h_apply_buffer_to_wf_tc`
       * :c:func:`diagonalize_ci_tc_bi_ortho`
       * :c:func:`increment_n_iter`
       * :c:func:`print_extrapolated_energy`
       * :c:func:`print_summary_tc`
       * :c:func:`provide_for_tc_pt2`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`pt2_stochastic`
       * :c:func:`save_energy`
       * :c:func:`write_cipsi_json`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`c0_weight`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_energy`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_energy`
       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`threshold_generators`
       * :c:data:`variance_match_weight`

 
.. c:function:: void select_connected(void)


    File : :file:`cipsi_tc_bi_ortho/selection.irp.f`

    .. code:: fortran

        subroutine select_connected(i_generator,E0,pt2_data,b,subset,csubset)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_generators`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection_task`
       * :c:func:`pt2_omp_task`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`build_fock_tmp_tc`
       * :c:func:`select_singles_and_doubles`

 
.. c:function:: void select_singles_and_doubles(void)


    File : :file:`cipsi_tc_bi_ortho/selection.irp.f`

    .. code:: fortran

        subroutine select_singles_and_doubles(i_generator, hole_mask, particle_mask, fock_diag_tmp, E0, pt2_data, buf, subset, csubset)


    WARNING /!\ : It is assumed that the generators and selectors are psi_det_sorted_tc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_coef_transp`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`select_connected`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_hole`
       * :c:func:`bitstring_to_list_ab`
       * :c:func:`fill_buffer_double`
       * :c:func:`get_excitation_degree_spin`
       * :c:func:`isort_noidx`
       * :c:func:`splash_pq`
       * :c:func:`spot_isinwf`

 
.. c:function:: void splash_pq(void)


    File : :file:`cipsi_tc_bi_ortho/selection.irp.f`

    .. code:: fortran

        subroutine splash_pq(mask, sp, det, i_gen, N_sel, bannedOrb, banned, mat, interesting, mat_l, mat_r)


    Computes the contributions A(r,s) by
    comparing the external determinant to all the internal determinants det(i).
    an applying two particles (r,s) to the mask.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_selectors_coef_transp`
       * :c:data:`psi_selectors_coef_transp`
       * :c:data:`transpose_two_e_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`select_singles_and_doubles`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_in_selection`
       * :c:func:`get_d0_new`
       * :c:func:`get_d0_transp`
       * :c:func:`get_d1_new`
       * :c:func:`get_d1_transp`
       * :c:func:`get_d2_new`
       * :c:func:`get_d2_new_transp`
       * :c:func:`get_mask_phase`
       * :c:func:`past_d1`
       * :c:func:`past_d2`

 
.. c:function:: void spot_isinwf(void)


    File : :file:`cipsi_tc_bi_ortho/selection.irp.f`

    .. code:: fortran

        subroutine spot_isinwf(mask, det, i_gen, N, banned, fullMatch, interesting)


    Identify the determinants in det which are in the internal space. These are
    the determinants that can be produced by creating two particles on the mask.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`select_singles_and_doubles`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_in_selection`

 
.. c:function:: void tc_pt2(void)


    File : :file:`cipsi_tc_bi_ortho/pt2.irp.f`

    .. code:: fortran

        subroutine tc_pt2


    Selected Full Configuration Interaction with Stochastic selection and PT2.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`distributed_davidson`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`mo_bi_ortho_tc_two_e_transp`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`psi_configuration`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy_with_nucl_rep`
       * :c:data:`pt2_relative_error`
       * :c:data:`tc_2e_3idx_coulomb_integrals_transp`
       * :c:data:`threshold_generators`
       * :c:data:`transpose_two_e_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_pt2_tc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`diagonalize_ci_tc_bi_ortho`
       * :c:func:`print_summary_tc`
       * :c:func:`provide_for_tc_pt2`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`pt2_stochastic`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`psi_coef`
       * :c:data:`psi_energy`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`psi_energy`
       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`threshold_generators`
       * :c:data:`variance_match_weight`

 
.. c:function:: void write_cipsi_json(void)


    File : :file:`cipsi_tc_bi_ortho/write_cipsi_json.irp.f`

    .. code:: fortran

        subroutine write_cipsi_json(pt2_data, pt2_data_err)


    Writes JSON data for CIPSI runs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`energy_iterations`
       * :c:data:`json_int_fmt`
       * :c:data:`json_unit`
       * :c:data:`n_csf`
       * :c:data:`n_det`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_configuration`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy_with_nucl_rep`
       * :c:data:`s2_eig`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_stochastic_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lock_io`
       * :c:func:`unlock_io`

