.. _module_cipsi_tc_bi_ortho: 
 
.. program:: cipsi_tc_bi_ortho 
 
.. default-role:: option 
 
cipsi_tc_bi_ortho
=================

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: build_fock_tmp_tc:


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

 
.. c:function:: get_d0_new:


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

 
.. c:function:: get_d1_new:


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

 
.. c:function:: get_d2_new:


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

 
.. c:function:: get_d3_h:


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

 
.. c:function:: get_d3_htc:


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

 
.. c:function:: past_d1:


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

 
.. c:function:: past_d2:


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

 
.. c:function:: pt2_tc_bi_ortho:


    File : :file:`cipsi_tc_bi_ortho/pt2.irp.f`

    .. code:: fortran

        subroutine pt2_tc_bi_ortho


    Selected Full Configuration Interaction with Stochastic selection and PT2.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`distributed_davidson`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_det`
       * :c:data:`n_det_max`
       * :c:data:`n_int`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_det_sorted_tc_gen`
       * :c:data:`pt2_relative_error`
       * :c:data:`s2_eig`
       * :c:data:`selection_factor`
       * :c:data:`threshold_generators`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`diagonalize_ci_tc_bi_ortho`
       * :c:func:`make_s2_eigenfunction`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`zmq_pt2`

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
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`threshold_generators`
       * :c:data:`variance_match_weight`

