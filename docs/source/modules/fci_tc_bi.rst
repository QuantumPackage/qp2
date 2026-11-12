.. _module_fci_tc_bi: 
 
.. program:: fci_tc_bi 
 
.. default-role:: option 
 
fci_tc_bi
=========

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
EZFIO parameters 
---------------- 
 
.. option:: energy
 
    Calculated Selected |FCI| energy
 
 
.. option:: energy_pt2
 
    Calculated |FCI| energy + |PT2|
 
 
.. option:: cipsi_tc
 
    TODO
 
    Default: h_tc
 
 
Programs 
-------- 
 
 * :ref:`fci_tc_bi_ortho` 
 * :ref:`pt2_tc` 
 
Providers 
--------- 
 
.. c:var:: bool do_ddci


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
.. c:var:: bool do_only_1h1p


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
.. c:var:: bool do_only_cas


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
.. c:var:: int n_det_selectors


    File : :file:`selectors.irp.f`

    .. code:: fortran

        integer	:: n_det_selectors	


    For Single reference wave functions, the number of selectors is 1 : the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`output_wall_time_0`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_coef_transp`

 
.. c:var:: int psi_det_sorted_bit_tc


    File : :file:`copy_wf.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_bit_tc	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_r_coef_sorted_bit	(N_det,N_states)
        double precision, allocatable	:: psi_l_coef_sorted_bit	(N_det,N_states)


    Determinants on which we apply :math:`\langle i|H|psi \rangle`  for perturbation.
    They are sorted by determinants interpreted as integers. Useful
    to accelerate the search of a random determinant in the wave
    function.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: double psi_l_coef_sorted_bit


    File : :file:`copy_wf.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_bit_tc	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_r_coef_sorted_bit	(N_det,N_states)
        double precision, allocatable	:: psi_l_coef_sorted_bit	(N_det,N_states)


    Determinants on which we apply :math:`\langle i|H|psi \rangle`  for perturbation.
    They are sorted by determinants interpreted as integers. Useful
    to accelerate the search of a random determinant in the wave
    function.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: double psi_r_coef_sorted_bit


    File : :file:`copy_wf.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_bit_tc	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_r_coef_sorted_bit	(N_det,N_states)
        double precision, allocatable	:: psi_l_coef_sorted_bit	(N_det,N_states)


    Determinants on which we apply :math:`\langle i|H|psi \rangle`  for perturbation.
    They are sorted by determinants interpreted as integers. Useful
    to accelerate the search of a random determinant in the wave
    function.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: int psi_selectors


    File : :file:`selectors.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_selectors	(N_int,2,psi_selectors_size)
        double precision, allocatable	:: psi_selectors_coef	(psi_selectors_size,N_states)
        double precision, allocatable	:: psi_selectors_coef_tc	(psi_selectors_size,2,N_states)


    Determinants on which we apply <i|H|psi> for perturbation.
    psi_selectors_coef_tc(iii,1,istate) = left  coefficient of the iii determinant
    psi_selectors_coef_tc(iii,2,istate) = right coefficient of the iii determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted_gen`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_selectors_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors_coef_transp`

 
.. c:var:: double psi_selectors_coef


    File : :file:`selectors.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_selectors	(N_int,2,psi_selectors_size)
        double precision, allocatable	:: psi_selectors_coef	(psi_selectors_size,N_states)
        double precision, allocatable	:: psi_selectors_coef_tc	(psi_selectors_size,2,N_states)


    Determinants on which we apply <i|H|psi> for perturbation.
    psi_selectors_coef_tc(iii,1,istate) = left  coefficient of the iii determinant
    psi_selectors_coef_tc(iii,2,istate) = right coefficient of the iii determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted_gen`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_selectors_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors_coef_transp`

 
.. c:var:: double psi_selectors_coef_tc


    File : :file:`selectors.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_selectors	(N_int,2,psi_selectors_size)
        double precision, allocatable	:: psi_selectors_coef	(psi_selectors_size,N_states)
        double precision, allocatable	:: psi_selectors_coef_tc	(psi_selectors_size,2,N_states)


    Determinants on which we apply <i|H|psi> for perturbation.
    psi_selectors_coef_tc(iii,1,istate) = left  coefficient of the iii determinant
    psi_selectors_coef_tc(iii,2,istate) = right coefficient of the iii determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted_gen`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_selectors_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors_coef_transp`

 
.. c:var:: double psi_selectors_coef_transp


    File : :file:`selectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_selectors_coef_transp	(N_states,psi_selectors_size)
        double precision, allocatable	:: psi_selectors_coef_transp_tc	(N_states,2,psi_selectors_size)


    Transposed psi_selectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_states`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_size`


 
.. c:var:: double psi_selectors_coef_transp_tc


    File : :file:`selectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_selectors_coef_transp	(N_states,psi_selectors_size)
        double precision, allocatable	:: psi_selectors_coef_transp_tc	(N_states,2,psi_selectors_size)


    Transposed psi_selectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_states`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_size`


 
.. c:var:: int psi_selectors_size


    File : :file:`selectors.irp.f`

    .. code:: fortran

        integer	:: psi_selectors_size	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_coef_transp`

 
.. c:var:: double threshold_selectors


    File : :file:`selectors.irp.f`

    .. code:: fortran

        double precision	:: threshold_selectors	


    Thresholds on selectors (fraction of the square of the norm)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`threshold_generators`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void copy_h_apply_buffer_to_wf_tc(void)


    File : :file:`copy_wf.irp.f`

    .. code:: fortran

        subroutine copy_H_apply_buffer_to_wf_tc


    Copies the H_apply buffer to psi_coef.
    After calling this subroutine, N_det, psi_det and psi_coef need to be touched

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nproc`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_stochastic_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bi_normalize`
       * :c:func:`remove_duplicates_in_psi_det_tc`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`c0_weight`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`

 
.. c:function:: void diagonalize_ci_tc_bi_ortho(void)


    File : :file:`diagonalize_ci.irp.f`

    .. code:: fortran

        subroutine diagonalize_CI_tc_bi_ortho(ndet, E_tc, norm )


    Replace the coefficients of the CI states by the coefficients of the
    eigenstates of the CI matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_coef`
       * :c:data:`psi_energy`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`psi_energy`
       * :c:data:`eigval_right_tc_bi_orth`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`
       * :c:func:`tc_pt2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_tc_bi_ortho_wavefunction`

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
       * :c:data:`eigval_right_tc_bi_orth`

 
.. c:function:: void remove_duplicates_in_psi_det_tc(void)


    File : :file:`copy_wf.irp.f`

    .. code:: fortran

        subroutine remove_duplicates_in_psi_det_tc(found_duplicates)


    Removes duplicate determinants in the wave function.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`c0_weight`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`copy_h_apply_buffer_to_wf_tc`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`c0_weight`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`

 
.. c:function:: void run_cipsi_tc(void)


    File : :file:`fci_tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine run_cipsi_tc()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`do_pt2`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`json_int_fmt`
       * :c:data:`json_unit`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`fci_tc_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`json_close`
       * :c:func:`provide_all_three_ints_bi_ortho`
       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

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

 
.. c:function:: void run_pt2_tc(void)


    File : :file:`pt2_tc.irp.f`

    .. code:: fortran

        subroutine run_pt2_tc()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_tc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`provide_all_three_ints_bi_ortho`
       * :c:func:`tc_pt2`

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

 
.. c:function:: void save_energy(void)


    File : :file:`save_energy.irp.f`

    .. code:: fortran

        subroutine save_energy(E,pt2)


    Saves the energy in |EZFIO|.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_stochastic_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_fci_tc_bi_energy`
       * :c:func:`ezfio_set_fci_tc_bi_energy_pt2`

