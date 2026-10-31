.. _module_fci_tc_bi: 
 
.. program:: fci_tc_bi 
 
.. default-role:: option 
 
fci_tc_bi
=========

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Providers 
--------- 
 
.. c:var:: do_ddci


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
.. c:var:: do_only_1h1p


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
.. c:var:: do_only_cas


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
.. c:var:: n_det_generators


    File : :file:`generators.irp.f`

    .. code:: fortran

        integer	:: n_det_generators	


    For Single reference wave functions, the number of generators is 1 : the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`output_wall_time_0`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`threshold_generators`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`global_selection_buffer`
       * :c:data:`pt2_f`
       * :c:data:`pt2_j`
       * :c:data:`pt2_n_tasks`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_u`
       * :c:data:`pt2_w`

 
.. c:var:: n_det_selectors


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

 
.. c:var:: psi_coef_generators


    File : :file:`generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_tc`


 
.. c:var:: psi_coef_sorted_tc_gen


    File : :file:`generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_tc_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: psi_det_generators


    File : :file:`generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_tc`


 
.. c:var:: psi_det_sorted_bit_tc


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


 
.. c:var:: psi_det_sorted_tc_gen


    File : :file:`generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_tc_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: psi_det_sorted_tc_gen_order


    File : :file:`generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_tc_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: psi_l_coef_sorted_bit


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


 
.. c:var:: psi_r_coef_sorted_bit


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


 
.. c:var:: psi_selectors


    File : :file:`selectors.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_selectors	(N_int,2,psi_selectors_size)
        double precision, allocatable	:: psi_selectors_coef	(psi_selectors_size,N_states)
        double precision, allocatable	:: psi_selectors_coef_tc	(psi_selectors_size,2,N_states)


    Determinants on which we apply <i|H|psi> for perturbation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_det_sorted_tc_gen`
       * :c:data:`psi_selectors_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors_coef_transp`

 
.. c:var:: psi_selectors_coef


    File : :file:`selectors.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_selectors	(N_int,2,psi_selectors_size)
        double precision, allocatable	:: psi_selectors_coef	(psi_selectors_size,N_states)
        double precision, allocatable	:: psi_selectors_coef_tc	(psi_selectors_size,2,N_states)


    Determinants on which we apply <i|H|psi> for perturbation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_det_sorted_tc_gen`
       * :c:data:`psi_selectors_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors_coef_transp`

 
.. c:var:: psi_selectors_coef_tc


    File : :file:`selectors.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_selectors	(N_int,2,psi_selectors_size)
        double precision, allocatable	:: psi_selectors_coef	(psi_selectors_size,N_states)
        double precision, allocatable	:: psi_selectors_coef_tc	(psi_selectors_size,2,N_states)


    Determinants on which we apply <i|H|psi> for perturbation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_det_sorted_tc_gen`
       * :c:data:`psi_selectors_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors_coef_transp`

 
.. c:var:: psi_selectors_coef_transp


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


 
.. c:var:: psi_selectors_coef_transp_tc


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


 
.. c:var:: psi_selectors_size


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

 
.. c:var:: threshold_selectors


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
 
.. c:function:: copy_h_apply_buffer_to_wf_tc:


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

 
.. c:function:: diagonalize_ci_tc_bi_ortho:


    File : :file:`diagonalize_ci.irp.f`

    .. code:: fortran

        subroutine diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)


    Replace the coefficients of the CI states by the coefficients of the
    eigenstates of the CI matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`psi_coef`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`eigval_right_tc_bi_orth`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_tc_bi_ortho`
       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

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
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`eigval_right_tc_bi_orth`

 
.. c:function:: fci:


    File : :file:`fci_tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine fci


    Selected Full Configuration Interaction with stochastic selection
    and PT2.
    
    This program performs a |CIPSI|-like selected |CI| using a
    stochastic scheme for both the selection of the important Slater
    determinants and the computation of the |PT2| correction. This
    |CIPSI|-like algorithm will be performed for the lowest states of
    the variational space (see :option:`determinants n_states`). The
    |FCI| program will stop when reaching at least one the two following
    conditions:
    
    * number of Slater determinants > :option:`determinants n_det_max`
    * abs(|PT2|) less than :option:`perturbation pt2_max`
    
    The following other options can be of interest:
    
    :option:`determinants read_wf`
      When set to |false|, the program starts with a ROHF-like Slater
      determinant as a guess wave function. When set to |true|, the
      program starts with the wave function(s) stored in the |EZFIO|
      directory as guess wave function(s).
    
    :option:`determinants s2_eig`
      When set to |true|, the selection will systematically add all the
      necessary Slater determinants in order to have a pure spin wave
      function with an |S^2| value corresponding to
      :option:`determinants expected_s2`.
    
    For excited states calculations, it is recommended to start with
    :ref:`cis` or :ref:`cisd` guess wave functions, eventually in
    a restricted set of |MOs|, and to set :option:`determinants s2_eig`
    to |true|.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`pruning`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi_tc`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`psi_configuration`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`c0_weight`
       * :c:data:`distributed_davidson`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`pruning`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`read_wf`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`threshold_generators`
       * :c:data:`variance_match_weight`

 
.. c:function:: print_ci_dressed:


    File : :file:`diagonalize_ci.irp.f`

    .. code:: fortran

        subroutine print_CI_dressed(ndet, E_tc,norm,pt2_data,print_pt2)


    Replace the coefficients of the CI states by the coefficients of the
    eigenstates of the CI matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`psi_coef`
       * :c:data:`eigval_right_tc_bi_orth`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`psi_coef`
       * :c:data:`eigval_right_tc_bi_orth`

 
.. c:function:: remove_duplicates_in_psi_det_tc:


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

 
.. c:function:: run_cipsi_tc:


    File : :file:`fci_tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine run_cipsi_tc



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`do_pt2`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`is_zmq_slave`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`pt2_min_parallel_tasks`
       * :c:data:`three_body_h_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`fci`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`provide_all_three_ints_bi_ortho`
       * :c:func:`run_cipsi`
       * :c:func:`run_slave_cipsi`
       * :c:func:`run_stochastic_cipsi`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`psi_configuration`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`c0_weight`
       * :c:data:`distributed_davidson`
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
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`psi_det_sorted_bit_tc`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`read_wf`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`threshold_generators`
       * :c:data:`variance_match_weight`

 
.. c:function:: save_energy:


    File : :file:`save_energy.irp.f`

    .. code:: fortran

        subroutine save_energy(E,pt2)


    Saves the energy in |EZFIO|.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_fci_tc_energy`
       * :c:func:`ezfio_set_fci_tc_energy_pt2`

 
.. c:function:: zmq_get_n_det_generators:


    File : :file:`zmq.irp.f_template_102`

    .. code:: fortran

        integer function zmq_get_N_det_generators(zmq_to_qp_run_socket, worker_id)


    Get N_det_generators from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det_generators`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_get_n_det_selectors:


    File : :file:`zmq.irp.f_template_102`

    .. code:: fortran

        integer function zmq_get_N_det_selectors(zmq_to_qp_run_socket, worker_id)


    Get N_det_selectors from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det_selectors`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_n_det_generators:


    File : :file:`zmq.irp.f_template_102`

    .. code:: fortran

        integer function zmq_put_N_det_generators(zmq_to_qp_run_socket,worker_id)


    Put N_det_generators on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_n_det_selectors:


    File : :file:`zmq.irp.f_template_102`

    .. code:: fortran

        integer function zmq_put_N_det_selectors(zmq_to_qp_run_socket,worker_id)


    Put N_det_selectors on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`zmq_state`

