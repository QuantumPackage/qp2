.. _module_two_rdm_routines: 
 
.. program:: two_rdm_routines 
 
.. default-role:: option 
 
two_rdm_routines
================

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void is_integer_in_string_local(void)


    File : :file:`two_rdm_routines/update_rdm.irp.f`

    .. code:: fortran

        logical function is_integer_in_string_local(orb,bitmask,Nint)



 
.. c:function:: void orb_range_2_rdm_openmp(void)


    File : :file:`two_rdm_routines/davidson_like_2rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_2_rdm_openmp(big_array,dim1,norb,list_orb,ispin,u_0,N_st,sze)


    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta + beta/alpha  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + ab + ba
    
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep are used in ZMQ parallelization.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`n_det`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dset_order`
       * :c:func:`dtranspose`
       * :c:func:`orb_range_2_rdm_openmp_work`

 
.. c:function:: void orb_range_2_rdm_openmp_work(void)


    File : :file:`two_rdm_routines/davidson_like_2rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_2_rdm_openmp_work(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes two-rdm
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`

 
.. c:function:: void orb_range_2_rdm_openmp_work_1(void)


    File : :file:`two_rdm_routines/davidson_like_2rdm.irp.f_template_572`

    .. code:: fortran

        subroutine orb_range_2_rdm_openmp_work_1(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_dm_buffer`
       * :c:func:`update_keys_values_n_states_local`

 
.. c:function:: void orb_range_2_rdm_openmp_work_2(void)


    File : :file:`two_rdm_routines/davidson_like_2rdm.irp.f_template_572`

    .. code:: fortran

        subroutine orb_range_2_rdm_openmp_work_2(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_dm_buffer`
       * :c:func:`update_keys_values_n_states_local`

 
.. c:function:: void orb_range_2_rdm_openmp_work_3(void)


    File : :file:`two_rdm_routines/davidson_like_2rdm.irp.f_template_572`

    .. code:: fortran

        subroutine orb_range_2_rdm_openmp_work_3(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_dm_buffer`
       * :c:func:`update_keys_values_n_states_local`

 
.. c:function:: void orb_range_2_rdm_openmp_work_4(void)


    File : :file:`two_rdm_routines/davidson_like_2rdm.irp.f_template_572`

    .. code:: fortran

        subroutine orb_range_2_rdm_openmp_work_4(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_dm_buffer`
       * :c:func:`update_keys_values_n_states_local`

 
.. c:function:: void orb_range_2_rdm_openmp_work_n_int(void)


    File : :file:`two_rdm_routines/davidson_like_2rdm.irp.f_template_572`

    .. code:: fortran

        subroutine orb_range_2_rdm_openmp_work_N_int(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_dm_buffer`
       * :c:func:`update_keys_values_n_states_local`

 
.. c:function:: void orb_range_2_rdm_state_av_openmp(void)


    File : :file:`two_rdm_routines/davidson_like_state_av_2rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_2_rdm_state_av_openmp(big_array,dim1,norb,list_orb,state_weights,ispin,u_0,N_st,sze)


    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep are used in ZMQ parallelization.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dset_order`
       * :c:func:`dtranspose`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work`

 
.. c:function:: void orb_range_2_rdm_state_av_openmp_work(void)


    File : :file:`two_rdm_routines/davidson_like_state_av_2rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_2_rdm_state_av_openmp_work(big_array,dim1,norb,list_orb,state_weights,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes two-rdm
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

 
.. c:function:: void orb_range_2_rdm_state_av_openmp_work_1(void)


    File : :file:`two_rdm_routines/davidson_like_state_av_2rdm.irp.f_template_557`

    .. code:: fortran

        subroutine orb_range_2_rdm_state_av_openmp_work_1(big_array,dim1,norb,list_orb,state_weights,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
    In any cases, the state average weights will be used with an array state_weights
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_bb_dm_buffer`
       * :c:func:`update_keys_values`

 
.. c:function:: void orb_range_2_rdm_state_av_openmp_work_2(void)


    File : :file:`two_rdm_routines/davidson_like_state_av_2rdm.irp.f_template_557`

    .. code:: fortran

        subroutine orb_range_2_rdm_state_av_openmp_work_2(big_array,dim1,norb,list_orb,state_weights,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
    In any cases, the state average weights will be used with an array state_weights
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_bb_dm_buffer`
       * :c:func:`update_keys_values`

 
.. c:function:: void orb_range_2_rdm_state_av_openmp_work_3(void)


    File : :file:`two_rdm_routines/davidson_like_state_av_2rdm.irp.f_template_557`

    .. code:: fortran

        subroutine orb_range_2_rdm_state_av_openmp_work_3(big_array,dim1,norb,list_orb,state_weights,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
    In any cases, the state average weights will be used with an array state_weights
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_bb_dm_buffer`
       * :c:func:`update_keys_values`

 
.. c:function:: void orb_range_2_rdm_state_av_openmp_work_4(void)


    File : :file:`two_rdm_routines/davidson_like_state_av_2rdm.irp.f_template_557`

    .. code:: fortran

        subroutine orb_range_2_rdm_state_av_openmp_work_4(big_array,dim1,norb,list_orb,state_weights,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
    In any cases, the state average weights will be used with an array state_weights
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_bb_dm_buffer`
       * :c:func:`update_keys_values`

 
.. c:function:: void orb_range_2_rdm_state_av_openmp_work_n_int(void)


    File : :file:`two_rdm_routines/davidson_like_state_av_2rdm.irp.f_template_557`

    .. code:: fortran

        subroutine orb_range_2_rdm_state_av_openmp_work_N_int(big_array,dim1,norb,list_orb,state_weights,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2rdm
             == 2 :: beta /beta  2rdm
             == 3 :: alpha/beta  2rdm
             == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
    The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
    In any cases, the state average weights will be used with an array state_weights
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_bb_dm_buffer`
       * :c:func:`update_keys_values`

 
.. c:function:: void orb_range_2_trans_rdm_openmp(void)


    File : :file:`two_rdm_routines/davidson_like_trans_2rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_2_trans_rdm_openmp(big_array,dim1,norb,list_orb,ispin,u_0,N_st,sze)


    if ispin == 1 :: alpha/alpha 2_rdm
             == 2 :: beta /beta  2_rdm
             == 3 :: alpha/beta + beta/alpha  2trans_rdm
             == 4 :: spin traced 2_rdm :: aa + bb + ab + ba
    
    notice that here it is the TRANSITION RDM THAT IS COMPUTED
    
    THE DIAGONAL PART IS THE USUAL ONE FOR A GIVEN STATE
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep are used in ZMQ parallelization.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`n_det`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_trans_spin_trace_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dset_order`
       * :c:func:`dtranspose`
       * :c:func:`orb_range_2_trans_rdm_openmp_work`

 
.. c:function:: void orb_range_2_trans_rdm_openmp_work(void)


    File : :file:`two_rdm_routines/davidson_like_trans_2rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_2_trans_rdm_openmp_work(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes two-trans_rdm
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work_1`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_2`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_3`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_4`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_n_int`

 
.. c:function:: void orb_range_2_trans_rdm_openmp_work_1(void)


    File : :file:`two_rdm_routines/davidson_like_trans_2rdm.irp.f_template_552`

    .. code:: fortran

        subroutine orb_range_2_trans_rdm_openmp_work_1(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two trans_rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2trans_rdm
             == 2 :: beta /beta  2trans_rdm
             == 3 :: alpha/beta  2trans_rdm
             == 4 :: spin traced 2trans_rdm :: aa + bb + 0.5 (ab + ba))
    The 2trans_rdm will be computed only on the list of orbitals list_orb, which contains norb
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_trans_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_trans_rdm_bb_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_trans_rdm_buffer`
       * :c:func:`update_keys_values_n_states_trans`

 
.. c:function:: void orb_range_2_trans_rdm_openmp_work_2(void)


    File : :file:`two_rdm_routines/davidson_like_trans_2rdm.irp.f_template_552`

    .. code:: fortran

        subroutine orb_range_2_trans_rdm_openmp_work_2(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two trans_rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2trans_rdm
             == 2 :: beta /beta  2trans_rdm
             == 3 :: alpha/beta  2trans_rdm
             == 4 :: spin traced 2trans_rdm :: aa + bb + 0.5 (ab + ba))
    The 2trans_rdm will be computed only on the list of orbitals list_orb, which contains norb
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_trans_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_trans_rdm_bb_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_trans_rdm_buffer`
       * :c:func:`update_keys_values_n_states_trans`

 
.. c:function:: void orb_range_2_trans_rdm_openmp_work_3(void)


    File : :file:`two_rdm_routines/davidson_like_trans_2rdm.irp.f_template_552`

    .. code:: fortran

        subroutine orb_range_2_trans_rdm_openmp_work_3(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two trans_rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2trans_rdm
             == 2 :: beta /beta  2trans_rdm
             == 3 :: alpha/beta  2trans_rdm
             == 4 :: spin traced 2trans_rdm :: aa + bb + 0.5 (ab + ba))
    The 2trans_rdm will be computed only on the list of orbitals list_orb, which contains norb
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_trans_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_trans_rdm_bb_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_trans_rdm_buffer`
       * :c:func:`update_keys_values_n_states_trans`

 
.. c:function:: void orb_range_2_trans_rdm_openmp_work_4(void)


    File : :file:`two_rdm_routines/davidson_like_trans_2rdm.irp.f_template_552`

    .. code:: fortran

        subroutine orb_range_2_trans_rdm_openmp_work_4(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two trans_rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2trans_rdm
             == 2 :: beta /beta  2trans_rdm
             == 3 :: alpha/beta  2trans_rdm
             == 4 :: spin traced 2trans_rdm :: aa + bb + 0.5 (ab + ba))
    The 2trans_rdm will be computed only on the list of orbitals list_orb, which contains norb
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_trans_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_trans_rdm_bb_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_trans_rdm_buffer`
       * :c:func:`update_keys_values_n_states_trans`

 
.. c:function:: void orb_range_2_trans_rdm_openmp_work_n_int(void)


    File : :file:`two_rdm_routines/davidson_like_trans_2rdm.irp.f_template_552`

    .. code:: fortran

        subroutine orb_range_2_trans_rdm_openmp_work_N_int(big_array,dim1,norb,list_orb,ispin,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes the two trans_rdm for the N_st vectors |u_t>
    if ispin == 1 :: alpha/alpha 2trans_rdm
             == 2 :: beta /beta  2trans_rdm
             == 3 :: alpha/beta  2trans_rdm
             == 4 :: spin traced 2trans_rdm :: aa + bb + 0.5 (ab + ba))
    The 2trans_rdm will be computed only on the list of orbitals list_orb, which contains norb
    Default should be 1,N_det,0,1 for istart,iend,ishift,istep

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`list_to_bitstring`
       * :c:func:`omp_init_lock`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_trans_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_trans_rdm_bb_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_trans_rdm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_trans_rdm_buffer`
       * :c:func:`update_keys_values_n_states_trans`

 
.. c:function:: void orb_range_diag_to_all_2_rdm_dm_buffer(void)


    File : :file:`two_rdm_routines/update_state_av_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_diag_to_all_2_rdm_dm_buffer(det_1,c_1,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the DIAGONAL PART of the two body rdms in a specific range of orbitals for a given determinant det_1
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-rdm

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: void orb_range_diag_to_all_states_2_rdm_dm_buffer(void)


    File : :file:`two_rdm_routines/update_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_diag_to_all_states_2_rdm_dm_buffer(det_1,c_1,N_st,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the DIAGONAL PART of the two body rdms in a specific range of orbitals for a given determinant det_1
    
    c_1 is the array of the contributions to the rdm for all states
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1                :: alpha/ alpha
    ispin   == 2                :: beta / beta
    ispin   == 3                :: alpha/ beta
    ispin   == 4                :: spin traced <=> total two-rdm

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: void orb_range_diag_to_all_states_2_rdm_trans_buffer(void)


    File : :file:`two_rdm_routines/update_trans_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_diag_to_all_states_2_rdm_trans_buffer(det_1,c_1,N_st,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the DIAGONAL PART of the two body trans_rdms in a specific range of orbitals for a given determinant det_1
    
    c_1 is the array of the contributions to the trans_rdm for all states
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-trans_rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-trans_rdm

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work_1`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_2`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_3`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_4`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: void orb_range_off_diag_double_to_2_rdm_aa_dm_buffer(void)


    File : :file:`two_rdm_routines/update_state_av_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_double_to_2_rdm_aa_dm_buffer(det_1,det_2,c_1,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a ALPHA/ALPHA DOUBLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-rdm
    
    here, only ispin == 1 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation_spin`

 
.. c:function:: void orb_range_off_diag_double_to_2_rdm_ab_dm_buffer(void)


    File : :file:`two_rdm_routines/update_state_av_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_double_to_2_rdm_ab_dm_buffer(det_1,det_2,c_1,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a alpha/beta DOUBLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-rdm
    
    here, only ispin == 3 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation`

 
.. c:function:: void orb_range_off_diag_double_to_2_rdm_bb_dm_buffer(void)


    File : :file:`two_rdm_routines/update_state_av_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_double_to_2_rdm_bb_dm_buffer(det_1,det_2,c_1,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a BETA /BETA  DOUBLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-rdm
    
    here, only ispin == 2 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation_spin`

 
.. c:function:: void orb_range_off_diag_double_to_all_states_aa_dm_buffer(void)


    File : :file:`two_rdm_routines/update_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_double_to_all_states_aa_dm_buffer(det_1,det_2,c_1,N_st,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a ALPHA/ALPHA DOUBLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1                 :: alpha/ alpha
    ispin   == 2                 :: beta / beta
    ispin   == 3                 :: alpha/ beta
    ispin   == 4                 :: spin traced <=> total two-rdm
    
    here, only ispin == 1 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation_spin`

 
.. c:function:: void orb_range_off_diag_double_to_all_states_aa_trans_rdm_buffer(void)


    File : :file:`two_rdm_routines/update_trans_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_double_to_all_states_aa_trans_rdm_buffer(det_1,det_2,c_1,N_st,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body trans_rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a ALPHA/ALPHA DOUBLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body trans_rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-trans_rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-trans_rdm
    
    here, only ispin == 1 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work_1`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_2`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_3`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_4`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation_spin`

 
.. c:function:: void orb_range_off_diag_double_to_all_states_ab_dm_buffer(void)


    File : :file:`two_rdm_routines/update_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_double_to_all_states_ab_dm_buffer(det_1,det_2,c_1,N_st,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that updates the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a alpha/beta DOUBLE excitation with respect to one another
    
     c_1 is the array of the contributions to the rdm for all states
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1                 :: alpha/ alpha
    ispin   == 2                 :: beta / beta
    ispin   == 3                 :: alpha/ beta
    ispin   == 4                 :: spin traced <=> total two-rdm
    
    here, only ispin == 3 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation`

 
.. c:function:: void orb_range_off_diag_double_to_all_states_ab_trans_rdm_buffer(void)


    File : :file:`two_rdm_routines/update_trans_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_double_to_all_states_ab_trans_rdm_buffer(det_1,det_2,c_1,N_st,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body trans_rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a alpha/beta DOUBLE excitation with respect to one another
    
     c_1 is the array of the contributions to the trans_rdm for all states
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-trans_rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-trans_rdm
    
    here, only ispin == 3 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work_1`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_2`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_3`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_4`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation`

 
.. c:function:: void orb_range_off_diag_double_to_all_states_bb_dm_buffer(void)


    File : :file:`two_rdm_routines/update_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_double_to_all_states_bb_dm_buffer(det_1,det_2,c_1,N_st,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a BETA /BETA  DOUBLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1                 :: alpha/ alpha
    ispin   == 2                 :: beta / beta
    ispin   == 3                 :: alpha/ beta
    ispin   == 4                 :: spin traced <=> total two-rdm
    
    here, only ispin == 2 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation_spin`

 
.. c:function:: void orb_range_off_diag_double_to_all_states_trans_rdm_bb_buffer(void)


    File : :file:`two_rdm_routines/update_trans_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_double_to_all_states_trans_rdm_bb_buffer(det_1,det_2,c_1,N_st,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body trans_rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a BETA /BETA  DOUBLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body trans_rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-trans_rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-trans_rdm
    
    here, only ispin == 2 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work_1`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_2`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_3`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_4`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation_spin`

 
.. c:function:: void orb_range_off_diag_single_to_2_rdm_aa_dm_buffer(void)


    File : :file:`two_rdm_routines/update_state_av_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_single_to_2_rdm_aa_dm_buffer(det_1,det_2,c_1,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a ALPHA SINGLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-rdm
    
    here, only ispin == 1 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_single_excitation`

 
.. c:function:: void orb_range_off_diag_single_to_2_rdm_ab_dm_buffer(void)


    File : :file:`two_rdm_routines/update_state_av_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_single_to_2_rdm_ab_dm_buffer(det_1,det_2,c_1,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a SINGLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-rdm
    
    here, only ispin == 3 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_single_excitation`

 
.. c:function:: void orb_range_off_diag_single_to_2_rdm_bb_dm_buffer(void)


    File : :file:`two_rdm_routines/update_state_av_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_single_to_2_rdm_bb_dm_buffer(det_1,det_2,c_1,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a BETA  SINGLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-rdm
    
    here, only ispin == 2 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_single_excitation`

 
.. c:function:: void orb_range_off_diag_single_to_all_states_aa_dm_buffer(void)


    File : :file:`two_rdm_routines/update_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_single_to_all_states_aa_dm_buffer(det_1,det_2,c_1,N_st,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a ALPHA SINGLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1                 :: alpha/ alpha
    ispin   == 2                 :: beta / beta
    ispin   == 3                 :: alpha/ beta
    ispin   == 4                 :: spin traced <=> total two-rdm
    
    here, only ispin == 1 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_single_excitation`

 
.. c:function:: void orb_range_off_diag_single_to_all_states_aa_trans_rdm_buffer(void)


    File : :file:`two_rdm_routines/update_trans_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_single_to_all_states_aa_trans_rdm_buffer(det_1,det_2,c_1,N_st,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body trans_rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a ALPHA SINGLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body trans_rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-trans_rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-trans_rdm
    
    here, only ispin == 1 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work_1`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_2`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_3`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_4`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_single_excitation`

 
.. c:function:: void orb_range_off_diag_single_to_all_states_ab_dm_buffer(void)


    File : :file:`two_rdm_routines/update_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_single_to_all_states_ab_dm_buffer(det_1,det_2,c_1,N_st,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that updates the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a SINGLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1                 :: alpha/ alpha
    ispin   == 2                 :: beta / beta
    ispin   == 3                 :: alpha/ beta
    ispin   == 4                 :: spin traced <=> total two-rdm
    
    here, only ispin == 3 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_single_excitation`

 
.. c:function:: void orb_range_off_diag_single_to_all_states_ab_trans_rdm_buffer(void)


    File : :file:`two_rdm_routines/update_trans_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_single_to_all_states_ab_trans_rdm_buffer(det_1,det_2,c_1,N_st,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body trans_rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a SINGLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body trans_rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-trans_rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-trans_rdm
    
    here, only ispin == 3 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work_1`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_2`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_3`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_4`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_single_excitation`

 
.. c:function:: void orb_range_off_diag_single_to_all_states_bb_dm_buffer(void)


    File : :file:`two_rdm_routines/update_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_single_to_all_states_bb_dm_buffer(det_1,det_2,c_1,N_st,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a BETA  SINGLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-rdm you will update
    
    ispin   == 1                 :: alpha/ alpha
    ispin   == 2                 :: beta / beta
    ispin   == 3                 :: alpha/ beta
    ispin   == 4                 :: spin traced <=> total two-rdm
    
    here, only ispin == 2 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_single_excitation`

 
.. c:function:: void orb_range_off_diag_single_to_all_states_bb_trans_rdm_buffer(void)


    File : :file:`two_rdm_routines/update_trans_rdm.irp.f`

    .. code:: fortran

        subroutine orb_range_off_diag_single_to_all_states_bb_trans_rdm_buffer(det_1,det_2,c_1,N_st,orb_bitmask,list_orb_reverse,ispin,sze_buff,nkeys,keys,values)


    routine that update the OFF DIAGONAL PART of the two body trans_rdms in a specific range of orbitals for
    
    a given couple of determinant det_1, det_2 being a BETA  SINGLE excitation with respect to one another
    
    c_1 is supposed to be a scalar quantity, such as state averaged coef of the determinant det_1
    
    big_array(dim1,dim1,dim1,dim1) is the two-body trans_rdm to be updated in physicist notation
    
    orb_bitmask(N_int) is the bitmask for the orbital range, list_orb_reverse(mo_num) is the inverse range of orbitals
    
    ispin determines which spin-spin component of the two-trans_rdm you will update
    
    ispin   == 1 :: alpha/ alpha
    ispin   == 2 :: beta / beta
    ispin   == 3 :: alpha/ beta
    ispin   == 4 :: spin traced <=> total two-trans_rdm
    
    here, only ispin == 2 or 4 will do something

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work_1`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_2`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_3`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_4`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_single_excitation`

 
.. c:function:: void update_keys_values(void)


    File : :file:`two_rdm_routines/davidson_like_state_av_2rdm.irp.f`

    .. code:: fortran

        subroutine update_keys_values(keys,values,nkeys,dim1,big_array,lock_2rdm)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:function:: void update_keys_values_n_states(void)


    File : :file:`two_rdm_routines/davidson_like_2rdm.irp.f`

    .. code:: fortran

        subroutine update_keys_values_n_states(keys,values,nkeys,dim1,n_st,big_array,lock_2rdm)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i8sort`
       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:function:: void update_keys_values_n_states_local(void)


    File : :file:`two_rdm_routines/davidson_like_2rdm.irp.f`

    .. code:: fortran

        subroutine update_keys_values_n_states_local(keys,values,nkeys,dim1,n_st,big_array_local)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`

 
.. c:function:: void update_keys_values_n_states_trans(void)


    File : :file:`two_rdm_routines/davidson_like_trans_2rdm.irp.f`

    .. code:: fortran

        subroutine update_keys_values_n_states_trans(keys,values,nkeys,dim1,n_st,big_array,lock_2rdm)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orb_range_2_trans_rdm_openmp_work_1`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_2`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_3`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_4`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

