.. _module_tc_scf: 
 
.. program:: tc_scf 
 
.. default-role:: option 
 
tc_scf
======

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Providers 
--------- 
 
.. c:var:: diag_three_elem_hf


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision	:: diag_three_elem_hf	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bi_ortho`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`fock_3_rho_beta`
       * :c:data:`fock_3_w_kk_sum`
       * :c:data:`fock_3_w_kl_mo_k_mo_l`
       * :c:data:`fock_3_w_kl_w_kl`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`three_body_h_tc`
       * :c:data:`x_w_ij_erf_rk`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_hf_energy`
       * :c:data:`vartc_hf_energy`

 
.. c:var:: eigval_fock_tc_mo


    File : :file:`diago_bi_ort_tcfock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_tc_reigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: fock_tc_leigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: eigval_fock_tc_mo	(mo_num)
        double precision, allocatable	:: overlap_fock_tc_eigvec_mo	(mo_num,mo_num)


    EIGENVECTORS OF FOCK MATRIX ON THE MO BASIS and their OVERLAP

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`im_thresh_tcscf`
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_num`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_tc_reigvec_ao`

 
.. c:var:: fock_3_mat


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_mat	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bi_ortho`
       * :c:data:`final_grid_points`
       * :c:data:`fock_3_rho_beta`
       * :c:data:`fock_3_trace_w_tilde`
       * :c:data:`fock_3_w_ki_mos_k`
       * :c:data:`fock_3_w_ki_wk_a`
       * :c:data:`fock_3_w_kk_sum`
       * :c:data:`fock_3_w_kl_mo_k_mo_l`
       * :c:data:`fock_3_w_kl_wla_phi_k`
       * :c:data:`mo_num`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`three_body_h_tc`
       * :c:data:`x_w_ij_erf_rk`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`

 
.. c:var:: fock_3_mat_a_op_sh


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_mat_a_op_sh	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ij_erf_rk`


 
.. c:var:: fock_3_mat_b_op_sh


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_mat_b_op_sh	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ij_erf_rk`


 
.. c:var:: fock_3_rho_beta


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_rho_beta	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mos_in_r_array`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diag_three_elem_hf`
       * :c:data:`fock_3_mat`

 
.. c:var:: fock_3_trace_w_tilde


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_trace_w_tilde	(n_points_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`fock_3_w_ki_wk_a`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3_mat`

 
.. c:var:: fock_3_w_ki_mos_k


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_w_ki_mos_k	(n_points_final_grid,3,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mos_in_r_array`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ij_erf_rk`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3_mat`

 
.. c:var:: fock_3_w_ki_wk_a


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_w_ki_wk_a	(n_points_final_grid,3,mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ij_erf_rk`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3_mat`
       * :c:data:`fock_3_trace_w_tilde`

 
.. c:var:: fock_3_w_kk_sum


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_w_kk_sum	(n_points_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ij_erf_rk`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diag_three_elem_hf`
       * :c:data:`fock_3_mat`

 
.. c:var:: fock_3_w_kl_mo_k_mo_l


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_w_kl_mo_k_mo_l	(n_points_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ij_erf_rk`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diag_three_elem_hf`
       * :c:data:`fock_3_mat`

 
.. c:var:: fock_3_w_kl_w_kl


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_w_kl_w_kl	(n_points_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ij_erf_rk`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diag_three_elem_hf`

 
.. c:var:: fock_3_w_kl_wla_phi_k


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3_w_kl_wla_phi_k	(n_points_final_grid,3,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ij_erf_rk`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3_mat`

 
.. c:var:: fock_3e_uhf_ao_a


    File : :file:`fock_3e_bi_ortho_uhf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3e_uhf_ao_a	(ao_num,ao_num)


    
    Equations (B6) and (B7)
    
    g <--> gamma
    d <--> delta
    e <--> eta
    k <--> kappa
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao_t`
       * :c:data:`n_points_final_grid`
       * :c:data:`tcscf_bi_ort_dm_ao_alpha`
       * :c:data:`tcscf_bi_ort_dm_ao_beta`


 
.. c:var:: fock_3e_uhf_ao_b


    File : :file:`fock_3e_bi_ortho_uhf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3e_uhf_ao_b	(ao_num,ao_num)


    
    Equations (B6) and (B7)
    
    g <--> gamma
    d <--> delta
    e <--> eta
    k <--> kappa
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao_t`
       * :c:data:`n_points_final_grid`
       * :c:data:`tcscf_bi_ort_dm_ao_alpha`
       * :c:data:`tcscf_bi_ort_dm_ao_beta`


 
.. c:var:: fock_3e_uhf_mo_a


    File : :file:`fock_3e_bi_ortho_uhf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3e_uhf_mo_a	(mo_num,mo_num)


    ALPHA part of the Fock matrix from three-electron terms
    
    WARNING :: non hermitian if bi-ortho MOS used

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`fock_3e_uhf_mo_cs`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_alpha`
       * :c:data:`fock_matrix_vartc_mo_alpha`

 
.. c:var:: fock_3e_uhf_mo_b


    File : :file:`fock_3e_bi_ortho_uhf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3e_uhf_mo_b	(mo_num,mo_num)


    BETA part of the Fock matrix from three-electron terms
    
    WARNING :: non hermitian if bi-ortho MOS used

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`fock_3e_uhf_mo_cs`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_beta`
       * :c:data:`fock_matrix_vartc_mo_beta`

 
.. c:var:: fock_3e_uhf_mo_cs


    File : :file:`fock_3e_bi_ortho_uhf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3e_uhf_mo_cs	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3e_uhf_mo_a`
       * :c:data:`fock_3e_uhf_mo_b`

 
.. c:var:: fock_a_tmp1_bi_ortho


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_a_tmp1_bi_ortho	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_a_tot_3e_bi_orth`

 
.. c:var:: fock_a_tmp2_bi_ortho


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_a_tmp2_bi_ortho	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_a_tot_3e_bi_orth`

 
.. c:var:: fock_a_tot_3e_bi_orth


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_a_tot_3e_bi_orth	(mo_num,mo_num)


    Alpha part of the Fock matrix from three-electron terms
    
    WARNING :: non hermitian if bi-ortho MOS used

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_a_tmp1_bi_ortho`
       * :c:data:`fock_a_tmp2_bi_ortho`
       * :c:data:`fock_cs_3e_bi_orth`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`fock_matrix_vartc_mo_tot`

 
.. c:var:: fock_b_tmp1_bi_ortho


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_b_tmp1_bi_ortho	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_b_tot_3e_bi_orth`

 
.. c:var:: fock_b_tmp2_bi_ortho


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_b_tmp2_bi_ortho	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_b_tot_3e_bi_orth`

 
.. c:var:: fock_b_tot_3e_bi_orth


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_b_tot_3e_bi_orth	(mo_num,mo_num)


    Beta part of the Fock matrix from three-electron terms
    
    WARNING :: non hermitian if bi-ortho MOS used

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_b_tmp1_bi_ortho`
       * :c:data:`fock_b_tmp2_bi_ortho`
       * :c:data:`fock_cs_3e_bi_orth`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`fock_matrix_vartc_mo_tot`

 
.. c:var:: fock_cs_3e_bi_orth


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_cs_3e_bi_orth	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`

 
.. c:var:: fock_matrix_tc_ao_alpha


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_ao_alpha	(ao_num,ao_num)


    Total alpha TC Fock matrix : h_c + Two-e^TC terms on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`two_e_tc_non_hermit_integral_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_alpha`

 
.. c:var:: fock_matrix_tc_ao_beta


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_ao_beta	(ao_num,ao_num)


    Total beta TC Fock matrix : h_c + Two-e^TC terms on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`two_e_tc_non_hermit_integral_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_beta`

 
.. c:var:: fock_matrix_tc_ao_tot


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_ao_tot	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fqs_sqf_ao`

 
.. c:var:: fock_matrix_tc_diag_mo_tot


    File : :file:`fock_tc_mo_tot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_mo_tot	(mo_num,mo_num)
        double precision, allocatable	:: fock_matrix_tc_diag_mo_tot	(mo_num)


    TC-Fock matrix on the MO basis. WARNING !!! NON HERMITIAN !!!
    For open shells, the ROHF Fock Matrix is ::
    
          |   F-K    |  F + K/2  |    F     |
          |---------------------------------|
          | F + K/2  |     F     |  F - K/2 |
          |---------------------------------|
          |    F     |  F - K/2  |  F + K   |
    
    
    F = 1/2 (Fa + Fb)
    
    K = Fb - Fa
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bi_ortho`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_3_mat`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`fock_matrix_tc_mo_alpha`
       * :c:data:`fock_matrix_tc_mo_beta`
       * :c:data:`frozen_orb_scf`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`no_oa_or_av_opt`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`fock_tc_reigvec_mo`
       * :c:data:`good_hermit_tc_fock_mat`
       * :c:data:`grad_non_hermit_left`
       * :c:data:`hermit_average_tc_fock_mat`

 
.. c:var:: fock_matrix_tc_mo_alpha


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_mo_alpha	(mo_num,mo_num)


    Total alpha TC Fock matrix : h_c + Two-e^TC terms on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`bi_ortho`
       * :c:data:`fock_3e_uhf_mo_a`
       * :c:data:`fock_matrix_tc_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`

 
.. c:var:: fock_matrix_tc_mo_beta


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_mo_beta	(mo_num,mo_num)


    Total beta TC Fock matrix : h_c + Two-e^TC terms on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`bi_ortho`
       * :c:data:`fock_3e_uhf_mo_b`
       * :c:data:`fock_matrix_tc_ao_beta`
       * :c:data:`mo_coef`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`

 
.. c:var:: fock_matrix_tc_mo_tot


    File : :file:`fock_tc_mo_tot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_mo_tot	(mo_num,mo_num)
        double precision, allocatable	:: fock_matrix_tc_diag_mo_tot	(mo_num)


    TC-Fock matrix on the MO basis. WARNING !!! NON HERMITIAN !!!
    For open shells, the ROHF Fock Matrix is ::
    
          |   F-K    |  F + K/2  |    F     |
          |---------------------------------|
          | F + K/2  |     F     |  F - K/2 |
          |---------------------------------|
          |    F     |  F - K/2  |  F + K   |
    
    
    F = 1/2 (Fa + Fb)
    
    K = Fb - Fa
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bi_ortho`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_3_mat`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`fock_matrix_tc_mo_alpha`
       * :c:data:`fock_matrix_tc_mo_beta`
       * :c:data:`frozen_orb_scf`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`no_oa_or_av_opt`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`fock_tc_reigvec_mo`
       * :c:data:`good_hermit_tc_fock_mat`
       * :c:data:`grad_non_hermit_left`
       * :c:data:`hermit_average_tc_fock_mat`

 
.. c:var:: fock_matrix_vartc_ao_alpha


    File : :file:`fock_vartc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_vartc_ao_alpha	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`two_e_vartc_integral_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_vartc_mo_alpha`

 
.. c:var:: fock_matrix_vartc_ao_beta


    File : :file:`fock_vartc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_vartc_ao_beta	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`two_e_vartc_integral_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_vartc_mo_beta`

 
.. c:var:: fock_matrix_vartc_ao_tot


    File : :file:`fock_vartc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_vartc_ao_tot	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_matrix_vartc_mo_tot`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fqs_sqf_ao`

 
.. c:var:: fock_matrix_vartc_diag_mo_tot


    File : :file:`fock_vartc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_vartc_mo_tot	(mo_num,mo_num)
        double precision, allocatable	:: fock_matrix_vartc_diag_mo_tot	(mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`fock_matrix_vartc_mo_alpha`
       * :c:data:`fock_matrix_vartc_mo_beta`
       * :c:data:`frozen_orb_scf`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`no_oa_or_av_opt`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_vartc_ao_tot`
       * :c:data:`fock_vartc_eigvec_mo`
       * :c:data:`grad_vartc`

 
.. c:var:: fock_matrix_vartc_mo_alpha


    File : :file:`fock_vartc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_vartc_mo_alpha	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_3e_uhf_mo_a`
       * :c:data:`fock_matrix_vartc_ao_alpha`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_vartc_mo_tot`

 
.. c:var:: fock_matrix_vartc_mo_beta


    File : :file:`fock_vartc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_vartc_mo_beta	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_3e_uhf_mo_b`
       * :c:data:`fock_matrix_vartc_ao_beta`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_vartc_mo_tot`

 
.. c:var:: fock_matrix_vartc_mo_tot


    File : :file:`fock_vartc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_vartc_mo_tot	(mo_num,mo_num)
        double precision, allocatable	:: fock_matrix_vartc_diag_mo_tot	(mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`fock_matrix_vartc_mo_alpha`
       * :c:data:`fock_matrix_vartc_mo_beta`
       * :c:data:`frozen_orb_scf`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`no_oa_or_av_opt`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_vartc_ao_tot`
       * :c:data:`fock_vartc_eigvec_mo`
       * :c:data:`grad_vartc`

 
.. c:var:: fock_tc_leigvec_ao


    File : :file:`diago_bi_ort_tcfock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_tc_reigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: fock_tc_leigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: overlap_fock_tc_eigvec_ao	(mo_num,mo_num)


    EIGENVECTORS OF FOCK MATRIX ON THE AO BASIS and their OVERLAP
    
    THE OVERLAP SHOULD BE THE SAME AS overlap_fock_tc_eigvec_mo

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_tc_reigvec_mo`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: fock_tc_leigvec_mo


    File : :file:`diago_bi_ort_tcfock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_tc_reigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: fock_tc_leigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: eigval_fock_tc_mo	(mo_num)
        double precision, allocatable	:: overlap_fock_tc_eigvec_mo	(mo_num,mo_num)


    EIGENVECTORS OF FOCK MATRIX ON THE MO BASIS and their OVERLAP

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`im_thresh_tcscf`
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_num`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_tc_reigvec_ao`

 
.. c:var:: fock_tc_reigvec_ao


    File : :file:`diago_bi_ort_tcfock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_tc_reigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: fock_tc_leigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: overlap_fock_tc_eigvec_ao	(mo_num,mo_num)


    EIGENVECTORS OF FOCK MATRIX ON THE AO BASIS and their OVERLAP
    
    THE OVERLAP SHOULD BE THE SAME AS overlap_fock_tc_eigvec_mo

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_tc_reigvec_mo`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: fock_tc_reigvec_mo


    File : :file:`diago_bi_ort_tcfock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_tc_reigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: fock_tc_leigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: eigval_fock_tc_mo	(mo_num)
        double precision, allocatable	:: overlap_fock_tc_eigvec_mo	(mo_num,mo_num)


    EIGENVECTORS OF FOCK MATRIX ON THE MO BASIS and their OVERLAP

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`im_thresh_tcscf`
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_num`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_tc_reigvec_ao`

 
.. c:var:: fock_vartc_eigvec_ao


    File : :file:`diago_vartcfock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_vartc_eigvec_ao	(ao_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_vartc_eigvec_mo`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: fock_vartc_eigvec_mo


    File : :file:`diago_vartcfock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_vartc_eigvec_mo	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_vartc_mo_tot`
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_vartc_eigvec_ao`

 
.. c:var:: fqs_sqf_ao


    File : :file:`diis_tcscf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fqs_sqf_ao	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`fock_matrix_vartc_ao_tot`
       * :c:data:`q_matrix`
       * :c:data:`var_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fqs_sqf_mo`

 
.. c:var:: fqs_sqf_mo


    File : :file:`diis_tcscf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fqs_sqf_mo	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fqs_sqf_ao`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: good_hermit_tc_fock_mat


    File : :file:`fock_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: good_hermit_tc_fock_mat	(mo_num,mo_num)


    good_hermit_tc_fock_mat = Hermitian Upper triangular Fock matrix
    
    The converged eigenvectors of such matrix yield to orthonormal vectors satisfying the left Brillouin theorem

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`grad_good_hermit_tc_fock_mat`

 
.. c:var:: grad_good_hermit_tc_fock_mat


    File : :file:`fock_hermit.irp.f`

    .. code:: fortran

        double precision	:: grad_good_hermit_tc_fock_mat	


    grad_good_hermit_tc_fock_mat = norm of gradients of the upper triangular TC fock

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`good_hermit_tc_fock_mat`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`grad_hermit`

 
.. c:var:: grad_hermit


    File : :file:`fock_hermit.irp.f`

    .. code:: fortran

        double precision	:: grad_hermit	


    square of gradient of the energy

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`grad_good_hermit_tc_fock_mat`
       * :c:data:`grad_hermit_average_tc_fock_mat`
       * :c:data:`symetric_fock_tc`


 
.. c:var:: grad_hermit_average_tc_fock_mat


    File : :file:`fock_hermit.irp.f`

    .. code:: fortran

        double precision	:: grad_hermit_average_tc_fock_mat	


    grad_hermit_average_tc_fock_mat = norm of gradients of the upper triangular TC fock

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`hermit_average_tc_fock_mat`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`grad_hermit`

 
.. c:var:: grad_non_hermit


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision	:: grad_non_hermit_left	
        double precision	:: grad_non_hermit_right	
        double precision	:: grad_non_hermit	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`


 
.. c:var:: grad_non_hermit_left


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision	:: grad_non_hermit_left	
        double precision	:: grad_non_hermit_right	
        double precision	:: grad_non_hermit	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`


 
.. c:var:: grad_non_hermit_right


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision	:: grad_non_hermit_left	
        double precision	:: grad_non_hermit_right	
        double precision	:: grad_non_hermit	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`


 
.. c:var:: grad_vartc


    File : :file:`fock_vartc.irp.f`

    .. code:: fortran

        double precision	:: grad_vartc	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_vartc_mo_tot`
       * :c:data:`mo_num`


 
.. c:var:: hermit_average_tc_fock_mat


    File : :file:`fock_hermit.irp.f`

    .. code:: fortran

        double precision, allocatable	:: hermit_average_tc_fock_mat	(mo_num,mo_num)


    hermit_average_tc_fock_mat = (F + F^\dagger)/2

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`grad_hermit_average_tc_fock_mat`

 
.. c:var:: overlap_fock_tc_eigvec_ao


    File : :file:`diago_bi_ort_tcfock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_tc_reigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: fock_tc_leigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: overlap_fock_tc_eigvec_ao	(mo_num,mo_num)


    EIGENVECTORS OF FOCK MATRIX ON THE AO BASIS and their OVERLAP
    
    THE OVERLAP SHOULD BE THE SAME AS overlap_fock_tc_eigvec_mo

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_tc_reigvec_mo`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: overlap_fock_tc_eigvec_mo


    File : :file:`diago_bi_ort_tcfock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_tc_reigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: fock_tc_leigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: eigval_fock_tc_mo	(mo_num)
        double precision, allocatable	:: overlap_fock_tc_eigvec_mo	(mo_num,mo_num)


    EIGENVECTORS OF FOCK MATRIX ON THE MO BASIS and their OVERLAP

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`im_thresh_tcscf`
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_num`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_tc_reigvec_ao`

 
.. c:var:: q_alpha


    File : :file:`diis_tcscf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: q_alpha	(ao_num,ao_num)


    
    Q_alpha = mo_r_coef x eta_occ_alpha x mo_l_coef.T
    
    [Q_alpha]_ij = \sum_{k=1}^{elec_alpha_num} [mo_r_coef]_ik [mo_l_coef]_jk
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`q_matrix`

 
.. c:var:: q_beta


    File : :file:`diis_tcscf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: q_beta	(ao_num,ao_num)


    
    Q_beta = mo_r_coef x eta_occ_beta x mo_l_coef.T
    
    [Q_beta]_ij = \sum_{k=1}^{elec_beta_num} [mo_r_coef]_ik [mo_l_coef]_jk
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`q_matrix`

 
.. c:var:: q_matrix


    File : :file:`diis_tcscf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: q_matrix	(ao_num,ao_num)


    
    Q_matrix = 2 mo_r_coef x eta_occ x mo_l_coef.T
    
    with:
                           | 1   if i = j = 1, ..., nb of occ orbitals
           [eta_occ]_ij =  |
                           | 0   otherwise
    
    the diis error is defines as:
                            e = F_ao x Q x ao_overlap - ao_overlap x Q x F_ao
    with:
          mo_l_coef.T x ao_overlap x mo_r_coef = I
          F_mo = mo_l_coef.T x F_ao x mo_r_coef
          F_ao = (ao_overlap x mo_r_coef) x F_mo x (ao_overlap x mo_l_coef).T
    
    ==> e = 2 ao_overlap x mo_r_coef x [ F_mo x eta_occ - eta_occ x F_mo ] x (ao_overlap x mo_l_coef).T
    
         at convergence:
                                         F_mo x eta_occ - eta_occ x F_mo = 0
                                     ==> [F_mo]_ij ([eta_occ]_ii - [eta_occ]_jj) = 0
                                     ==> [F_mo]_ia = [F_mo]_ai = 0 where: i = occ and a = vir
                                     ==> Brillouin conditions
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`q_alpha`
       * :c:data:`q_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fqs_sqf_ao`

 
.. c:var:: tc_hf_energy


    File : :file:`tc_scf_energy.irp.f`

    .. code:: fortran

        double precision	:: tc_hf_energy	
        double precision	:: tc_hf_one_e_energy	
        double precision	:: tc_hf_two_e_energy	


    TC Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`diag_three_elem_hf`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`nuclear_repulsion`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`two_e_tc_non_hermit_integral_alpha`


 
.. c:var:: tc_hf_one_e_energy


    File : :file:`tc_scf_energy.irp.f`

    .. code:: fortran

        double precision	:: tc_hf_energy	
        double precision	:: tc_hf_one_e_energy	
        double precision	:: tc_hf_two_e_energy	


    TC Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`diag_three_elem_hf`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`nuclear_repulsion`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`two_e_tc_non_hermit_integral_alpha`


 
.. c:var:: tc_hf_two_e_energy


    File : :file:`tc_scf_energy.irp.f`

    .. code:: fortran

        double precision	:: tc_hf_energy	
        double precision	:: tc_hf_one_e_energy	
        double precision	:: tc_hf_two_e_energy	


    TC Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`diag_three_elem_hf`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`nuclear_repulsion`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`two_e_tc_non_hermit_integral_alpha`


 
.. c:var:: tc_scf_dm_in_r


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_scf_dm_in_r	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mos_l_in_r_array`
       * :c:data:`mos_r_in_r_array`
       * :c:data:`n_points_final_grid`


 
.. c:var:: tcscf_density_matrix_ao_alpha


    File : :file:`tc_scf_dm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tcscf_density_matrix_ao_alpha	(ao_num,ao_num)


    TC-SCF transition density matrix on the AO basis for ALPHA electrons

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`bi_ortho`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`tcscf_bi_ort_dm_ao_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_hf_energy`
       * :c:data:`tcscf_density_matrix_ao_tot`
       * :c:data:`two_e_tc_non_hermit_integral_alpha`
       * :c:data:`two_e_tc_non_hermit_integral_seq_alpha`
       * :c:data:`two_e_vartc_integral_alpha`
       * :c:data:`vartc_hf_energy`

 
.. c:var:: tcscf_density_matrix_ao_beta


    File : :file:`tc_scf_dm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tcscf_density_matrix_ao_beta	(ao_num,ao_num)


    TC-SCF transition density matrix on the AO basis for BETA electrons

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`bi_ortho`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`scf_density_matrix_ao_beta`
       * :c:data:`tcscf_bi_ort_dm_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_hf_energy`
       * :c:data:`tcscf_density_matrix_ao_tot`
       * :c:data:`two_e_tc_non_hermit_integral_alpha`
       * :c:data:`two_e_tc_non_hermit_integral_seq_alpha`
       * :c:data:`two_e_vartc_integral_alpha`
       * :c:data:`vartc_hf_energy`

 
.. c:var:: tcscf_density_matrix_ao_tot


    File : :file:`tc_scf_dm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tcscf_density_matrix_ao_tot	(ao_num,ao_num)


    TC-SCF transition density matrix on the AO basis for ALPHA+BETA electrons

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`


 
.. c:var:: two_e_tc_non_hermit_integral_alpha


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_tc_non_hermit_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: two_e_tc_non_hermit_integral_beta	(ao_num,ao_num)


    
    two_e_tc_non_hermit_integral_alpha(k,i) = <k| F^tc_alpha |i> ON THE AO BASIS
    
    where F^tc is the TWO-BODY part of the TC Fock matrix and k,i are AO basis functions
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_ao_alpha`
       * :c:data:`fock_matrix_tc_ao_beta`
       * :c:data:`tc_hf_energy`

 
.. c:var:: two_e_tc_non_hermit_integral_beta


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_tc_non_hermit_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: two_e_tc_non_hermit_integral_beta	(ao_num,ao_num)


    
    two_e_tc_non_hermit_integral_alpha(k,i) = <k| F^tc_alpha |i> ON THE AO BASIS
    
    where F^tc is the TWO-BODY part of the TC Fock matrix and k,i are AO basis functions
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_ao_alpha`
       * :c:data:`fock_matrix_tc_ao_beta`
       * :c:data:`tc_hf_energy`

 
.. c:var:: two_e_tc_non_hermit_integral_seq_alpha


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_tc_non_hermit_integral_seq_alpha	(ao_num,ao_num)
        double precision, allocatable	:: two_e_tc_non_hermit_integral_seq_beta	(ao_num,ao_num)


    
    two_e_tc_non_hermit_integral_seq_alpha(k,i) = <k| F^tc_alpha |i> ON THE AO BASIS
    
    where F^tc is the TWO-BODY part of the TC Fock matrix and k,i are AO basis functions
    
    works in SEQUENTIAL

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`


 
.. c:var:: two_e_tc_non_hermit_integral_seq_beta


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_tc_non_hermit_integral_seq_alpha	(ao_num,ao_num)
        double precision, allocatable	:: two_e_tc_non_hermit_integral_seq_beta	(ao_num,ao_num)


    
    two_e_tc_non_hermit_integral_seq_alpha(k,i) = <k| F^tc_alpha |i> ON THE AO BASIS
    
    where F^tc is the TWO-BODY part of the TC Fock matrix and k,i are AO basis functions
    
    works in SEQUENTIAL

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`


 
.. c:var:: two_e_vartc_integral_alpha


    File : :file:`fock_vartc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_vartc_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: two_e_vartc_integral_beta	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_vartc_tot`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_vartc_ao_alpha`
       * :c:data:`fock_matrix_vartc_ao_beta`
       * :c:data:`vartc_hf_energy`

 
.. c:var:: two_e_vartc_integral_beta


    File : :file:`fock_vartc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_vartc_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: two_e_vartc_integral_beta	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_vartc_tot`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_vartc_ao_alpha`
       * :c:data:`fock_matrix_vartc_ao_beta`
       * :c:data:`vartc_hf_energy`

 
.. c:var:: v_sum_in_r


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        double precision, allocatable	:: v_sum_in_r	(n_points_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ki_bi_ortho_erf_rk`


 
.. c:var:: vartc_hf_energy


    File : :file:`tc_scf_energy.irp.f`

    .. code:: fortran

        double precision	:: vartc_hf_energy	
        double precision	:: vartc_hf_one_e_energy	
        double precision	:: vartc_hf_two_e_energy	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`diag_three_elem_hf`
       * :c:data:`mo_r_coef`
       * :c:data:`nuclear_repulsion`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`two_e_vartc_integral_alpha`


 
.. c:var:: vartc_hf_one_e_energy


    File : :file:`tc_scf_energy.irp.f`

    .. code:: fortran

        double precision	:: vartc_hf_energy	
        double precision	:: vartc_hf_one_e_energy	
        double precision	:: vartc_hf_two_e_energy	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`diag_three_elem_hf`
       * :c:data:`mo_r_coef`
       * :c:data:`nuclear_repulsion`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`two_e_vartc_integral_alpha`


 
.. c:var:: vartc_hf_two_e_energy


    File : :file:`tc_scf_energy.irp.f`

    .. code:: fortran

        double precision	:: vartc_hf_energy	
        double precision	:: vartc_hf_one_e_energy	
        double precision	:: vartc_hf_two_e_energy	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`diag_three_elem_hf`
       * :c:data:`mo_r_coef`
       * :c:data:`nuclear_repulsion`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`two_e_vartc_integral_alpha`


 
.. c:var:: w1_diag_in_r


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        double precision, allocatable	:: w1_diag_in_r	(n_points_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ki_bi_ortho_erf_rk_diag`


 
.. c:var:: w1_in_r


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        double precision, allocatable	:: w1_in_r	(n_points_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`w1_l_in_r`


 
.. c:var:: w1_l_in_r


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        double precision, allocatable	:: w1_l_in_r	(n_points_final_grid,3,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ki_bi_ortho_erf_rk`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`w1_in_r`
       * :c:data:`w1_w1_l_in_r`

 
.. c:var:: w1_r_in_r


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        double precision, allocatable	:: w1_r_in_r	(n_points_final_grid,3,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ki_bi_ortho_erf_rk`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`w1_w1_r_in_r`

 
.. c:var:: w1_w1_l_in_r


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        double precision, allocatable	:: w1_w1_l_in_r	(n_points_final_grid,3,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`w1_l_in_r`
       * :c:data:`x_w_ki_bi_ortho_erf_rk`


 
.. c:var:: w1_w1_r_in_r


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        double precision, allocatable	:: w1_w1_r_in_r	(n_points_final_grid,3,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`w1_r_in_r`
       * :c:data:`x_w_ki_bi_ortho_erf_rk`


 
.. c:var:: w_sum_in_r


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        double precision, allocatable	:: w_sum_in_r	(n_points_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ki_bi_ortho_erf_rk_diag`


 
.. c:var:: ww_sum_in_r


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ww_sum_in_r	(n_points_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ki_bi_ortho_erf_rk_diag`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: build_s_matrix:


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine build_s_matrix(m, n, C1, C2, overlap, smat)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`orthog_functions`
       * :c:func:`routine_save_rotated_mos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: comb_orbitals:


    File : :file:`combine_lr_tcscf.irp.f`

    .. code:: fortran

        subroutine comb_orbitals()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`elec_alpha_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`combine_lr_tcscf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_weighted_biorthog`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_r_coef`
       * :c:func:`impose_weighted_orthog_svd`

 
.. c:function:: combine_lr_tcscf:


    File : :file:`combine_lr_tcscf.irp.f`

    .. code:: fortran

        subroutine combine_lr_tcscf


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bi_ortho`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`comb_orbitals`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`bi_ortho`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: contrib_3e_diag_soo:


    File : :file:`three_e_energy_bi_ortho.irp.f`

    .. code:: fortran

        subroutine contrib_3e_diag_soo(i,j,k,integral)


    returns the pure same spin contribution to diagonal matrix element of 3e term

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_aab_contrib`
       * :c:func:`give_aab_contrib_bis`
       * :c:func:`give_abb_contrib`
       * :c:func:`give_abb_contrib_bis`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`

 
.. c:function:: contrib_3e_diag_sss:


    File : :file:`three_e_energy_bi_ortho.irp.f`

    .. code:: fortran

        subroutine contrib_3e_diag_sss(i,j,k,integral)


    returns the pure same spin contribution to diagonal matrix element of 3e term

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_aaa_contrib`
       * :c:func:`give_aaa_contrib_bis`
       * :c:func:`give_bbb_contrib`
       * :c:func:`give_bbb_contrib_bis`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`

 
.. c:function:: contrib_3e_soo:


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        subroutine contrib_3e_soo(a,i,j,k,integral)


    returns the same spin / opposite spin / opposite spin contribution to F(a,i) from two orbitals j,k

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_b_tmp2_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`

 
.. c:function:: contrib_3e_sos:


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        subroutine contrib_3e_sos(a, i, j, k, integral)


    returns the same spin / opposite spin / same spin contribution to F(a,i) from two orbitals j,k

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`

 
.. c:function:: contrib_3e_sss:


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        subroutine contrib_3e_sss(a, i, j, k, integral)


    returns the pure same spin contribution to F(a,i) from two orbitals j,k

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_a_tmp2_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`

 
.. c:function:: create_guess:


    File : :file:`tc_scf.irp.f`

    .. code:: fortran

        subroutine create_guess()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_ortho_lowdin_coef`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_coef`
       * :c:data:`mo_guess_type`
       * :c:data:`mo_label`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`huckel_guess`
       * :c:func:`mo_as_eigvectors_of_mo_matrix`
       * :c:func:`restore_symmetry`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`

 
.. c:function:: cyclic_term_jim_bi_ortho:


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        subroutine cyclic_term_jim_bi_ortho(a, i, integral)


    computes sum_(j,m = 1, elec_beta_num) < a m j | j i m > with bi ortho mos

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`w1_l_in_r`
       * :c:data:`w1_r_in_r`
       * :c:data:`w1_w1_l_in_r`
       * :c:data:`w1_w1_r_in_r`

 
.. c:function:: cyclic_term_mji_bi_ortho:


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        subroutine cyclic_term_mji_bi_ortho(a, i, integral)


    computes sum_(j,m = 1, elec_beta_num) < a m j | m j i > with bi ortho mos

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`w1_l_in_r`
       * :c:data:`w1_r_in_r`
       * :c:data:`w1_w1_l_in_r`
       * :c:data:`w1_w1_r_in_r`

 
.. c:function:: direct_term_ijj_bi_ortho:


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        subroutine direct_term_ijj_bi_ortho(a, i, integral)


    computes sum_(j = 1, elec_beta_num) < a j j | i j j > with bi ortho mos

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`w1_diag_in_r`
       * :c:data:`ww_sum_in_r`
       * :c:data:`x_w_ki_bi_ortho_erf_rk`

 
.. c:function:: direct_term_imj_bi_ortho:


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        subroutine direct_term_imj_bi_ortho(a, i, integral)


    computes sum_(j,m = 1, elec_beta_num) < a m j | i m j > with bi ortho mos

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`tc_scf_dm_in_r`
       * :c:data:`w_sum_in_r`
       * :c:data:`x_w_ki_bi_ortho_erf_rk`

 
.. c:function:: exch_term_ijm_bi_ortho:


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        subroutine exch_term_ijm_bi_ortho(a, i, integral)


    computes sum_(j,m = 1, elec_beta_num) < a m j | i j m > with bi ortho mos

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_sum_in_r`
       * :c:data:`w1_in_r`
       * :c:data:`x_w_ki_bi_ortho_erf_rk`

 
.. c:function:: exch_term_jmi_bi_ortho:


    File : :file:`integrals_in_r_stuff.irp.f`

    .. code:: fortran

        subroutine exch_term_jmi_bi_ortho(a, i, integral)


    computes sum_(j,m = 1, elec_beta_num) < a m j | j m i > with bi ortho mos

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`tc_scf_dm_in_r`
       * :c:data:`w1_l_in_r`
       * :c:data:`w1_r_in_r`
       * :c:data:`w_sum_in_r`
       * :c:data:`x_w_ki_bi_ortho_erf_rk`

 
.. c:function:: extrapolate_tc_fock_matrix:


    File : :file:`rh_tcscf_diis.irp.f`

    .. code:: fortran

        subroutine extrapolate_TC_Fock_matrix(E_DIIS, F_DIIS, F_ao, size_F_ao, it, dim_DIIS)


    
    Compute the extrapolated Fock matrix using the DIIS procedure
    
    e = \sum_i c_i e_i and \sum_i c_i = 1
    ==> lagrange multiplier with L = |e|^2 - \lambda (\sum_i c_i = 1)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`max_dim_diis_tcscf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`rh_tcscf_diis`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgecon`
       * :c:func:`dgesv`
       * :c:func:`dgetrf`

 
.. c:function:: fix_right_to_one:


    File : :file:`rotate_tcscf_orbitals.irp.f`

    .. code:: fortran

        subroutine fix_right_to_one()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_weighted_biorthog_binormalize`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_l_coef`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_r_coef`
       * :c:func:`impose_weighted_biorthog_qr`
       * :c:func:`impose_weighted_orthog_svd`

 
.. c:function:: give_aaa_contrib:


    File : :file:`three_e_energy_bi_ortho.irp.f`

    .. code:: fortran

        subroutine give_aaa_contrib(integral_aaa)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`diag_three_elem_hf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`contrib_3e_diag_sss`

 
.. c:function:: give_aaa_contrib_bis:


    File : :file:`three_e_energy_bi_ortho.irp.f`

    .. code:: fortran

        subroutine give_aaa_contrib_bis(integral_aaa)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`contrib_3e_diag_sss`

 
.. c:function:: give_aab_contrib:


    File : :file:`three_e_energy_bi_ortho.irp.f`

    .. code:: fortran

        subroutine give_aab_contrib(integral_aab)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`diag_three_elem_hf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`contrib_3e_diag_soo`

 
.. c:function:: give_aab_contrib_bis:


    File : :file:`three_e_energy_bi_ortho.irp.f`

    .. code:: fortran

        subroutine give_aab_contrib_bis(integral_aab)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`contrib_3e_diag_soo`

 
.. c:function:: give_abb_contrib:


    File : :file:`three_e_energy_bi_ortho.irp.f`

    .. code:: fortran

        subroutine give_abb_contrib(integral_abb)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`diag_three_elem_hf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`contrib_3e_diag_soo`

 
.. c:function:: give_abb_contrib_bis:


    File : :file:`three_e_energy_bi_ortho.irp.f`

    .. code:: fortran

        subroutine give_abb_contrib_bis(integral_abb)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`contrib_3e_diag_soo`

 
.. c:function:: give_bbb_contrib:


    File : :file:`three_e_energy_bi_ortho.irp.f`

    .. code:: fortran

        subroutine give_bbb_contrib(integral_bbb)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`diag_three_elem_hf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`contrib_3e_diag_sss`

 
.. c:function:: give_bbb_contrib_bis:


    File : :file:`three_e_energy_bi_ortho.irp.f`

    .. code:: fortran

        subroutine give_bbb_contrib_bis(integral_bbb)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`contrib_3e_diag_sss`

 
.. c:function:: give_fock_ia_three_e_total:


    File : :file:`fock_three_hermit.irp.f`

    .. code:: fortran

        subroutine give_fock_ia_three_e_total(i,a,contrib)


    contrib is the TOTAL (same spins / opposite spins) contribution from the three body term to the Fock operator
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`fock_3_rho_beta`
       * :c:data:`fock_3_trace_w_tilde`
       * :c:data:`fock_3_w_ki_mos_k`
       * :c:data:`fock_3_w_ki_wk_a`
       * :c:data:`fock_3_w_kk_sum`
       * :c:data:`fock_3_w_kl_mo_k_mo_l`
       * :c:data:`fock_3_w_kl_wla_phi_k`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_w_ij_erf_rk`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3_mat`

 
.. c:function:: ltxsxr:


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine LTxSxR(n, m, L, S, R, C)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`maximize_overlap`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: main:


    File : :file:`tc_petermann_factor.irp.f`

    .. code:: fortran

        subroutine main()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_fit_param`
       * :c:func:`tc_petermann_factor`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: maximize_overlap:


    File : :file:`rotate_tcscf_orbitals.irp.f`

    .. code:: fortran

        subroutine maximize_overlap()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_coef`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`rotate_tcscf_orbitals`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_bi_ortho_mos_mo_l_coef`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_r_coef`
       * :c:func:`ltxsxr`
       * :c:func:`rotate_degen_eigvec_to_maximize_overlap`

 
.. c:function:: minimize_tc_orb_angles:


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine minimize_tc_orb_angles()


    routine that minimizes the angle between left- and right-orbitals when degeneracies are found

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thr_degen_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_angles`
       * :c:func:`tc_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_angles_tc`
       * :c:func:`print_energy_and_mos`
       * :c:func:`routine_save_rotated_mos`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

 
.. c:function:: molden:


    File : :file:`molden_lr_mos.irp.f`

    .. code:: fortran

        subroutine molden


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`molden_lr`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: molden_lr:


    File : :file:`molden_lr_mos.irp.f`

    .. code:: fortran

        subroutine molden_lr


    Produces a Molden file

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_l`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`element_name`
       * :c:data:`ezfio_filename`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_occ`
       * :c:data:`mo_r_coef`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_list_shell_aos`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`molden`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`isort`

 
.. c:function:: orthog_functions:


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine orthog_functions(m, n, coef, overlap)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`routine_save_rotated_mos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`build_s_matrix`
       * :c:func:`impose_orthog_svd_overlap`

 
.. c:function:: print_angles:


    File : :file:`minimize_tc_angles.irp.f`

    .. code:: fortran

        subroutine print_angles


    program that minimizes the angle between left- and right-orbitals when degeneracies are found in the TC-Fock matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`minimize_tc_orb_angles`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: print_angles_tc:


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine print_angles_tc()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`angle_left_right`
       * :c:data:`mo_num`
       * :c:data:`overlap_mo_r`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`minimize_tc_orb_angles`

 
.. c:function:: print_energy_and_mos:


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine print_energy_and_mos()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`angle_left_right`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`grad_non_hermit_left`
       * :c:data:`mo_num`
       * :c:data:`overlap_mo_r`
       * :c:data:`tc_hf_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`minimize_tc_orb_angles`
       * :c:func:`rh_tcscf_diis`
       * :c:func:`rh_tcscf_simple`
       * :c:func:`tc_scf`

 
.. c:function:: print_fit_param:


    File : :file:`print_fit_param.irp.f`

    .. code:: fortran

        subroutine print_fit_param


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`main`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: rh_tcscf_diis:


    File : :file:`rh_tcscf_diis.irp.f`

    .. code:: fortran

        subroutine rh_tcscf_diis()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`diag_three_elem_hf`
       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`fock_tc_reigvec_ao`
       * :c:data:`fqs_sqf_ao`
       * :c:data:`fqs_sqf_mo`
       * :c:data:`level_shift_tcscf`
       * :c:data:`max_dim_diis_tcscf`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_it_tcscf_max`
       * :c:data:`tc_hf_energy`
       * :c:data:`three_body_h_tc`
       * :c:data:`thresh_tcscf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`tc_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_to_mo_bi_ortho`
       * :c:func:`extrapolate_tc_fock_matrix`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_l_coef`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_r_coef`
       * :c:func:`ezfio_set_tc_scf_bitc_energy`
       * :c:func:`print_energy_and_mos`
       * :c:func:`wall_time`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

 
.. c:function:: rh_tcscf_simple:


    File : :file:`rh_tcscf_simple.irp.f`

    .. code:: fortran

        subroutine rh_tcscf_simple()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`bi_ortho`
       * :c:data:`diag_three_elem_hf`
       * :c:data:`fock_tc_reigvec_ao`
       * :c:data:`fqs_sqf_mo`
       * :c:data:`grad_hermit`
       * :c:data:`grad_non_hermit_left`
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_coef`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`n_it_tcscf_max`
       * :c:data:`tc_hf_energy`
       * :c:data:`three_body_h_tc`
       * :c:data:`thresh_tcscf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`tc_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_bi_ortho_mos_mo_l_coef`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_r_coef`
       * :c:func:`ezfio_set_tc_scf_bitc_energy`
       * :c:func:`print_energy_and_mos`
       * :c:func:`save_good_hermit_tc_eigvectors`
       * :c:func:`save_mos`
       * :c:func:`wall_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_coef`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

 
.. c:function:: rh_vartcscf_simple:


    File : :file:`rh_vartcscf_simple.irp.f`

    .. code:: fortran

        subroutine rh_vartcscf_simple()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`diag_three_elem_hf`
       * :c:data:`fock_vartc_eigvec_ao`
       * :c:data:`fqs_sqf_mo`
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`n_it_tcscf_max`
       * :c:data:`three_body_h_tc`
       * :c:data:`thresh_tcscf`
       * :c:data:`vartc_hf_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`tc_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_bi_ortho_mos_mo_l_coef`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_r_coef`
       * :c:func:`ezfio_set_tc_scf_bitc_energy`
       * :c:func:`wall_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

 
.. c:function:: rotate_degen_eigvec_to_maximize_overlap:


    File : :file:`rotate_tcscf_orbitals.irp.f`

    .. code:: fortran

        subroutine rotate_degen_eigvec_to_maximize_overlap(n, m, e0, C0, W0, L0, R0)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thr_degen_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`maximize_overlap`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`maxovl`

 
.. c:function:: rotate_tcscf_orbitals:


    File : :file:`rotate_tcscf_orbitals.irp.f`

    .. code:: fortran

        subroutine rotate_tcscf_orbitals


    TODO : Rotate the bi-orthonormal orbitals in order to minimize left-right angles when degenerate

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bi_ortho`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`maximize_overlap`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`bi_ortho`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: routine_int2_grad1u2_grad2u2_j1b2:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine routine_int2_grad1u2_grad2u2_j1b2



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`elec_alpha_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1u2_grad2u2_j1b2_test`
       * :c:data:`mo_coef`
       * :c:data:`n_points_final_grid`

 
.. c:function:: routine_int2_u2_j1b2:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine routine_int2_u2_j1b2



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u2_j1b2`
       * :c:data:`int2_u2_j1b2_test`
       * :c:data:`n_points_final_grid`

 
.. c:function:: routine_int2_u_grad1u_j1b2:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine routine_int2_u_grad1u_j1b2



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u_grad1u_j1b2`
       * :c:data:`int2_u_grad1u_j1b2_test`
       * :c:data:`n_points_final_grid`

 
.. c:function:: routine_int2_u_grad1u_x_j1b2:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine routine_int2_u_grad1u_x_j1b2



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u_grad1u_x_j1b2`
       * :c:data:`int2_u_grad1u_x_j1b2_test`
       * :c:data:`n_points_final_grid`

 
.. c:function:: routine_save_rotated_mos:


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine routine_save_rotated_mos(thr_deg, good_angles)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`angle_left_right`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`grad_non_hermit_left`
       * :c:data:`maxovl_tc`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`overlap_mo_r`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`minimize_tc_orb_angles`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`build_s_matrix`
       * :c:func:`dgemm`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_l_coef`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_r_coef`
       * :c:func:`give_degen_full_list`
       * :c:func:`impose_weighted_biorthog_svd`
       * :c:func:`maxovl`
       * :c:func:`orthog_functions`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

 
.. c:function:: routine_test_j1b:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine routine_test_j1b



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`list_all_comb_b3_coef`
       * :c:data:`list_all_comb_b3_size`
       * :c:data:`list_comb_thr_b3_coef`
       * :c:data:`list_comb_thr_b3_size`

 
.. c:function:: routine_v_ij_erf_rk_cst_mu_j1b:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine routine_v_ij_erf_rk_cst_mu_j1b



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_erf_rk_cst_mu_j1b`
       * :c:data:`v_ij_erf_rk_cst_mu_j1b_test`

 
.. c:function:: routine_v_ij_u_cst_mu_j1b:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine routine_v_ij_u_cst_mu_j1b



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_u_cst_mu_j1b`
       * :c:data:`v_ij_u_cst_mu_j1b_test`

 
.. c:function:: routine_v_ij_u_cst_mu_j1b_test:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine routine_v_ij_u_cst_mu_j1b_test



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_u_cst_mu_j1b`
       * :c:data:`v_ij_u_cst_mu_j1b_test`

 
.. c:function:: routine_x_v_ij_erf_rk_cst_mu_j1b:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine routine_x_v_ij_erf_rk_cst_mu_j1b



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_v_ij_erf_rk_cst_mu_j1b`
       * :c:data:`x_v_ij_erf_rk_cst_mu_j1b_test`

 
.. c:function:: save_good_hermit_tc_eigvectors:


    File : :file:`fock_hermit.irp.f`

    .. code:: fortran

        subroutine save_good_hermit_tc_eigvectors()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`good_hermit_tc_fock_mat`
       * :c:data:`hermit_average_tc_fock_mat`
       * :c:data:`mo_num`
       * :c:data:`symetric_fock_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`rh_tcscf_simple`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`mo_as_eigvectors_of_mo_matrix`

 
.. c:function:: sort_by_tc_fock:


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine sort_by_tc_fock



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsort`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

 
.. c:function:: tc_petermann_factor:


    File : :file:`tc_petermann_factor.irp.f`

    .. code:: fortran

        subroutine tc_petermann_factor


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`main`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: tc_scf:


    File : :file:`tc_scf.irp.f`

    .. code:: fortran

        subroutine tc_scf


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`j1b_pen`
       * :c:data:`j1b_type`
       * :c:data:`mu_erf`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`tcscf_algorithm`
       * :c:data:`var_tc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`minimize_tc_orb_angles`
       * :c:func:`print_energy_and_mos`
       * :c:func:`rh_tcscf_diis`
       * :c:func:`rh_tcscf_simple`
       * :c:func:`rh_vartcscf_simple`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_coef`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: test_ao_tc_int_chemist:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_ao_tc_int_chemist



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_tc_int_chemist`

 
.. c:function:: test_fock_3e_uhf_ao:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_fock_3e_uhf_ao()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3e_uhf_ao_a`
       * :c:data:`fock_3e_uhf_ao_b`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_to_mo_bi_ortho`

 
.. c:function:: test_fock_3e_uhf_mo:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_fock_3e_uhf_mo()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3e_uhf_mo_a`
       * :c:data:`fock_3e_uhf_mo_b`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`mo_num`

 
.. c:function:: test_grid_points_ao:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_grid_points_ao



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_u_cst_mu_j1b_ng_1_test`
       * :c:data:`v_ij_u_cst_mu_j1b_test`

 
.. c:function:: test_int_gauss:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_int_gauss



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs_grid`
       * :c:data:`ao_prod_center`
       * :c:data:`ao_prod_sigma`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

 
.. c:function:: test_ints:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_ints


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_extra_grid_becke`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_extra_grid`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_extra_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`tc_hf_energy`
       * :c:data:`vartc_hf_energy`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`test_tc_grad_square_ao_test`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_extra_grid_becke`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_extra_grid`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_extra_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: test_tc_grad_and_lapl_ao:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_tc_grad_and_lapl_ao()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`tc_grad_and_lapl_ao`
       * :c:data:`tc_grad_and_lapl_ao_loop`

 
.. c:function:: test_tc_grad_square_ao:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_tc_grad_square_ao()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`tc_grad_square_ao`
       * :c:data:`tc_grad_square_ao_loop`

 
.. c:function:: test_tc_grad_square_ao_test:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_tc_grad_square_ao_test()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`tc_grad_square_ao_test`
       * :c:data:`tc_grad_square_ao_test_ref`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_ints`

 
.. c:function:: test_tc_scf:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_tc_scf



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`x_v_ij_erf_rk_cst_mu_j1b_test`

 
.. c:function:: test_total_grad_lapl:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_total_grad_lapl



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`tc_grad_and_lapl_ao`
       * :c:data:`tc_grad_and_lapl_ao_test`

 
.. c:function:: test_total_grad_square:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_total_grad_square



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`tc_grad_square_ao`
       * :c:data:`tc_grad_square_ao_test`

 
.. c:function:: test_two_e_tc_non_hermit_integral:


    File : :file:`test_int.irp.f`

    .. code:: fortran

        subroutine test_two_e_tc_non_hermit_integral()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`two_e_tc_non_hermit_integral_alpha`
       * :c:data:`two_e_tc_non_hermit_integral_alpha`
       * :c:data:`two_e_tc_non_hermit_integral_seq_alpha`
       * :c:data:`two_e_tc_non_hermit_integral_seq_alpha`

