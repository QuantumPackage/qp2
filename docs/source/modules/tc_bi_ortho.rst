.. _module_tc_bi_ortho: 
 
.. program:: tc_bi_ortho 
 
.. default-role:: option 
 
tc_bi_ortho
===========

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Providers 
--------- 
 
.. c:var:: coef_pt1_bi_ortho


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: coef_pt1_bi_ortho	(N_det)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`e_tilde_00`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`


 
.. c:var:: e_corr_bi_orth


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`


 
.. c:var:: e_corr_bi_orth_proj


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`


 
.. c:var:: e_corr_double_bi_orth


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`


 
.. c:var:: e_corr_single_bi_orth


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`


 
.. c:var:: e_pt2_tc_bi_orth


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_pt2_tc_bi_orth	
        double precision	:: e_pt2_tc_bi_orth_single	
        double precision	:: e_pt2_tc_bi_orth_double	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`e_tilde_00`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`


 
.. c:var:: e_pt2_tc_bi_orth_double


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_pt2_tc_bi_orth	
        double precision	:: e_pt2_tc_bi_orth_single	
        double precision	:: e_pt2_tc_bi_orth_double	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`e_tilde_00`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`


 
.. c:var:: e_pt2_tc_bi_orth_single


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_pt2_tc_bi_orth	
        double precision	:: e_pt2_tc_bi_orth_single	
        double precision	:: e_pt2_tc_bi_orth_double	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`e_tilde_00`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`


 
.. c:var:: e_tc_left_right


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_tc_left_right	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`


 
.. c:var:: e_tilde_00


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_tilde_00	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`hf_bitmask`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_pt1_bi_ortho`
       * :c:data:`e_pt2_tc_bi_orth`

 
.. c:var:: e_tilde_bi_orth_00


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_tilde_bi_orth_00	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`hf_bitmask`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`

 
.. c:var:: eff_2_e_from_3_e_aa


    File : :file:`slater_tc_opt_double.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eff_2_e_from_3_e_aa	(mo_num,mo_num,mo_num,mo_num)


    eff_2_e_from_3_e_ab(p2,p1,h2,h1) = Effective Two-electron operator for alpha/alpha double excitations
    
    from contractionelec_alpha_num with HF density = a^{dagger}_p1_alpha a^{dagger}_p2_alpha a_h2_alpha a_h1_alpha
    
    WARNING :: to be coherent with the phase convention used in the Hamiltonian matrix elements, you must fulfill
    
    ||||    h2>h1, p2>p1   ||||

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_double_parrallel_spin_prov`


 
.. c:var:: eff_2_e_from_3_e_ab


    File : :file:`slater_tc_opt_double.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eff_2_e_from_3_e_ab	(mo_num,mo_num,mo_num,mo_num)


    eff_2_e_from_3_e_ab(p2,p1,h2,h1) = Effective Two-electron operator for alpha/beta double excitations
    
    from contraction with HF density = a^{dagger}_p1_alpha a^{dagger}_p2_beta a_h2_beta a_h1_alpha

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch13_bi_ort`
       * :c:data:`three_e_5_idx_exch23_bi_ort`


 
.. c:var:: eff_2_e_from_3_e_bb


    File : :file:`slater_tc_opt_double.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eff_2_e_from_3_e_bb	(mo_num,mo_num,mo_num,mo_num)


    eff_2_e_from_3_e_ab(p2,p1,h2,h1) = Effective Two-electron operator for beta/beta double excitations
    
    from contractionelec_beta_num with HF density = a^{dagger}_p1_beta a^{dagger}_p2_beta a_h2_beta a_h1_beta
    
    WARNING :: to be coherent with the phase convention used in the Hamiltonian matrix elements, you must fulfill
    
    ||||    h2>h1, p2>p1   ||||

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_double_parrallel_spin_prov`


 
.. c:var:: eigval_left_tc_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision	:: norm_ground_left_right_bi_orth	


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`qp_max_mem`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: eigval_right_tc_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision	:: norm_ground_left_right_bi_orth	


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`qp_max_mem`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: eigval_tc_cisd_sc2_bi_ortho


    File : :file:`tc_cisd_sc2_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: reigvec_tc_cisd_sc2_bi_ortho	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_cisd_sc2_bi_ortho	(N_det,N_states)
        double precision, allocatable	:: eigval_tc_cisd_sc2_bi_ortho	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`disk_based_davidson`
       * :c:data:`double_normal_ord`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`qp_max_mem`
       * :c:data:`three_body_h_tc`
       * :c:data:`threshold_nonsym_davidson`


 
.. c:var:: fock_diag_sorted_l_natorb


    File : :file:`tc_natorb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_diag_sorted_r_natorb	(mo_num,mo_num)
        double precision, allocatable	:: fock_diag_sorted_l_natorb	(mo_num,mo_num)
        double precision, allocatable	:: fock_diag_sorted_v_natorb	(mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`
       * :c:data:`natorb_tc_reigvec_mo`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`natorb_tc_reigvec_ao`

 
.. c:var:: fock_diag_sorted_r_natorb


    File : :file:`tc_natorb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_diag_sorted_r_natorb	(mo_num,mo_num)
        double precision, allocatable	:: fock_diag_sorted_l_natorb	(mo_num,mo_num)
        double precision, allocatable	:: fock_diag_sorted_v_natorb	(mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`
       * :c:data:`natorb_tc_reigvec_mo`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`natorb_tc_reigvec_ao`

 
.. c:var:: fock_diag_sorted_v_natorb


    File : :file:`tc_natorb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_diag_sorted_r_natorb	(mo_num,mo_num)
        double precision, allocatable	:: fock_diag_sorted_l_natorb	(mo_num,mo_num)
        double precision, allocatable	:: fock_diag_sorted_v_natorb	(mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`
       * :c:data:`natorb_tc_reigvec_mo`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`natorb_tc_reigvec_ao`

 
.. c:var:: fock_op_2_e_tc_closed_shell


    File : :file:`slater_tc_opt_single.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_op_2_e_tc_closed_shell	(mo_num,mo_num)


    Closed-shell part of the Fock operator for the TC operator

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`full_ijkl_bitmask`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`ref_closed_shell_bitmask`
       * :c:data:`tc_2e_3idx_coulomb_integrals`


 
.. c:var:: htilde_matrix_elmt_bi_ortho


    File : :file:`tc_hmat.irp.f`

    .. code:: fortran

        double precision, allocatable	:: htilde_matrix_elmt_bi_ortho	(N_det,N_det)


    htilde_matrix_elmt_bi_ortho(j,i) = <J| H^tilde |I>
    
    WARNING !!!!!!!!! IT IS NOT HERMITIAN !!!!!!!!!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`htilde_matrix_elmt_bi_ortho_tranp`

 
.. c:var:: htilde_matrix_elmt_bi_ortho_tranp


    File : :file:`tc_hmat.irp.f`

    .. code:: fortran

        double precision, allocatable	:: htilde_matrix_elmt_bi_ortho_tranp	(N_det,N_det)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`n_det`


 
.. c:var:: index_hf_psi_det


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        integer	:: index_hf_psi_det	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`

 
.. c:var:: leigvec_tc_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision	:: norm_ground_left_right_bi_orth	


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`qp_max_mem`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: leigvec_tc_cisd_sc2_bi_ortho


    File : :file:`tc_cisd_sc2_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: reigvec_tc_cisd_sc2_bi_ortho	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_cisd_sc2_bi_ortho	(N_det,N_states)
        double precision, allocatable	:: eigval_tc_cisd_sc2_bi_ortho	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`disk_based_davidson`
       * :c:data:`double_normal_ord`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`qp_max_mem`
       * :c:data:`three_body_h_tc`
       * :c:data:`threshold_nonsym_davidson`


 
.. c:var:: mo_r_coef_normalized


    File : :file:`test_spin_dens.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_r_coef_normalized	(ao_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_spin_dens_right_only`

 
.. c:var:: natorb_tc_eigval


    File : :file:`tc_natorb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: natorb_tc_reigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: natorb_tc_leigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: natorb_tc_eigval	(mo_num)


    
    natorb_tc_reigvec_mo : RIGHT eigenvectors of the ground state transition matrix (equivalent of natural orbitals)
    natorb_tc_leigvec_mo : LEFT  eigenvectors of the ground state transition matrix (equivalent of natural orbitals)
    natorb_tc_eigval     : eigenvalues of the ground state transition matrix (equivalent of the occupation numbers). WARNINING :: can be negative !!
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`
       * :c:data:`tc_transition_matrix_mo`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_diag_sorted_r_natorb`

 
.. c:var:: natorb_tc_leigvec_ao


    File : :file:`tc_natorb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: natorb_tc_reigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: natorb_tc_leigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: overlap_natorb_tc_eigvec_ao	(mo_num,mo_num)


    EIGENVECTORS OF FOCK MATRIX ON THE AO BASIS and their OVERLAP
    
    THE OVERLAP SHOULD BE THE SAME AS overlap_natorb_tc_eigvec_mo

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_diag_sorted_r_natorb`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: natorb_tc_leigvec_mo


    File : :file:`tc_natorb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: natorb_tc_reigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: natorb_tc_leigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: natorb_tc_eigval	(mo_num)


    
    natorb_tc_reigvec_mo : RIGHT eigenvectors of the ground state transition matrix (equivalent of natural orbitals)
    natorb_tc_leigvec_mo : LEFT  eigenvectors of the ground state transition matrix (equivalent of natural orbitals)
    natorb_tc_eigval     : eigenvalues of the ground state transition matrix (equivalent of the occupation numbers). WARNINING :: can be negative !!
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`
       * :c:data:`tc_transition_matrix_mo`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_diag_sorted_r_natorb`

 
.. c:var:: natorb_tc_reigvec_ao


    File : :file:`tc_natorb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: natorb_tc_reigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: natorb_tc_leigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: overlap_natorb_tc_eigvec_ao	(mo_num,mo_num)


    EIGENVECTORS OF FOCK MATRIX ON THE AO BASIS and their OVERLAP
    
    THE OVERLAP SHOULD BE THE SAME AS overlap_natorb_tc_eigvec_mo

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_diag_sorted_r_natorb`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: natorb_tc_reigvec_mo


    File : :file:`tc_natorb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: natorb_tc_reigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: natorb_tc_leigvec_mo	(mo_num,mo_num)
        double precision, allocatable	:: natorb_tc_eigval	(mo_num)


    
    natorb_tc_reigvec_mo : RIGHT eigenvectors of the ground state transition matrix (equivalent of natural orbitals)
    natorb_tc_leigvec_mo : LEFT  eigenvectors of the ground state transition matrix (equivalent of natural orbitals)
    natorb_tc_eigval     : eigenvalues of the ground state transition matrix (equivalent of the occupation numbers). WARNINING :: can be negative !!
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`mo_num`
       * :c:data:`tc_transition_matrix_mo`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_diag_sorted_r_natorb`

 
.. c:var:: norm_ground_left_right_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision	:: norm_ground_left_right_bi_orth	


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`qp_max_mem`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: normal_two_body_bi_orth


    File : :file:`normal_ordered.irp.f`

    .. code:: fortran

        double precision, allocatable	:: normal_two_body_bi_orth	(mo_num,mo_num,mo_num,mo_num)


    Normal ordering of the three body interaction on the HF density

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`


 
.. c:var:: normal_two_body_bi_orth_aa_bb


    File : :file:`normal_ordered.irp.f`

    .. code:: fortran

        double precision, allocatable	:: normal_two_body_bi_orth_aa_bb	(n_act_orb,n_act_orb,n_act_orb,n_act_orb)


    Normal ordered two-body sector of the three-body terms for same spin double excitations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`


 
.. c:var:: normal_two_body_bi_orth_ab


    File : :file:`normal_ordered.irp.f`

    .. code:: fortran

        double precision, allocatable	:: normal_two_body_bi_orth_ab	(mo_num,mo_num,mo_num,mo_num)


    Normal ordered two-body sector of the three-body terms for opposite spin double excitations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`


 
.. c:var:: overlap_natorb_tc_eigvec_ao


    File : :file:`tc_natorb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: natorb_tc_reigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: natorb_tc_leigvec_ao	(ao_num,mo_num)
        double precision, allocatable	:: overlap_natorb_tc_eigvec_ao	(mo_num,mo_num)


    EIGENVECTORS OF FOCK MATRIX ON THE AO BASIS and their OVERLAP
    
    THE OVERLAP SHOULD BE THE SAME AS overlap_natorb_tc_eigvec_mo

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_diag_sorted_r_natorb`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: psi_average_norm_contrib_sorted_tc


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc	(psi_det_size,N_states)
        double precision, allocatable	:: psi_average_norm_contrib_sorted_tc	(psi_det_size)
        integer, allocatable	:: psi_det_sorted_tc_order	(psi_det_size)
        double precision, allocatable	:: psi_r_coef_sorted_bi_ortho	(psi_det_size,N_states)
        double precision, allocatable	:: psi_l_coef_sorted_bi_ortho	(psi_det_size,N_states)


    Wave function sorted by determinants contribution to the norm (state-averaged)
    
    psi_det_sorted_tc_order(i) -> k : index in psi_det
    
    psi_r_coef_sorted_bi_ortho : right coefficients corresponding to psi_det_sorted_tc
    
    psi_l_coef_sorted_bi_ortho : left  coefficients corresponding to psi_det_sorted_tc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_average_norm_contrib_tc`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: psi_average_norm_contrib_tc


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_average_norm_contrib_tc	(psi_det_size)


    Contribution of determinants to the state-averaged density.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_sorted_tc`

 
.. c:var:: psi_bitcleft_bilinear_matrix_values


    File : :file:`psi_left_qmc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_bitcleft_bilinear_matrix_values	(N_det,N_states)


    Sparse coefficient matrix if the wave function is expressed in a bilinear form :
     :math:`D_\alpha^\dagger.C.D_\beta` 
    
    Rows are :math:`\alpha`  determinants and columns are :math:`\beta` .
    
    Order refers to psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_l_coef_bi_ortho`


 
.. c:var:: psi_coef_sorted_tc


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc	(psi_det_size,N_states)
        double precision, allocatable	:: psi_average_norm_contrib_sorted_tc	(psi_det_size)
        integer, allocatable	:: psi_det_sorted_tc_order	(psi_det_size)
        double precision, allocatable	:: psi_r_coef_sorted_bi_ortho	(psi_det_size,N_states)
        double precision, allocatable	:: psi_l_coef_sorted_bi_ortho	(psi_det_size,N_states)


    Wave function sorted by determinants contribution to the norm (state-averaged)
    
    psi_det_sorted_tc_order(i) -> k : index in psi_det
    
    psi_r_coef_sorted_bi_ortho : right coefficients corresponding to psi_det_sorted_tc
    
    psi_l_coef_sorted_bi_ortho : left  coefficients corresponding to psi_det_sorted_tc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_average_norm_contrib_tc`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: psi_coef_sorted_tc_bit


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc_bit	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc_bit	(psi_det_size,N_states)


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
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`


 
.. c:var:: psi_det_sorted_tc


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc	(psi_det_size,N_states)
        double precision, allocatable	:: psi_average_norm_contrib_sorted_tc	(psi_det_size)
        integer, allocatable	:: psi_det_sorted_tc_order	(psi_det_size)
        double precision, allocatable	:: psi_r_coef_sorted_bi_ortho	(psi_det_size,N_states)
        double precision, allocatable	:: psi_l_coef_sorted_bi_ortho	(psi_det_size,N_states)


    Wave function sorted by determinants contribution to the norm (state-averaged)
    
    psi_det_sorted_tc_order(i) -> k : index in psi_det
    
    psi_r_coef_sorted_bi_ortho : right coefficients corresponding to psi_det_sorted_tc
    
    psi_l_coef_sorted_bi_ortho : left  coefficients corresponding to psi_det_sorted_tc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_average_norm_contrib_tc`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: psi_det_sorted_tc_bit


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc_bit	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc_bit	(psi_det_size,N_states)


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
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`


 
.. c:var:: psi_det_sorted_tc_left


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc_left	(N_int,2,N_det)
        double precision, allocatable	:: psi_l_coef_sorted_bi_ortho_left	(N_det)


    psi_det_sorted_tc_left : Slater determinants sorted by decreasing value of |LEFTt- coefficients|
    
    psi_r_coef_sorted_bi_ortho_left : LEFT wave function according to psi_det_sorted_tc_left

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: psi_det_sorted_tc_order


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc	(psi_det_size,N_states)
        double precision, allocatable	:: psi_average_norm_contrib_sorted_tc	(psi_det_size)
        integer, allocatable	:: psi_det_sorted_tc_order	(psi_det_size)
        double precision, allocatable	:: psi_r_coef_sorted_bi_ortho	(psi_det_size,N_states)
        double precision, allocatable	:: psi_l_coef_sorted_bi_ortho	(psi_det_size,N_states)


    Wave function sorted by determinants contribution to the norm (state-averaged)
    
    psi_det_sorted_tc_order(i) -> k : index in psi_det
    
    psi_r_coef_sorted_bi_ortho : right coefficients corresponding to psi_det_sorted_tc
    
    psi_l_coef_sorted_bi_ortho : left  coefficients corresponding to psi_det_sorted_tc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_average_norm_contrib_tc`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: psi_det_sorted_tc_right


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc_right	(N_int,2,N_det)
        double precision, allocatable	:: psi_r_coef_sorted_bi_ortho_right	(N_det)


    psi_det_sorted_tc_right : Slater determinants sorted by decreasing value of |right- coefficients|
    
    psi_r_coef_sorted_bi_ortho_right : right wave function according to psi_det_sorted_tc_right

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: psi_l_coef_bi_ortho


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_l_coef_bi_ortho	(psi_det_size,N_states)


    The wave function coefficients. Initialized with Hartree-Fock if the |EZFIO| file
    is empty.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_label`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`nproc`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_size`
       * :c:data:`read_wf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`psi_average_norm_contrib_tc`
       * :c:data:`psi_bitcleft_bilinear_matrix_values`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_det_sorted_tc_left`
       * :c:data:`tc_transition_matrix_mo_beta`

 
.. c:var:: psi_l_coef_sorted_bi_ortho


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc	(psi_det_size,N_states)
        double precision, allocatable	:: psi_average_norm_contrib_sorted_tc	(psi_det_size)
        integer, allocatable	:: psi_det_sorted_tc_order	(psi_det_size)
        double precision, allocatable	:: psi_r_coef_sorted_bi_ortho	(psi_det_size,N_states)
        double precision, allocatable	:: psi_l_coef_sorted_bi_ortho	(psi_det_size,N_states)


    Wave function sorted by determinants contribution to the norm (state-averaged)
    
    psi_det_sorted_tc_order(i) -> k : index in psi_det
    
    psi_r_coef_sorted_bi_ortho : right coefficients corresponding to psi_det_sorted_tc
    
    psi_l_coef_sorted_bi_ortho : left  coefficients corresponding to psi_det_sorted_tc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_average_norm_contrib_tc`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: psi_l_coef_sorted_bi_ortho_left


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc_left	(N_int,2,N_det)
        double precision, allocatable	:: psi_l_coef_sorted_bi_ortho_left	(N_det)


    psi_det_sorted_tc_left : Slater determinants sorted by decreasing value of |LEFTt- coefficients|
    
    psi_r_coef_sorted_bi_ortho_left : LEFT wave function according to psi_det_sorted_tc_left

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: psi_r_coef_bi_ortho


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_r_coef_bi_ortho	(psi_det_size,N_states)


    The wave function coefficients. Initialized with Hartree-Fock if the |EZFIO| file
    is empty.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_label`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`nproc`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_size`
       * :c:data:`read_wf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`psi_average_norm_contrib_tc`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_det_sorted_tc_left`
       * :c:data:`psi_det_sorted_tc_right`
       * :c:data:`reigvec_tc_cisd_sc2_bi_ortho`
       * :c:data:`tc_transition_matrix_mo_beta`

 
.. c:var:: psi_r_coef_sorted_bi_ortho


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_tc	(psi_det_size,N_states)
        double precision, allocatable	:: psi_average_norm_contrib_sorted_tc	(psi_det_size)
        integer, allocatable	:: psi_det_sorted_tc_order	(psi_det_size)
        double precision, allocatable	:: psi_r_coef_sorted_bi_ortho	(psi_det_size,N_states)
        double precision, allocatable	:: psi_l_coef_sorted_bi_ortho	(psi_det_size,N_states)


    Wave function sorted by determinants contribution to the norm (state-averaged)
    
    psi_det_sorted_tc_order(i) -> k : index in psi_det
    
    psi_r_coef_sorted_bi_ortho : right coefficients corresponding to psi_det_sorted_tc
    
    psi_l_coef_sorted_bi_ortho : left  coefficients corresponding to psi_det_sorted_tc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_average_norm_contrib_tc`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: psi_r_coef_sorted_bi_ortho_right


    File : :file:`psi_det_tc_sorted.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_tc_right	(N_int,2,N_det)
        double precision, allocatable	:: psi_r_coef_sorted_bi_ortho_right	(N_det)


    psi_det_sorted_tc_right : Slater determinants sorted by decreasing value of |right- coefficients|
    
    psi_r_coef_sorted_bi_ortho_right : right wave function according to psi_det_sorted_tc_right

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_r_coef_bi_ortho`


 
.. c:var:: ref_tc_energy_1e


    File : :file:`slater_tc_opt_diag.irp.f`

    .. code:: fortran

        double precision	:: ref_tc_energy_tot	
        double precision	:: ref_tc_energy_1e	
        double precision	:: ref_tc_energy_2e	
        double precision	:: ref_tc_energy_3e	


    Various component of the TC energy for the reference "HF" Slater determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`hf_bitmask`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`three_body_h_tc`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`


 
.. c:var:: ref_tc_energy_2e


    File : :file:`slater_tc_opt_diag.irp.f`

    .. code:: fortran

        double precision	:: ref_tc_energy_tot	
        double precision	:: ref_tc_energy_1e	
        double precision	:: ref_tc_energy_2e	
        double precision	:: ref_tc_energy_3e	


    Various component of the TC energy for the reference "HF" Slater determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`hf_bitmask`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`three_body_h_tc`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`


 
.. c:var:: ref_tc_energy_3e


    File : :file:`slater_tc_opt_diag.irp.f`

    .. code:: fortran

        double precision	:: ref_tc_energy_tot	
        double precision	:: ref_tc_energy_1e	
        double precision	:: ref_tc_energy_2e	
        double precision	:: ref_tc_energy_3e	


    Various component of the TC energy for the reference "HF" Slater determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`hf_bitmask`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`three_body_h_tc`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`


 
.. c:var:: ref_tc_energy_tot


    File : :file:`slater_tc_opt_diag.irp.f`

    .. code:: fortran

        double precision	:: ref_tc_energy_tot	
        double precision	:: ref_tc_energy_1e	
        double precision	:: ref_tc_energy_2e	
        double precision	:: ref_tc_energy_3e	


    Various component of the TC energy for the reference "HF" Slater determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`hf_bitmask`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`three_body_h_tc`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`


 
.. c:var:: reigvec_tc_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision	:: norm_ground_left_right_bi_orth	


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`qp_max_mem`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: reigvec_tc_cisd_sc2_bi_ortho


    File : :file:`tc_cisd_sc2_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: reigvec_tc_cisd_sc2_bi_ortho	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_cisd_sc2_bi_ortho	(N_det,N_states)
        double precision, allocatable	:: eigval_tc_cisd_sc2_bi_ortho	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`disk_based_davidson`
       * :c:data:`double_normal_ord`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`qp_max_mem`
       * :c:data:`three_body_h_tc`
       * :c:data:`threshold_nonsym_davidson`


 
.. c:var:: tc_bi_ortho_dipole


    File : :file:`tc_prop.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_bi_ortho_dipole	(3,N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_orth_bipole_x`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`tc_transition_matrix_mo`


 
.. c:var:: tc_mulliken_spin_densities


    File : :file:`spin_mulliken.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_mulliken_spin_densities	(nucl_num,N_states)


    ATOMIC SPIN POPULATION (ALPHA MINUS BETA)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`nucl_num`
       * :c:data:`tc_spin_gross_orbital_product`


 
.. c:var:: tc_spin_dens_right_only


    File : :file:`test_spin_dens.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_spin_dens_right_only	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_r_coef_normalized`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_spin_population`

 
.. c:var:: tc_spin_gross_orbital_product


    File : :file:`spin_mulliken.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_spin_gross_orbital_product	(ao_num,N_states)


    gross orbital product for the spin population

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`tc_spin_population`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_mulliken_spin_densities`
       * :c:data:`tc_spin_population_angular_momentum`

 
.. c:var:: tc_spin_population


    File : :file:`spin_mulliken.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_spin_population	(ao_num,ao_num,N_states)


    spin population on the ao basis :
    tc_spin_population(i,j) = rho_AO(alpha)(i,j) - rho_AO(beta)(i,j) * <AO_i|AO_j>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`n_states`
       * :c:data:`only_spin_tc_right`
       * :c:data:`tc_spin_dens_right_only`
       * :c:data:`tc_spin_transition_matrix_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_spin_gross_orbital_product`

 
.. c:var:: tc_spin_population_angular_momentum


    File : :file:`spin_mulliken.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_spin_population_angular_momentum	(0:ao_l_max,N_states)
        double precision, allocatable	:: tc_spin_population_angular_momentum_per_atom	(0:ao_l_max,nucl_num,N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`nucl_num`
       * :c:data:`tc_spin_gross_orbital_product`


 
.. c:var:: tc_spin_population_angular_momentum_per_atom


    File : :file:`spin_mulliken.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_spin_population_angular_momentum	(0:ao_l_max,N_states)
        double precision, allocatable	:: tc_spin_population_angular_momentum_per_atom	(0:ao_l_max,nucl_num,N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`nucl_num`
       * :c:data:`tc_spin_gross_orbital_product`


 
.. c:var:: tc_spin_transition_matrix_ao


    File : :file:`tc_prop.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_spin_transition_matrix_ao	(ao_num,ao_num,N_states,N_states)


    tc_spin_transition_matrix_ao(p,h,istate,jstate) in the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_states`
       * :c:data:`tc_spin_transition_matrix_mo`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_spin_population`

 
.. c:var:: tc_spin_transition_matrix_mo


    File : :file:`tc_prop.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_spin_transition_matrix_mo	(mo_num,mo_num,N_states,N_states)


    tc_spin_transition_matrix_mo = tc_transition_matrix_mo_alpha - tc_transition_matrix_mo_beta
    
    where <Chi_istate| and |Phi_jstate> are the left/right eigenvectors on a bi-ortho basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`tc_transition_matrix_mo_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_spin_transition_matrix_ao`

 
.. c:var:: tc_transition_matrix_ao


    File : :file:`tc_prop.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_transition_matrix_ao	(ao_num,ao_num,N_states,N_states)


    tc_transition_matrix(p,h,istate,jstate) in the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_states`
       * :c:data:`tc_transition_matrix_mo`


 
.. c:var:: tc_transition_matrix_mo


    File : :file:`tc_prop.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_transition_matrix_mo	(mo_num,mo_num,N_states,N_states)


    tc_transition_matrix_mo(p,h,istate,jstate) = \sum_{sigma=alpha,beta} <Chi_istate| a^\dagger_p,sigma a_h,sigma |Phi_jstate>
    
    where <Chi_istate| and |Phi_jstate> are the left/right eigenvectors on a bi-ortho basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`tc_transition_matrix_mo_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`natorb_tc_reigvec_mo`
       * :c:data:`tc_bi_ortho_dipole`
       * :c:data:`tc_transition_matrix_ao`

 
.. c:var:: tc_transition_matrix_mo_alpha


    File : :file:`tc_prop.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_transition_matrix_mo_beta	(mo_num,mo_num,N_states,N_states)
        double precision, allocatable	:: tc_transition_matrix_mo_alpha	(mo_num,mo_num,N_states,N_states)


    tc_transition_matrix_mo_alpha(p,h,istate,jstate) = <Chi_istate| a^\dagger_p,alpha a_h,alpha |Phi_jstate>
    
    tc_transition_matrix_mo_beta(p,h,istate,jstate) = <Chi_istate| a^\dagger_p,beta a_h,beta |Phi_jstate>
    
    where <Chi_istate| and |Phi_jstate> are the left/right eigenvectors on a bi-ortho basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_spin_transition_matrix_mo`
       * :c:data:`tc_transition_matrix_mo`

 
.. c:var:: tc_transition_matrix_mo_beta


    File : :file:`tc_prop.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_transition_matrix_mo_beta	(mo_num,mo_num,N_states,N_states)
        double precision, allocatable	:: tc_transition_matrix_mo_alpha	(mo_num,mo_num,N_states,N_states)


    tc_transition_matrix_mo_alpha(p,h,istate,jstate) = <Chi_istate| a^\dagger_p,alpha a_h,alpha |Phi_jstate>
    
    tc_transition_matrix_mo_beta(p,h,istate,jstate) = <Chi_istate| a^\dagger_p,beta a_h,beta |Phi_jstate>
    
    where <Chi_istate| and |Phi_jstate> are the left/right eigenvectors on a bi-ortho basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_spin_transition_matrix_mo`
       * :c:data:`tc_transition_matrix_mo`

 
.. c:var:: three_e_diag_parrallel_spin_prov


    File : :file:`symmetrized_3_e_int_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_diag_parrallel_spin_prov	(mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator ON A BI ORTHONORMAL BASIS
    
    three_e_diag_parrallel_spin_prov(m,j,i) = All combinations of the form <mji|-L|mji> for same spin matrix elements
    
    notice the -1 sign: in this way three_e_diag_parrallel_spin_prov can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`three_e_3_idx_cycle_1_bi_ort`
       * :c:data:`three_e_3_idx_cycle_2_bi_ort`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`
       * :c:data:`three_e_3_idx_exch13_bi_ort`
       * :c:data:`three_e_3_idx_exch23_bi_ort`


 
.. c:var:: three_e_double_parrallel_spin_prov


    File : :file:`symmetrized_3_e_int_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_double_parrallel_spin_prov	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_double_parrallel_spin_prov(m,l,j,k,i) = <mlk|-L|mji> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`three_e_5_idx_cycle_1_bi_ort`
       * :c:data:`three_e_5_idx_cycle_2_bi_ort`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_5_idx_exch13_bi_ort`
       * :c:data:`three_e_5_idx_exch23_bi_ort`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eff_2_e_from_3_e_aa`
       * :c:data:`eff_2_e_from_3_e_bb`

 
.. c:var:: three_e_single_parrallel_spin_prov


    File : :file:`symmetrized_3_e_int_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_single_parrallel_spin_prov	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_single_parrallel_spin_prov(m,j,k,i) = All combination of <mjk|-L|mji> for same spin matrix elements
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`three_e_4_idx_cycle_1_bi_ort`
       * :c:data:`three_e_4_idx_cycle_2_bi_ort`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_exch12_bi_ort`
       * :c:data:`three_e_4_idx_exch13_bi_ort`
       * :c:data:`three_e_4_idx_exch23_bi_ort`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: a_tc_operator:


    File : :file:`slater_tc_opt_diag.irp.f`

    .. code:: fortran

        subroutine a_tc_operator(iorb,ispin,key,hmono,htwoe,hthree,Nint,na,nb)


    Routine that computes one- and two-body energy corresponding
    
    to the REMOVAL of an electron in an orbital 'iorb' of spin 'ispin'
    
    onto a determinant 'key'.
    
    in output, the determinant key is changed by the REMOVAL of that electron
    
    and the quantities hmono,htwoe,hthree are INCREMENTED

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e_jj`
       * :c:data:`three_body_h_tc`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`
       * :c:data:`three_e_3_idx_exch23_bi_ort`
       * :c:data:`three_e_diag_parrallel_spin_prov`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: a_tc_operator_no_3e:


    File : :file:`slater_tc_opt_diag.irp.f`

    .. code:: fortran

        subroutine a_tc_operator_no_3e(iorb,ispin,key,hmono,htwoe,Nint,na,nb)


    Routine that computes one- and two-body energy corresponding
    
    to the REMOVAL of an electron in an orbital 'iorb' of spin 'ispin'
    
    onto a determinant 'key'.
    
    in output, the determinant key is changed by the REMOVAL of that electron
    
    and the quantities hmono,htwoe are INCREMENTED

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e_jj`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho_no_3e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: ac_tc_operator:


    File : :file:`slater_tc_opt_diag.irp.f`

    .. code:: fortran

        subroutine ac_tc_operator(iorb,ispin,key,hmono,htwoe,hthree,Nint,na,nb)


    Routine that computes one- and two-body energy corresponding
    
    to the ADDITION of an electron in an orbital 'iorb' of spin 'ispin'
    
    onto a determinant 'key'.
    
    in output, the determinant key is changed by the ADDITION of that electron
    
    and the quantities hmono,htwoe,hthree are INCREMENTED

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e_jj`
       * :c:data:`mo_num`
       * :c:data:`three_body_h_tc`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`
       * :c:data:`three_e_3_idx_exch23_bi_ort`
       * :c:data:`three_e_diag_parrallel_spin_prov`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: ac_tc_operator_no_3e:


    File : :file:`slater_tc_opt_diag.irp.f`

    .. code:: fortran

        subroutine ac_tc_operator_no_3e(iorb,ispin,key,hmono,htwoe,Nint,na,nb)


    Routine that computes one- and two-body energy corresponding
    
    to the ADDITION of an electron in an orbital 'iorb' of spin 'ispin'
    
    onto a determinant 'key'.
    
    in output, the determinant key is changed by the ADDITION of that electron
    
    and the quantities hmono,htwoe are INCREMENTED

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e_jj`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho_no_3e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: bi_normalize:


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        subroutine bi_normalize(u_l,u_r,n,ld,nstates)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`

 
.. c:function:: compute_deltamu_right:


    File : :file:`compute_deltamu_right.irp.f`

    .. code:: fortran

        subroutine compute_deltamu_right



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`n_int`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`delta_right`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

 
.. c:function:: delta_right:


    File : :file:`compute_deltamu_right.irp.f`

    .. code:: fortran

        subroutine delta_right()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`j1b_pen`
       * :c:data:`j1b_type`
       * :c:data:`mu_erf`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det`
       * :c:data:`psi_r_coef_bi_ortho`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_deltamu_right`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dset_order`
       * :c:func:`get_delta_bitc_right`

 
.. c:function:: diag_hmat_bi_ortho:


    File : :file:`h_biortho.irp.f`

    .. code:: fortran

        subroutine diag_hmat_bi_ortho(Nint, key_i, hmono, htwoe)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_coul_e`
       * :c:data:`mo_bi_ortho_one_e`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`hmat_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: diag_htilde_mu_mat_bi_ortho:


    File : :file:`slater_tc.irp.f`

    .. code:: fortran

        subroutine diag_htilde_mu_mat_bi_ortho(Nint, key_i, hmono, htwoe, htot)


    diagonal element of htilde ONLY FOR ONE- AND TWO-BODY TERMS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho`
       * :c:data:`ref_tc_energy_tot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: diag_htilde_mu_mat_fock_bi_ortho:


    File : :file:`slater_tc_opt_diag.irp.f`

    .. code:: fortran

        subroutine diag_htilde_mu_mat_fock_bi_ortho(Nint, det_in, hmono, htwoe, hthree, htot)


    Computes $\langle i|H|i \rangle$.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`nuclear_repulsion`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_tc_energy_tot`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp_work_1`
       * :c:func:`h_tc_u_0_nstates_openmp_work_2`
       * :c:func:`h_tc_u_0_nstates_openmp_work_3`
       * :c:func:`h_tc_u_0_nstates_openmp_work_4`
       * :c:func:`h_tc_u_0_nstates_openmp_work_n_int`
       * :c:func:`htilde_mu_mat_opt_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`a_tc_operator`
       * :c:func:`ac_tc_operator`
       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: diag_htilde_mu_mat_fock_bi_ortho_no_3e:


    File : :file:`slater_tc_opt_diag.irp.f`

    .. code:: fortran

        subroutine diag_htilde_mu_mat_fock_bi_ortho_no_3e(Nint, det_in,htot)


    Computes $\langle i|H|i \rangle$. WITHOUT ANY CONTRIBUTIONS FROM 3E TERMS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_tc_energy_tot`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_opt_bi_ortho_no_3e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`a_tc_operator_no_3e`
       * :c:func:`ac_tc_operator_no_3e`
       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: diag_htilde_three_body_ints_bi_ort:


    File : :file:`slater_tc_3e.irp.f`

    .. code:: fortran

        subroutine diag_htilde_three_body_ints_bi_ort(Nint, key_i, hthree)


    diagonal element of htilde ONLY FOR THREE-BODY TERMS WITH BI ORTHONORMAL ORBITALS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho`
       * :c:data:`ref_tc_energy_tot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: diagonalize_ci_tc:


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        subroutine diagonalize_CI_tc


    Replace the coefficients of the |CI| states by the coefficients of the
    eigenstates of the |CI| matrix.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`

 
.. c:function:: direct_sym_3_e_int:


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        double precision function direct_sym_3_e_int(n,l,k,m,j,i)


    returns all good combinations of permutations of integrals with the good signs
    
    for a given (k^dagger l^dagger n^dagger m j i)  <nlk|L|mji> when all indices have the same spins

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`

 
.. c:function:: double_hmat_bi_ortho:


    File : :file:`h_biortho.irp.f`

    .. code:: fortran

        subroutine double_hmat_bi_ortho(Nint, key_j, key_i, hmono, htwoe)


    
    < key_j | H | key_i> for double excitation
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_coul_e`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`hmat_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`decode_exc`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`

 
.. c:function:: double_htilde_mu_mat_bi_ortho:


    File : :file:`slater_tc.irp.f`

    .. code:: fortran

        subroutine double_htilde_mu_mat_bi_ortho(Nint, key_j, key_i, hmono, htwoe, htot)


     <key_j | H_tilde | key_i> for double excitation  ONLY FOR ONE- AND TWO-BODY TERMS
    !
    ! WARNING !!
    
     Non hermitian !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`normal_two_body_bi_orth`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`decode_exc`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`

 
.. c:function:: double_htilde_mu_mat_fock_bi_ortho:


    File : :file:`slater_tc_opt_double.irp.f`

    .. code:: fortran

        subroutine double_htilde_mu_mat_fock_bi_ortho(Nint, key_j, key_i, hmono, htwoe, hthree, htot)


     <key_j | H_tilde | key_i> for double excitation  ONLY FOR ONE- AND TWO-BODY TERMS
    !
    ! WARNING !!
    
     Non hermitian !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`elec_num`
       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`n_int`
       * :c:data:`normal_two_body_bi_orth`
       * :c:data:`ref_bitmask`
       * :c:data:`three_body_h_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_opt_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`decode_exc`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`
       * :c:func:`three_comp_two_e_elem`

 
.. c:function:: double_htilde_mu_mat_fock_bi_ortho_no_3e:


    File : :file:`slater_tc_opt_double.irp.f`

    .. code:: fortran

        subroutine double_htilde_mu_mat_fock_bi_ortho_no_3e(Nint, key_j, key_i, htot)


     <key_j | H_tilde | key_i> for double excitation  ONLY FOR ONE- AND TWO-BODY TERMS
    !
    ! WARNING !!
    
     Non hermitian !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_opt_bi_ortho_no_3e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`decode_exc`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`

 
.. c:function:: double_htilde_three_body_ints_bi_ort:


    File : :file:`slater_tc_3e.irp.f`

    .. code:: fortran

        subroutine double_htilde_three_body_ints_bi_ort(Nint, key_j, key_i, hthree)


     <key_j | H_tilde | key_i> for double excitation ONLY FOR THREE-BODY TERMS  WITH BI ORTHONORMAL ORBITALS
    !
    ! WARNING !!
    
     Non hermitian !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_5_idx_exch13_bi_ort`
       * :c:data:`three_e_5_idx_exch23_bi_ort`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`decode_exc`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`

 
.. c:function:: fock_a_tc_operator:


    File : :file:`slater_tc_opt_single.irp.f`

    .. code:: fortran

        subroutine fock_a_tc_operator(iorb,ispin,key, h_fock,p_fock, ispin_fock,hthree,Nint,na,nb)


    Routine that computes the contribution to the three-electron part of the Fock operator
    
    a^dagger_{p_fock} a_{h_fock} of spin ispin_fock
    
    on top of a determinant 'key' on which you REMOVE an electron of spin ispin in orbital iorb
    
    in output, the determinant key is changed by the REMOVAL of that electron
    
    the output hthree is INCREMENTED

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_exch12_bi_ort`
       * :c:data:`three_e_4_idx_exch13_bi_ort`
       * :c:data:`three_e_4_idx_exch23_bi_ort`
       * :c:data:`three_e_single_parrallel_spin_prov`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`three_comp_fock_elem`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: fock_ac_tc_operator:


    File : :file:`slater_tc_opt_single.irp.f`

    .. code:: fortran

        subroutine fock_ac_tc_operator(iorb,ispin,key, h_fock,p_fock, ispin_fock,hthree,Nint,na,nb)


    Routine that computes the contribution to the three-electron part of the Fock operator
    
    a^dagger_{p_fock} a_{h_fock} of spin ispin_fock
    
    on top of a determinant 'key' on which you ADD an electron of spin ispin in orbital iorb
    
    in output, the determinant key is changed by the ADDITION of that electron
    
    the output hthree is INCREMENTED

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_exch12_bi_ort`
       * :c:data:`three_e_4_idx_exch13_bi_ort`
       * :c:data:`three_e_4_idx_exch23_bi_ort`
       * :c:data:`three_e_single_parrallel_spin_prov`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`three_comp_fock_elem`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: get_cisd_sc2_dressing:


    File : :file:`tc_cisd_sc2_utils.irp.f`

    .. code:: fortran

        subroutine get_cisd_sc2_dressing(dets,e_corr_dets,ndet,dressing_dets)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`hf_bitmask`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`reigvec_tc_cisd_sc2_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`decode_exc`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`

 
.. c:function:: get_delta_bitc_right:


    File : :file:`dressing_vectors_lr.irp.f`

    .. code:: fortran

        subroutine get_delta_bitc_right(psidet, psicoef, ndet, Nint, delta)


    
    delta(I) = < I_left | H_TC - H | Psi_right >
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`delta_right`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`hmat_bi_ortho`
       * :c:func:`htilde_mu_mat_bi_ortho`

 
.. c:function:: get_h_bitc_right:


    File : :file:`dressing_vectors_lr.irp.f`

    .. code:: fortran

        subroutine get_h_bitc_right(psidet, psicoef, ndet, Nint, delta)


    
    delta(I) = < I_left | H | Psi_right >
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`hmat_bi_ortho`

 
.. c:function:: get_htc_bitc_right:


    File : :file:`dressing_vectors_lr.irp.f`

    .. code:: fortran

        subroutine get_htc_bitc_right(psidet, psicoef, ndet, Nint, delta)


    
    delta(I) = < I_left | H_TC | Psi_right >
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho`

 
.. c:function:: get_single_excitation_from_fock_tc:


    File : :file:`slater_tc_opt_single.irp.f`

    .. code:: fortran

        subroutine get_single_excitation_from_fock_tc(key_i,key_j,h,p,spin,phase,hmono,htwoe,hthree,htot)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_op_2_e_tc_closed_shell`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`ref_closed_shell_bitmask`
       * :c:data:`tc_2e_3idx_coulomb_integrals`
       * :c:data:`three_body_h_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`single_htilde_mu_mat_fock_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`three_comp_fock_elem`

 
.. c:function:: get_single_excitation_from_fock_tc_no_3e:


    File : :file:`slater_tc_opt_single.irp.f`

    .. code:: fortran

        subroutine get_single_excitation_from_fock_tc_no_3e(key_i,key_j,h,p,spin,phase,hmono,htwoe,htot)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_op_2_e_tc_closed_shell`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`ref_closed_shell_bitmask`
       * :c:data:`tc_2e_3idx_coulomb_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`single_htilde_mu_mat_fock_bi_ortho_no_3e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: give_aaa_contraction:


    File : :file:`normal_ordered.irp.f`

    .. code:: fortran

        subroutine give_aaa_contraction(Nint, h1, h2, p1, p2, Ne, occ, hthree)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`normal_two_body_bi_orth`
       * :c:data:`normal_two_body_bi_orth_aa_bb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`

 
.. c:function:: give_aab_contraction:


    File : :file:`normal_ordered.irp.f`

    .. code:: fortran

        subroutine give_aab_contraction(Nint, h1, h2, p1, p2, Ne, occ, hthree)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`normal_two_body_bi_orth`
       * :c:data:`normal_two_body_bi_orth_aa_bb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`

 
.. c:function:: give_aba_contraction:


    File : :file:`normal_ordered.irp.f`

    .. code:: fortran

        subroutine give_aba_contraction(Nint, h1, h2, p1, p2, Ne, occ, hthree)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`normal_two_body_bi_orth`
       * :c:data:`normal_two_body_bi_orth_ab`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`

 
.. c:function:: give_all_perm_for_three_e:


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        subroutine give_all_perm_for_three_e(n,l,k,m,j,i,idx_list,phase)


    returns all the list of permutting indices for the antimmetrization of
    
    (k^dagger l^dagger n^dagger m j i)  <nlk|L|mji>   when all indices have the same spins
    
    idx_list(:,i) == list of the 6 indices corresponding the permutation "i"
    
    phase(i) == phase of the permutation "i"
    
    there are in total 6 permutations with different indices

 
.. c:function:: give_contrib_for_aaaa:


    File : :file:`slater_tc_opt_double.irp.f`

    .. code:: fortran

        subroutine give_contrib_for_aaaa(h1,h2,p1,p2,occ,Ne,contrib)


    gives the contribution for a double excitation (h1,p1)_alpha (h2,p2)_alpha
    
    on top of a determinant whose occupied orbitals is in (occ, Ne)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_double_parrallel_spin_prov`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`eff_2_e_from_3_e_aa`

 
.. c:function:: give_contrib_for_abab:


    File : :file:`slater_tc_opt_double.irp.f`

    .. code:: fortran

        subroutine give_contrib_for_abab(h1,h2,p1,p2,occ,Ne,contrib)


    gives the contribution for a double excitation (h1,p1)_alpha (h2,p2)_beta
    
    on top of a determinant whose occupied orbitals is in (occ, Ne)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch13_bi_ort`
       * :c:data:`three_e_5_idx_exch23_bi_ort`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`eff_2_e_from_3_e_ab`

 
.. c:function:: give_contrib_for_bbbb:


    File : :file:`slater_tc_opt_double.irp.f`

    .. code:: fortran

        subroutine give_contrib_for_bbbb(h1,h2,p1,p2,occ,Ne,contrib)


    gives the contribution for a double excitation (h1,p1)_beta (h2,p2)_beta
    
    on top of a determinant whose occupied orbitals is in (occ, Ne)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_double_parrallel_spin_prov`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`eff_2_e_from_3_e_bb`

 
.. c:function:: h_tc_dagger_u_0_opt:


    File : :file:`u0_h_u0.irp.f`

    .. code:: fortran

        subroutine H_tc_dagger_u_0_opt(v_0,u_0,N_st,sze)


    Computes $v_0 = H | u_0\rangle$.
    
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep are used in ZMQ parallelization.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp`

 
.. c:function:: h_tc_u_0_nstates_openmp:


    File : :file:`u0_h_u0.irp.f`

    .. code:: fortran

        subroutine H_tc_u_0_nstates_openmp(v_0,u_0,N_st,sze, do_right)


    Computes $v_0 = H | u_0\rangle$.
    
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep are used in ZMQ parallelization.
    
    if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_dagger_u_0_opt`
       * :c:func:`h_tc_u_0_opt`
       * :c:func:`test_h_u0`
       * :c:func:`u_0_h_tc_u_0`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dset_order`
       * :c:func:`dtranspose`
       * :c:func:`h_tc_u_0_nstates_openmp_work`

 
.. c:function:: h_tc_u_0_nstates_openmp_work:


    File : :file:`u0_h_u0.irp.f`

    .. code:: fortran

        subroutine H_tc_u_0_nstates_openmp_work(v_t,u_t,N_st,sze,istart,iend,ishift,istep, do_right)


    Computes $v_t = H | u_t\rangle$
    
    Default should be 1,N_det,0,1
    
    if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp_work_1`
       * :c:func:`h_tc_u_0_nstates_openmp_work_2`
       * :c:func:`h_tc_u_0_nstates_openmp_work_3`
       * :c:func:`h_tc_u_0_nstates_openmp_work_4`
       * :c:func:`h_tc_u_0_nstates_openmp_work_n_int`

 
.. c:function:: h_tc_u_0_nstates_openmp_work_1:


    File : :file:`u0_h_u0.irp.f_template_768`

    .. code:: fortran

        subroutine H_tc_u_0_nstates_openmp_work_1(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)


    Computes $v_t = H | u_t \rangle$
    
    Default should be 1,N_det,0,1
    
    if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho`
       * :c:func:`get_all_spin_singles_1`
       * :c:func:`get_all_spin_singles_and_doubles_1`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_tot`

 
.. c:function:: h_tc_u_0_nstates_openmp_work_2:


    File : :file:`u0_h_u0.irp.f_template_768`

    .. code:: fortran

        subroutine H_tc_u_0_nstates_openmp_work_2(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)


    Computes $v_t = H | u_t \rangle$
    
    Default should be 1,N_det,0,1
    
    if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho`
       * :c:func:`get_all_spin_singles_2`
       * :c:func:`get_all_spin_singles_and_doubles_2`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_tot`

 
.. c:function:: h_tc_u_0_nstates_openmp_work_3:


    File : :file:`u0_h_u0.irp.f_template_768`

    .. code:: fortran

        subroutine H_tc_u_0_nstates_openmp_work_3(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)


    Computes $v_t = H | u_t \rangle$
    
    Default should be 1,N_det,0,1
    
    if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho`
       * :c:func:`get_all_spin_singles_3`
       * :c:func:`get_all_spin_singles_and_doubles_3`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_tot`

 
.. c:function:: h_tc_u_0_nstates_openmp_work_4:


    File : :file:`u0_h_u0.irp.f_template_768`

    .. code:: fortran

        subroutine H_tc_u_0_nstates_openmp_work_4(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)


    Computes $v_t = H | u_t \rangle$
    
    Default should be 1,N_det,0,1
    
    if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho`
       * :c:func:`get_all_spin_singles_4`
       * :c:func:`get_all_spin_singles_and_doubles_4`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_tot`

 
.. c:function:: h_tc_u_0_nstates_openmp_work_n_int:


    File : :file:`u0_h_u0.irp.f_template_768`

    .. code:: fortran

        subroutine H_tc_u_0_nstates_openmp_work_N_int(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)


    Computes $v_t = H | u_t \rangle$
    
    Default should be 1,N_det,0,1
    
    if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho`
       * :c:func:`get_all_spin_singles_and_doubles_n_int`
       * :c:func:`get_all_spin_singles_n_int`
       * :c:func:`htilde_mu_mat_opt_bi_ortho_tot`

 
.. c:function:: h_tc_u_0_opt:


    File : :file:`u0_h_u0.irp.f`

    .. code:: fortran

        subroutine H_tc_u_0_opt(v_0,u_0,N_st,sze)


    Computes $v_0 = H | u_0\rangle$.
    
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep are used in ZMQ parallelization.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp`

 
.. c:function:: hmat_bi_ortho:


    File : :file:`h_biortho.irp.f`

    .. code:: fortran

        subroutine hmat_bi_ortho(key_j, key_i, Nint, hmono, htwoe, htot)


    
    < key_j | H | key_i > where | key_j > is developed on the LEFT basis and | key_i > is developed on the RIGHT basis
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nuclear_repulsion`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_delta_bitc_right`
       * :c:func:`get_h_bitc_right`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diag_hmat_bi_ortho`
       * :c:func:`double_hmat_bi_ortho`
       * :c:func:`get_excitation_degree`
       * :c:func:`single_hmat_bi_ortho`

 
.. c:function:: htc_bi_ortho_calc_tdav:


    File : :file:`h_tc_bi_ortho_psi.irp.f`

    .. code:: fortran

        subroutine htc_bi_ortho_calc_tdav(v, u, N_st, sze)


    Application of H_TC on a vector
    
    v(i,istate) = \sum_j u(j,istate) H_TC(i,j), with:
      H_TC(i,j) = < Di | H_TC | Dj >
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_h_u0`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho_tot`

 
.. c:function:: htcdag_bi_ortho_calc_tdav:


    File : :file:`h_tc_bi_ortho_psi.irp.f`

    .. code:: fortran

        subroutine htcdag_bi_ortho_calc_tdav(v, u, N_st, sze)


    Application of (H_TC)^dagger on a vector
    
    v(i,istate) = \sum_j u(j,istate) H_TC(j,i), with:
      H_TC(i,j) = < Di | H_TC | Dj >
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_h_u0`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho_tot`

 
.. c:function:: htilde_mu_mat_bi_ortho:


    File : :file:`slater_tc.irp.f`

    .. code:: fortran

        subroutine htilde_mu_mat_bi_ortho(key_j, key_i, Nint, hmono, htwoe, hthree, htot)


    
     <key_j | H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
    !
     Returns the detail of the matrix element in terms of single, two and three electron contribution.
    ! WARNING !!
    
     Non hermitian !!
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`nuclear_repulsion`
       * :c:data:`three_body_h_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_pt1_bi_ortho`
       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_pt2_tc_bi_orth`
       * :c:data:`e_tc_left_right`
       * :c:data:`e_tilde_00`
       * :c:data:`e_tilde_bi_orth_00`
       * :c:func:`get_delta_bitc_right`
       * :c:func:`get_htc_bitc_right`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:func:`htilde_mu_mat_bi_ortho_tot`
       * :c:func:`main`
       * :c:data:`reigvec_tc_cisd_sc2_bi_ortho`
       * :c:func:`routine_3`
       * :c:func:`routine_tot`
       * :c:func:`test`
       * :c:func:`test_3e`
       * :c:func:`test_slater_tc_opt`
       * :c:func:`timing_diag`
       * :c:func:`timing_double`
       * :c:func:`timing_single`
       * :c:func:`timing_tot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_bi_ortho`
       * :c:func:`diag_htilde_three_body_ints_bi_ort`
       * :c:func:`double_htilde_mu_mat_bi_ortho`
       * :c:func:`double_htilde_three_body_ints_bi_ort`
       * :c:func:`get_excitation_degree`
       * :c:func:`single_htilde_mu_mat_bi_ortho`
       * :c:func:`single_htilde_three_body_ints_bi_ort`

 
.. c:function:: htilde_mu_mat_bi_ortho_tot:


    File : :file:`slater_tc.irp.f`

    .. code:: fortran

        subroutine htilde_mu_mat_bi_ortho_tot(key_j, key_i, Nint, htot)


     <key_j | H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
    !
    ! WARNING !!
    
     Non hermitian !!

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:func:`htc_bi_ortho_calc_tdav`
       * :c:func:`htcdag_bi_ortho_calc_tdav`
       * :c:data:`reigvec_tc_cisd_sc2_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_degree`
       * :c:func:`htilde_mu_mat_bi_ortho`

 
.. c:function:: htilde_mu_mat_opt_bi_ortho:


    File : :file:`slater_tc_opt.irp.f`

    .. code:: fortran

        subroutine htilde_mu_mat_opt_bi_ortho(key_j, key_i, Nint, hmono, htwoe, hthree, htot)


    
     <key_j | H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
    !
     Returns the detail of the matrix element in terms of single, two and three electron contribution.
    ! WARNING !!
    
     Non hermitian !!
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nuclear_repulsion`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_opt_bi_ortho_tot`
       * :c:func:`test_slater_tc_opt`
       * :c:func:`timing_diag`
       * :c:func:`timing_double`
       * :c:func:`timing_single`
       * :c:func:`timing_tot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho`
       * :c:func:`double_htilde_mu_mat_fock_bi_ortho`
       * :c:func:`get_excitation_degree`
       * :c:func:`single_htilde_mu_mat_fock_bi_ortho`

 
.. c:function:: htilde_mu_mat_opt_bi_ortho_no_3e:


    File : :file:`slater_tc_opt.irp.f`

    .. code:: fortran

        subroutine htilde_mu_mat_opt_bi_ortho_no_3e(key_j, key_i, Nint, htot)


    
     <key_j | H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
    !
     Returns the detail of the matrix element WITHOUT ANY CONTRIBUTION FROM THE THREE ELECTRON TERMS
    ! WARNING !!
    
     Non hermitian !!
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nuclear_repulsion`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diag_htilde_mu_mat_fock_bi_ortho_no_3e`
       * :c:func:`double_htilde_mu_mat_fock_bi_ortho_no_3e`
       * :c:func:`get_excitation_degree`
       * :c:func:`single_htilde_mu_mat_fock_bi_ortho_no_3e`

 
.. c:function:: htilde_mu_mat_opt_bi_ortho_tot:


    File : :file:`slater_tc_opt.irp.f`

    .. code:: fortran

        subroutine htilde_mu_mat_opt_bi_ortho_tot(key_j, key_i, Nint, htot)


    
     <key_j | H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
    !
     Returns the total matrix element
    ! WARNING !!
    
     Non hermitian !!
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp_work_1`
       * :c:func:`h_tc_u_0_nstates_openmp_work_2`
       * :c:func:`h_tc_u_0_nstates_openmp_work_3`
       * :c:func:`h_tc_u_0_nstates_openmp_work_4`
       * :c:func:`h_tc_u_0_nstates_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_opt_bi_ortho`

 
.. c:function:: main:


    File : :file:`tc_som.irp.f`

    .. code:: fortran

        subroutine main()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_fit_param`
       * :c:func:`tc_petermann_factor`
       * :c:func:`tc_som`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_degree`
       * :c:func:`htilde_mu_mat_bi_ortho`

 
.. c:function:: print_energy1:


    File : :file:`print_he_tc_energy.irp.f`

    .. code:: fortran

        subroutine print_energy1()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_non_hermit_term_chemist`
       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_overlap`
       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_he_tc_energy`

 
.. c:function:: print_he_tc_energy:


    File : :file:`print_he_tc_energy.irp.f`

    .. code:: fortran

        subroutine print_he_tc_energy



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_energy1`
       * :c:func:`print_overlap`

 
.. c:function:: print_overlap:


    File : :file:`print_he_tc_energy.irp.f`

    .. code:: fortran

        subroutine print_overlap()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_he_tc_energy`

 
.. c:function:: print_tc_bi_ortho:


    File : :file:`print_tc_wf.irp.f`

    .. code:: fortran

        subroutine print_tc_bi_ortho


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`write_l_r_wf`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

 
.. c:function:: provide_all_three_ints_bi_ortho:


    File : :file:`slater_tc_3e.irp.f`

    .. code:: fortran

        subroutine provide_all_three_ints_bi_ortho


    routine that provides all necessary three-electron integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`double_normal_ord`
       * :c:data:`normal_two_body_bi_orth`
       * :c:data:`three_body_h_tc`
       * :c:data:`three_e_3_idx_cycle_1_bi_ort`
       * :c:data:`three_e_3_idx_cycle_2_bi_ort`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`
       * :c:data:`three_e_3_idx_exch13_bi_ort`
       * :c:data:`three_e_3_idx_exch23_bi_ort`
       * :c:data:`three_e_4_idx_cycle_1_bi_ort`
       * :c:data:`three_e_4_idx_cycle_2_bi_ort`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_exch12_bi_ort`
       * :c:data:`three_e_4_idx_exch13_bi_ort`
       * :c:data:`three_e_4_idx_exch23_bi_ort`
       * :c:data:`three_e_5_idx_cycle_1_bi_ort`
       * :c:data:`three_e_5_idx_cycle_2_bi_ort`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_5_idx_exch13_bi_ort`
       * :c:data:`three_e_5_idx_exch23_bi_ort`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_normal_order`

 
.. c:function:: routine:


    File : :file:`test_natorb.irp.f`

    .. code:: fortran

        subroutine routine



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_mo`
       * :c:data:`mo_num`
       * :c:data:`one_e_dm_mo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_dav`
       * :c:func:`test_natorb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diag_mat_per_fock_degen`
       * :c:func:`non_hrmt_bieig`

 
.. c:function:: routine_3:


    File : :file:`test_tc_fock.irp.f`

    .. code:: fortran

        subroutine routine_3()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`
       * :c:func:`do_single_excitation`
       * :c:func:`htilde_mu_mat_bi_ortho`

 
.. c:function:: routine_diag:


    File : :file:`tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine routine_diag



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_pt2_tc_bi_orth`
       * :c:data:`e_tc_left_right`
       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`n_det`

 
.. c:function:: routine_save_left_right_bi_ortho:


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        subroutine routine_save_left_right_bi_ortho



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_wavefunction_general_unormalized`

 
.. c:function:: routine_save_right_bi_ortho:


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        subroutine routine_save_right_bi_ortho



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`save_tc_bi_ortho_wavefunction`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_wavefunction_general_unormalized`

 
.. c:function:: routine_test:


    File : :file:`select_dets_bi_ortho.irp.f`

    .. code:: fortran

        subroutine routine_test



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`coef_pt1_bi_ortho`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsort`
       * :c:func:`get_excitation_degree`
       * :c:func:`save_wavefunction_general`

 
.. c:function:: routine_tot:


    File : :file:`test_tc_fock.irp.f`

    .. code:: fortran

        subroutine routine_tot()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_tc_mo_alpha`
       * :c:data:`fock_matrix_tc_mo_beta`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_tc_fock`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`
       * :c:func:`do_single_excitation`
       * :c:func:`htilde_mu_mat_bi_ortho`

 
.. c:function:: save_ref_determinant_nstates_1:


    File : :file:`save_tc_bi_ortho_nat.irp.f`

    .. code:: fortran

        subroutine save_ref_determinant_nstates_1



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`ref_bitmask`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`save_tc_natorb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_wavefunction_general`

 
.. c:function:: save_tc_bi_ortho_wavefunction:


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        subroutine save_tc_bi_ortho_wavefunction



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`save_sorted_tc_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`routine_save_right_bi_ortho`
       * :c:func:`save_tc_wavefunction_general`

 
.. c:function:: save_tc_natorb:


    File : :file:`save_tc_bi_ortho_nat.irp.f`

    .. code:: fortran

        subroutine save_tc_natorb



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`natorb_tc_reigvec_ao`
       * :c:data:`natorb_tc_reigvec_ao`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`tc_natorb_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_bi_ortho_mos_mo_l_coef`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_r_coef`
       * :c:func:`ezfio_set_determinants_read_wf`
       * :c:func:`save_ref_determinant_nstates_1`

 
.. c:function:: save_tc_wavefunction_general:


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        subroutine save_tc_wavefunction_general(ndet,nstates,psidet,sze,dim_psicoef,psilcoef,psircoef)


    Save the wave function into the |EZFIO| file

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_label`
       * :c:data:`mpi_master`
       * :c:data:`n_det_qp_edit`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`save_tc_bi_ortho_wavefunction`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_determinants_bit_kind`
       * :c:func:`ezfio_set_determinants_mo_label`
       * :c:func:`ezfio_set_determinants_n_det`
       * :c:func:`ezfio_set_determinants_n_det_qp_edit`
       * :c:func:`ezfio_set_determinants_n_int`
       * :c:func:`ezfio_set_determinants_n_states`
       * :c:func:`ezfio_set_determinants_psi_det`
       * :c:func:`ezfio_set_determinants_psi_det_qp_edit`
       * :c:func:`ezfio_set_tc_bi_ortho_psi_l_coef_bi_ortho`
       * :c:func:`ezfio_set_tc_bi_ortho_psi_r_coef_bi_ortho`
       * :c:func:`write_int`

 
.. c:function:: single_hmat_bi_ortho:


    File : :file:`h_biortho.irp.f`

    .. code:: fortran

        subroutine single_hmat_bi_ortho(Nint, key_j, key_i, hmono, htwoe)


    
    < key_j | H | key_i > for single excitation
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_coul_e`
       * :c:data:`mo_bi_ortho_one_e`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`hmat_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`decode_exc`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`

 
.. c:function:: single_htilde_mu_mat_bi_ortho:


    File : :file:`slater_tc.irp.f`

    .. code:: fortran

        subroutine single_htilde_mu_mat_bi_ortho(Nint, key_j, key_i, hmono, htwoe, htot)


     <key_j | H_tilde | key_i> for single excitation ONLY FOR ONE- AND TWO-BODY TERMS
    !
    ! WARNING !!
    
     Non hermitian !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_bi_ortho_tc_two_e`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`decode_exc`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`

 
.. c:function:: single_htilde_mu_mat_fock_bi_ortho:


    File : :file:`slater_tc_opt_single.irp.f`

    .. code:: fortran

        subroutine  single_htilde_mu_mat_fock_bi_ortho (Nint, key_j, key_i, hmono, htwoe, hthree, htot)


     <key_j | H_tilde | key_i> for single excitation ONLY FOR ONE- AND TWO-BODY TERMS
    !
    ! WARNING !!
    
     Non hermitian !!

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_opt_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`decode_exc`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`
       * :c:func:`get_single_excitation_from_fock_tc`

 
.. c:function:: single_htilde_mu_mat_fock_bi_ortho_no_3e:


    File : :file:`slater_tc_opt_single.irp.f`

    .. code:: fortran

        subroutine  single_htilde_mu_mat_fock_bi_ortho_no_3e(Nint, key_j, key_i, htot)


     <key_j | H_tilde | key_i> for single excitation ONLY FOR ONE- AND TWO-BODY TERMS
    !
    ! WARNING !!
    
     Non hermitian !!

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_opt_bi_ortho_no_3e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`decode_exc`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`
       * :c:func:`get_single_excitation_from_fock_tc_no_3e`

 
.. c:function:: single_htilde_three_body_ints_bi_ort:


    File : :file:`slater_tc_3e.irp.f`

    .. code:: fortran

        subroutine single_htilde_three_body_ints_bi_ort(Nint, key_j, key_i, hthree)


     <key_j | H_tilde | key_i> for single excitation ONLY FOR THREE-BODY TERMS WITH BI ORTHONORMAL ORBITALS
    !
    ! WARNING !!
    
     Non hermitian !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_exch13_bi_ort`
       * :c:data:`three_e_4_idx_exch23_bi_ort`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`decode_exc`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`

 
.. c:function:: sym_3_e_int_from_6_idx_tensor:


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        double precision function sym_3_e_int_from_6_idx_tensor(n,l,k,m,j,i)


    returns all good combinations of permutations of integrals with the good signs
    
    for a given (k^dagger l^dagger n^dagger m j i)  <nlk|L|mji> when all indices have the same spins

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`three_body_ints_bi_ort`

 
.. c:function:: tc_bi_ortho:


    File : :file:`test_tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine tc_bi_ortho


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`test_h_u0`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

 
.. c:function:: tc_bi_ortho_prop:


    File : :file:`tc_bi_ortho_prop.irp.f`

    .. code:: fortran

        subroutine tc_bi_ortho_prop


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`test`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

 
.. c:function:: tc_natorb_bi_ortho:


    File : :file:`save_tc_bi_ortho_nat.irp.f`

    .. code:: fortran

        subroutine tc_natorb_bi_ortho


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_energy_and_mos`
       * :c:func:`save_tc_natorb`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

 
.. c:function:: tc_print_mulliken_sd:


    File : :file:`spin_mulliken.irp.f`

    .. code:: fortran

        subroutine tc_print_mulliken_sd



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`element_name`
       * :c:data:`l_to_character`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_num`
       * :c:data:`tc_mulliken_spin_densities`
       * :c:data:`tc_spin_gross_orbital_product`
       * :c:data:`tc_spin_population_angular_momentum`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_spin_dens`

 
.. c:function:: tc_som:


    File : :file:`tc_som.irp.f`

    .. code:: fortran

        subroutine tc_som


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
       * :c:data:`read_wf`

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
       * :c:data:`read_wf`

 
.. c:function:: test:


    File : :file:`test_normal_order.irp.f`

    .. code:: fortran

        subroutine test



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`eff_2_e_from_3_e_aa`
       * :c:data:`eff_2_e_from_3_e_bb`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`hf_bitmask`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`tc_bi_ortho_prop`
       * :c:func:`test_normal_order`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`decode_exc`
       * :c:func:`do_single_excitation`
       * :c:func:`get_excitation`
       * :c:func:`get_excitation_degree`
       * :c:func:`htilde_mu_mat_bi_ortho`
       * :c:func:`three_comp_two_e_elem`

 
.. c:function:: test_3e:


    File : :file:`test_tc_fock.irp.f`

    .. code:: fortran

        subroutine test_3e



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`diag_three_elem_hf`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`
       * :c:data:`tc_hf_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`bi_ort_ints`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_aaa_contrib`
       * :c:func:`give_aab_contrib`
       * :c:func:`give_abb_contrib`
       * :c:func:`give_bbb_contrib`
       * :c:func:`htilde_mu_mat_bi_ortho`

 
.. c:function:: test_h_u0:


    File : :file:`test_tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine test_h_u0



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_r_coef_bi_ortho`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`tc_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp`
       * :c:func:`htc_bi_ortho_calc_tdav`
       * :c:func:`htcdag_bi_ortho_calc_tdav`

 
.. c:function:: test_natorb:


    File : :file:`test_natorb.irp.f`

    .. code:: fortran

        subroutine test_natorb


    TODO : Reads psi_det in the EZFIO folder and prints out the left- and right-eigenvectors together with the energy. Saves the left-right wave functions at the end.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`routine`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

 
.. c:function:: test_normal_order:


    File : :file:`test_normal_order.irp.f`

    .. code:: fortran

        subroutine test_normal_order


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`provide_all_three_ints_bi_ortho`
       * :c:func:`test`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

 
.. c:function:: test_slater_tc_opt:


    File : :file:`test_tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine test_slater_tc_opt



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`
       * :c:func:`get_excitation_degree`
       * :c:func:`htilde_mu_mat_bi_ortho`
       * :c:func:`htilde_mu_mat_opt_bi_ortho`

 
.. c:function:: test_spin_dens:


    File : :file:`print_tc_spin_dens.irp.f`

    .. code:: fortran

        subroutine test_spin_dens


    TODO : Reads psi_det in the EZFIO folder and prints out the left- and right-eigenvectors together with the energy. Saves the left-right wave functions at the end.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`tc_print_mulliken_sd`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

 
.. c:function:: test_tc_fock:


    File : :file:`test_tc_fock.irp.f`

    .. code:: fortran

        subroutine test_tc_fock


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`routine_tot`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`
       * :c:data:`read_wf`

 
.. c:function:: three_comp_fock_elem:


    File : :file:`slater_tc_opt_single.irp.f`

    .. code:: fortran

        subroutine three_comp_fock_elem(key_i,h_fock,p_fock,ispin_fock,hthree)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_single_excitation_from_fock_tc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`fock_a_tc_operator`
       * :c:func:`fock_ac_tc_operator`

 
.. c:function:: three_comp_two_e_elem:


    File : :file:`slater_tc_opt_double.irp.f`

    .. code:: fortran

        subroutine three_comp_two_e_elem(key_i,h1,h2,p1,p2,s1,s2,hthree)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`eff_2_e_from_3_e_aa`
       * :c:data:`eff_2_e_from_3_e_ab`
       * :c:data:`eff_2_e_from_3_e_bb`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_5_idx_exch13_bi_ort`
       * :c:data:`three_e_5_idx_exch23_bi_ort`
       * :c:data:`three_e_double_parrallel_spin_prov`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`double_htilde_mu_mat_fock_bi_ortho`
       * :c:func:`test`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: three_e_diag_parrallel_spin:


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        double precision function three_e_diag_parrallel_spin(m,j,i)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`three_e_3_idx_cycle_1_bi_ort`
       * :c:data:`three_e_3_idx_cycle_2_bi_ort`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`
       * :c:data:`three_e_3_idx_exch13_bi_ort`
       * :c:data:`three_e_3_idx_exch23_bi_ort`

 
.. c:function:: three_e_double_parrallel_spin:


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        double precision function three_e_double_parrallel_spin(m,l,j,k,i)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`three_e_5_idx_cycle_1_bi_ort`
       * :c:data:`three_e_5_idx_cycle_2_bi_ort`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_5_idx_exch13_bi_ort`
       * :c:data:`three_e_5_idx_exch23_bi_ort`

 
.. c:function:: three_e_single_parrallel_spin:


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        double precision function three_e_single_parrallel_spin(m,j,k,i)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`three_e_4_idx_cycle_1_bi_ort`
       * :c:data:`three_e_4_idx_cycle_2_bi_ort`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_exch12_bi_ort`
       * :c:data:`three_e_4_idx_exch13_bi_ort`
       * :c:data:`three_e_4_idx_exch23_bi_ort`

 
.. c:function:: timing_diag:


    File : :file:`test_tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine timing_diag



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho`
       * :c:func:`htilde_mu_mat_opt_bi_ortho`
       * :c:func:`wall_time`

 
.. c:function:: timing_double:


    File : :file:`test_tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine timing_double



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_degree`
       * :c:func:`htilde_mu_mat_bi_ortho`
       * :c:func:`htilde_mu_mat_opt_bi_ortho`
       * :c:func:`wall_time`

 
.. c:function:: timing_single:


    File : :file:`test_tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine timing_single



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_degree`
       * :c:func:`htilde_mu_mat_bi_ortho`
       * :c:func:`htilde_mu_mat_opt_bi_ortho`
       * :c:func:`wall_time`

 
.. c:function:: timing_tot:


    File : :file:`test_tc_bi_ortho.irp.f`

    .. code:: fortran

        subroutine timing_tot



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_bi_ortho`
       * :c:func:`htilde_mu_mat_opt_bi_ortho`
       * :c:func:`i_h_j`
       * :c:func:`wall_time`

 
.. c:function:: u_0_h_tc_u_0:


    File : :file:`u0_h_u0.irp.f`

    .. code:: fortran

        subroutine u_0_H_tc_u_0(e_0,u_0,n,keys_tmp,Nint,N_st,sze, do_right)


    Computes $E_0 = \frac{\langle u_0 | H_TC | u_0 \rangle}{\langle u_0 | u_0 \rangle}$
    
    n : number of determinants
    
    if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp`

 
.. c:function:: write_l_r_wf:


    File : :file:`print_tc_wf.irp.f`

    .. code:: fortran

        subroutine write_l_r_wf



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`n_det`
       * :c:data:`psi_det_sorted_tc_left`
       * :c:data:`psi_det_sorted_tc_right`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_tc_bi_ortho`

