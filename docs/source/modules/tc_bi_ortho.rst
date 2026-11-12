.. _module_tc_bi_ortho: 
 
.. program:: tc_bi_ortho 
 
.. default-role:: option 
 
tc_bi_ortho
===========

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
EZFIO parameters 
---------------- 
 
.. option:: psi_l_coef_bi_ortho
 
    Coefficients for the left wave function
 
 
.. option:: psi_r_coef_bi_ortho
 
    Coefficients for the right wave function
 
 
.. option:: tc_gs_energy
 
    TC GS Energy
 
 
.. option:: tc_gs_var
 
    TC GS VAR
 
 
 
Programs 
-------- 
 
 * :ref:`diagonalize_tc_h` 
 
Providers 
--------- 
 
.. c:var:: double coef_pt1_bi_ortho


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: coef_pt1_bi_ortho	(N_det)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_00`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_corr_bi_orth


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	
        double precision	:: e_corr_bi_orth_proj_abs	
        double precision	:: e_corr_single_bi_orth_abs	
        double precision	:: e_corr_double_bi_orth_abs	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_corr_bi_orth_proj


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	
        double precision	:: e_corr_bi_orth_proj_abs	
        double precision	:: e_corr_single_bi_orth_abs	
        double precision	:: e_corr_double_bi_orth_abs	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_corr_bi_orth_proj_abs


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	
        double precision	:: e_corr_bi_orth_proj_abs	
        double precision	:: e_corr_single_bi_orth_abs	
        double precision	:: e_corr_double_bi_orth_abs	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_corr_double_bi_orth


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	
        double precision	:: e_corr_bi_orth_proj_abs	
        double precision	:: e_corr_single_bi_orth_abs	
        double precision	:: e_corr_double_bi_orth_abs	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_corr_double_bi_orth_abs


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	
        double precision	:: e_corr_bi_orth_proj_abs	
        double precision	:: e_corr_single_bi_orth_abs	
        double precision	:: e_corr_double_bi_orth_abs	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_corr_single_bi_orth


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	
        double precision	:: e_corr_bi_orth_proj_abs	
        double precision	:: e_corr_single_bi_orth_abs	
        double precision	:: e_corr_double_bi_orth_abs	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_corr_single_bi_orth_abs


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_corr_bi_orth	
        double precision	:: e_corr_bi_orth_proj	
        double precision	:: e_corr_single_bi_orth	
        double precision	:: e_corr_double_bi_orth	
        double precision	:: e_corr_bi_orth_proj_abs	
        double precision	:: e_corr_single_bi_orth_abs	
        double precision	:: e_corr_double_bi_orth_abs	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_bi_orth_00`
       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_pt2_tc_bi_orth


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_pt2_tc_bi_orth	
        double precision	:: e_pt2_tc_bi_orth_single	
        double precision	:: e_pt2_tc_bi_orth_double	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_00`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_pt2_tc_bi_orth_double


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_pt2_tc_bi_orth	
        double precision	:: e_pt2_tc_bi_orth_single	
        double precision	:: e_pt2_tc_bi_orth_double	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_00`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_pt2_tc_bi_orth_single


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_pt2_tc_bi_orth	
        double precision	:: e_pt2_tc_bi_orth_single	
        double precision	:: e_pt2_tc_bi_orth_double	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_tilde_00`
       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_tc_left_right


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_tc_left_right	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_det`
       * :c:data:`pure_three_body_h_tc`


 
.. c:var:: double e_tilde_00


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_tilde_00	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_int`
       * :c:data:`pure_three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_pt1_bi_ortho`
       * :c:data:`e_pt2_tc_bi_orth`

 
.. c:var:: double e_tilde_bi_orth_00


    File : :file:`e_corr_bi_ortho.irp.f`

    .. code:: fortran

        double precision	:: e_tilde_bi_orth_00	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`
       * :c:data:`pure_three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`

 
.. c:var:: double eigval_left_tc_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: s2_eigvec_tc_bi_orth	(N_states)
        double precision, allocatable	:: norm_ground_left_right_bi_orth	(N_states)


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`expected_s2`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`pure_three_body_h_tc`
       * :c:data:`qp_max_mem`
       * :c:data:`s2_eig`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: double eigval_right_tc_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: s2_eigvec_tc_bi_orth	(N_states)
        double precision, allocatable	:: norm_ground_left_right_bi_orth	(N_states)


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`expected_s2`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`pure_three_body_h_tc`
       * :c:data:`qp_max_mem`
       * :c:data:`s2_eig`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: double fock_diag_sorted_l_natorb


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

 
.. c:var:: double fock_diag_sorted_r_natorb


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

 
.. c:var:: double fock_diag_sorted_v_natorb


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

 
.. c:var:: int index_hf_psi_det


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

 
.. c:var:: double leigvec_tc_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: s2_eigvec_tc_bi_orth	(N_states)
        double precision, allocatable	:: norm_ground_left_right_bi_orth	(N_states)


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`expected_s2`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`pure_three_body_h_tc`
       * :c:data:`qp_max_mem`
       * :c:data:`s2_eig`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: double mo_r_coef_normalized


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

 
.. c:var:: double natorb_tc_eigval


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

       * :c:data:`ao_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`im_thresh_tc`
       * :c:data:`mo_num`
       * :c:data:`tc_transition_matrix_ao`
       * :c:data:`tc_transition_matrix_mo`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_diag_sorted_r_natorb`

 
.. c:var:: double natorb_tc_leigvec_ao


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


 
.. c:var:: double natorb_tc_leigvec_mo


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

       * :c:data:`ao_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`im_thresh_tc`
       * :c:data:`mo_num`
       * :c:data:`tc_transition_matrix_ao`
       * :c:data:`tc_transition_matrix_mo`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_diag_sorted_r_natorb`

 
.. c:var:: double natorb_tc_reigvec_ao


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


 
.. c:var:: double natorb_tc_reigvec_mo


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

       * :c:data:`ao_num`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`im_thresh_tc`
       * :c:data:`mo_num`
       * :c:data:`tc_transition_matrix_ao`
       * :c:data:`tc_transition_matrix_mo`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_diag_sorted_r_natorb`

 
.. c:var:: double norm_ground_left_right_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: s2_eigvec_tc_bi_orth	(N_states)
        double precision, allocatable	:: norm_ground_left_right_bi_orth	(N_states)


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`expected_s2`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`pure_three_body_h_tc`
       * :c:data:`qp_max_mem`
       * :c:data:`s2_eig`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: double overlap_natorb_tc_eigvec_ao


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


 
.. c:var:: double psi_average_norm_contrib_sorted_tc


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


 
.. c:var:: double psi_average_norm_contrib_tc


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

 
.. c:var:: double psi_bitcleft_bilinear_matrix_values


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


 
.. c:var:: double psi_coef_sorted_tc


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


 
.. c:var:: double psi_coef_sorted_tc_bit


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


 
.. c:var:: int psi_det_sorted_tc


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


 
.. c:var:: int psi_det_sorted_tc_bit


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


 
.. c:var:: int psi_det_sorted_tc_left


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


 
.. c:var:: int psi_det_sorted_tc_order


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


 
.. c:var:: int psi_det_sorted_tc_right


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


 
.. c:var:: double psi_l_coef_bi_ortho


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_l_coef_bi_ortho	(N_det,N_states)


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
       * :c:data:`tc_two_rdm_chemist`

 
.. c:var:: double psi_l_coef_sorted_bi_ortho


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


 
.. c:var:: double psi_l_coef_sorted_bi_ortho_left


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


 
.. c:var:: double psi_r_coef_bi_ortho


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_r_coef_bi_ortho	(N_det,N_states)


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
       * :c:data:`read_wf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`psi_average_norm_contrib_tc`
       * :c:data:`psi_det_sorted_tc`
       * :c:data:`psi_det_sorted_tc_left`
       * :c:data:`psi_det_sorted_tc_right`
       * :c:data:`tc_transition_matrix_mo_beta`
       * :c:data:`tc_two_rdm_chemist`

 
.. c:var:: double psi_r_coef_sorted_bi_ortho


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


 
.. c:var:: double psi_r_coef_sorted_bi_ortho_right


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


 
.. c:var:: double reigvec_tc_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: s2_eigvec_tc_bi_orth	(N_states)
        double precision, allocatable	:: norm_ground_left_right_bi_orth	(N_states)


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`expected_s2`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`pure_three_body_h_tc`
       * :c:data:`qp_max_mem`
       * :c:data:`s2_eig`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: double s2_eigvec_tc_bi_orth


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigval_right_tc_bi_orth	(N_states)
        double precision, allocatable	:: eigval_left_tc_bi_orth	(N_states)
        double precision, allocatable	:: reigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: leigvec_tc_bi_orth	(N_det,N_states)
        double precision, allocatable	:: s2_eigvec_tc_bi_orth	(N_states)
        double precision, allocatable	:: norm_ground_left_right_bi_orth	(N_states)


    eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`expected_s2`
       * :c:data:`htilde_matrix_elmt_bi_ortho`
       * :c:data:`index_hf_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`
       * :c:data:`pure_three_body_h_tc`
       * :c:data:`qp_max_mem`
       * :c:data:`s2_eig`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`state_following_tc`
       * :c:data:`threshold_nonsym_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_corr_bi_orth`
       * :c:data:`e_tc_left_right`

 
.. c:var:: double tc_bi_ortho_dipole


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


 
.. c:var:: double tc_mulliken_spin_densities


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


 
.. c:var:: double tc_spin_dens_right_only


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

 
.. c:var:: double tc_spin_gross_orbital_product


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

 
.. c:var:: double tc_spin_population


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

 
.. c:var:: double tc_spin_population_angular_momentum


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


 
.. c:var:: double tc_spin_population_angular_momentum_per_atom


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


 
.. c:var:: double tc_spin_transition_matrix_ao


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

 
.. c:var:: double tc_spin_transition_matrix_mo


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

 
.. c:var:: double tc_transition_matrix_ao


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`natorb_tc_reigvec_mo`

 
.. c:var:: double tc_transition_matrix_mo


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

 
.. c:var:: double tc_transition_matrix_mo_alpha


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

 
.. c:var:: double tc_transition_matrix_mo_beta


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

 
.. c:var:: double tc_two_rdm


    File : :file:`two_rdm_naive.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_two_rdm	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: tc_two_rdm_s1s2	(mo_num,mo_num,mo_num,mo_num,2,2)


    tc_two_rdm(p,q,s,r) = <Phi| a^dagger_p a^dagger_q q_r a_s |Phi> = PHYSICIST NOTATION

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`tc_two_rdm_chemist`


 
.. c:var:: double tc_two_rdm_chemist


    File : :file:`two_rdm_naive.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_two_rdm_chemist	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: tc_two_rdm_chemist_s1s2	(mo_num,mo_num,mo_num,mo_num,2,2)


    tc_two_rdm_chemist(p,s,q,r) = <Phi| a^dagger_p a^dagger_q q_r a_s |Phi> = CHEMIST NOTATION

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
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_two_rdm`

 
.. c:var:: double tc_two_rdm_chemist_s1s2


    File : :file:`two_rdm_naive.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_two_rdm_chemist	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: tc_two_rdm_chemist_s1s2	(mo_num,mo_num,mo_num,mo_num,2,2)


    tc_two_rdm_chemist(p,s,q,r) = <Phi| a^dagger_p a^dagger_q q_r a_s |Phi> = CHEMIST NOTATION

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
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_two_rdm`

 
.. c:var:: double tc_two_rdm_s1s2


    File : :file:`two_rdm_naive.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_two_rdm	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: tc_two_rdm_s1s2	(mo_num,mo_num,mo_num,mo_num,2,2)


    tc_two_rdm(p,q,s,r) = <Phi| a^dagger_p a^dagger_q q_r a_s |Phi> = PHYSICIST NOTATION

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`tc_two_rdm_chemist`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void bi_normalize(void)


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        subroutine bi_normalize(u_l, u_r, n, ld, nstates)


    !!! Normalization of the scalar product of the left/right eigenvectors

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`

 
.. c:function:: void davidson_hs2_nonsym_b1space(void)


    File : :file:`dav_h_tc_s2.irp.f`

    .. code:: fortran

        subroutine davidson_hs2_nonsym_b1space(u_in, H_jj, s2_out,energies, sze, N_st, N_st_diag_in, n_it_max_dav, converged, hcalc)


    Generic modified-Davidson diagonalization
    
    H_jj : specific diagonal H matrix elements to diagonalize de Davidson
    
    u_in : guess coefficients on the various states. Overwritten on exit by right eigenvectors
    
    sze : Number of determinants
    
    N_st : Number of eigenstates
    
    N_st_diag_in : Number of states in which H is diagonalized. Assumed > N_st
    
    Initial guess vectors are not necessarily orthonormal
    
    hcalc subroutine to compute W = H U (see routine hcalc_template for template of input/output)
    
    !!! WARNING !!! IT SEEMS THAT IF THE NUMBER OF MACRO ITERATIONS EXCEEDS n_it_max_dav,
    
                    THE RECONTRACTION IS WRONG. YOU SHOULD CONSIDER CALLING MULTIPLE TIME THE ROUTINE
    
                    SEE FOR INSTANCE IN tc_bi_ortho/tc_h_eigvectors.irp.f

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`qp_max_mem`
       * :c:data:`threshold_nonsym_davidson`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`dgemm`
       * :c:func:`diag_nonsym_right`
       * :c:func:`hcalc`
       * :c:func:`normalize`
       * :c:func:`ortho_qr`
       * :c:func:`random_number`
       * :c:func:`resident_memory`
       * :c:func:`write_double`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:function:: void diagonalize_ci_tc(void)


    File : :file:`tc_h_eigvectors.irp.f`

    .. code:: fortran

        subroutine diagonalize_CI_tc()


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

 
.. c:function:: void direct_sym_3_e_int(void)


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        double precision function direct_sym_3_e_int(n,l,k,m,j,i)


    returns all good combinations of permutations of integrals with the good signs
    
    for a given (k^dagger l^dagger n^dagger m j i)  <nlk|L|mji> when all indices have the same spins

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`

 
.. c:function:: void give_all_perm_for_three_e(void)


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        subroutine give_all_perm_for_three_e(n,l,k,m,j,i,idx_list,phase)


    returns all the list of permutting indices for the antimmetrization of
    
    (k^dagger l^dagger n^dagger m j i)  <nlk|L|mji>   when all indices have the same spins
    
    idx_list(:,i) == list of the 6 indices corresponding the permutation "i"
    
    phase(i) == phase of the permutation "i"
    
    there are in total 6 permutations with different indices

 
.. c:function:: void h_tc_dagger_u_0_opt(void)


    File : :file:`h_tc_u0.irp.f`

    .. code:: fortran

        subroutine H_tc_dagger_u_0_opt(v_0,u_0,N_st,sze)


    Computes $v_0 = H | u_0\rangle$.
    
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep split the work between the OpenMP threads.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp`

 
.. c:function:: void h_tc_u_0_nstates_openmp(void)


    File : :file:`h_tc_u0.irp.f`

    .. code:: fortran

        subroutine H_tc_u_0_nstates_openmp(v_0,u_0,N_st,sze, do_right)


    Computes $v_0 = H | u_0\rangle$.
    
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep split the work between the OpenMP threads.
    
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
       * :c:func:`u_0_h_tc_u_0`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dset_order`
       * :c:func:`dtranspose`
       * :c:func:`h_tc_u_0_nstates_openmp_work`

 
.. c:function:: void h_tc_u_0_nstates_openmp_work(void)


    File : :file:`h_tc_u0.irp.f`

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

 
.. c:function:: void h_tc_u_0_nstates_openmp_work_1(void)


    File : :file:`h_tc_u0.irp.f_template_765`

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

 
.. c:function:: void h_tc_u_0_nstates_openmp_work_2(void)


    File : :file:`h_tc_u0.irp.f_template_765`

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

 
.. c:function:: void h_tc_u_0_nstates_openmp_work_3(void)


    File : :file:`h_tc_u0.irp.f_template_765`

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

 
.. c:function:: void h_tc_u_0_nstates_openmp_work_4(void)


    File : :file:`h_tc_u0.irp.f_template_765`

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

 
.. c:function:: void h_tc_u_0_nstates_openmp_work_n_int(void)


    File : :file:`h_tc_u0.irp.f_template_765`

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

 
.. c:function:: void h_tc_u_0_opt(void)


    File : :file:`h_tc_u0.irp.f`

    .. code:: fortran

        subroutine H_tc_u_0_opt(v_0,u_0,N_st,sze)


    Computes $v_0 = H | u_0\rangle$.
    
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep split the work between the OpenMP threads.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp`

 
.. c:function:: void routine_diag(void)


    File : :file:`diagonalize_tc_h.irp.f`

    .. code:: fortran

        subroutine routine_diag()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`eigval_right_tc_bi_orth`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`nol_0e`
       * :c:data:`nol_1e`
       * :c:data:`nol_2e`
       * :c:data:`nol_standard`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diagonalize_tc_h`

 
.. c:function:: void routine_save_left_right_sorted_bi_ortho(void)


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        subroutine routine_save_left_right_sorted_bi_ortho()



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

 
.. c:function:: void routine_save_right_bi_ortho(void)


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        subroutine routine_save_right_bi_ortho()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_det`
       * :c:data:`psi_r_coef_bi_ortho`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`save_tc_bi_ortho_wavefunction`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_wavefunction_general_unormalized`

 
.. c:function:: void routine_save_right_sorted_bi_ortho(void)


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        subroutine routine_save_right_sorted_bi_ortho()



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

 
.. c:function:: void save_tc_bi_ortho_wavefunction(void)


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        subroutine save_tc_bi_ortho_wavefunction()



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
       * :c:func:`routine_save_right_sorted_bi_ortho`
       * :c:func:`save_tc_wavefunction_general`

 
.. c:function:: void save_tc_wavefunction_general(void)


    File : :file:`psi_r_l_prov.irp.f`

    .. code:: fortran

        subroutine save_tc_wavefunction_general(ndet, nstates, psidet, sze, dim_psicoef, psilcoef, psircoef)


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

 
.. c:function:: void sym_3_e_int_from_6_idx_tensor(void)


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        double precision function sym_3_e_int_from_6_idx_tensor(n, l, k, m, j, i)


    returns all good combinations of permutations of integrals with the good signs
    
    for a given (k^dagger l^dagger n^dagger m j i)  <nlk|L|mji> when all indices have the same spins

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`three_body_ints_bi_ort`

 
.. c:function:: void tc_print_mulliken_sd(void)


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

 
.. c:function:: void three_e_diag_parrallel_spin(void)


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        double precision function three_e_diag_parrallel_spin(m, j, i)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`three_e_3_idx_cycle_1_bi_ort`
       * :c:data:`three_e_3_idx_cycle_2_bi_ort`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`
       * :c:data:`three_e_3_idx_exch13_bi_ort`
       * :c:data:`three_e_3_idx_exch23_bi_ort`

 
.. c:function:: void three_e_double_parrallel_spin(void)


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        double precision function three_e_double_parrallel_spin(m,l,j,k,i)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`three_e_5_idx_direct_bi_ort`

 
.. c:function:: void three_e_single_parrallel_spin(void)


    File : :file:`symmetrized_3_e_int.irp.f`

    .. code:: fortran

        double precision function three_e_single_parrallel_spin(m,j,k,i)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`three_e_4_idx_direct_bi_ort`

 
.. c:function:: void u_0_h_tc_u_0(void)


    File : :file:`h_tc_u0.irp.f`

    .. code:: fortran

        subroutine u_0_H_tc_u_0(e_0,u_0,n,keys_tmp,Nint,N_st,sze, do_right)


    Computes $E_0 = \frac{\langle u_0 | H_TC | u_0 \rangle}{\langle u_0 | u_0 \rangle}$
    
    n : number of determinants
    
    if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_tc_u_0_nstates_openmp`

 
.. c:function:: void update_tc_rdm(void)


    File : :file:`two_rdm_naive.irp.f`

    .. code:: fortran

        subroutine update_tc_rdm(h1,p1,h2,p2,s1,s2,array,sze,contrib)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_two_rdm_chemist`

 
.. c:function:: void write_tc_energy(void)


    File : :file:`tc_utils.irp.f`

    .. code:: fortran

        subroutine write_tc_energy()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det`
       * :c:data:`psi_l_coef_bi_ortho`
       * :c:data:`psi_r_coef_bi_ortho`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_tc_bi_ortho_tc_gs_energy`
       * :c:func:`htilde_mu_mat_opt_bi_ortho`

 
.. c:function:: void write_tc_gs_var_hf(void)


    File : :file:`tc_utils.irp.f`

    .. code:: fortran

        subroutine write_tc_gs_var_HF()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_tc_bi_ortho_tc_gs_var`
       * :c:func:`htilde_mu_mat_opt_bi_ortho`

 
.. c:function:: void write_tc_var(void)


    File : :file:`tc_utils.irp.f`

    .. code:: fortran

        subroutine write_tc_var()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`htilde_mu_mat_opt_bi_ortho`

