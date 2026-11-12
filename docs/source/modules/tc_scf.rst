.. _module_tc_scf: 
 
.. program:: tc_scf 
 
.. default-role:: option 
 
tc_scf
======

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
EZFIO parameters 
---------------- 
 
.. option:: tcscf_energy
 
    TC-SCF ENERGY
 
 
.. option:: converged_tcscf
 
    If |true|, tc-scf has converged
 
    Default: False
 
.. option:: max_dim_diis_tcscf
 
    Maximum size of the DIIS extrapolation procedure
 
    Default: 15
 
.. option:: level_shift_tcscf
 
    Energy shift on the virtual MOs to improve TCSCF convergence
 
    Default: 0.
 
.. option:: thresh_tcscf
 
    Threshold on the convergence of the Hartree Fock energy.
 
    Default: 1.e-8
 
.. option:: n_it_tcscf_max
 
    Maximum number of SCF iterations
 
    Default: 50
 
.. option:: tc_Brillouin_Right
 
    If |true|, impose only right-Brillouin condition
 
    Default: False
 
 
Programs 
-------- 
 
 * :ref:`minimize_tc_angles` 
 * :ref:`molden_lr_mos` 
 * :ref:`print_overlap_bi_ort` 
 * :ref:`tc_scf` 
 * :ref:`write_ao_2e_tc_integ` 
 
Providers 
--------- 
 
.. c:var:: double ao_two_e_core_integral


    File : :file:`fock_tc_core.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_core_integral	(ao_num,ao_num)


    Fock matrices in AO basis set for only the core orbitals
    
    WARNING : use the TCSCF_bi_ort_core_dm_ao matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`tcscf_bi_ort_core_dm_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_core_eri`
       * :c:data:`tc_hf_energy`

 
.. c:var:: double ao_two_e_core_integral_chol


    File : :file:`fock_tc_core.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_core_integral_chol	(ao_num,ao_num)


    Fock matrices in AO basis set for only the core orbitals
    
    WARNING : use the TCSCF_bi_ort_core_dm_ao matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`cholesky_ao_num`
       * :c:data:`qp_max_mem`
       * :c:data:`tcscf_bi_ort_core_dm_ao`


 
.. c:var:: double ao_two_e_valence_integral_alpha


    File : :file:`fock_tc_core.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_valence_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_valence_integral_beta	(ao_num,ao_num)


    Fock matrices in AO basis set for only the core orbitals
    
    WARNING : use the TCSCF_bi_ort_core_dm_ao matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`tcscf_bi_ort_valence_dm_ao_alpha`
       * :c:data:`tcscf_bi_ort_valence_dm_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_alpha_tc_mo_valence_eri`
       * :c:data:`fock_matrix_beta_tc_mo_valence_eri`
       * :c:data:`tc_hf_energy`

 
.. c:var:: double ao_two_e_valence_integral_alpha_chol


    File : :file:`fock_tc_core.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_valence_integral_alpha_chol	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_valence_integral_beta_chol	(ao_num,ao_num)


    Alpha and Beta Fock matrices in AO basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`cholesky_ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`qp_max_mem`
       * :c:data:`tcscf_bi_ort_valence_dm_ao_alpha`
       * :c:data:`tcscf_bi_ort_valence_dm_ao_beta`
       * :c:data:`tcscf_bi_ort_valence_dm_ao_full`


 
.. c:var:: double ao_two_e_valence_integral_beta


    File : :file:`fock_tc_core.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_valence_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_valence_integral_beta	(ao_num,ao_num)


    Fock matrices in AO basis set for only the core orbitals
    
    WARNING : use the TCSCF_bi_ort_core_dm_ao matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`tcscf_bi_ort_valence_dm_ao_alpha`
       * :c:data:`tcscf_bi_ort_valence_dm_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_alpha_tc_mo_valence_eri`
       * :c:data:`fock_matrix_beta_tc_mo_valence_eri`
       * :c:data:`tc_hf_energy`

 
.. c:var:: double ao_two_e_valence_integral_beta_chol


    File : :file:`fock_tc_core.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_valence_integral_alpha_chol	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_valence_integral_beta_chol	(ao_num,ao_num)


    Alpha and Beta Fock matrices in AO basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`cholesky_ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`qp_max_mem`
       * :c:data:`tcscf_bi_ort_valence_dm_ao_alpha`
       * :c:data:`tcscf_bi_ort_valence_dm_ao_beta`
       * :c:data:`tcscf_bi_ort_valence_dm_ao_full`


 
.. c:var:: double eigval_fock_tc_mo


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
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_num`
       * :c:data:`overlap_mo_r`
       * :c:data:`thr_degen_tc`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_tc_reigvec_ao`

 
.. c:var:: double fock_3e_mo_a


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3e_mo_a	(mo_num,mo_num)


    
    Fock matrix alpha from three-electron terms
    
    WARNING :: non hermitian if bi-ortho MOS used
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_3e_mo_a_os`
       * :c:data:`fock_3e_mo_cs`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_alpha`

 
.. c:var:: double fock_3e_mo_a_os


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3e_mo_a_os	(mo_num,mo_num)
        double precision, allocatable	:: fock_3e_mo_b_os	(mo_num,mo_num)


    
    Open Shell part of the Fock matrix from three-electron terms
    
    WARNING :: non hermitian if bi-ortho MOS used
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_tc_op`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3e_mo_a`
       * :c:data:`fock_3e_mo_b`

 
.. c:var:: double fock_3e_mo_b


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3e_mo_b	(mo_num,mo_num)


    
    Fock matrix beta from three-electron terms
    
    WARNING :: non hermitian if bi-ortho MOS used
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_3e_mo_a_os`
       * :c:data:`fock_3e_mo_a_os`
       * :c:data:`fock_3e_mo_cs`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_beta`

 
.. c:var:: double fock_3e_mo_b_os


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3e_mo_a_os	(mo_num,mo_num)
        double precision, allocatable	:: fock_3e_mo_b_os	(mo_num,mo_num)


    
    Open Shell part of the Fock matrix from three-electron terms
    
    WARNING :: non hermitian if bi-ortho MOS used
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_tc_op`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3e_mo_a`
       * :c:data:`fock_3e_mo_b`

 
.. c:var:: double fock_3e_mo_cs


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_3e_mo_cs	(mo_num,mo_num)


    
    Closed Shell part of the Fock matrix from three-electron terms
    
    WARNING: if core_tc_op then only the valence electrons are considered
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_tc_op`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`tc_save_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3e_mo_a`
       * :c:data:`fock_3e_mo_b`

 
.. c:var:: double fock_a_tmp1_bi_ortho


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
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_a_tot_3e_bi_orth`

 
.. c:var:: double fock_a_tmp2_bi_ortho


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

 
.. c:var:: double fock_a_tot_3e_bi_orth


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_a_tot_3e_bi_orth	(mo_num,mo_num)


    
    Alpha part of the Fock matrix from three-electron terms
    
    WARNING :: non hermitian if bi-ortho MOS used
    
    This calculation becomes the dominant part one the integrals are provided
    

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

 
.. c:var:: double fock_b_tmp1_bi_ortho


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
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_b_tot_3e_bi_orth`

 
.. c:var:: double fock_b_tmp2_bi_ortho


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

 
.. c:var:: double fock_b_tot_3e_bi_orth


    File : :file:`fock_three_bi_ortho.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_b_tot_3e_bi_orth	(mo_num,mo_num)


    
    Beta part of the Fock matrix from three-electron terms
    
    WARNING :: non hermitian if bi-ortho MOS used
    
    This calculation becomes the dominant part one the integrals are provided
    

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

 
.. c:var:: double fock_cs_3e_bi_orth


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
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`

 
.. c:var:: double fock_matrix_alpha_tc_mo_valence_eri


    File : :file:`fock_tc_core.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_alpha_tc_mo_valence_eri	(mo_num,mo_num)


    Two-electron part with BARE COULOMB interaction of the Fock matrix
    
    with contributions coming only from the VALENCE ELECTRONS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_valence_integral_alpha`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_eri_mo_valence`

 
.. c:var:: double fock_matrix_beta_tc_mo_valence_eri


    File : :file:`fock_tc_core.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_beta_tc_mo_valence_eri	(mo_num,mo_num)


    Two-electron part with BARE COULOMB interaction of the Fock matrix
    
    with contributions coming only from the VALENCE ELECTRONS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_valence_integral_alpha`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_eri_mo_valence`

 
.. c:var:: double fock_matrix_tc_ao_alpha


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_ao_alpha	(ao_num,ao_num)


    Total alpha TC Fock matrix : h_c + Two-e^TC terms on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`two_e_tc_integral_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_alpha`

 
.. c:var:: double fock_matrix_tc_ao_beta


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_ao_beta	(ao_num,ao_num)


    Total beta TC Fock matrix : h_c + Two-e^TC terms on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`two_e_tc_integral_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_beta`

 
.. c:var:: double fock_matrix_tc_ao_tot


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

 
.. c:var:: double fock_matrix_tc_diag_mo_tot


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

       * :c:data:`all_shells_closed`
       * :c:data:`core_tc_op`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`fock_matrix_param`
       * :c:data:`fock_matrix_tc_eri_mo_valence`
       * :c:data:`fock_matrix_tc_mo_alpha`
       * :c:data:`fock_matrix_tc_mo_beta`
       * :c:data:`fock_matrix_tc_mo_core_eri`
       * :c:data:`frozen_orb_scf`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`no_oa_or_av_opt`
       * :c:data:`tc_brillouin_right`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`fock_tc_reigvec_mo`
       * :c:data:`grad_non_hermit_left`

 
.. c:var:: double fock_matrix_tc_eri_mo_valence


    File : :file:`fock_tc_core.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_eri_mo_valence	(mo_num,mo_num)


    Fock matrix containing only the 1/r12 interaction and only the VALENCE density matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_shells_closed`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_alpha_tc_mo_valence_eri`
       * :c:data:`fock_matrix_beta_tc_mo_valence_eri`
       * :c:data:`fock_matrix_param`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`

 
.. c:var:: double fock_matrix_tc_mo_alpha


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_mo_alpha	(mo_num,mo_num)


    Total alpha TC Fock matrix : h_c + Two-e^TC terms on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_3e_mo_a`
       * :c:data:`fock_matrix_tc_ao_alpha`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`

 
.. c:var:: double fock_matrix_tc_mo_beta


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_mo_beta	(mo_num,mo_num)


    Total beta TC Fock matrix : h_c + Two-e^TC terms on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_3e_mo_b`
       * :c:data:`fock_matrix_tc_ao_beta`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`

 
.. c:var:: double fock_matrix_tc_mo_core_eri


    File : :file:`fock_tc_core.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_tc_mo_core_eri	(mo_num,mo_num)


    Two-electron part with BARE COULOMB interaction of the Fock matrix
    
    with contributions coming only from the CORE ELECTRONS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_core_integral`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_tot`

 
.. c:var:: double fock_matrix_tc_mo_tot


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

       * :c:data:`all_shells_closed`
       * :c:data:`core_tc_op`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`fock_matrix_param`
       * :c:data:`fock_matrix_tc_eri_mo_valence`
       * :c:data:`fock_matrix_tc_mo_alpha`
       * :c:data:`fock_matrix_tc_mo_beta`
       * :c:data:`fock_matrix_tc_mo_core_eri`
       * :c:data:`frozen_orb_scf`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`no_oa_or_av_opt`
       * :c:data:`tc_brillouin_right`
       * :c:data:`three_body_h_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`fock_tc_reigvec_mo`
       * :c:data:`grad_non_hermit_left`

 
.. c:var:: double fock_tc_leigvec_ao


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


 
.. c:var:: double fock_tc_leigvec_mo


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
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_num`
       * :c:data:`overlap_mo_r`
       * :c:data:`thr_degen_tc`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_tc_reigvec_ao`

 
.. c:var:: double fock_tc_reigvec_ao


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


 
.. c:var:: double fock_tc_reigvec_mo


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
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_num`
       * :c:data:`overlap_mo_r`
       * :c:data:`thr_degen_tc`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_tc_reigvec_ao`

 
.. c:var:: double fqs_sqf_ao


    File : :file:`diis_tcscf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fqs_sqf_ao	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`q_matrix`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fqs_sqf_mo`

 
.. c:var:: double fqs_sqf_mo


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


 
.. c:var:: double grad_non_hermit


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


 
.. c:var:: double grad_non_hermit_left


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


 
.. c:var:: double grad_non_hermit_right


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


 
.. c:var:: double overlap_fock_tc_eigvec_ao


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


 
.. c:var:: double overlap_fock_tc_eigvec_mo


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
       * :c:data:`level_shift_tcscf`
       * :c:data:`mo_num`
       * :c:data:`overlap_mo_r`
       * :c:data:`thr_degen_tc`
       * :c:data:`thresh_biorthog_diag`
       * :c:data:`thresh_biorthog_nondiag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_tc_reigvec_ao`

 
.. c:var:: double q_alpha


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

 
.. c:var:: double q_beta


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

 
.. c:var:: double q_matrix


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

 
.. c:var:: double tc_hf_energy


    File : :file:`tc_scf_energy.irp.f`

    .. code:: fortran

        double precision	:: tc_hf_energy	
        double precision	:: tc_hf_one_e_energy	
        double precision	:: tc_hf_two_e_energy	
        double precision	:: tc_hf_three_e_energy	


    TC Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`ao_two_e_core_integral`
       * :c:data:`ao_two_e_valence_integral_alpha`
       * :c:data:`core_tc_op`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`nol_0e`
       * :c:data:`nol_standard`
       * :c:data:`nuclear_repulsion`
       * :c:data:`tcscf_bi_ort_core_dm_ao`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`three_body_h_tc`
       * :c:data:`two_e_tc_integral_alpha`
       * :c:data:`two_e_tc_integral_alpha`


 
.. c:var:: double tc_hf_one_e_energy


    File : :file:`tc_scf_energy.irp.f`

    .. code:: fortran

        double precision	:: tc_hf_energy	
        double precision	:: tc_hf_one_e_energy	
        double precision	:: tc_hf_two_e_energy	
        double precision	:: tc_hf_three_e_energy	


    TC Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`ao_two_e_core_integral`
       * :c:data:`ao_two_e_valence_integral_alpha`
       * :c:data:`core_tc_op`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`nol_0e`
       * :c:data:`nol_standard`
       * :c:data:`nuclear_repulsion`
       * :c:data:`tcscf_bi_ort_core_dm_ao`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`three_body_h_tc`
       * :c:data:`two_e_tc_integral_alpha`
       * :c:data:`two_e_tc_integral_alpha`


 
.. c:var:: double tc_hf_three_e_energy


    File : :file:`tc_scf_energy.irp.f`

    .. code:: fortran

        double precision	:: tc_hf_energy	
        double precision	:: tc_hf_one_e_energy	
        double precision	:: tc_hf_two_e_energy	
        double precision	:: tc_hf_three_e_energy	


    TC Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`ao_two_e_core_integral`
       * :c:data:`ao_two_e_valence_integral_alpha`
       * :c:data:`core_tc_op`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`nol_0e`
       * :c:data:`nol_standard`
       * :c:data:`nuclear_repulsion`
       * :c:data:`tcscf_bi_ort_core_dm_ao`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`three_body_h_tc`
       * :c:data:`two_e_tc_integral_alpha`
       * :c:data:`two_e_tc_integral_alpha`


 
.. c:var:: double tc_hf_two_e_energy


    File : :file:`tc_scf_energy.irp.f`

    .. code:: fortran

        double precision	:: tc_hf_energy	
        double precision	:: tc_hf_one_e_energy	
        double precision	:: tc_hf_two_e_energy	
        double precision	:: tc_hf_three_e_energy	


    TC Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`ao_two_e_core_integral`
       * :c:data:`ao_two_e_valence_integral_alpha`
       * :c:data:`core_tc_op`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`nol_0e`
       * :c:data:`nol_standard`
       * :c:data:`nuclear_repulsion`
       * :c:data:`tcscf_bi_ort_core_dm_ao`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`three_body_h_tc`
       * :c:data:`two_e_tc_integral_alpha`
       * :c:data:`two_e_tc_integral_alpha`


 
.. c:var:: double tcscf_density_matrix_ao_alpha


    File : :file:`tc_scf_dm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tcscf_density_matrix_ao_alpha	(ao_num,ao_num)


    
    TC-SCF transition density matrix on the AO basis for ALPHA electrons
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`tcscf_bi_ort_dm_ao_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_hf_energy`
       * :c:data:`tcscf_density_matrix_ao_tot`
       * :c:data:`two_e_tc_integral_alpha`

 
.. c:var:: double tcscf_density_matrix_ao_beta


    File : :file:`tc_scf_dm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tcscf_density_matrix_ao_beta	(ao_num,ao_num)


    
    TC-SCF transition density matrix on the AO basis for BETA electrons
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`tcscf_bi_ort_dm_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_hf_energy`
       * :c:data:`tcscf_density_matrix_ao_tot`
       * :c:data:`two_e_tc_integral_alpha`

 
.. c:var:: double tcscf_density_matrix_ao_tot


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


 
.. c:var:: double two_e_tc_integral_alpha


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_tc_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: two_e_tc_integral_beta	(ao_num,ao_num)


    
    two_e_tc_integral_alpha(k,i) = <k| F^tc_2e_alpha |i> ON THE AO BASIS
    
    where F^tc_2e is the TWO-BODY part of the TC Fock matrix and k,i are AO basis functions
    
    WARNING: if core_tc_op then only the valence electrons are considered
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_ao_alpha`
       * :c:data:`fock_matrix_tc_ao_beta`
       * :c:data:`tc_hf_energy`

 
.. c:var:: double two_e_tc_integral_beta


    File : :file:`fock_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_tc_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: two_e_tc_integral_beta	(ao_num,ao_num)


    
    two_e_tc_integral_alpha(k,i) = <k| F^tc_2e_alpha |i> ON THE AO BASIS
    
    where F^tc_2e is the TWO-BODY part of the TC Fock matrix and k,i are AO basis functions
    
    WARNING: if core_tc_op then only the valence electrons are considered
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_ao_alpha`
       * :c:data:`fock_matrix_tc_ao_beta`
       * :c:data:`tc_hf_energy`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void build_s_matrix(void)


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

 
.. c:function:: void contrib_3e_soo(void)


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

 
.. c:function:: void contrib_3e_sos(void)


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

 
.. c:function:: void contrib_3e_sss(void)


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

 
.. c:function:: void create_guess(void)


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

 
.. c:function:: void extrapolate_tc_fock_matrix(void)


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

 
.. c:function:: void ltxsxr(void)


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine LTxSxR(n, m, L, S, R, C)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void main(void)


    File : :file:`write_ao_2e_tc_integ.irp.f`

    .. code:: fortran

        subroutine main()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`io_tc_integ`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`tc_scf`
       * :c:func:`write_ao_2e_tc_integ`

 
.. c:function:: void minimize_tc_orb_angles(void)


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

       * :c:func:`minimize_tc_angles`

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

 
.. c:function:: void molden_l(void)


    File : :file:`molden_lr_mos.irp.f`

    .. code:: fortran

        subroutine molden_l()


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
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_list_shell_aos`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`molden_lr_mos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`isort`

 
.. c:function:: void molden_lr(void)


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

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`isort`

 
.. c:function:: void molden_r(void)


    File : :file:`molden_lr_mos.irp.f`

    .. code:: fortran

        subroutine molden_r()


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

       * :c:func:`molden_lr_mos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`isort`

 
.. c:function:: void orthog_functions(void)


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

 
.. c:function:: void print_angles_tc(void)


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
       * :c:func:`print_energy_and_mos`

 
.. c:function:: void print_energy_and_mos(void)


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine print_energy_and_mos(good_angles)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`angle_left_right`
       * :c:data:`grad_non_hermit_left`
       * :c:data:`tc_hf_energy`
       * :c:data:`thresh_lr_angle`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`minimize_tc_orb_angles`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_angles_tc`

 
.. c:function:: void print_strong_overlap(void)


    File : :file:`routines_rotates.irp.f`

    .. code:: fortran

        subroutine print_strong_overlap(i_list, j_list)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`overlap_mo_r`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`routine_save_rotated_mos`

 
.. c:function:: void rh_tcscf_diis(void)


    File : :file:`rh_tcscf_diis.irp.f`

    .. code:: fortran

        subroutine rh_tcscf_diis()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`fock_matrix_tc_mo_tot`
       * :c:data:`fock_tc_reigvec_ao`
       * :c:data:`fqs_sqf_ao`
       * :c:data:`fqs_sqf_mo`
       * :c:data:`json_int_fmt`
       * :c:data:`json_unit`
       * :c:data:`level_shift_tcscf`
       * :c:data:`max_dim_diis_tcscf`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_it_tcscf_max`
       * :c:data:`tc_hf_energy`
       * :c:data:`thresh_tcscf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_to_mo_bi_ortho`
       * :c:func:`extrapolate_tc_fock_matrix`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_l_coef`
       * :c:func:`ezfio_set_bi_ortho_mos_mo_r_coef`
       * :c:func:`ezfio_set_tc_scf_converged_tcscf`
       * :c:func:`ezfio_set_tc_scf_tcscf_energy`
       * :c:func:`lock_io`
       * :c:func:`unlock_io`
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

 
.. c:function:: void routine_save_rotated_mos(void)


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
       * :c:data:`list_core`
       * :c:data:`maxovl_tc`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_core_orb`
       * :c:data:`overlap_mo_r`
       * :c:data:`tc_hf_energy`
       * :c:data:`thresh_de_tc_angles`
       * :c:data:`thresh_lr_angle`

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
       * :c:func:`give_degen_full_listcore`
       * :c:func:`impose_weighted_biorthog_svd`
       * :c:func:`maxovl`
       * :c:func:`orthog_functions`
       * :c:func:`print_strong_overlap`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

 
.. c:function:: void sort_by_tc_fock(void)


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

