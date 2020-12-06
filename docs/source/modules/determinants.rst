.. _module_determinants: 
 
.. program:: determinants 
 
.. default-role:: option 
 
============
determinants
============

Contains everything for the computation of the Hamiltonian matrix elements in the basis of orthogonal Slater determinants built on a restricted spin-orbitals basis.

The main providers for this module are:

* :option:`determinants n_states`: number of states to be computed
* :c:data:`psi_det`: list of determinants in the wave function used in many routines/providers of the |QP|.
* :c:data:`psi_coef`: list of coefficients, for all :option:`determinants n_states` states, and all determinants.

The main routines for this module are:

* :c:func:`i_H_j`: computes the Hamiltonian matrix element between two arbitrary Slater determinants.
* :c:func:`i_H_j_s2`: computes the Hamiltonian and (|S^2|) matrix element between two arbitrary Slater determinants.
* :c:func:`i_H_j_verbose`: returns the decomposition in terms of one- and two-body components of the Hamiltonian matrix elements between two arbitrary Slater determinants. Also return the fermionic phase factor.
* :c:func:`i_H_psi`: computes the Hamiltonian matrix element between an arbitrary Slater determinant and a wave function composed of a sum of arbitrary Slater determinants.


For an example of how to use these routines and providers, take a look at :file:`example.irp.f`.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: n_det_max
 
    Maximum number of determinants in the wave function
 
    Default: 1000000
 
.. option:: n_det_print_wf
 
    Maximum number of determinants to be printed with the program print_wf
 
    Default: 10000
 
.. option:: n_states
 
    Number of states to consider
 
    Default: 1
 
.. option:: read_wf
 
    If |true|, read the wave function from the |EZFIO| file
 
    Default: False
 
.. option:: pruning
 
    If p>0., remove p*Ndet determinants at every iteration
 
    Default: 0.
 
.. option:: s2_eig
 
    Force the wave function to be an eigenfunction of |S^2|
 
    Default: True
 
.. option:: weight_one_e_dm
 
    Weight used in the calculation of the one-electron density matrix. 0: 1./(c_0^2), 1: 1/N_states, 2: input state-average weight, 3: 1/(Norm_L3(Psi))
 
    Default: 2
 
.. option:: weight_selection
 
    Weight used in the selection. 0: input state-average weight, 1: 1./(c_0^2), 2: rPT2 matching, 3: variance matching, 4: variance and rPT2 matching, 5: variance minimization and matching, 6: CI coefficients
 
    Default: 2
 
.. option:: threshold_generators
 
    Thresholds on generators (fraction of the square of the norm)
 
    Default: 0.999
 
.. option:: n_int
 
    Number of integers required to represent bitstrings (set in module :ref:`module_bitmask`)
 
 
.. option:: bit_kind
 
    (set in module :ref:`module_bitmask`)
 
 
.. option:: mo_label
 
    Label of the |MOs| on which the determinants are expressed
 
 
.. option:: n_det
 
    Number of determinants in the current wave function
 
 
.. option:: n_det_qp_edit
 
    Number of determinants to print in qp_edit
 
 
.. option:: psi_coef
 
    Coefficients of the wave function
 
 
.. option:: psi_det
 
    Determinants of the variational space
 
 
.. option:: psi_coef_qp_edit
 
    Coefficients of the wave function
 
 
.. option:: psi_det_qp_edit
 
    Determinants of the variational space
 
 
.. option:: expected_s2
 
    Expected value of |S^2|
 
 
.. option:: target_energy
 
    Energy that should be obtained when truncating the wave function (optional)
 
    Default: 0.
 
.. option:: state_average_weight
 
    Weight of the states in state-average calculations.
 
 
.. option:: selection_factor
 
    f such that the number of determinants to add is f * N_det * sqrt(N_states)
 
    Default: 1.
 
.. option:: thresh_sym
 
    Thresholds to check if a determinant is connected with HF
 
    Default: 1.e-15
 
.. option:: pseudo_sym
 
    If |true|, discard any Slater determinants with an interaction smaller than thresh_sym with HF.
 
    Default: False
 
 
Providers 
--------- 
 
.. c:var:: abs_psi_coef_max


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_coef_max	(N_states)
        double precision, allocatable	:: psi_coef_min	(N_states)
        double precision, allocatable	:: abs_psi_coef_max	(N_states)
        double precision, allocatable	:: abs_psi_coef_min	(N_states)


    Max and min values of the coefficients

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_states`
       * :c:data:`psi_coef`


 
.. c:var:: abs_psi_coef_min


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_coef_max	(N_states)
        double precision, allocatable	:: psi_coef_min	(N_states)
        double precision, allocatable	:: abs_psi_coef_max	(N_states)
        double precision, allocatable	:: abs_psi_coef_min	(N_states)


    Max and min values of the coefficients

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_states`
       * :c:data:`psi_coef`


 
.. c:var:: barycentric_electronic_energy


    File : :file:`determinants/energy.irp.f`

    .. code:: fortran

        double precision, allocatable	:: barycentric_electronic_energy	(N_states)


    :math:`E_n = \sum_i {c_i^{(n)}}^2 H_{ii}` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`


 
.. c:var:: c0_weight


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: c0_weight	(N_states)


    Weight of the states in the selection : :math:`\frac{1}{c_0^2}` .

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`psi_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`state_average_weight`

 
.. c:var:: det_alpha_norm


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: det_alpha_norm	(N_det_alpha_unique)
        double precision, allocatable	:: det_beta_norm	(N_det_beta_unique)


    Norm of the :math:`\alpha`  and :math:`\beta`  spin determinants in the wave function:
    
    :math:`||D_\alpha||_i = \sum_j C_{ij}^2` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`state_average_weight`


 
.. c:var:: det_beta_norm


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: det_alpha_norm	(N_det_alpha_unique)
        double precision, allocatable	:: det_beta_norm	(N_det_beta_unique)


    Norm of the :math:`\alpha`  and :math:`\beta`  spin determinants in the wave function:
    
    :math:`||D_\alpha||_i = \sum_j C_{ij}^2` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`state_average_weight`


 
.. c:var:: det_to_occ_pattern


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        integer, allocatable	:: det_to_occ_pattern	(N_det)


    Returns the index of the occupation pattern for each determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_occ_pattern`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pruned`
       * :c:data:`psi_occ_pattern_hii`
       * :c:data:`weight_occ_pattern`
       * :c:data:`weight_occ_pattern_average`

 
.. c:var:: diagonal_h_matrix_on_psi_det


    File : :file:`determinants/energy.irp.f`

    .. code:: fortran

        double precision, allocatable	:: diagonal_h_matrix_on_psi_det	(N_det)


    Diagonal of the Hamiltonian ordered as psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`barycentric_electronic_energy`

 
.. c:var:: double_exc_bitmask


    File : :file:`determinants/determinants_bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: double_exc_bitmask	(N_int,4,N_double_exc_bitmasks)


    double_exc_bitmask(:,1,i) is the bitmask for holes of excitation 1
    
    double_exc_bitmask(:,2,i) is the bitmask for particles of excitation 1
    
    double_exc_bitmask(:,3,i) is the bitmask for holes of excitation 2
    
    double_exc_bitmask(:,4,i) is the bitmask for particles of excitation 2
    
    for a given couple of hole/particle excitations i.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_double_exc_bitmasks`
       * :c:data:`n_int`


 
.. c:var:: expected_s2


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        double precision	:: expected_s2	


    Expected value of |S^2| : S*(S+1)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`

 
.. c:var:: fock_operator_closed_shell_ref_bitmask


    File : :file:`determinants/single_excitations.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_operator_closed_shell_ref_bitmask	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`full_ijkl_bitmask`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_int`
       * :c:data:`ref_closed_shell_bitmask`


 
.. c:var:: fock_wee_closed_shell


    File : :file:`determinants/single_excitation_two_e.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_wee_closed_shell	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`full_ijkl_bitmask`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_int`
       * :c:data:`ref_closed_shell_bitmask`


 
.. c:var:: h_apply_buffer_allocated


    File : :file:`determinants/h_apply.irp.f`

    .. code:: fortran

        logical	:: h_apply_buffer_allocated	
        integer(omp_lock_kind), allocatable	:: h_apply_buffer_lock	(64,0:nproc-1)


    Buffer of determinants/coefficients/perturbative energy for H_apply.
    Uninitialized. Filled by H_apply subroutines.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nproc`


 
.. c:var:: h_apply_buffer_lock


    File : :file:`determinants/h_apply.irp.f`

    .. code:: fortran

        logical	:: h_apply_buffer_allocated	
        integer(omp_lock_kind), allocatable	:: h_apply_buffer_lock	(64,0:nproc-1)


    Buffer of determinants/coefficients/perturbative energy for H_apply.
    Uninitialized. Filled by H_apply subroutines.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nproc`


 
.. c:var:: h_matrix_all_dets


    File : :file:`determinants/utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: h_matrix_all_dets	(N_det,N_det)


    |H| matrix on the basis of the Slater determinants defined by psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`psi_energy`

 
.. c:var:: h_matrix_cas


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        double precision, allocatable	:: h_matrix_cas	(N_det_cas,N_det_cas)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_int`
       * :c:data:`psi_cas`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_cas_energy`
       * :c:data:`psi_coef_cas_diagonalized`

 
.. c:var:: idx_cas


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_cas	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_cas_coef	(psi_det_size,n_states)
        integer, allocatable	:: idx_cas	(psi_det_size)
        integer	:: n_det_cas	


    |CAS| wave function, defined from the application of the |CAS| bitmask on the
    determinants. idx_cas gives the indice of the |CAS| determinant in psi_det.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`hf_bitmask`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_cas`
       * :c:data:`psi_cas_energy`
       * :c:data:`psi_cas_sorted_bit`
       * :c:data:`psi_coef_cas_diagonalized`
       * :c:data:`psi_non_cas`
       * :c:data:`psi_non_cas_sorted_bit`

 
.. c:var:: idx_non_cas


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_non_cas	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_non_cas_coef	(psi_det_size,n_states)
        integer, allocatable	:: idx_non_cas	(psi_det_size)
        integer	:: n_det_non_cas	


    Set of determinants which are not part of the |CAS|, defined from the application
    of the |CAS| bitmask on the determinants.
    idx_non_cas gives the indice of the determinant in psi_det.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_cas`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_non_cas_sorted_bit`

 
.. c:var:: max_degree_exc


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer	:: max_degree_exc	


    Maximum degree of excitation in the wave function with respect to the Hartree-Fock
    determinant.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`


 
.. c:var:: n_det


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer	:: n_det	


    Number of determinants in the wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_label`
       * :c:data:`mpi_master`
       * :c:data:`nproc`
       * :c:data:`read_wf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`barycentric_electronic_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
       * :c:data:`det_alpha_norm`
       * :c:data:`det_to_occ_pattern`
       * :c:data:`diag_algorithm`
       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`dressed_column_idx`
       * :c:data:`dressing_column_h`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`max_degree_exc`
       * :c:data:`n_det_qp_edit`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`pruned`
       * :c:data:`psi_average_norm_contrib`
       * :c:data:`psi_bilinear_matrix`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_cas`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_hii`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy_two_e`
       * :c:data:`psi_non_cas`
       * :c:data:`psi_occ_pattern`
       * :c:data:`psi_occ_pattern_hii`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`s2_values`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`
       * :c:data:`weight_occ_pattern`
       * :c:data:`weight_occ_pattern_average`

 
.. c:var:: n_det_alpha_unique


    File : :file:`determinants/spindeterminants.irp.f_template_144`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_alpha_unique	(N_int,psi_det_size)
        integer	:: n_det_alpha_unique	


    Unique :math:`\alpha`  determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det_alpha`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`det_alpha_norm`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`singles_alpha_csc`
       * :c:data:`singles_alpha_csc_idx`

 
.. c:var:: n_det_beta_unique


    File : :file:`determinants/spindeterminants.irp.f_template_144`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_beta_unique	(N_int,psi_det_size)
        integer	:: n_det_beta_unique	


    Unique :math:`\beta`  determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det_beta`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`det_alpha_norm`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

 
.. c:var:: n_det_cas


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_cas	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_cas_coef	(psi_det_size,n_states)
        integer, allocatable	:: idx_cas	(psi_det_size)
        integer	:: n_det_cas	


    |CAS| wave function, defined from the application of the |CAS| bitmask on the
    determinants. idx_cas gives the indice of the |CAS| determinant in psi_det.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`hf_bitmask`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_cas`
       * :c:data:`psi_cas_energy`
       * :c:data:`psi_cas_sorted_bit`
       * :c:data:`psi_coef_cas_diagonalized`
       * :c:data:`psi_non_cas`
       * :c:data:`psi_non_cas_sorted_bit`

 
.. c:var:: n_det_non_cas


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_non_cas	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_non_cas_coef	(psi_det_size,n_states)
        integer, allocatable	:: idx_non_cas	(psi_det_size)
        integer	:: n_det_non_cas	


    Set of determinants which are not part of the |CAS|, defined from the application
    of the |CAS| bitmask on the determinants.
    idx_non_cas gives the indice of the determinant in psi_det.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_cas`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_non_cas_sorted_bit`

 
.. c:var:: n_det_qp_edit


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer	:: n_det_qp_edit	


    Number of determinants to print in qp_edit

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`


 
.. c:var:: n_double_exc_bitmasks


    File : :file:`determinants/determinants_bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_double_exc_bitmasks	


    Number of double excitation bitmasks

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`double_exc_bitmask`

 
.. c:var:: n_occ_pattern


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_occ_pattern	(N_int,2,psi_det_size)
        integer	:: n_occ_pattern	


    Array of the occ_patterns present in the wave function.
    
    psi_occ_pattern(:,1,j) = j-th occ_pattern of the wave function : represents all the single occupations
    
    psi_occ_pattern(:,2,j) = j-th occ_pattern of the wave function : represents all the double occupations
    
    The occ patterns are sorted by :c:func:`occ_pattern_search_key`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_occ_pattern`
       * :c:data:`pruned`
       * :c:data:`psi_occ_pattern_hii`
       * :c:data:`psi_occ_pattern_sorted`
       * :c:data:`weight_occ_pattern`
       * :c:data:`weight_occ_pattern_average`

 
.. c:var:: n_single_exc_bitmasks


    File : :file:`determinants/determinants_bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_single_exc_bitmasks	


    Number of single excitation bitmasks

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`single_exc_bitmask`

 
.. c:var:: one_e_dm_ao_alpha


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_ao_alpha	(ao_num,ao_num)
        double precision, allocatable	:: one_e_dm_ao_beta	(ao_num,ao_num)


    One body density matrix on the |AO| basis : :math:`\rho_{AO}(\alpha), \rho_{AO}(\beta)` .

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`one_e_dm_mo_alpha_average`


 
.. c:var:: one_e_dm_ao_beta


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_ao_alpha	(ao_num,ao_num)
        double precision, allocatable	:: one_e_dm_ao_beta	(ao_num,ao_num)


    One body density matrix on the |AO| basis : :math:`\rho_{AO}(\alpha), \rho_{AO}(\beta)` .

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`one_e_dm_mo_alpha_average`


 
.. c:var:: one_e_dm_dagger_mo_spin_index


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_dagger_mo_spin_index	(mo_num,mo_num,N_states,2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha`


 
.. c:var:: one_e_dm_mo


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo	(mo_num,mo_num)


    One-body density matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`one_e_dm_mo_alpha_average`


 
.. c:var:: one_e_dm_mo_alpha


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_alpha	(mo_num,mo_num,N_states)
        double precision, allocatable	:: one_e_dm_mo_beta	(mo_num,mo_num,N_states)


    :math:`\alpha`  and :math:`\beta`  one-body density matrix for each state

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`one_e_dm_dagger_mo_spin_index`
       * :c:data:`one_e_dm_mo_alpha_average`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`one_e_dm_mo_diff`
       * :c:data:`one_e_dm_mo_spin_index`
       * :c:data:`psi_energy_h_core`

 
.. c:var:: one_e_dm_mo_alpha_average


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_alpha_average	(mo_num,mo_num)
        double precision, allocatable	:: one_e_dm_mo_beta_average	(mo_num,mo_num)


    :math:`\alpha`  and :math:`\beta`  one-body density matrix for each state

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_ao_alpha`
       * :c:data:`one_e_dm_mo`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`one_e_spin_density_mo`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`

 
.. c:var:: one_e_dm_mo_beta


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_alpha	(mo_num,mo_num,N_states)
        double precision, allocatable	:: one_e_dm_mo_beta	(mo_num,mo_num,N_states)


    :math:`\alpha`  and :math:`\beta`  one-body density matrix for each state

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`one_e_dm_dagger_mo_spin_index`
       * :c:data:`one_e_dm_mo_alpha_average`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`one_e_dm_mo_diff`
       * :c:data:`one_e_dm_mo_spin_index`
       * :c:data:`psi_energy_h_core`

 
.. c:var:: one_e_dm_mo_beta_average


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_alpha_average	(mo_num,mo_num)
        double precision, allocatable	:: one_e_dm_mo_beta_average	(mo_num,mo_num)


    :math:`\alpha`  and :math:`\beta`  one-body density matrix for each state

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_ao_alpha`
       * :c:data:`one_e_dm_mo`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`one_e_spin_density_mo`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`

 
.. c:var:: one_e_dm_mo_diff


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_diff	(mo_num,mo_num,2:N_states)


    Difference of the one-body density matrix with respect to the ground state

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha`


 
.. c:var:: one_e_dm_mo_spin_index


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_spin_index	(mo_num,mo_num,N_states,2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha`


 
.. c:var:: one_e_spin_density_ao


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_spin_density_ao	(ao_num,ao_num)


    One body spin density matrix on the |AO| basis : :math:`\rho_{AO}(\alpha) - \rho_{AO}(\beta)` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`one_e_spin_density_mo`


 
.. c:var:: one_e_spin_density_mo


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_spin_density_mo	(mo_num,mo_num)


    :math:`\rho(\alpha) - \rho(\beta)` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`one_e_dm_mo_alpha_average`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_spin_density_ao`

 
.. c:var:: pruned


    File : :file:`determinants/prune_wf.irp.f`

    .. code:: fortran

        logical, allocatable	:: pruned	(N_det)


    True if determinant is removed by pruning

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_occ_pattern`
       * :c:data:`n_det`
       * :c:data:`pruning`
       * :c:data:`psi_average_norm_contrib`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_occ_pattern`
       * :c:data:`psi_occ_pattern_sorted`
       * :c:data:`s2_eig`


 
.. c:var:: psi_average_norm_contrib


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_average_norm_contrib	(psi_det_size)


    Contribution of determinants to the state-averaged density.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_size`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pruned`
       * :c:data:`psi_det_sorted`

 
.. c:var:: psi_average_norm_contrib_sorted


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted	(psi_det_size,N_states)
        double precision, allocatable	:: psi_average_norm_contrib_sorted	(psi_det_size)
        integer, allocatable	:: psi_det_sorted_order	(psi_det_size)


    Wave function sorted by determinants contribution to the norm (state-averaged)
    
    psi_det_sorted_order(i) -> k : index in psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_average_norm_contrib`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pruned`

 
.. c:var:: psi_bilinear_matrix


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_bilinear_matrix	(N_det_alpha_unique,N_det_beta_unique,N_states)


    Coefficient matrix if the wave function is expressed in a bilinear form :
    
    :math:`D_\alpha^\dagger.C.D_\beta` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`


 
.. c:var:: psi_bilinear_matrix_columns


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_bilinear_matrix_values	(N_det,N_states)
        integer, allocatable	:: psi_bilinear_matrix_rows	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_columns	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_order	(N_det)


    Sparse coefficient matrix if the wave function is expressed in a bilinear form :
     :math:`D_\alpha^\dagger.C.D_\beta` 
    
    Rows are :math:`\alpha`  determinants and columns are :math:`\beta` .
    
    Order refers to psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`det_alpha_norm`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`

 
.. c:var:: psi_bilinear_matrix_columns_loc


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer, allocatable	:: psi_bilinear_matrix_columns_loc	(N_det_beta_unique+1)


    Sparse coefficient matrix if the wave function is expressed in a bilinear form :
    
    :math:`D_\alpha^\dagger.C.D_\beta` 
    
    Rows are :math:`\alpha`  determinants and columns are :math:`\beta` .
    
    Order refers to :c:data:`psi_det`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_beta_unique`


 
.. c:var:: psi_bilinear_matrix_order


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_bilinear_matrix_values	(N_det,N_states)
        integer, allocatable	:: psi_bilinear_matrix_rows	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_columns	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_order	(N_det)


    Sparse coefficient matrix if the wave function is expressed in a bilinear form :
     :math:`D_\alpha^\dagger.C.D_\beta` 
    
    Rows are :math:`\alpha`  determinants and columns are :math:`\beta` .
    
    Order refers to psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`det_alpha_norm`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`

 
.. c:var:: psi_bilinear_matrix_order_reverse


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer, allocatable	:: psi_bilinear_matrix_order_reverse	(N_det)


    Order which allows to go from :c:data:`psi_bilinear_matrix` to :c:data:`psi_det`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_bilinear_matrix_values`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`

 
.. c:var:: psi_bilinear_matrix_order_transp_reverse


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer, allocatable	:: psi_bilinear_matrix_order_transp_reverse	(N_det)


    Order which allows to go from :c:data:`psi_bilinear_matrix_order_transp` to
    :c:data:`psi_bilinear_matrix`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_bilinear_matrix_transp_values`


 
.. c:var:: psi_bilinear_matrix_rows


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_bilinear_matrix_values	(N_det,N_states)
        integer, allocatable	:: psi_bilinear_matrix_rows	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_columns	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_order	(N_det)


    Sparse coefficient matrix if the wave function is expressed in a bilinear form :
     :math:`D_\alpha^\dagger.C.D_\beta` 
    
    Rows are :math:`\alpha`  determinants and columns are :math:`\beta` .
    
    Order refers to psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`det_alpha_norm`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`

 
.. c:var:: psi_bilinear_matrix_transp_columns


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_bilinear_matrix_transp_values	(N_det,N_states)
        integer, allocatable	:: psi_bilinear_matrix_transp_rows	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_transp_columns	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_transp_order	(N_det)


    Transpose of :c:data:`psi_bilinear_matrix`
    
    :math:`D_\beta^\dagger.C^\dagger.D_\alpha` 
    
    Rows are :math:`\alpha`  determinants and columns are :math:`\beta` , but the matrix is stored in row major
    format.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`

 
.. c:var:: psi_bilinear_matrix_transp_order


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_bilinear_matrix_transp_values	(N_det,N_states)
        integer, allocatable	:: psi_bilinear_matrix_transp_rows	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_transp_columns	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_transp_order	(N_det)


    Transpose of :c:data:`psi_bilinear_matrix`
    
    :math:`D_\beta^\dagger.C^\dagger.D_\alpha` 
    
    Rows are :math:`\alpha`  determinants and columns are :math:`\beta` , but the matrix is stored in row major
    format.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`

 
.. c:var:: psi_bilinear_matrix_transp_rows


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_bilinear_matrix_transp_values	(N_det,N_states)
        integer, allocatable	:: psi_bilinear_matrix_transp_rows	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_transp_columns	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_transp_order	(N_det)


    Transpose of :c:data:`psi_bilinear_matrix`
    
    :math:`D_\beta^\dagger.C^\dagger.D_\alpha` 
    
    Rows are :math:`\alpha`  determinants and columns are :math:`\beta` , but the matrix is stored in row major
    format.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`

 
.. c:var:: psi_bilinear_matrix_transp_rows_loc


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer, allocatable	:: psi_bilinear_matrix_transp_rows_loc	(N_det_alpha_unique+1)


    Location of the columns in the :c:data:`psi_bilinear_matrix`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_det_alpha_unique`


 
.. c:var:: psi_bilinear_matrix_transp_values


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_bilinear_matrix_transp_values	(N_det,N_states)
        integer, allocatable	:: psi_bilinear_matrix_transp_rows	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_transp_columns	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_transp_order	(N_det)


    Transpose of :c:data:`psi_bilinear_matrix`
    
    :math:`D_\beta^\dagger.C^\dagger.D_\alpha` 
    
    Rows are :math:`\alpha`  determinants and columns are :math:`\beta` , but the matrix is stored in row major
    format.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`

 
.. c:var:: psi_bilinear_matrix_values


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_bilinear_matrix_values	(N_det,N_states)
        integer, allocatable	:: psi_bilinear_matrix_rows	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_columns	(N_det)
        integer, allocatable	:: psi_bilinear_matrix_order	(N_det)


    Sparse coefficient matrix if the wave function is expressed in a bilinear form :
     :math:`D_\alpha^\dagger.C.D_\beta` 
    
    Rows are :math:`\alpha`  determinants and columns are :math:`\beta` .
    
    Order refers to psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`det_alpha_norm`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`

 
.. c:var:: psi_cas


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_cas	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_cas_coef	(psi_det_size,n_states)
        integer, allocatable	:: idx_cas	(psi_det_size)
        integer	:: n_det_cas	


    |CAS| wave function, defined from the application of the |CAS| bitmask on the
    determinants. idx_cas gives the indice of the |CAS| determinant in psi_det.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`hf_bitmask`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_cas`
       * :c:data:`psi_cas_energy`
       * :c:data:`psi_cas_sorted_bit`
       * :c:data:`psi_coef_cas_diagonalized`
       * :c:data:`psi_non_cas`
       * :c:data:`psi_non_cas_sorted_bit`

 
.. c:var:: psi_cas_coef


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_cas	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_cas_coef	(psi_det_size,n_states)
        integer, allocatable	:: idx_cas	(psi_det_size)
        integer	:: n_det_cas	


    |CAS| wave function, defined from the application of the |CAS| bitmask on the
    determinants. idx_cas gives the indice of the |CAS| determinant in psi_det.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`hf_bitmask`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_cas`
       * :c:data:`psi_cas_energy`
       * :c:data:`psi_cas_sorted_bit`
       * :c:data:`psi_coef_cas_diagonalized`
       * :c:data:`psi_non_cas`
       * :c:data:`psi_non_cas_sorted_bit`

 
.. c:var:: psi_cas_coef_sorted_bit


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_cas_sorted_bit	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_cas_coef_sorted_bit	(psi_det_size,N_states)


    |CAS| determinants sorted to accelerate the search of a random determinant in the wave
    function.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_cas`
       * :c:data:`psi_det_size`


 
.. c:var:: psi_cas_energy


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_cas_energy	(N_states)


    Variational energy of :math:`\Psi_{CAS}` , where :math:`\Psi_{CAS} =  \sum_{I \in CAS} \I \rangle \langle I | \Psi \rangle` .

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_cas`
       * :c:data:`n_states`
       * :c:data:`psi_cas`


 
.. c:var:: psi_cas_energy_diagonalized


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_coef_cas_diagonalized	(N_det_cas,N_states)
        double precision, allocatable	:: psi_cas_energy_diagonalized	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_cas`
       * :c:data:`n_states`
       * :c:data:`psi_cas`


 
.. c:var:: psi_cas_sorted_bit


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_cas_sorted_bit	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_cas_coef_sorted_bit	(psi_det_size,N_states)


    |CAS| determinants sorted to accelerate the search of a random determinant in the wave
    function.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_cas`
       * :c:data:`psi_det_size`


 
.. c:var:: psi_coef


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_coef	(psi_det_size,N_states)


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
       * :c:data:`psi_det_size`
       * :c:data:`read_wf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`barycentric_electronic_energy`
       * :c:data:`c0_weight`
       * :c:data:`ci_electronic_energy`
       * :c:data:`dressed_column_idx`
       * :c:data:`psi_average_norm_contrib`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_cas`
       * :c:data:`psi_coef_max`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy_two_e`
       * :c:data:`psi_non_cas`
       * :c:data:`s2_values`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`
       * :c:data:`weight_occ_pattern`
       * :c:data:`weight_occ_pattern_average`

 
.. c:var:: psi_coef_cas_diagonalized


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_coef_cas_diagonalized	(N_det_cas,N_states)
        double precision, allocatable	:: psi_cas_energy_diagonalized	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_cas`
       * :c:data:`n_states`
       * :c:data:`psi_cas`


 
.. c:var:: psi_coef_max


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_coef_max	(N_states)
        double precision, allocatable	:: psi_coef_min	(N_states)
        double precision, allocatable	:: abs_psi_coef_max	(N_states)
        double precision, allocatable	:: abs_psi_coef_min	(N_states)


    Max and min values of the coefficients

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_states`
       * :c:data:`psi_coef`


 
.. c:var:: psi_coef_min


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_coef_max	(N_states)
        double precision, allocatable	:: psi_coef_min	(N_states)
        double precision, allocatable	:: abs_psi_coef_max	(N_states)
        double precision, allocatable	:: abs_psi_coef_min	(N_states)


    Max and min values of the coefficients

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_states`
       * :c:data:`psi_coef`


 
.. c:var:: psi_coef_sorted


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted	(psi_det_size,N_states)
        double precision, allocatable	:: psi_average_norm_contrib_sorted	(psi_det_size)
        integer, allocatable	:: psi_det_sorted_order	(psi_det_size)


    Wave function sorted by determinants contribution to the norm (state-averaged)
    
    psi_det_sorted_order(i) -> k : index in psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_average_norm_contrib`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pruned`

 
.. c:var:: psi_coef_sorted_bit


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_bit	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_bit	(psi_det_size,N_states)


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`

 
.. c:var:: psi_det


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det	(N_int,2,psi_det_size)


    The determinants of the wave function. Initialized with Hartree-Fock if the |EZFIO| file
    is empty.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`hf_bitmask`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det_size`
       * :c:data:`read_wf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`det_to_occ_pattern`
       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`max_degree_exc`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_cas`
       * :c:data:`psi_det_alpha`
       * :c:data:`psi_det_beta`
       * :c:data:`psi_det_hii`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy_two_e`
       * :c:data:`psi_non_cas`
       * :c:data:`psi_occ_pattern`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`s2_values`

 
.. c:var:: psi_det_alpha


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_alpha	(N_int,psi_det_size)


    List of :math:`\alpha`  determinants of psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_alpha_unique`

 
.. c:var:: psi_det_alpha_unique


    File : :file:`determinants/spindeterminants.irp.f_template_144`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_alpha_unique	(N_int,psi_det_size)
        integer	:: n_det_alpha_unique	


    Unique :math:`\alpha`  determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det_alpha`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`det_alpha_norm`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`singles_alpha_csc`
       * :c:data:`singles_alpha_csc_idx`

 
.. c:var:: psi_det_beta


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_beta	(N_int,psi_det_size)


    List of :math:`\beta`  determinants of psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_beta_unique`

 
.. c:var:: psi_det_beta_unique


    File : :file:`determinants/spindeterminants.irp.f_template_144`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_beta_unique	(N_int,psi_det_size)
        integer	:: n_det_beta_unique	


    Unique :math:`\beta`  determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det_beta`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`det_alpha_norm`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

 
.. c:var:: psi_det_hii


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_det_hii	(N_det)


    :math:`\langle i|h|i \rangle`  for all determinants.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_occ_pattern_hii`

 
.. c:var:: psi_det_size


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer	:: psi_det_size	


    Size of the psi_det and psi_coef arrays

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_average_norm_contrib`
       * :c:data:`psi_cas`
       * :c:data:`psi_cas_sorted_bit`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy_two_e`
       * :c:data:`psi_non_cas`
       * :c:data:`psi_non_cas_sorted_bit`
       * :c:data:`psi_occ_pattern`
       * :c:data:`s2_values`

 
.. c:var:: psi_det_sorted


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted	(psi_det_size,N_states)
        double precision, allocatable	:: psi_average_norm_contrib_sorted	(psi_det_size)
        integer, allocatable	:: psi_det_sorted_order	(psi_det_size)


    Wave function sorted by determinants contribution to the norm (state-averaged)
    
    psi_det_sorted_order(i) -> k : index in psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_average_norm_contrib`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pruned`

 
.. c:var:: psi_det_sorted_bit


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_bit	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_bit	(psi_det_size,N_states)


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`

 
.. c:var:: psi_det_sorted_order


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted	(psi_det_size,N_states)
        double precision, allocatable	:: psi_average_norm_contrib_sorted	(psi_det_size)
        integer, allocatable	:: psi_det_sorted_order	(psi_det_size)


    Wave function sorted by determinants contribution to the norm (state-averaged)
    
    psi_det_sorted_order(i) -> k : index in psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_average_norm_contrib`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pruned`

 
.. c:var:: psi_energy_h_core


    File : :file:`determinants/psi_energy_mono_elec.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_energy_h_core	(N_states)


    psi_energy_h_core = :math:`\langle \Psi | h_{core} |\Psi \rangle` 
    
    computed using the :c:data:`one_e_dm_mo_alpha` +
    :c:data:`one_e_dm_mo_beta` and :c:data:`mo_one_e_integrals`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha`


 
.. c:var:: psi_non_cas


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_non_cas	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_non_cas_coef	(psi_det_size,n_states)
        integer, allocatable	:: idx_non_cas	(psi_det_size)
        integer	:: n_det_non_cas	


    Set of determinants which are not part of the |CAS|, defined from the application
    of the |CAS| bitmask on the determinants.
    idx_non_cas gives the indice of the determinant in psi_det.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_cas`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_non_cas_sorted_bit`

 
.. c:var:: psi_non_cas_coef


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_non_cas	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_non_cas_coef	(psi_det_size,n_states)
        integer, allocatable	:: idx_non_cas	(psi_det_size)
        integer	:: n_det_non_cas	


    Set of determinants which are not part of the |CAS|, defined from the application
    of the |CAS| bitmask on the determinants.
    idx_non_cas gives the indice of the determinant in psi_det.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_cas`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_non_cas_sorted_bit`

 
.. c:var:: psi_non_cas_coef_sorted_bit


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_non_cas_sorted_bit	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_non_cas_coef_sorted_bit	(psi_det_size,N_states)


    |CAS| determinants sorted to accelerate the search of a random determinant in the wave
    function.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_cas`
       * :c:data:`psi_det_size`
       * :c:data:`psi_non_cas`


 
.. c:var:: psi_non_cas_sorted_bit


    File : :file:`determinants/psi_cas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_non_cas_sorted_bit	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_non_cas_coef_sorted_bit	(psi_det_size,N_states)


    |CAS| determinants sorted to accelerate the search of a random determinant in the wave
    function.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_cas`
       * :c:data:`psi_det_size`
       * :c:data:`psi_non_cas`


 
.. c:var:: psi_occ_pattern


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_occ_pattern	(N_int,2,psi_det_size)
        integer	:: n_occ_pattern	


    Array of the occ_patterns present in the wave function.
    
    psi_occ_pattern(:,1,j) = j-th occ_pattern of the wave function : represents all the single occupations
    
    psi_occ_pattern(:,2,j) = j-th occ_pattern of the wave function : represents all the double occupations
    
    The occ patterns are sorted by :c:func:`occ_pattern_search_key`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_occ_pattern`
       * :c:data:`pruned`
       * :c:data:`psi_occ_pattern_hii`
       * :c:data:`psi_occ_pattern_sorted`
       * :c:data:`weight_occ_pattern`
       * :c:data:`weight_occ_pattern_average`

 
.. c:var:: psi_occ_pattern_hii


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_occ_pattern_hii	(N_occ_pattern)


    :math:`\langle I|H|I \rangle`  where :math:`|I\rangle`  is an occupation pattern.
    This is the minimum :math:`H_{ii}` , where the :math:`|i\rangle`  are the
    determinants of :math:`|I\rangle` .

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_occ_pattern`
       * :c:data:`n_det`
       * :c:data:`psi_det_hii`
       * :c:data:`psi_occ_pattern`


 
.. c:var:: psi_occ_pattern_sorted


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_occ_pattern_sorted	(N_int,2,N_occ_pattern)
        double precision, allocatable	:: weight_occ_pattern_average_sorted	(N_occ_pattern)
        integer, allocatable	:: psi_occ_pattern_sorted_order	(N_occ_pattern)
        integer, allocatable	:: psi_occ_pattern_sorted_order_reverse	(N_occ_pattern)


    Occupation patterns sorted by weight

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_occ_pattern`
       * :c:data:`weight_occ_pattern_average`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pruned`

 
.. c:var:: psi_occ_pattern_sorted_order


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_occ_pattern_sorted	(N_int,2,N_occ_pattern)
        double precision, allocatable	:: weight_occ_pattern_average_sorted	(N_occ_pattern)
        integer, allocatable	:: psi_occ_pattern_sorted_order	(N_occ_pattern)
        integer, allocatable	:: psi_occ_pattern_sorted_order_reverse	(N_occ_pattern)


    Occupation patterns sorted by weight

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_occ_pattern`
       * :c:data:`weight_occ_pattern_average`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pruned`

 
.. c:var:: psi_occ_pattern_sorted_order_reverse


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_occ_pattern_sorted	(N_int,2,N_occ_pattern)
        double precision, allocatable	:: weight_occ_pattern_average_sorted	(N_occ_pattern)
        integer, allocatable	:: psi_occ_pattern_sorted_order	(N_occ_pattern)
        integer, allocatable	:: psi_occ_pattern_sorted_order_reverse	(N_occ_pattern)


    Occupation patterns sorted by weight

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_occ_pattern`
       * :c:data:`weight_occ_pattern_average`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pruned`

 
.. c:var:: ref_bitmask_energy


    File : :file:`determinants/ref_bitmask.irp.f`

    .. code:: fortran

        double precision	:: ref_bitmask_energy	
        double precision	:: ref_bitmask_one_e_energy	
        double precision	:: ref_bitmask_kinetic_energy	
        double precision	:: ref_bitmask_n_e_energy	
        double precision	:: ref_bitmask_two_e_energy	
        double precision	:: ref_bitmask_energy_ab	
        double precision	:: ref_bitmask_energy_bb	
        double precision	:: ref_bitmask_energy_aa	


    Energy of the reference bitmask used in Slater rules

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`psi_det_hii`

 
.. c:var:: ref_bitmask_energy_aa


    File : :file:`determinants/ref_bitmask.irp.f`

    .. code:: fortran

        double precision	:: ref_bitmask_energy	
        double precision	:: ref_bitmask_one_e_energy	
        double precision	:: ref_bitmask_kinetic_energy	
        double precision	:: ref_bitmask_n_e_energy	
        double precision	:: ref_bitmask_two_e_energy	
        double precision	:: ref_bitmask_energy_ab	
        double precision	:: ref_bitmask_energy_bb	
        double precision	:: ref_bitmask_energy_aa	


    Energy of the reference bitmask used in Slater rules

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`psi_det_hii`

 
.. c:var:: ref_bitmask_energy_ab


    File : :file:`determinants/ref_bitmask.irp.f`

    .. code:: fortran

        double precision	:: ref_bitmask_energy	
        double precision	:: ref_bitmask_one_e_energy	
        double precision	:: ref_bitmask_kinetic_energy	
        double precision	:: ref_bitmask_n_e_energy	
        double precision	:: ref_bitmask_two_e_energy	
        double precision	:: ref_bitmask_energy_ab	
        double precision	:: ref_bitmask_energy_bb	
        double precision	:: ref_bitmask_energy_aa	


    Energy of the reference bitmask used in Slater rules

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`psi_det_hii`

 
.. c:var:: ref_bitmask_energy_bb


    File : :file:`determinants/ref_bitmask.irp.f`

    .. code:: fortran

        double precision	:: ref_bitmask_energy	
        double precision	:: ref_bitmask_one_e_energy	
        double precision	:: ref_bitmask_kinetic_energy	
        double precision	:: ref_bitmask_n_e_energy	
        double precision	:: ref_bitmask_two_e_energy	
        double precision	:: ref_bitmask_energy_ab	
        double precision	:: ref_bitmask_energy_bb	
        double precision	:: ref_bitmask_energy_aa	


    Energy of the reference bitmask used in Slater rules

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`psi_det_hii`

 
.. c:var:: ref_bitmask_kinetic_energy


    File : :file:`determinants/ref_bitmask.irp.f`

    .. code:: fortran

        double precision	:: ref_bitmask_energy	
        double precision	:: ref_bitmask_one_e_energy	
        double precision	:: ref_bitmask_kinetic_energy	
        double precision	:: ref_bitmask_n_e_energy	
        double precision	:: ref_bitmask_two_e_energy	
        double precision	:: ref_bitmask_energy_ab	
        double precision	:: ref_bitmask_energy_bb	
        double precision	:: ref_bitmask_energy_aa	


    Energy of the reference bitmask used in Slater rules

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`psi_det_hii`

 
.. c:var:: ref_bitmask_n_e_energy


    File : :file:`determinants/ref_bitmask.irp.f`

    .. code:: fortran

        double precision	:: ref_bitmask_energy	
        double precision	:: ref_bitmask_one_e_energy	
        double precision	:: ref_bitmask_kinetic_energy	
        double precision	:: ref_bitmask_n_e_energy	
        double precision	:: ref_bitmask_two_e_energy	
        double precision	:: ref_bitmask_energy_ab	
        double precision	:: ref_bitmask_energy_bb	
        double precision	:: ref_bitmask_energy_aa	


    Energy of the reference bitmask used in Slater rules

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`psi_det_hii`

 
.. c:var:: ref_bitmask_one_e_energy


    File : :file:`determinants/ref_bitmask.irp.f`

    .. code:: fortran

        double precision	:: ref_bitmask_energy	
        double precision	:: ref_bitmask_one_e_energy	
        double precision	:: ref_bitmask_kinetic_energy	
        double precision	:: ref_bitmask_n_e_energy	
        double precision	:: ref_bitmask_two_e_energy	
        double precision	:: ref_bitmask_energy_ab	
        double precision	:: ref_bitmask_energy_bb	
        double precision	:: ref_bitmask_energy_aa	


    Energy of the reference bitmask used in Slater rules

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`psi_det_hii`

 
.. c:var:: ref_bitmask_two_e_energy


    File : :file:`determinants/ref_bitmask.irp.f`

    .. code:: fortran

        double precision	:: ref_bitmask_energy	
        double precision	:: ref_bitmask_one_e_energy	
        double precision	:: ref_bitmask_kinetic_energy	
        double precision	:: ref_bitmask_n_e_energy	
        double precision	:: ref_bitmask_two_e_energy	
        double precision	:: ref_bitmask_energy_ab	
        double precision	:: ref_bitmask_energy_bb	
        double precision	:: ref_bitmask_energy_aa	


    Energy of the reference bitmask used in Slater rules

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`psi_det_hii`

 
.. c:var:: ref_closed_shell_bitmask


    File : :file:`determinants/single_excitations.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: ref_closed_shell_bitmask	(N_int,2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`

 
.. c:var:: s2_matrix_all_dets


    File : :file:`determinants/utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: s2_matrix_all_dets	(N_det,N_det)


    |S^2| matrix on the basis of the Slater determinants defined by psi_det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`psi_energy`

 
.. c:var:: s2_values


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        double precision, allocatable	:: s2_values	(N_states)


    array of the averaged values of the S^2 operator on the various states

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`


 
.. c:var:: s_z


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        double precision	:: s_z	
        double precision	:: s_z2_sz	


    z component of the Spin

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`


 
.. c:var:: s_z2_sz


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        double precision	:: s_z	
        double precision	:: s_z2_sz	


    z component of the Spin

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`


 
.. c:var:: single_exc_bitmask


    File : :file:`determinants/determinants_bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: single_exc_bitmask	(N_int,2,N_single_exc_bitmasks)


    single_exc_bitmask(:,1,i) is the bitmask for holes
    
    single_exc_bitmask(:,2,i) is the bitmask for particles
    
    for a given couple of hole/particle excitations i.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_int`
       * :c:data:`n_single_exc_bitmasks`


 
.. c:var:: singles_alpha_csc


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer, allocatable	:: singles_alpha_csc	(singles_alpha_csc_size)


    Indices of all single excitations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`singles_alpha_csc_idx`


 
.. c:var:: singles_alpha_csc_idx


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer*8, allocatable	:: singles_alpha_csc_idx	(N_det_alpha_unique+1)
        integer*8	:: singles_alpha_csc_size	


    singles_alpha_csc_size : Dimension of the :c:data:`singles_alpha_csc` array
    
    singles_alpha_csc_idx  : Index where the single excitations of determinant i start

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det_alpha_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`singles_alpha_csc`

 
.. c:var:: singles_alpha_csc_size


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer*8, allocatable	:: singles_alpha_csc_idx	(N_det_alpha_unique+1)
        integer*8	:: singles_alpha_csc_size	


    singles_alpha_csc_size : Dimension of the :c:data:`singles_alpha_csc` array
    
    singles_alpha_csc_idx  : Index where the single excitations of determinant i start

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det_alpha_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`singles_alpha_csc`

 
.. c:var:: singles_beta_csc


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer, allocatable	:: singles_beta_csc	(singles_beta_csc_size)


    Indices of all single excitations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc_idx`


 
.. c:var:: singles_beta_csc_idx


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer*8, allocatable	:: singles_beta_csc_idx	(N_det_beta_unique+1)
        integer*8	:: singles_beta_csc_size	


    singles_beta_csc_size : Dimension of the :c:data:`singles_beta_csc` array
    
    singles_beta_csc_idx  : Index where the single excitations of determinant i start

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`singles_beta_csc`

 
.. c:var:: singles_beta_csc_size


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer*8, allocatable	:: singles_beta_csc_idx	(N_det_beta_unique+1)
        integer*8	:: singles_beta_csc_size	


    singles_beta_csc_size : Dimension of the :c:data:`singles_beta_csc` array
    
    singles_beta_csc_idx  : Index where the single excitations of determinant i start

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det_beta_unique`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`singles_beta_csc`

 
.. c:var:: state_average_weight


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: state_average_weight	(N_states)


    Weights in the state-average calculation of the density matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`c0_weight`
       * :c:data:`n_states`
       * :c:data:`weight_one_e_dm`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`det_alpha_norm`
       * :c:data:`one_e_dm_average_alpha_mo_for_dft`
       * :c:data:`one_e_dm_average_beta_mo_for_dft`
       * :c:data:`one_e_dm_mo_alpha_average`
       * :c:data:`psi_average_norm_contrib`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`
       * :c:data:`weight_occ_pattern_average`

 
.. c:var:: weight_occ_pattern


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        double precision, allocatable	:: weight_occ_pattern	(N_occ_pattern,N_states)


    Weight of the occupation patterns in the wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_occ_pattern`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_occ_pattern`


 
.. c:var:: weight_occ_pattern_average


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        double precision, allocatable	:: weight_occ_pattern_average	(N_occ_pattern)


    State-average weight of the occupation patterns in the wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_occ_pattern`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_occ_pattern`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_occ_pattern_sorted`

 
.. c:var:: weight_occ_pattern_average_sorted


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_occ_pattern_sorted	(N_int,2,N_occ_pattern)
        double precision, allocatable	:: weight_occ_pattern_average_sorted	(N_occ_pattern)
        integer, allocatable	:: psi_occ_pattern_sorted_order	(N_occ_pattern)
        integer, allocatable	:: psi_occ_pattern_sorted_order_reverse	(N_occ_pattern)


    Occupation patterns sorted by weight

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_occ_pattern`
       * :c:data:`weight_occ_pattern_average`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pruned`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: a_operator:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine a_operator(iorb,ispin,key,hjj,Nint,na,nb)


    Needed for :c:func:`diag_H_mat_elem`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_h_mat_elem`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: a_operator_two_e:


    File : :file:`determinants/slater_rules_wee_mono.irp.f`

    .. code:: fortran

        subroutine a_operator_two_e(iorb,ispin,key,hjj,Nint,na,nb)


    Needed for :c:func:`diag_Wee_mat_elem`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_jj`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_wee_mat_elem`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: ac_operator:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine ac_operator(iorb,ispin,key,hjj,Nint,na,nb)


    Needed for :c:func:`diag_H_mat_elem`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_h_mat_elem`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: ac_operator_two_e:


    File : :file:`determinants/slater_rules_wee_mono.irp.f`

    .. code:: fortran

        subroutine ac_operator_two_e(iorb,ispin,key,hjj,Nint,na,nb)


    Needed for :c:func:`diag_Wee_mat_elem`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_two_e_integrals_jj`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_wee_mat_elem`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: apply_excitation:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine apply_excitation(det, exc, res, ok, Nint)



 
.. c:function:: apply_hole:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine apply_hole(det, s1, h1, res, ok, Nint)



 
.. c:function:: apply_holes:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine apply_holes(det, s1, h1, s2, h2, res, ok, Nint)



 
.. c:function:: apply_particle:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine apply_particle(det, s1, p1, res, ok, Nint)



 
.. c:function:: apply_particles:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine apply_particles(det, s1, p1, s2, p2, res, ok, Nint)



 
.. c:function:: bitstring_to_list_ab:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine bitstring_to_list_ab( string, list, n_elements, Nint)


    Gives the inidices(+1) of the bits set to 1 in the bit string
    For alpha/beta determinants.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`a_operator`
       * :c:func:`a_operator_two_e`
       * :c:func:`ac_operator`
       * :c:func:`ac_operator_two_e`
       * :c:func:`build_fock_tmp`
       * :c:func:`diag_h_mat_elem`
       * :c:func:`diag_h_mat_elem_one_e`
       * :c:func:`diag_wee_mat_elem`
       * :c:func:`example_determinants`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:func:`get_occupation_from_dets`
       * :c:func:`get_single_excitation_from_fock`
       * :c:func:`i_h_j`
       * :c:func:`i_h_j_s2`
       * :c:func:`i_h_j_two_e`
       * :c:func:`i_h_j_verbose`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:func:`orb_range_diag_to_all_2_rdm_dm_buffer`
       * :c:func:`orb_range_diag_to_all_states_2_rdm_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_dm_buffer`
       * :c:data:`ref_closed_shell_bitmask`
       * :c:func:`single_excitation_wee`

 
.. c:function:: build_fock_tmp:


    File : :file:`determinants/fock_diag.irp.f`

    .. code:: fortran

        subroutine build_fock_tmp(fock_diag_tmp,det_ref,Nint)


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

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`debug_det`

 
.. c:function:: build_singly_excited_wavefunction:


    File : :file:`determinants/create_excitations.irp.f`

    .. code:: fortran

        subroutine build_singly_excited_wavefunction(i_hole,i_particle,ispin,det_out,coef_out)


    Applies the single excitation operator : a^{dager}_(i_particle) a_(i_hole) of
    spin = ispin to the current wave function (psi_det, psi_coef)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`do_single_excitation`
       * :c:func:`get_phase`

 
.. c:function:: connected_to_hf:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine connected_to_hf(key_i,yes_no)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`
       * :c:data:`thresh_sym`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`
       * :c:func:`i_h_j`

 
.. c:function:: connected_to_ref:


    File : :file:`determinants/connected_to_ref.irp.f`

    .. code:: fortran

        integer function connected_to_ref(key,keys,Nint,N_past_in,Ndet)


    input  : key : a given Slater determinant
    
           : keys: a list of Slater determinants
    
           : Ndet: the number of Slater determinants in keys
    
           : N_past_in the number of Slater determinants for the connectivity research
    
    output :   0 : key not connected to the N_past_in first Slater determinants in keys
    
               i : key is connected to determinant i of keys
    
              -i : key is the ith determinant of the reference wf keys

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: connected_to_ref_by_single:


    File : :file:`determinants/connected_to_ref.irp.f`

    .. code:: fortran

        integer function connected_to_ref_by_single(key,keys,Nint,N_past_in,Ndet)


    Returns |true| is ``key`` is connected to the reference by a single excitation.
    input  : key : a given Slater determinant
    
           : keys: a list of Slater determinants
    
           : Ndet: the number of Slater determinants in keys
    
           : N_past_in the number of Slater determinants for the connectivity research
    
    output :   0 : key not connected by a MONO EXCITATION to the N_past_in first Slater determinants in keys
    
               i : key is connected by a MONO EXCITATION to determinant i of keys
    
              -i : key is the ith determinant of the reference wf keys

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: copy_h_apply_buffer_to_wf:


    File : :file:`determinants/h_apply.irp.f`

    .. code:: fortran

        subroutine copy_H_apply_buffer_to_wf


    Copies the H_apply buffer to psi_coef.
    After calling this subroutine, N_det, psi_det and psi_coef need to be touched

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nproc`
       * :c:data:`pruned`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`generate_all_alpha_beta_det_products`
       * :c:func:`make_s2_eigenfunction`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`normalize`
       * :c:func:`remove_duplicates_in_psi_det`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`c0_weight`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit`

 
.. c:function:: copy_psi_bilinear_to_psi:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine copy_psi_bilinear_to_psi(psi, isize)


    Overwrites :c:data:`psi_det` and :c:data:`psi_coef` with the wave function
    in bilinear order

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

 
.. c:function:: create_microlist:


    File : :file:`determinants/filter_connected.irp.f`

    .. code:: fortran

        subroutine create_microlist(minilist, N_minilist, key_mask, microlist, idx_microlist, N_microlist, ptr_microlist, Nint)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list`

 
.. c:function:: create_minilist:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine create_minilist(key_mask, fullList, miniList, idx_miniList, N_fullList, N_miniList, Nint)



 
.. c:function:: create_minilist_find_previous:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine create_minilist_find_previous(key_mask, fullList, miniList, N_fullList, N_miniList, fullMatch, Nint)



 
.. c:function:: create_wf_of_psi_bilinear_matrix:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine create_wf_of_psi_bilinear_matrix(truncate)


    Generates a wave function containing all possible products
    of $\alpha$ and $\beta$ determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det_sorted_bit`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`generate_all_alpha_beta_det_products`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`c0_weight`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit`

 
.. c:function:: decode_exc:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)


    Decodes the exc arrays returned by get_excitation.
    h1,h2 : Holes
    p1,p2 : Particles
    s1,s2 : Spins (1:alpha, 2:beta)
    degree : Degree of excitation

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_h_mat_elem_fock`
       * :c:func:`example_determinants`

 
.. c:function:: decode_exc_spin:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine decode_exc_spin(exc,h1,p1,h2,p2)


    Decodes the exc arrays returned by get_excitation.
    
    h1,h2 : Holes
    
    p1,p2 : Particles

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_mo_alpha`

 
.. c:function:: det_inf:


    File : :file:`determinants/sort_dets_ab.irp.f`

    .. code:: fortran

        logical function det_inf(key1, key2, Nint)


    Ordering function for determinants.

 
.. c:function:: det_search_key:


    File : :file:`determinants/connected_to_ref.irp.f`

    .. code:: fortran

        integer*8 function det_search_key(det,Nint)


    Return an integer*8 corresponding to a determinant index for searching

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`

 
.. c:function:: detcmp:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        integer function detCmp(a,b,Nint)



 
.. c:function:: deteq:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        logical function detEq(a,b,Nint)



 
.. c:function:: diag_h_mat_elem:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        double precision function diag_H_mat_elem(det_in,Nint)


    Computes $\langle i|H|i \rangle$.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`a_operator`
       * :c:func:`ac_operator`
       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: diag_h_mat_elem_fock:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        double precision function diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)


    Computes $\langle i|H|i \rangle$ when $i$ is at most a double excitation from
    a reference.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_jj`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`decode_exc`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`

 
.. c:function:: diag_h_mat_elem_one_e:


    File : :file:`determinants/slater_rules_wee_mono.irp.f`

    .. code:: fortran

        double precision function diag_H_mat_elem_one_e(det_in,Nint)


    Computes $\langle i|H|i \rangle$.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_one_e_integrals`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: diag_s_mat_elem:


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        double precision function diag_S_mat_elem(key_i,Nint)


    Returns <i|S^2|i>

 
.. c:function:: diag_wee_mat_elem:


    File : :file:`determinants/slater_rules_wee_mono.irp.f`

    .. code:: fortran

        double precision function diag_wee_mat_elem(det_in,Nint)


    Computes $\langle i|H|i \rangle$.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`a_operator_two_e`
       * :c:func:`ac_operator_two_e`
       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: do_single_excitation:


    File : :file:`determinants/create_excitations.irp.f`

    .. code:: fortran

        subroutine do_single_excitation(key_in,i_hole,i_particle,ispin,i_ok)


    Apply the single excitation operator : a^{dager}_(i_particle) a_(i_hole) of spin = ispin
    on key_in
    ispin = 1  == alpha
    ispin = 2  == beta
    i_ok = 1  == the excitation is possible
    i_ok = -1 == the excitation is not possible

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`build_singly_excited_wavefunction`
       * :c:func:`example_determinants`

 
.. c:function:: example_determinants:


    File : :file:`determinants/example.irp.f`

    .. code:: fortran

        subroutine example_determinants


    subroutine that illustrates the main features available in determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`debug_det`
       * :c:func:`decode_exc`
       * :c:func:`do_single_excitation`
       * :c:func:`get_excitation`
       * :c:func:`get_excitation_degree`
       * :c:func:`i_h_j`
       * :c:func:`print_det`

 
.. c:function:: example_determinants_psi_det:


    File : :file:`determinants/example.irp.f`

    .. code:: fortran

        subroutine example_determinants_psi_det


    subroutine that illustrates the main features available in determinants using the psi_det/psi_coef

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`routine_example_psi_det`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`read_wf`

 
.. c:function:: fill_h_apply_buffer_no_selection:


    File : :file:`determinants/h_apply.irp.f`

    .. code:: fortran

        subroutine fill_H_apply_buffer_no_selection(n_selected,det_buffer,Nint,iproc)


    Fill the H_apply buffer with determiants for |CISD|

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`generate_all_alpha_beta_det_products`
       * :c:func:`make_s2_eigenfunction`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`
       * :c:func:`resize_h_apply_buffer`

 
.. c:function:: filter_connected:


    File : :file:`determinants/filter_connected.irp.f`

    .. code:: fortran

        subroutine filter_connected(key1,key2,Nint,sze,idx)


    Filters out the determinants that are not connected by H
    
    returns the array idx which contains the index of the
    
    determinants in the array key1 that interact
    
    via the H operator with key2.
    
    idx(0) is the number of determinants that interact with key1

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_uj_s2_ui`

 
.. c:function:: filter_connected_i_h_psi0:


    File : :file:`determinants/filter_connected.irp.f`

    .. code:: fortran

        subroutine filter_connected_i_H_psi0(key1,key2,Nint,sze,idx)


    Returns the array idx which contains the index of the
    
    determinants in the array key1 that interact
    
    via the H operator with key2.
    
    idx(0) is the number of determinants that interact with key1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_h_psi`
       * :c:func:`i_h_psi_minilist`
       * :c:func:`i_s2_psi_minilist`

 
.. c:function:: filter_not_connected:


    File : :file:`determinants/filter_connected.irp.f`

    .. code:: fortran

        subroutine filter_not_connected(key1,key2,Nint,sze,idx)


    Returns the array idx which contains the index of the
    
    determinants in the array key1 that DO NOT interact
    
    via the H operator with key2.
    
    idx(0) is the number of determinants that DO NOT interact with key1

 
.. c:function:: generate_all_alpha_beta_det_products:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine generate_all_alpha_beta_det_products


    Creates a wave function from all possible $\alpha \times \beta$ determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`create_wf_of_psi_bilinear_matrix`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`copy_h_apply_buffer_to_wf`
       * :c:func:`fill_h_apply_buffer_no_selection`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`c0_weight`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit`

 
.. c:function:: get_all_spin_doubles:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine get_all_spin_doubles(buffer, idx, spindet, Nint, size_buffer, doubles, n_doubles)


    
    Returns the indices of all the double excitations in the list of
    unique $\alpha$ determinants.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_doubles_1`
       * :c:func:`get_all_spin_doubles_2`
       * :c:func:`get_all_spin_doubles_3`
       * :c:func:`get_all_spin_doubles_4`
       * :c:func:`get_all_spin_doubles_n_int`

 
.. c:function:: get_all_spin_doubles_1:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine get_all_spin_doubles_1(buffer, idx, spindet, size_buffer, doubles, n_doubles)


    
    Returns the indices of all the double excitations in the list of
    unique $\alpha$ determinants.
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_doubles`

 
.. c:function:: get_all_spin_doubles_2:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_doubles_2(buffer, idx, spindet, size_buffer, doubles, n_doubles)


    
    Returns the indices of all the double excitations in the list of
    unique $lpha$ determinants.
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_doubles`

 
.. c:function:: get_all_spin_doubles_3:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_doubles_3(buffer, idx, spindet, size_buffer, doubles, n_doubles)


    
    Returns the indices of all the double excitations in the list of
    unique $lpha$ determinants.
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_doubles`

 
.. c:function:: get_all_spin_doubles_4:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_doubles_4(buffer, idx, spindet, size_buffer, doubles, n_doubles)


    
    Returns the indices of all the double excitations in the list of
    unique $lpha$ determinants.
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_doubles`

 
.. c:function:: get_all_spin_doubles_n_int:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_doubles_N_int(buffer, idx, spindet, size_buffer, doubles, n_doubles)


    
    Returns the indices of all the double excitations in the list of
    unique $lpha$ determinants.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_doubles`

 
.. c:function:: get_all_spin_singles:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine get_all_spin_singles(buffer, idx, spindet, Nint, size_buffer, singles, n_singles)


    
    Returns the indices of all the single excitations in the list of
    unique $\alpha$ determinants.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`singles_alpha_csc`
       * :c:data:`singles_alpha_csc_idx`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_1`
       * :c:func:`get_all_spin_singles_2`
       * :c:func:`get_all_spin_singles_3`
       * :c:func:`get_all_spin_singles_4`
       * :c:func:`get_all_spin_singles_n_int`

 
.. c:function:: get_all_spin_singles_1:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine get_all_spin_singles_1(buffer, idx, spindet, size_buffer, singles, n_singles)


    
    Returns the indices of all the single excitations in the list of
    unique $\alpha$ determinants.
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`h_s2_u_0_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`

 
.. c:function:: get_all_spin_singles_2:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_singles_2(buffer, idx, spindet, size_buffer, singles, n_singles)


    
    Returns the indices of all the single excitations in the list of
    unique $lpha$ determinants.
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`h_s2_u_0_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`

 
.. c:function:: get_all_spin_singles_3:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_singles_3(buffer, idx, spindet, size_buffer, singles, n_singles)


    
    Returns the indices of all the single excitations in the list of
    unique $lpha$ determinants.
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`h_s2_u_0_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`

 
.. c:function:: get_all_spin_singles_4:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_singles_4(buffer, idx, spindet, size_buffer, singles, n_singles)


    
    Returns the indices of all the single excitations in the list of
    unique $lpha$ determinants.
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`h_s2_u_0_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`

 
.. c:function:: get_all_spin_singles_and_doubles:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine get_all_spin_singles_and_doubles(buffer, idx, spindet, Nint, size_buffer, singles, doubles, n_singles, n_doubles)


    
    Returns the indices of all the single and double excitations in the list of
    unique $\alpha$ determinants.
    
    Warning: The buffer is transposed.
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_and_doubles_1`
       * :c:func:`get_all_spin_singles_and_doubles_2`
       * :c:func:`get_all_spin_singles_and_doubles_3`
       * :c:func:`get_all_spin_singles_and_doubles_4`
       * :c:func:`get_all_spin_singles_and_doubles_n_int`

 
.. c:function:: get_all_spin_singles_and_doubles_1:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine get_all_spin_singles_and_doubles_1(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)


    
    Returns the indices of all the single and double excitations in the list of
    unique $\alpha$ determinants.
    
    /!\ : The buffer is transposed !
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`h_s2_u_0_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`

 
.. c:function:: get_all_spin_singles_and_doubles_2:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_singles_and_doubles_2(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)


    
    Returns the indices of all the single and double excitations in the list of
    unique $lpha$ determinants.
    
    /!\ : The buffer is transposed !
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`h_s2_u_0_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`

 
.. c:function:: get_all_spin_singles_and_doubles_3:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_singles_and_doubles_3(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)


    
    Returns the indices of all the single and double excitations in the list of
    unique $lpha$ determinants.
    
    /!\ : The buffer is transposed !
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`h_s2_u_0_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`

 
.. c:function:: get_all_spin_singles_and_doubles_4:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_singles_and_doubles_4(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)


    
    Returns the indices of all the single and double excitations in the list of
    unique $lpha$ determinants.
    
    /!\ : The buffer is transposed !
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`h_s2_u_0_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`

 
.. c:function:: get_all_spin_singles_and_doubles_n_int:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_singles_and_doubles_N_int(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)


    
    Returns the indices of all the single and double excitations in the list of
    unique $lpha$ determinants.
    
    /!\ : The buffer is transposed !
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_and_doubles`
       * :c:func:`h_s2_u_0_nstates_openmp_work_n_int`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_n_int`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

 
.. c:function:: get_all_spin_singles_n_int:


    File : :file:`determinants/spindeterminants.irp.f_template_1316`

    .. code:: fortran

        subroutine get_all_spin_singles_N_int(buffer, idx, spindet, size_buffer, singles, n_singles)


    
    Returns the indices of all the single excitations in the list of
    unique $lpha$ determinants.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles`
       * :c:func:`h_s2_u_0_nstates_openmp_work_n_int`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_n_int`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`

 
.. c:function:: get_double_excitation:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_double_excitation(det1,det2,exc,phase,Nint)


    Returns the two excitation operators between two doubly excited determinants and the phase.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_h_mat_elem_fock`
       * :c:func:`get_excitation`
       * :c:func:`get_s2`
       * :c:func:`i_h_j`
       * :c:func:`i_h_j_s2`
       * :c:func:`i_h_j_two_e`
       * :c:func:`i_h_j_verbose`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_ab_dm_buffer`

 
.. c:function:: get_double_excitation_spin:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_double_excitation_spin(det1,det2,exc,phase,Nint)


    Returns the two excitation operators between two doubly excited spin-determinants
    and the phase.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_spin`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_2_rdm_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_double_to_all_states_bb_dm_buffer`

 
.. c:function:: get_excitation:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_excitation(det1,det2,exc,degree,phase,Nint)


    Returns the excitation operators between two determinants and the phase.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`example_determinants`
       * :c:func:`get_phase`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`

 
.. c:function:: get_excitation_degree:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_excitation_degree(key1,key2,degree,Nint)


    Returns the excitation degree between two determinants.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`connected_to_hf`
       * :c:func:`diag_h_mat_elem_fock`
       * :c:func:`example_determinants`
       * :c:func:`get_excitation`
       * :c:func:`get_s2`
       * :c:func:`i_h_j`
       * :c:func:`i_h_j_one_e`
       * :c:func:`i_h_j_s2`
       * :c:func:`i_h_j_two_e`
       * :c:func:`i_h_j_verbose`
       * :c:data:`max_degree_exc`
       * :c:data:`psi_non_cas`

 
.. c:function:: get_excitation_degree_spin:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_excitation_degree_spin(key1,key2,degree,Nint)


    Returns the excitation degree between two determinants.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_spin`
       * :c:data:`one_e_dm_mo_alpha`

 
.. c:function:: get_excitation_degree_vector:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_excitation_degree_vector(key1,key2,degree,Nint,sze,idx)


    Applies get_excitation_degree to an array of determinants.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`routine_example_psi_det`

 
.. c:function:: get_excitation_degree_vector_double_alpha_beta:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_excitation_degree_vector_double_alpha_beta(key1,key2,degree,Nint,sze,idx)


    Applies get_excitation_degree to an array of determinants and return only the
    single excitations and the connections through exchange integrals.

 
.. c:function:: get_excitation_degree_vector_single:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_excitation_degree_vector_single(key1,key2,degree,Nint,sze,idx)


    Applies get_excitation_degree to an array of determinants and returns only
    the single excitations.

 
.. c:function:: get_excitation_degree_vector_single_or_exchange:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_excitation_degree_vector_single_or_exchange(key1,key2,degree,Nint,sze,idx)


    Applies get_excitation_degree to an array of determinants and return only the
    single excitations and the connections through exchange integrals.

 
.. c:function:: get_excitation_degree_vector_single_or_exchange_verbose:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_excitation_degree_vector_single_or_exchange_verbose(key1,key2,degree,Nint,sze,idx)


    Applies get_excitation_degree to an array of determinants and return only the single
    excitations and the connections through exchange integrals.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`

 
.. c:function:: get_excitation_spin:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_excitation_spin(det1,det2,exc,degree,phase,Nint)


    Returns the excitation operators between two determinants and the phase.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation_spin`
       * :c:func:`get_excitation_degree_spin`
       * :c:func:`get_single_excitation_spin`

 
.. c:function:: get_index_in_psi_det_alpha_unique:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer function get_index_in_psi_det_alpha_unique(key,Nint)


    Returns the index of the determinant in the :c:data:`psi_det_alpha_unique` array

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_alpha_unique`

 
.. c:function:: get_index_in_psi_det_beta_unique:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer function get_index_in_psi_det_beta_unique(key,Nint)


    Returns the index of the determinant in the :c:data:`psi_det_beta_unique` array

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_beta_unique`

 
.. c:function:: get_index_in_psi_det_sorted_bit:


    File : :file:`determinants/connected_to_ref.irp.f`

    .. code:: fortran

        integer function get_index_in_psi_det_sorted_bit(key,Nint)


    Returns the index of the determinant in the ``psi_det_sorted_bit`` array

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_det_sorted_bit`

 
.. c:function:: get_occupation_from_dets:


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        subroutine get_occupation_from_dets(istate,occupation)


    Returns the average occupation of the MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: get_phase:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_phase(key1,key2,phase,Nint)


    Returns the phase between key1 and key2.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`build_singly_excited_wavefunction`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation`

 
.. c:function:: get_phasemask_bit:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_phasemask_bit(det1, pm, Nint)



 
.. c:function:: get_s2:


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        subroutine get_s2(key_i,key_j,Nint,s2)


    Returns $\langle S^2 \rangle - S_z^2 S_z$

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_uj_s2_ui`
       * :c:func:`h_s2_u_0_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_nstates_openmp_work_n_int`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_n_int`
       * :c:func:`i_s2_psi_minilist`
       * :c:data:`s2_matrix_all_dets`
       * :c:func:`s2_u_0_nstates`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`

 
.. c:function:: get_single_excitation:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_single_excitation(det1,det2,exc,phase,Nint)


    Returns the excitation operator between two singly excited determinants and the phase.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`connected_to_hf`
       * :c:func:`diag_h_mat_elem_fock`
       * :c:func:`get_excitation`
       * :c:func:`i_h_j`
       * :c:func:`i_h_j_one_e`
       * :c:func:`i_h_j_s2`
       * :c:func:`i_h_j_two_e`
       * :c:func:`i_h_j_verbose`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_2_rdm_bb_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_aa_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_ab_dm_buffer`
       * :c:func:`orb_range_off_diag_single_to_all_states_bb_dm_buffer`

 
.. c:function:: get_single_excitation_from_fock:


    File : :file:`determinants/single_excitations.irp.f`

    .. code:: fortran

        subroutine get_single_excitation_from_fock(det_1,det_2,h,p,spin,phase,hij)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`ref_closed_shell_bitmask`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_h_j`
       * :c:func:`i_h_j_s2`
       * :c:func:`i_h_j_single_spin`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: get_single_excitation_spin:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine get_single_excitation_spin(det1,det2,exc,phase,Nint)


    Returns the excitation operator between two singly excited determinants and the phase.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_spin`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_mono_spin_one_e`
       * :c:func:`i_h_j_single_spin`
       * :c:func:`i_wee_j_single`
       * :c:data:`one_e_dm_mo_alpha`

 
.. c:function:: get_uj_s2_ui:


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        subroutine get_uJ_s2_uI(psi_keys_tmp,psi_coefs_tmp,n,nmax_coefs,nmax_keys,s2,nstates)


    returns the matrix elements of S^2 "s2(i,j)" between the "nstates" states
    psi_coefs_tmp(:,i) and psi_coefs_tmp(:,j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`filter_connected`
       * :c:func:`get_s2`

 
.. c:function:: getmobiles:


    File : :file:`determinants/filter_connected.irp.f`

    .. code:: fortran

        subroutine getMobiles(key,key_mask, mobiles,Nint)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list`

 
.. c:function:: i_h_j:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine i_H_j(key_i,key_j,Nint,hij)


    Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`connected_to_hf`
       * :c:func:`example_determinants`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:func:`i_h_psi`
       * :c:func:`i_h_psi_minilist`
       * :c:func:`routine_example_psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`
       * :c:func:`get_single_excitation_from_fock`

 
.. c:function:: i_h_j_double_alpha_beta:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine i_H_j_double_alpha_beta(key_i,key_j,Nint,hij)


    Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants differing by
    an opposite-spin double excitation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_nstates_openmp_work_n_int`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_single_excitation_spin`

 
.. c:function:: i_h_j_double_spin:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine i_H_j_double_spin(key_i,key_j,Nint,hij)


    Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants differing by
    a same-spin double excitation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_nstates_openmp_work_n_int`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_double_excitation_spin`

 
.. c:function:: i_h_j_mono_spin_one_e:


    File : :file:`determinants/slater_rules_wee_mono.irp.f`

    .. code:: fortran

        subroutine i_H_j_mono_spin_one_e(key_i,key_j,Nint,spin,hij)


    Returns $\langle i|H|j \rangle$  where $i$ and $j$ are determinants differing by
    a single excitation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_one_e_integrals`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_single_excitation_spin`

 
.. c:function:: i_h_j_one_e:


    File : :file:`determinants/slater_rules_wee_mono.irp.f`

    .. code:: fortran

        subroutine i_H_j_one_e(key_i,key_j,Nint,hij)


    Returns $\langle i|H|j \rangle$  where $i$ and $j$ are determinants.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`

 
.. c:function:: i_h_j_s2:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine i_H_j_s2(key_i,key_j,Nint,hij,s2)


    Returns $\langle i|H|j \rangle$ and $\langle i|S^2|j \rangle$
    where $i$ and $j$ are determinants.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`
       * :c:func:`get_single_excitation_from_fock`

 
.. c:function:: i_h_j_single_spin:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine i_H_j_single_spin(key_i,key_j,Nint,spin,hij)


    Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants differing by
    a single excitation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`mo_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_nstates_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_single_excitation_from_fock`
       * :c:func:`get_single_excitation_spin`

 
.. c:function:: i_h_j_two_e:


    File : :file:`determinants/slater_rules_wee_mono.irp.f`

    .. code:: fortran

        subroutine i_H_j_two_e(key_i,key_j,Nint,hij)


    Returns $\langle i|H|j \rangle$  where $i$ and $j$ are determinants.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask_energy`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`
       * :c:func:`single_excitation_wee`

 
.. c:function:: i_h_j_verbose:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine i_H_j_verbose(key_i,key_j,Nint,hij,hmono,hdouble,phase)


    Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`get_double_excitation`
       * :c:func:`get_excitation_degree`
       * :c:func:`get_single_excitation`

 
.. c:function:: i_h_psi:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine i_H_psi(key,keys,coef,Nint,Ndet,Ndet_max,Nstate,i_H_psi_array)


    Computes $\langle i|H|Psi \rangle  = \sum_J c_J \langle i | H | J \rangle$.
    
    Uses filter_connected_i_H_psi0 to get all the $|J \rangle$ to which $|i \rangle$
    is connected.
    The i_H_psi_minilist is much faster but requires to build the
    minilists.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`filter_connected_i_h_psi0`
       * :c:func:`i_h_j`

 
.. c:function:: i_h_psi_minilist:


    File : :file:`determinants/slater_rules.irp.f`

    .. code:: fortran

        subroutine i_H_psi_minilist(key,keys,idx_key,N_minilist,coef,Nint,Ndet,Ndet_max,Nstate,i_H_psi_array)


    Computes $\langle i|H|\Psi \rangle = \sum_J c_J \langle i|H|J\rangle$.
    
    Uses filter_connected_i_H_psi0 to get all the $|J \rangle$ to which $|i \rangle$
    is connected. The $|J\rangle$ are searched in short pre-computed lists.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`filter_connected_i_h_psi0`
       * :c:func:`i_h_j`

 
.. c:function:: i_s2_psi_minilist:


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        subroutine i_S2_psi_minilist(key,keys,idx_key,N_minilist,coef,Nint,Ndet,Ndet_max,Nstate,i_S2_psi_array)


    Computes $\langle i|S^2|\Psi \rangle = \sum_J c_J \langle i|S^2|J \rangle$.
    
    Uses filter_connected_i_H_psi0 to get all the $|J\rangle$ to which $|i\rangle$
    is connected. The $|J\rangle$ are searched in short pre-computed lists.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`filter_connected_i_h_psi0`
       * :c:func:`get_s2`

 
.. c:function:: i_wee_j_single:


    File : :file:`determinants/slater_rules_wee_mono.irp.f`

    .. code:: fortran

        subroutine i_Wee_j_single(key_i,key_j,Nint,spin,hij)


    Returns $\langle i|H|j \rangle$  where $i$ and $j$ are determinants differing by a
    single excitation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`mo_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_single_excitation_spin`
       * :c:func:`single_excitation_wee`

 
.. c:function:: is_connected_to:


    File : :file:`determinants/connected_to_ref.irp.f`

    .. code:: fortran

        logical function is_connected_to(key,keys,Nint,Ndet)


    Returns |true| if determinant ``key`` is connected to ``keys``

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_connected_to_by_single:


    File : :file:`determinants/connected_to_ref.irp.f`

    .. code:: fortran

        logical function is_connected_to_by_single(key,keys,Nint,Ndet)


    Returns |true| is ``key`` is connected to ``keys`` by a single excitation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_in_wavefunction:


    File : :file:`determinants/connected_to_ref.irp.f`

    .. code:: fortran

        logical function is_in_wavefunction(key,Nint)


    |true| if the determinant ``det`` is in the wave function

 
.. c:function:: is_spin_flip_possible:


    File : :file:`determinants/create_excitations.irp.f`

    .. code:: fortran

        logical function is_spin_flip_possible(key_in,i_flip,ispin)


    returns |true| if the spin-flip of spin ispin in the orbital i_flip is possible
    on key_in

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: make_s2_eigenfunction:


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        subroutine make_s2_eigenfunction



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_occ_pattern`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_occ_pattern`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`copy_h_apply_buffer_to_wf`
       * :c:func:`fill_h_apply_buffer_no_selection`
       * :c:func:`occ_pattern_to_dets`
       * :c:func:`occ_pattern_to_dets_size`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_occ_pattern`
       * :c:data:`c0_weight`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_occ_pattern`

 
.. c:function:: occ_pattern_of_det:


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        subroutine occ_pattern_of_det(d,o,Nint)


    Transforms a determinant to an occupation pattern
    
    occ(:,1) : Single occupations
    
    occ(:,2) : Double occupations
    

 
.. c:function:: occ_pattern_search_key:


    File : :file:`determinants/connected_to_ref.irp.f`

    .. code:: fortran

        integer*8 function occ_pattern_search_key(det,Nint)


    Return an integer*8 corresponding to a determinant index for searching

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`

 
.. c:function:: occ_pattern_to_dets:


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        subroutine occ_pattern_to_dets(o,d,sze,n_alpha,Nint)


    Generate all possible determinants for a given occ_pattern
    
    Input :
       o   : occupation pattern : (doubly occupied, singly occupied)
       sze : Number of produced determinants, computed by `occ_pattern_to_dets_size`
       n_alpha : Number of $\alpha$ electrons
       Nint    : N_int
    
    Output:
       d : determinants
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`make_s2_eigenfunction`

 
.. c:function:: occ_pattern_to_dets_size:


    File : :file:`determinants/occ_pattern.irp.f`

    .. code:: fortran

        subroutine occ_pattern_to_dets_size(o,sze,n_alpha,Nint)


    Number of possible determinants for a given occ_pattern

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`make_s2_eigenfunction`

 
.. c:function:: read_dets:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine read_dets(det,Nint,Ndet)


    Reads the determinants from the |EZFIO| file

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_get_determinants_bit_kind`
       * :c:func:`ezfio_get_determinants_n_int`
       * :c:func:`ezfio_get_determinants_psi_det`

 
.. c:function:: read_spindeterminants:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine read_spindeterminants



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_get_spindeterminants_n_det`
       * :c:func:`ezfio_get_spindeterminants_n_det_alpha`
       * :c:func:`ezfio_get_spindeterminants_n_det_beta`
       * :c:func:`ezfio_get_spindeterminants_n_states`
       * :c:func:`ezfio_get_spindeterminants_psi_coef_matrix_columns`
       * :c:func:`ezfio_get_spindeterminants_psi_coef_matrix_rows`
       * :c:func:`ezfio_get_spindeterminants_psi_coef_matrix_values`
       * :c:func:`ezfio_get_spindeterminants_psi_det_alpha`
       * :c:func:`ezfio_get_spindeterminants_psi_det_beta`
       * :c:func:`wf_of_psi_bilinear_matrix`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

 
.. c:function:: remove_duplicates_in_psi_det:


    File : :file:`determinants/h_apply.irp.f`

    .. code:: fortran

        subroutine remove_duplicates_in_psi_det(found_duplicates)


    Removes duplicate determinants in the wave function.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`c0_weight`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det_sorted_bit`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`copy_h_apply_buffer_to_wf`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`c0_weight`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

 
.. c:function:: resize_h_apply_buffer:


    File : :file:`determinants/h_apply.irp.f`

    .. code:: fortran

        subroutine resize_H_apply_buffer(new_size,iproc)


    Resizes the H_apply buffer of proc iproc. The buffer lock should
    be set before calling this function.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nproc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`fill_h_apply_buffer_no_selection`

 
.. c:function:: routine_example_psi_det:


    File : :file:`determinants/example.irp.f`

    .. code:: fortran

        subroutine routine_example_psi_det


    subroutine that illustrates the main features available in determinants using many determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`example_determinants_psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`
       * :c:func:`get_excitation_degree_vector`
       * :c:func:`i_h_j`

 
.. c:function:: s2_u_0:


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        subroutine S2_u_0(v_0,u_0,n,keys_tmp,Nint)


    Computes v_0 = S^2|u_0>
    
    n : number of determinants
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`s2_u_0_nstates`

 
.. c:function:: s2_u_0_nstates:


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        subroutine S2_u_0_nstates(v_0,u_0,n,keys_tmp,Nint,N_st,sze_8)


    Computes v_0  = S^2|u_0>
    
    n : number of determinants
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`ref_bitmask_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`s2_u_0`
       * :c:func:`u_0_s2_u_0`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_s2`
       * :c:func:`sort_dets_ab_v`
       * :c:func:`sort_dets_ba_v`

 
.. c:function:: save_natural_mos:


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        subroutine save_natural_mos


    Save natural orbitals, obtained by diagonalization of the one-body density matrix in
    the |MO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`nullify_small_elements`
       * :c:func:`orthonormalize_mos`
       * :c:func:`save_mos`
       * :c:func:`set_natural_mos`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_coef`
       * :c:data:`mo_occ`

 
.. c:function:: save_ref_determinant:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine save_ref_determinant



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`ref_bitmask`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_wavefunction_general`

 
.. c:function:: save_wavefunction:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine save_wavefunction


    Save the wave function into the |EZFIO| file

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_wavefunction_general`

 
.. c:function:: save_wavefunction_general:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine save_wavefunction_general(ndet,nstates,psidet,dim_psicoef,psicoef)


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

       * :c:func:`save_ref_determinant`
       * :c:func:`save_wavefunction`
       * :c:func:`save_wavefunction_truncated`
       * :c:func:`save_wavefunction_unsorted`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_determinants_bit_kind`
       * :c:func:`ezfio_set_determinants_mo_label`
       * :c:func:`ezfio_set_determinants_n_det`
       * :c:func:`ezfio_set_determinants_n_det_qp_edit`
       * :c:func:`ezfio_set_determinants_n_int`
       * :c:func:`ezfio_set_determinants_n_states`
       * :c:func:`ezfio_set_determinants_psi_coef`
       * :c:func:`ezfio_set_determinants_psi_coef_qp_edit`
       * :c:func:`ezfio_set_determinants_psi_det`
       * :c:func:`ezfio_set_determinants_psi_det_qp_edit`
       * :c:func:`normalize`
       * :c:func:`write_int`

 
.. c:function:: save_wavefunction_specified:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine save_wavefunction_specified(ndet,nstates,psidet,psicoef,ndetsave,index_det_save)


    Save the wave function into the |EZFIO| file

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_label`
       * :c:data:`mpi_master`
       * :c:data:`n_det_qp_edit`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_determinants_bit_kind`
       * :c:func:`ezfio_set_determinants_mo_label`
       * :c:func:`ezfio_set_determinants_n_det`
       * :c:func:`ezfio_set_determinants_n_det_qp_edit`
       * :c:func:`ezfio_set_determinants_n_int`
       * :c:func:`ezfio_set_determinants_n_states`
       * :c:func:`ezfio_set_determinants_psi_coef`
       * :c:func:`ezfio_set_determinants_psi_det`
       * :c:func:`ezfio_set_determinants_psi_det_qp_edit`
       * :c:func:`write_int`

 
.. c:function:: save_wavefunction_truncated:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine save_wavefunction_truncated(thr)


    Save the wave function into the |EZFIO| file

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`nullify_small_elements`
       * :c:func:`save_wavefunction_general`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`psi_coef`

 
.. c:function:: save_wavefunction_unsorted:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine save_wavefunction_unsorted


    Save the wave function into the |EZFIO| file

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_wavefunction_general`

 
.. c:function:: set_natural_mos:


    File : :file:`determinants/density_matrix.irp.f`

    .. code:: fortran

        subroutine set_natural_mos


    Set natural orbitals, obtained by diagonalization of the one-body density matrix
    in the |MO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core_inact_act`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`mo_occ`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`one_e_dm_mo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`save_natural_mos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`mo_as_svd_vectors_of_mo_matrix_eig`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_occ`

 
.. c:function:: single_excitation_wee:


    File : :file:`determinants/single_excitation_two_e.irp.f`

    .. code:: fortran

        subroutine single_excitation_wee(det_1,det_2,h,p,spin,phase,hij)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`n_int`
       * :c:data:`ref_closed_shell_bitmask`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_h_j_two_e`
       * :c:func:`i_wee_j_single`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: sort_dets_ab:


    File : :file:`determinants/sort_dets_ab.irp.f`

    .. code:: fortran

        subroutine sort_dets_ab(key, idx, shortcut, N_key, Nint)


    Deprecated routine

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`tamiser`

 
.. c:function:: sort_dets_ab_v:


    File : :file:`determinants/sort_dets_ab.irp.f`

    .. code:: fortran

        subroutine sort_dets_ab_v(key_in, key_out, idx, shortcut, version, N_key, Nint)


    Deprecated routine

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`s2_u_0_nstates`
       * :c:func:`sort_dets_ba_v`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`tamiser`

 
.. c:function:: sort_dets_ba_v:


    File : :file:`determinants/sort_dets_ab.irp.f`

    .. code:: fortran

        subroutine sort_dets_ba_v(key_in, key_out, idx, shortcut, version, N_key, Nint)


    Deprecated routine

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`s2_u_0_nstates`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`sort_dets_ab_v`

 
.. c:function:: sort_dets_by_det_search_key:


    File : :file:`determinants/determinants.irp.f`

    .. code:: fortran

        subroutine sort_dets_by_det_search_key(Ndet, det_in, coef_in, sze, det_out, coef_out, N_st)


    Determinants are sorted according to their :c:func:`det_search_key`.
    Useful to accelerate the search of a random determinant in the wave
    function.
    
    /!\ The first dimension of coef_out and coef_in need to be psi_det_size
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_cas_sorted_bit`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_non_cas_sorted_bit`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i8sort`

 
.. c:function:: spin_det_search_key:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        integer*8 function spin_det_search_key(det,Nint)


    Returns an integer(8) corresponding to a determinant index for searching

 
.. c:function:: tamiser:


    File : :file:`determinants/sort_dets_ab.irp.f`

    .. code:: fortran

        subroutine tamiser(key, idx, no, n, Nint, N_key)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`sort_dets_ab`
       * :c:func:`sort_dets_ab_v`

 
.. c:function:: u_0_s2_u_0:


    File : :file:`determinants/s2.irp.f`

    .. code:: fortran

        subroutine u_0_S2_u_0(e_0,u_0,n,keys_tmp,Nint,N_st,sze_8)


    Computes e_0 = <u_0|S2|u_0>/<u_0|u_0>
    
    n : number of determinants
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`s2_values`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`s2_u_0_nstates`

 
.. c:function:: wf_of_psi_bilinear_matrix:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine wf_of_psi_bilinear_matrix(truncate)


    Generate a wave function containing all possible products
    of $\alpha$ and $\beta$ determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_sorted`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`read_spindeterminants`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

 
.. c:function:: write_spindeterminants:


    File : :file:`determinants/spindeterminants.irp.f`

    .. code:: fortran

        subroutine write_spindeterminants



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_spindeterminants_bit_kind`
       * :c:func:`ezfio_set_spindeterminants_n_det`
       * :c:func:`ezfio_set_spindeterminants_n_det_alpha`
       * :c:func:`ezfio_set_spindeterminants_n_det_beta`
       * :c:func:`ezfio_set_spindeterminants_n_int`
       * :c:func:`ezfio_set_spindeterminants_n_states`
       * :c:func:`ezfio_set_spindeterminants_psi_coef_matrix_columns`
       * :c:func:`ezfio_set_spindeterminants_psi_coef_matrix_rows`
       * :c:func:`ezfio_set_spindeterminants_psi_coef_matrix_values`
       * :c:func:`ezfio_set_spindeterminants_psi_det_alpha`
       * :c:func:`ezfio_set_spindeterminants_psi_det_beta`

 
.. c:function:: zmq_get_n_det:


    File : :file:`determinants/zmq.irp.f_template_379`

    .. code:: fortran

        integer function zmq_get_N_det(zmq_to_qp_run_socket, worker_id)


    Get N_det from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_get_n_det_alpha_unique:


    File : :file:`determinants/zmq.irp.f_template_379`

    .. code:: fortran

        integer function zmq_get_N_det_alpha_unique(zmq_to_qp_run_socket, worker_id)


    Get N_det_alpha_unique from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_get_n_det_beta_unique:


    File : :file:`determinants/zmq.irp.f_template_379`

    .. code:: fortran

        integer function zmq_get_N_det_beta_unique(zmq_to_qp_run_socket, worker_id)


    Get N_det_beta_unique from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_get_n_states:


    File : :file:`determinants/zmq.irp.f_template_379`

    .. code:: fortran

        integer function zmq_get_N_states(zmq_to_qp_run_socket, worker_id)


    Get N_states from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_states`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_get_psi:


    File : :file:`determinants/zmq.irp.f`

    .. code:: fortran

        integer function zmq_get_psi(zmq_to_qp_run_socket, worker_id)


    Get the wave function from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

 
.. c:function:: zmq_get_psi_bilinear:


    File : :file:`determinants/zmq.irp.f`

    .. code:: fortran

        integer function zmq_get_psi_bilinear(zmq_to_qp_run_socket, worker_id)


    Get the wave function from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_size`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_size`

 
.. c:function:: zmq_get_psi_bilinear_matrix_columns:


    File : :file:`determinants/zmq.irp.f_template_500`

    .. code:: fortran

        integer*8 function zmq_get_psi_bilinear_matrix_columns(zmq_to_qp_run_socket,worker_id)


    Get psi_bilinear_matrix_columns on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_values`

 
.. c:function:: zmq_get_psi_bilinear_matrix_order:


    File : :file:`determinants/zmq.irp.f_template_500`

    .. code:: fortran

        integer*8 function zmq_get_psi_bilinear_matrix_order(zmq_to_qp_run_socket,worker_id)


    Get psi_bilinear_matrix_order on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_values`

 
.. c:function:: zmq_get_psi_bilinear_matrix_rows:


    File : :file:`determinants/zmq.irp.f_template_500`

    .. code:: fortran

        integer*8 function zmq_get_psi_bilinear_matrix_rows(zmq_to_qp_run_socket,worker_id)


    Get psi_bilinear_matrix_rows on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_values`

 
.. c:function:: zmq_get_psi_bilinear_matrix_values:


    File : :file:`determinants/zmq.irp.f_template_564`

    .. code:: fortran

        integer*8 function zmq_get_psi_bilinear_matrix_values(zmq_to_qp_run_socket,worker_id)


    get psi_bilinear_matrix_values on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_values`

 
.. c:function:: zmq_get_psi_coef:


    File : :file:`determinants/zmq.irp.f_template_564`

    .. code:: fortran

        integer*8 function zmq_get_psi_coef(zmq_to_qp_run_socket,worker_id)


    get psi_coef on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_coef`

 
.. c:function:: zmq_get_psi_det:


    File : :file:`determinants/zmq.irp.f_template_440`

    .. code:: fortran

        integer*8 function zmq_get_psi_det(zmq_to_qp_run_socket,worker_id)


    Get psi_det on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det`

 
.. c:function:: zmq_get_psi_det_alpha_unique:


    File : :file:`determinants/zmq.irp.f_template_440`

    .. code:: fortran

        integer*8 function zmq_get_psi_det_alpha_unique(zmq_to_qp_run_socket,worker_id)


    Get psi_det_alpha_unique on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_alpha_unique`

 
.. c:function:: zmq_get_psi_det_beta_unique:


    File : :file:`determinants/zmq.irp.f_template_440`

    .. code:: fortran

        integer*8 function zmq_get_psi_det_beta_unique(zmq_to_qp_run_socket,worker_id)


    Get psi_det_beta_unique on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_beta_unique`

 
.. c:function:: zmq_get_psi_det_size:


    File : :file:`determinants/zmq.irp.f_template_379`

    .. code:: fortran

        integer function zmq_get_psi_det_size(zmq_to_qp_run_socket, worker_id)


    Get psi_det_size from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`psi_det_size`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_get_psi_notouch:


    File : :file:`determinants/zmq.irp.f`

    .. code:: fortran

        integer function zmq_get_psi_notouch(zmq_to_qp_run_socket, worker_id)


    Get the wave function from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

 
.. c:function:: zmq_put_n_det:


    File : :file:`determinants/zmq.irp.f_template_379`

    .. code:: fortran

        integer function zmq_put_N_det(zmq_to_qp_run_socket,worker_id)


    Put N_det on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_n_det_alpha_unique:


    File : :file:`determinants/zmq.irp.f_template_379`

    .. code:: fortran

        integer function zmq_put_N_det_alpha_unique(zmq_to_qp_run_socket,worker_id)


    Put N_det_alpha_unique on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_alpha_unique`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_n_det_beta_unique:


    File : :file:`determinants/zmq.irp.f_template_379`

    .. code:: fortran

        integer function zmq_put_N_det_beta_unique(zmq_to_qp_run_socket,worker_id)


    Put N_det_beta_unique on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_beta_unique`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_n_states:


    File : :file:`determinants/zmq.irp.f_template_379`

    .. code:: fortran

        integer function zmq_put_N_states(zmq_to_qp_run_socket,worker_id)


    Put N_states on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_psi:


    File : :file:`determinants/zmq.irp.f`

    .. code:: fortran

        integer function zmq_put_psi(zmq_to_qp_run_socket,worker_id)


    Put the wave function on the qp_run scheduler

 
.. c:function:: zmq_put_psi_bilinear:


    File : :file:`determinants/zmq.irp.f`

    .. code:: fortran

        integer function zmq_put_psi_bilinear(zmq_to_qp_run_socket,worker_id)


    Put the wave function on the qp_run scheduler

 
.. c:function:: zmq_put_psi_bilinear_matrix_columns:


    File : :file:`determinants/zmq.irp.f_template_500`

    .. code:: fortran

        integer*8 function zmq_put_psi_bilinear_matrix_columns(zmq_to_qp_run_socket,worker_id)


    Put psi_bilinear_matrix_columns on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_values`

 
.. c:function:: zmq_put_psi_bilinear_matrix_order:


    File : :file:`determinants/zmq.irp.f_template_500`

    .. code:: fortran

        integer*8 function zmq_put_psi_bilinear_matrix_order(zmq_to_qp_run_socket,worker_id)


    Put psi_bilinear_matrix_order on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_values`

 
.. c:function:: zmq_put_psi_bilinear_matrix_rows:


    File : :file:`determinants/zmq.irp.f_template_500`

    .. code:: fortran

        integer*8 function zmq_put_psi_bilinear_matrix_rows(zmq_to_qp_run_socket,worker_id)


    Put psi_bilinear_matrix_rows on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_values`

 
.. c:function:: zmq_put_psi_bilinear_matrix_values:


    File : :file:`determinants/zmq.irp.f_template_564`

    .. code:: fortran

        integer*8 function zmq_put_psi_bilinear_matrix_values(zmq_to_qp_run_socket,worker_id)


    Put psi_bilinear_matrix_values on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_values`

 
.. c:function:: zmq_put_psi_coef:


    File : :file:`determinants/zmq.irp.f_template_564`

    .. code:: fortran

        integer*8 function zmq_put_psi_coef(zmq_to_qp_run_socket,worker_id)


    Put psi_coef on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_coef`

 
.. c:function:: zmq_put_psi_det:


    File : :file:`determinants/zmq.irp.f_template_440`

    .. code:: fortran

        integer*8 function zmq_put_psi_det(zmq_to_qp_run_socket,worker_id)


    Put psi_det on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det`

 
.. c:function:: zmq_put_psi_det_alpha_unique:


    File : :file:`determinants/zmq.irp.f_template_440`

    .. code:: fortran

        integer*8 function zmq_put_psi_det_alpha_unique(zmq_to_qp_run_socket,worker_id)


    Put psi_det_alpha_unique on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_alpha_unique`

 
.. c:function:: zmq_put_psi_det_beta_unique:


    File : :file:`determinants/zmq.irp.f_template_440`

    .. code:: fortran

        integer*8 function zmq_put_psi_det_beta_unique(zmq_to_qp_run_socket,worker_id)


    Put psi_det_beta_unique on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_beta_unique`

 
.. c:function:: zmq_put_psi_det_size:


    File : :file:`determinants/zmq.irp.f_template_379`

    .. code:: fortran

        integer function zmq_put_psi_det_size(zmq_to_qp_run_socket,worker_id)


    Put psi_det_size on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_size`
       * :c:data:`zmq_state`

