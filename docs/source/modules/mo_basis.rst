.. _module_mo_basis: 
 
.. program:: mo_basis 
 
.. default-role:: option 
 
========
mo_basis
========

Molecular orbitals are expressed as

.. math::

  \phi_k({\bf r}) = \sum_i C_{ik} \chi_k({\bf r})


where :math:`\chi_k` are *normalized* atomic basis functions.

The current set of |MOs| has a label `mo_label`.
When the orbitals are modified, the label should also be updated to keep
everything consistent.

When saving the |MOs|, the :file:`mo_basis` directory of the |EZFIO| database
is copied in the :file:`save` directory, named by the current `mo_label`. All
this is done with the script named :file:`save_current_mos.sh` in the
:file:`$QP_ROOT/scripts` directory.



 
 
 
EZFIO parameters 
---------------- 
 
.. option:: mo_num
 
    Total number of |MOs|
 
 
.. option:: mo_coef
 
    Coefficient of the i-th |AO| on the j-th |MO|
 
 
.. option:: mo_coef_aux
 
    AUX Coefficient of the i-th |AO| on the j-th |MO|
 
 
.. option:: mo_coef_imag
 
    Imaginary part of the MO coefficient of the i-th |AO| on the j-th |MO|
 
 
.. option:: mo_label
 
    Label characterizing the MOS (Local, Canonical, Natural, *etc*)
 
 
.. option:: mo_occ
 
    |MO| occupation numbers
 
 
.. option:: mo_symmetry
 
    MOs with the same integer belong to the same irrep.
 
 
.. option:: mo_class
 
    [ Core | Inactive | Active | Virtual | Deleted ], as defined by :ref:`qp_set_mo_class`
 
 
.. option:: ao_md5
 
    MD5 checksum characterizing the |AO| basis set.
 
 
 
Providers 
--------- 
 
.. c:var:: mo_class


    File : :file:`mo_basis/mo_class.irp.f`

    .. code:: fortran

        character*(32), allocatable	:: mo_class	(mo_num)


    [ Core | Inactive | Active | Virtual | Deleted ], as defined by :ref:`qp_set_mo_class`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_all_but_del_orb`
       * :c:data:`list_core`
       * :c:data:`list_del`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_act_orb`
       * :c:data:`n_all_but_del_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_del_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`

 
.. c:var:: mo_coef


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_coef	(ao_num,mo_num)


    Molecular orbital coefficients on |AO| basis set
    
    mo_coef(i,j) = coefficient of the i-th |AO| on the jth |MO|
    
    mo_label : Label characterizing the |MOs| (local, canonical, natural, etc)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`attachment_orbitals`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`d0tu_alpha_ao`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo_alpha`
       * :c:data:`fock_matrix_mo_beta`
       * :c:data:`fps_spf_matrix_mo`
       * :c:data:`mcscf_fock_alpha_mo`
       * :c:data:`mo_coef_in_ao_ortho_basis`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_deriv_1_x`
       * :c:data:`mo_dipole_x`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_integrals_n_e_per_atom`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_overlap`
       * :c:data:`mo_pseudo_integrals`
       * :c:data:`mo_pseudo_integrals_local`
       * :c:data:`mo_pseudo_integrals_non_local`
       * :c:data:`mo_spread_x`
       * :c:data:`mo_two_e_int_erf_jj_from_ao`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`natorbsfci`
       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`one_e_dm_alpha_ao_for_dft_no_core`
       * :c:data:`one_e_dm_ao_alpha`
       * :c:data:`one_e_dm_ao_alpha_nstates`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`one_e_spin_density_ao`
       * :c:data:`psi_det`
       * :c:data:`s_mo_coef`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

 
.. c:var:: mo_coef_aux


    File : :file:`mo_basis/mos_aux.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_coef_aux	(ao_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`


 
.. c:var:: mo_coef_imag


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_coef_imag	(ao_num,mo_num)


    Molecular orbital coefficients on |AO| basis set
    
    mo_coef_imag(i,j) = coefficient of the i-th |AO| on the jth |MO|
    
    mo_label : Label characterizing the |MOs| (local, canonical, natural, etc)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`


 
.. c:var:: mo_coef_in_ao_ortho_basis


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_coef_in_ao_ortho_basis	(ao_num,mo_num)


    |MO| coefficients in orthogonalized |AO| basis
    
    :math:`C^{-1}.C_{mo}` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_ortho_canonical_coef_inv`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: mo_coef_transp


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_coef_transp	(mo_num,ao_num)


    |MO| coefficients on |AO| basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:data:`mo_two_e_int_erf_jj_from_ao`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`

 
.. c:var:: mo_label


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        character*(64)	:: mo_label	


    |MO| coefficients on |AO| basis set
    
    mo_coef(i,j) = coefficient of the i-th |AO| on the j-th |MO|
    
    mo_label : Label characterizing the |MOs| (local, canonical, natural, etc)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

 
.. c:var:: mo_num


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        integer	:: mo_num	


    Number of MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_one_e_integrals_from_mo`
       * :c:data:`ao_ortho_canonical_nucl_elec_integrals`
       * :c:data:`ao_ortho_lowdin_nucl_elec_integrals`
       * :c:data:`attachment_numbers_sorted`
       * :c:data:`attachment_orbitals`
       * :c:data:`banned_excitation`
       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pqxx_no_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`bielec_pxxq_no_array`
       * :c:data:`bielecci`
       * :c:data:`bielecci_no`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`d0tu_alpha_ao`
       * :c:data:`data_one_e_dm_alpha_mo`
       * :c:data:`data_one_e_dm_beta_mo`
       * :c:data:`difference_dm`
       * :c:data:`difference_dm_eigvect`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`fock_matrix_ao`
       * :c:data:`fock_matrix_mo`
       * :c:data:`fock_matrix_mo_alpha`
       * :c:data:`fock_matrix_mo_beta`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`fps_spf_matrix_mo`
       * :c:data:`full_ijkl_bitmask`
       * :c:data:`h_core_ri`
       * :c:data:`int_erf_3_index`
       * :c:data:`list_act`
       * :c:data:`list_all_but_del_orb`
       * :c:data:`list_core`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_del`
       * :c:data:`list_inact`
       * :c:data:`list_inact_act`
       * :c:data:`list_virt`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mcscf_fock_alpha_mo`
       * :c:data:`mcscf_fock_mo`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_aux`
       * :c:data:`mo_coef_begin_iteration`
       * :c:data:`mo_coef_imag`
       * :c:data:`mo_coef_in_ao_ortho_basis`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_deriv_1_x`
       * :c:data:`mo_dipole_x`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_erf_cache_min`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_integrals_n_e_per_atom`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_occ`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_overlap`
       * :c:data:`mo_pseudo_integrals`
       * :c:data:`mo_pseudo_integrals_local`
       * :c:data:`mo_pseudo_integrals_non_local`
       * :c:data:`mo_spread_centered_x`
       * :c:data:`mo_spread_x`
       * :c:data:`mo_two_e_int_erf_jj`
       * :c:data:`mo_two_e_int_erf_jj_from_ao`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`multi_s_deriv_1`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_act_orb`
       * :c:data:`n_all_but_del_orb`
       * :c:data:`n_attachment`
       * :c:data:`n_core_orb`
       * :c:data:`n_del_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_int`
       * :c:data:`n_virt_orb`
       * :c:data:`natorbsci_mos`
       * :c:data:`natorbsfci`
       * :c:data:`neworbs`
       * :c:data:`occnum`
       * :c:data:`one_body_dm_mo_alpha_one_det`
       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`one_e_dm_alpha_ao_for_dft_no_core`
       * :c:data:`one_e_dm_ao_alpha`
       * :c:data:`one_e_dm_ao_alpha_nstates`
       * :c:data:`one_e_dm_average_alpha_mo_for_dft`
       * :c:data:`one_e_dm_average_beta_mo_for_dft`
       * :c:data:`one_e_dm_average_mo_for_dft`
       * :c:data:`one_e_dm_dagger_mo_spin_index`
       * :c:data:`one_e_dm_mo`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`one_e_dm_mo_alpha_average`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_alpha_for_dft_no_core`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft_no_core`
       * :c:data:`one_e_dm_mo_diff`
       * :c:data:`one_e_dm_mo_for_dft`
       * :c:data:`one_e_dm_mo_spin_index`
       * :c:data:`one_e_spin_density_ao`
       * :c:data:`one_e_spin_density_mo`
       * :c:data:`one_e_tr_dm_mo`
       * :c:data:`one_e_tr_dm_mo_alpha`
       * :c:data:`one_ints_no`
       * :c:data:`ormas_n_orb`
       * :c:data:`psi_energy_h_core`
       * :c:data:`s_mo_coef`
       * :c:data:`super_ci_dm`
       * :c:data:`superci_natorb`
       * :c:data:`switch_mo_coef`
       * :c:data:`two_e_dm_mo`
       * :c:data:`umat`
       * :c:data:`v_ne_psi_energy`
       * :c:data:`z_dipole_moment`

 
.. c:var:: mo_occ


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_occ	(mo_num)


    |MO| occupation numbers

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: ao_ortho_cano_to_ao:


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        subroutine ao_ortho_cano_to_ao(A_ao,LDA_ao,A,LDA)


    Transform A from the |AO| basis to the orthogonal |AO| basis
    
    $C^{-1}.A_{ao}.C^{\dagger-1}$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_ortho_canonical_coef_inv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: ao_to_mo:


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        subroutine ao_to_mo(A_ao,LDA_ao,A_mo,LDA_mo)


    Transform A from the |AO| basis to the |MO| basis
    
    $C^\dagger.A_{ao}.C$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_mo_alpha`
       * :c:data:`fock_matrix_mo_beta`
       * :c:data:`fps_spf_matrix_mo`
       * :c:data:`mcscf_fock_alpha_mo`
       * :c:data:`mo_deriv_1_x`
       * :c:data:`mo_dipole_x`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_integrals_n_e_per_atom`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_pseudo_integrals`
       * :c:data:`mo_pseudo_integrals_local`
       * :c:data:`mo_pseudo_integrals_non_local`
       * :c:data:`mo_spread_x`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`restore_symmetry`

 
.. c:function:: give_all_mos_and_grad_and_lapl_at_r:


    File : :file:`mo_basis/mos_in_r.irp.f`

    .. code:: fortran

        subroutine give_all_mos_and_grad_and_lapl_at_r(r,mos_array,mos_grad_array,mos_lapl_array)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_aos_and_grad_and_lapl_at_r`

 
.. c:function:: give_all_mos_and_grad_at_r:


    File : :file:`mo_basis/mos_in_r.irp.f`

    .. code:: fortran

        subroutine give_all_mos_and_grad_at_r(r,mos_array,mos_grad_array)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_aos_and_grad_at_r`

 
.. c:function:: give_all_mos_at_r:


    File : :file:`mo_basis/mos_in_r.irp.f`

    .. code:: fortran

        subroutine give_all_mos_at_r(r,mos_array)


    mos_array(i) = ith MO function evaluated at "r"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`give_all_aos_at_r`

 
.. c:function:: mix_mo_jk:


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        subroutine mix_mo_jk(j,k)


    Rotates the j-th |MO| with the k-th |MO| to give two new |MOs| that are
    
    * $+ = \frac{1}{\sqrt{2}} ( | j\rangle +  | k\rangle)$
    
    * $- = \frac{1}{\sqrt{2}} ( | j\rangle -  | k\rangle)$
    
    by convention, the '+' |MO| is in the lowest  index (min(j,k))
    by convention, the '-' |MO| is in the highest index (max(j,k))

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`

 
.. c:function:: mo_as_eigvectors_of_mo_matrix:


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine mo_as_eigvectors_of_mo_matrix(matrix,n,m,label,sign,output)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`create_guess`
       * :c:func:`damping_scf`
       * :c:func:`hcore_guess`
       * :c:func:`roothaan_hall_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`lapack_diag`
       * :c:func:`write_time`

 
.. c:function:: mo_as_svd_vectors_of_mo_matrix:


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine mo_as_svd_vectors_of_mo_matrix(matrix,lda,m,n,label)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`
       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`svd`
       * :c:func:`write_time`

 
.. c:function:: mo_as_svd_vectors_of_mo_matrix_eig:


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine mo_as_svd_vectors_of_mo_matrix_eig(matrix,lda,m,n,eig,label)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`set_natorb_no_ov_rot`
       * :c:func:`set_natural_mos`
       * :c:func:`set_natural_mos_canon_label`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`svd`
       * :c:func:`write_time`

 
.. c:function:: mo_coef_new_as_svd_vectors_of_mo_matrix_eig:


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine mo_coef_new_as_svd_vectors_of_mo_matrix_eig(matrix,lda,m,n,mo_coef_before,eig,mo_coef_new)


    You enter with matrix in the MO basis defined with the mo_coef_before.
    
    You SVD the matrix and set the eigenvectors as mo_coef_new ordered by increasing singular values

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`svd`
       * :c:func:`write_time`

 
.. c:function:: save_mos:


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine save_mos



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_md5`
       * :c:data:`ao_num`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`
       * :c:data:`mo_num`
       * :c:data:`mo_occ`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`damping_scf`
       * :c:func:`hcore_guess`
       * :c:func:`huckel_guess`
       * :c:func:`roothaan_hall_scf`
       * :c:func:`run_orb_opt_trust_v2`
       * :c:func:`save_natural_mos`
       * :c:func:`save_natural_mos_canon_label`
       * :c:func:`save_natural_mos_no_ov_rot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_mo_basis_ao_md5`
       * :c:func:`ezfio_set_mo_basis_mo_class`
       * :c:func:`ezfio_set_mo_basis_mo_coef`
       * :c:func:`ezfio_set_mo_basis_mo_label`
       * :c:func:`ezfio_set_mo_basis_mo_num`
       * :c:func:`ezfio_set_mo_basis_mo_occ`

 
.. c:function:: save_mos_no_occ:


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine save_mos_no_occ



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_mo_basis_mo_coef`

 
.. c:function:: save_mos_truncated:


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine save_mos_truncated(n)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_md5`
       * :c:data:`ao_num`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`
       * :c:data:`mo_occ`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_mo_basis_ao_md5`
       * :c:func:`ezfio_set_mo_basis_mo_class`
       * :c:func:`ezfio_set_mo_basis_mo_coef`
       * :c:func:`ezfio_set_mo_basis_mo_label`
       * :c:func:`ezfio_set_mo_basis_mo_num`
       * :c:func:`ezfio_set_mo_basis_mo_occ`

