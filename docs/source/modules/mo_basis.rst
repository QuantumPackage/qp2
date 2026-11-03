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
 
.. c:var:: char mo_class


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

 
.. c:var:: double mo_coef


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_coef	(ao_num,mo_num)


    Molecular orbital coefficients on |AO| basis set
    
    mo_coef(i,j) = coefficient of the i-th |AO| on the jth |MO|

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
       * :c:data:`mo_coef_in_ao_ortho_basis`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_deriv_1_x`
       * :c:data:`mo_dipole_x`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_integrals_n_e_per_atom`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_overlap`
       * :c:data:`mo_pseudo_integrals`
       * :c:data:`mo_pseudo_integrals_local`
       * :c:data:`mo_pseudo_integrals_non_local`
       * :c:data:`mo_sphe_coef`
       * :c:data:`mo_spread_x`
       * :c:data:`mo_two_e_int_erf_jj_from_ao`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`one_e_dm_alpha_ao_for_dft_no_core`
       * :c:data:`one_e_dm_ao_alpha`
       * :c:data:`one_e_dm_ao_alpha_nstates`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`one_e_spin_density_ao`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_xg`
       * :c:data:`psi_det`
       * :c:data:`s_mo_coef`
       * :c:data:`s_mo_sphe_coef`

 
.. c:var:: double mo_coef_aux


    File : :file:`mo_basis/mos_aux.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_coef_aux	(ao_num,mo_num)


    Auxiliary MO coefficients, same (ao_num,mo_num) layout as
    mo_coef: coefficient of the i-th AO on the j-th aux MO.
    Read from the EZFIO database if the entry exists (master
    reads, then MPI broadcast), otherwise set to the :math:`S^{-1/2}` 
    canonical orthonormalized AOs (ao_ortho_canonical_coef).
    Auxiliary basis, probably to go beyond the basis set for
    integral evaluation (author's reading, 2026-10-31); introduced
    by AbdAmmar (commit 868d5a162, 2023-05-06, "jast 201 added").

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`


 
.. c:var:: double mo_coef_imag


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_coef_imag	(ao_num,mo_num)


    Molecular orbital coefficients on |AO| basis set
    
    mo_coef_imag(i,j) = coefficient of the i-th |AO| on the jth |MO|

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`


 
.. c:var:: double mo_coef_in_ao_ortho_basis


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


 
.. c:var:: double mo_coef_transp


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
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`mos_in_r_array`
       * :c:data:`mos_in_r_array_extra_omp`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`mos_lapl_in_r_array`

 
.. c:var:: char mo_label


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        character*(64)	:: mo_label	


    |MO| coefficients on |AO| basis set
    
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

 
.. c:var:: int mo_num


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
       * :c:data:`ao_sphe_one_e_integrals_from_mo`
       * :c:data:`attachment_numbers_sorted`
       * :c:data:`attachment_orbitals`
       * :c:data:`banned_excitation`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo`
       * :c:data:`cholesky_mo_sp`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`data_one_e_dm_alpha_mo`
       * :c:data:`data_one_e_dm_beta_mo`
       * :c:data:`difference_dm`
       * :c:data:`difference_dm_eigvect`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`full_ijkl_bitmask`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`h_core_ri`
       * :c:data:`int_erf_3_index`
       * :c:data:`kinetic_density_generalized`
       * :c:data:`list_act`
       * :c:data:`list_all_but_del_orb`
       * :c:data:`list_core`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_del`
       * :c:data:`list_inact`
       * :c:data:`list_inact_act`
       * :c:data:`list_virt`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_aux`
       * :c:data:`mo_coef_begin_iteration`
       * :c:data:`mo_coef_imag`
       * :c:data:`mo_coef_in_ao_ortho_basis`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_deriv_1_x`
       * :c:data:`mo_dipole_x`
       * :c:data:`mo_dipole_x_becke`
       * :c:data:`mo_grad_ints`
       * :c:data:`mo_grad_ints_transp`
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
       * :c:data:`mo_prop_dipole_x`
       * :c:data:`mo_pseudo_integrals`
       * :c:data:`mo_pseudo_integrals_local`
       * :c:data:`mo_pseudo_integrals_non_local`
       * :c:data:`mo_sphe_coef`
       * :c:data:`mo_spread_centered_x`
       * :c:data:`mo_spread_x`
       * :c:data:`mo_symmetry`
       * :c:data:`mo_two_e_int_erf_jj`
       * :c:data:`mo_two_e_int_erf_jj_from_ao`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`mos_grad_in_r_array_tranp`
       * :c:data:`mos_grad_in_r_array_transp_3`
       * :c:data:`mos_grad_in_r_array_transp_bis`
       * :c:data:`mos_in_r_array`
       * :c:data:`mos_in_r_array_extra_omp`
       * :c:data:`mos_in_r_array_extra_transp`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`mos_lapl_in_r_array`
       * :c:data:`mos_lapl_in_r_array_tranp`
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
       * :c:data:`ormas_n_orb`
       * :c:data:`partial_occupation`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_mo_mixed`
       * :c:data:`pg_xg`
       * :c:data:`psi_energy_h_core`
       * :c:data:`s_mo_coef`
       * :c:data:`s_mo_sphe_coef`
       * :c:data:`v_ne_psi_energy`
       * :c:data:`z_dipole_moment`

 
.. c:var:: double mo_occ


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


 
.. c:var:: double mo_sphe_coef


    File : :file:`mo_basis/mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_sphe_coef	(ao_sphe_num,mo_num)


    MO coefficients in the basis of spherical harmonics AOs.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`
       * :c:data:`ao_num`
       * :c:data:`ao_sphe_num`
       * :c:data:`ao_sphe_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`s_mo_sphe_coef`

 
.. c:var:: int mo_symmetry


    File : :file:`mo_basis/mo_class.irp.f`

    .. code:: fortran

        integer, allocatable	:: mo_symmetry	(mo_num)


    MOs with the same integer belong to the same irrep.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`


 
.. c:var:: double partial_occupation


    File : :file:`mo_basis/partial_occupation.irp.f`

    .. code:: fortran

        double precision, allocatable	:: partial_occupation	(mo_num)


    Partial occupation numbers for optimal damping algorithm

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void ao_ortho_cano_to_ao(void)


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

 
.. c:function:: void ao_to_mo(void)


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

       * :c:data:`mo_deriv_1_x`
       * :c:data:`mo_dipole_x`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_integrals_n_e_per_atom`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_one_e_integrals`
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

 
.. c:function:: void give_all_mos_and_grad_and_lapl_at_r(void)


    File : :file:`mo_basis/mos_in_r.irp.f`

    .. code:: fortran

        subroutine give_all_mos_and_grad_and_lapl_at_r(r,mos_array,mos_grad_array,mos_lapl_array)


    mos_array(j) = jth MO function evaluated at r
    mos_grad_array(i,j) = derivative of the jth MO along axis i (x, y, z) at r
    mos_lapl_array(i,j) = second derivative of the jth MO along axis i at r:
    d2/dx2 for i=1, d2/dy2 for i=2, d2/dz2 for i=3; the three pure
    second derivatives are stored separately, not their sum

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

 
.. c:function:: void give_all_mos_and_grad_at_r(void)


    File : :file:`mo_basis/mos_in_r.irp.f`

    .. code:: fortran

        subroutine give_all_mos_and_grad_at_r(r,mos_array,mos_grad_array)


    mos_array(j) = jth MO function evaluated at r
    mos_grad_array(i,j) = derivative of the jth MO along axis i (x, y, z) at r

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

 
.. c:function:: void give_all_mos_at_r(void)


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

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`mos_in_r_array`
       * :c:data:`mos_in_r_array_extra_omp`
       * :c:data:`mos_in_r_array_omp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`give_all_aos_at_r`

 
.. c:function:: void mix_mo_jk(void)


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

 
.. c:function:: void mo_as_eigvectors_of_mo_matrix(void)


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine mo_as_eigvectors_of_mo_matrix(matrix,n,m,label,sign,output)


    Set the MOs to the eigenvectors of matrix in the current MO
    basis (mo_coef = mo_coef * eigvectors, via lapack_diag);
    mo_label is set to label. m must equal mo_num. If
    sign == -1, -matrix is diagonalized and the printed
    eigenvalues are negated back; if output, they are printed.
    The sign argument exists because diagonalization sorts the
    eigenvalues; changing the sign is a trick to reverse the
    sorting order when needed (author, 2026-10-31).

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
       * :c:func:`lapack_diag`
       * :c:func:`write_time`

 
.. c:function:: void mo_as_svd_vectors_of_mo_matrix(void)


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine mo_as_svd_vectors_of_mo_matrix(matrix,lda,m,n,label)


    Set the MOs to the left singular vectors of matrix in the
    current MO basis: mo_coef is replaced by mo_coef times U,
    U being the left singular vector matrix of matrix (svd,
    LAPACK dgesdd). m must equal mo_num. A table of the
    singular values and their cumulative sum is printed on
    unit 6. mo_label is set to label.

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

 
.. c:function:: void mo_as_svd_vectors_of_mo_matrix_eig(void)


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine mo_as_svd_vectors_of_mo_matrix_eig(matrix,lda,m,n,eig,label)


    Same as mo_as_svd_vectors_of_mo_matrix, plus eig(i)
    receives the i-th singular value of matrix; singular
    values at or below 1.d-20 are written as 0.d0.
    1.d-20 is a threshold for numerical noise (author-confirmed
    2026-10-31).

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

 
.. c:function:: void mo_coef_new_as_svd_vectors_of_mo_matrix_eig(void)


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

 
.. c:function:: void save_mos(void)


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine save_mos


    Write the current MO basis to the EZFIO database: mo_num,
    mo_label, ao_md5, mo_coef, mo_occ and mo_class.

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

 
.. c:function:: void save_mos_no_occ(void)


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine save_mos_no_occ


    Write only mo_coef to the EZFIO database. mo_num, mo_label,
    ao_md5, mo_occ and mo_class are not written (the
    corresponding set calls are commented out).

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

 
.. c:function:: void save_mos_truncated(void)


    File : :file:`mo_basis/utils.irp.f`

    .. code:: fortran

        subroutine save_mos_truncated(n)


    Write the first n MOs to the EZFIO database: mo_num is set
    to n, and mo_label, ao_md5, mo_occ, mo_class and the first
    n columns of mo_coef are stored.
    mo_occ and mo_class are providers of the full mo_num size and
    are written with that full size even when n < mo_num; this is
    intended and does not matter (author, 2026-10-31).

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

