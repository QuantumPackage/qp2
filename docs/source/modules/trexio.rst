.. _module_trexio: 
 
.. program:: trexio 
 
.. default-role:: option 
 
======
trexio
======

Reading and writing of TREXIO files
(see https://github.com/trex-coe/trexio).

The :ref:`export_trexio` program exports the wave function in TREXIO
format, to the file named by :c:data:`trexio_file`, with the back-end
selected by :c:data:`backend` (0 for HDF5, 1 for text). The
`export_*` providers gather the basis, the |MOs|, the one- and
two-electron integrals (including the Cholesky decompositions), the
determinants and the density matrices of the wave function for the
export.

The program :ref:`import_trexio_integrals` imports the one- and
two-electron integrals from a TREXIO file into the |EZFIO| database,
and :ref:`import_trexio_determinants` imports determinants from a
TREXIO file.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: backend
 
    Back-end used in TREXIO. 0: HDF5, 1:Text
 
    Default: 0
 
.. option:: trexio_file
 
    Name of the exported TREXIO file
 
    Default: None
 
.. option:: export_basis
 
    If True, export basis set and AOs
 
    Default: True
 
.. option:: export_mos
 
    If True, export MO coefficients
 
    Default: True
 
.. option:: export_cartesian
 
    If False, export everything in the spherical AO basis
 
    Default: True
 
.. option:: export_ao_one_e_ints
 
    If True, export one-electron integrals in AO basis
 
    Default: False
 
.. option:: export_ao_two_e_ints
 
    If True, export two-electron integrals in AO basis
 
    Default: False
 
.. option:: export_ao_two_e_ints_cholesky
 
    If True, export Cholesky-decomposed two-electron integrals in AO basis
 
    Default: False
 
.. option:: export_mo_one_e_ints
 
    If True, export one-electron integrals in MO basis
 
    Default: False
 
.. option:: export_mo_two_e_ints
 
    If True, export two-electron integrals in MO basis
 
    Default: False
 
.. option:: export_mo_two_e_ints_cholesky
 
    If True, export Cholesky-decomposed two-electron integrals in MO basis
 
    Default: False
 
.. option:: export_rdm
 
    If True, export two-body reduced density matrix
 
    Default: False
 
.. option:: export_csf
 
    If True, export wave function in CSF instead of determinants
 
    Default: False
 
 
Programs 
-------- 
 
 * :ref:`export_trexio` 
 * :ref:`import_trexio_determinants` 
 * :ref:`import_trexio_integrals` 
 
Providers 
--------- 
 
.. c:var:: trexio_filename


    File : :file:`trexio_file.irp.f`

    .. code:: fortran

        character*(1024)	:: trexio_filename	


    Name of the TREXIO file

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`backend`
       * :c:data:`ezfio_filename`
       * :c:data:`ezfio_work_dir`
       * :c:data:`trexio_file`


 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: export_trexio_subroutine:


    File : :file:`export_trexio_routines.irp.f`

    .. code:: fortran

        subroutine export_trexio_subroutine(update,full_path)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_kinetic_integrals`
       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_overlap`
       * :c:data:`ao_pseudo_integrals`
       * :c:data:`ao_shell`
       * :c:data:`ao_sphe_coef_normalization_factor`
       * :c:data:`ao_sphe_integrals_n_e`
       * :c:data:`ao_sphe_kinetic_integrals`
       * :c:data:`ao_sphe_num`
       * :c:data:`ao_sphe_one_e_integrals`
       * :c:data:`ao_sphe_overlap`
       * :c:data:`ao_sphe_pseudo_integrals`
       * :c:data:`ao_sphe_shell`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`backend`
       * :c:data:`basis_nucleus_index`
       * :c:data:`basis_type`
       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`do_pseudo`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`export_ao_one_e_ints`
       * :c:data:`export_ao_two_e_ints`
       * :c:data:`export_ao_two_e_ints_cholesky`
       * :c:data:`export_basis`
       * :c:data:`export_cartesian`
       * :c:data:`export_csf`
       * :c:data:`export_mo_one_e_ints`
       * :c:data:`export_mo_two_e_ints`
       * :c:data:`export_mo_two_e_ints_cholesky`
       * :c:data:`export_mos`
       * :c:data:`export_rdm`
       * :c:data:`ezfio_convention`
       * :c:data:`fock_matrix_diag_mo`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_label`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_pseudo_integrals_local`
       * :c:data:`mo_sphe_coef`
       * :c:data:`n_csf`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_charge_remove`
       * :c:data:`nucl_coord_transp`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`
       * :c:data:`nuclear_repulsion`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`prim_coef`
       * :c:data:`prim_expo`
       * :c:data:`prim_normalization_factor`
       * :c:data:`prim_num`
       * :c:data:`primitives_normalized`
       * :c:data:`pseudo_dz_k`
       * :c:data:`pseudo_dz_kl`
       * :c:data:`pseudo_klocmax`
       * :c:data:`pseudo_kmax`
       * :c:data:`pseudo_lmax`
       * :c:data:`pseudo_v_k_transp`
       * :c:data:`pseudo_v_kl_transp`
       * :c:data:`psi_coef`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_csf_coef`
       * :c:data:`psi_det`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_index`
       * :c:data:`shell_num`
       * :c:data:`slater_r_power`
       * :c:data:`trexio_filename`
       * :c:data:`two_e_dm_mo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`export_trexio`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_trexio_trexio_file`
       * :c:func:`get_det_csf_transformation`
       * :c:func:`getenv`
       * :c:func:`system`
       * :c:func:`trexio_assert`

 
.. c:macro:: run:


    File : :file:`import_trexio_integrals.irp.f`

    .. code:: fortran

        subroutine run(f)


    Program to import integrals from TREXIO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`
       * :c:data:`ezfio_filename`
       * :c:data:`ezfio_work_dir`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`import_trexio_determinants`
       * :c:func:`import_trexio_integrals`
       * :c:func:`scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_ao_one_e_ints_ao_integrals_kinetic`
       * :c:func:`ezfio_set_ao_one_e_ints_ao_integrals_n_e`
       * :c:func:`ezfio_set_ao_one_e_ints_ao_integrals_overlap`
       * :c:func:`ezfio_set_ao_one_e_ints_ao_integrals_pseudo`
       * :c:func:`ezfio_set_ao_one_e_ints_ao_one_e_integrals`
       * :c:func:`ezfio_set_ao_one_e_ints_io_ao_integrals_kinetic`
       * :c:func:`ezfio_set_ao_one_e_ints_io_ao_integrals_n_e`
       * :c:func:`ezfio_set_ao_one_e_ints_io_ao_integrals_overlap`
       * :c:func:`ezfio_set_ao_one_e_ints_io_ao_integrals_pseudo`
       * :c:func:`ezfio_set_ao_one_e_ints_io_ao_one_e_integrals`
       * :c:func:`ezfio_set_ao_two_e_ints_do_ao_cholesky`
       * :c:func:`ezfio_set_ao_two_e_ints_io_ao_cholesky`
       * :c:func:`ezfio_set_ao_two_e_ints_io_ao_two_e_integrals`
       * :c:func:`ezfio_set_mo_one_e_ints_io_mo_one_e_integrals`
       * :c:func:`ezfio_set_mo_one_e_ints_mo_one_e_integrals`
       * :c:func:`ezfio_set_mo_two_e_ints_io_mo_cholesky`
       * :c:func:`ezfio_set_mo_two_e_ints_io_mo_two_e_integrals`
       * :c:func:`ezfio_set_nuclei_io_nuclear_repulsion`
       * :c:func:`ezfio_set_nuclei_nuclear_repulsion`
       * :c:func:`ezfio_set_pseudo_do_pseudo`
       * :c:func:`insert_into_ao_integrals_map`
       * :c:func:`map_append`
       * :c:func:`map_save_to_disk`
       * :c:func:`map_sort`
       * :c:func:`map_unique`
       * :c:func:`trexio_assert`
       * :c:func:`two_e_integrals_index`

