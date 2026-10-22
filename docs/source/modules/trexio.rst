.. _module_trexio: 
 
.. program:: trexio 
 
.. default-role:: option 
 
======
trexio
======

Module for handling TREXIO files.
See https://github.com/trex-coe/trexio
 
 
 
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

