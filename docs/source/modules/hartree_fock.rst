.. _module_hartree_fock: 
 
.. program:: hartree_fock 
 
.. default-role:: option 
 
============
hartree_fock
============


Quick description
-----------------

The :ref:`scf` program performs *Restricted* Hartree-Fock
calculations (the spatial part of the |MOs| is common for alpha and beta
spinorbitals).

.. seealso:: 
   To see the keywords/options associated to the |SCF| algorithm itself,  
   see the documentation of the :ref:`module_scf_utils` module. 


More advanced description
-------------------------

The Hartree-Fock algorithm is a |SCF| and therefore is based on the
:ref:`module_scf_utils` module. 

The Fock matrix is defined in :file:`fock_matrix_hf.irp.f`.


.. seealso:: 
   For a more detailed description of the |SCF| structure, 
   see the documentation of the :ref:`module_scf_utils` module. 


 
 
 
EZFIO parameters 
---------------- 
 
.. option:: energy
 
    Energy HF
 
 
 
Programs 
-------- 
 
 * :ref:`scf` 
 
Providers 
--------- 
 
.. c:var:: ao_two_e_integral_alpha


    File : :file:`hartree_fock/fock_matrix_hf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_integral_beta	(ao_num,ao_num)


    Alpha and Beta Fock matrices in AO basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_two_e_integral_alpha_chol`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_ao_cholesky`
       * :c:data:`do_direct_integrals`
       * :c:data:`is_periodic`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`
       * :c:data:`read_ao_two_e_integrals`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`
       * :c:data:`use_cgtos`
       * :c:data:`use_only_lr`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`hf_energy`

 
.. c:var:: ao_two_e_integral_alpha_chol


    File : :file:`hartree_fock/fock_matrix_hf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_integral_alpha_chol	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_integral_beta_chol	(ao_num,ao_num)


    Alpha and Beta Fock matrices in AO basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`cholesky_ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`qp_max_mem`
       * :c:data:`scf_density_matrix_ao`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha`

 
.. c:var:: ao_two_e_integral_beta


    File : :file:`hartree_fock/fock_matrix_hf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_integral_beta	(ao_num,ao_num)


    Alpha and Beta Fock matrices in AO basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_two_e_integral_alpha_chol`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_ao_cholesky`
       * :c:data:`do_direct_integrals`
       * :c:data:`is_periodic`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`
       * :c:data:`read_ao_two_e_integrals`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`
       * :c:data:`use_cgtos`
       * :c:data:`use_only_lr`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`hf_energy`

 
.. c:var:: ao_two_e_integral_beta_chol


    File : :file:`hartree_fock/fock_matrix_hf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_integral_alpha_chol	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_integral_beta_chol	(ao_num,ao_num)


    Alpha and Beta Fock matrices in AO basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`cholesky_ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`qp_max_mem`
       * :c:data:`scf_density_matrix_ao`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha`

 
.. c:var:: extra_e_contrib_density


    File : :file:`hartree_fock/hf_energy.irp.f`

    .. code:: fortran

        double precision	:: extra_e_contrib_density	


    Extra contribution to the SCF energy coming from the density.
    
    For a Hartree-Fock calculation: extra_e_contrib_density = 0
    
    For a Kohn-Sham or Range-separated Kohn-Sham: the exchange/correlation - trace of the V_xc potential

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`scf_energy`

 
.. c:var:: fock_matrix_ao_alpha


    File : :file:`hartree_fock/fock_matrix_hf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_ao_alpha	(ao_num,ao_num)
        double precision, allocatable	:: fock_matrix_ao_beta	(ao_num,ao_num)


    Alpha Fock matrix in AO basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_two_e_integral_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao`
       * :c:data:`fock_matrix_mo_alpha`
       * :c:data:`fock_matrix_mo_beta`
       * :c:data:`mcscf_fock_alpha_ao`
       * :c:data:`scf_energy`

 
.. c:var:: fock_matrix_ao_beta


    File : :file:`hartree_fock/fock_matrix_hf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_ao_alpha	(ao_num,ao_num)
        double precision, allocatable	:: fock_matrix_ao_beta	(ao_num,ao_num)


    Alpha Fock matrix in AO basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_two_e_integral_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao`
       * :c:data:`fock_matrix_mo_alpha`
       * :c:data:`fock_matrix_mo_beta`
       * :c:data:`mcscf_fock_alpha_ao`
       * :c:data:`scf_energy`

 
.. c:var:: hf_energy


    File : :file:`hartree_fock/hf_energy.irp.f`

    .. code:: fortran

        double precision	:: hf_energy	
        double precision	:: hf_two_electron_energy	
        double precision	:: hf_one_electron_energy	


    Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`nuclear_repulsion`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`


 
.. c:var:: hf_kinetic_energy


    File : :file:`hartree_fock/hf_energy.irp.f`

    .. code:: fortran

        double precision	:: hf_kinetic_energy	
        double precision	:: hf_n_e_energy	


    Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_kinetic_integrals`
       * :c:data:`ao_num`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`


 
.. c:var:: hf_n_e_energy


    File : :file:`hartree_fock/hf_energy.irp.f`

    .. code:: fortran

        double precision	:: hf_kinetic_energy	
        double precision	:: hf_n_e_energy	


    Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_kinetic_integrals`
       * :c:data:`ao_num`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`


 
.. c:var:: hf_one_electron_energy


    File : :file:`hartree_fock/hf_energy.irp.f`

    .. code:: fortran

        double precision	:: hf_energy	
        double precision	:: hf_two_electron_energy	
        double precision	:: hf_one_electron_energy	


    Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`nuclear_repulsion`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`


 
.. c:var:: hf_two_electron_energy


    File : :file:`hartree_fock/hf_energy.irp.f`

    .. code:: fortran

        double precision	:: hf_energy	
        double precision	:: hf_two_electron_energy	
        double precision	:: hf_one_electron_energy	


    Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`nuclear_repulsion`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: create_guess:


    File : :file:`hartree_fock/scf.irp.f`

    .. code:: fortran

        subroutine create_guess


    Create a MO guess if no MOs are present in the EZFIO directory

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
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

       * :c:func:`ezfio_has_mo_basis_mo_coef`
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

 
.. c:function:: main:


    File : :file:`hartree_fock/print_scf_int.irp.f`

    .. code:: fortran

        subroutine main()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`fock_matrix_ao`
       * :c:data:`fock_matrix_mo_alpha`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_scf_int`

 
.. c:function:: print_fock_diag:


    File : :file:`hartree_fock/print_fock_diag.irp.f`

    .. code:: fortran

        subroutine print_fock_diag



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_mo`
       * :c:data:`mo_num`

 
.. c:function:: print_pseudo_overlap:


    File : :file:`hartree_fock/print_pseudo_overlap.irp.f`

    .. code:: fortran

        subroutine print_pseudo_overlap



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`list_core`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`

 
.. c:function:: print_scf_int:


    File : :file:`hartree_fock/print_scf_int.irp.f`

    .. code:: fortran

        subroutine print_scf_int



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`main`

 
.. c:function:: run:


    File : :file:`hartree_fock/scf.irp.f`

    .. code:: fortran

        subroutine run


    Run SCF calculation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`json_int_fmt`
       * :c:data:`json_unit`
       * :c:data:`mo_label`
       * :c:data:`scf_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`casscf`
       * :c:func:`scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_hartree_fock_energy`
       * :c:func:`json_close`
       * :c:func:`roothaan_hall_scf`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`level_shift`
       * :c:data:`mo_coef`

