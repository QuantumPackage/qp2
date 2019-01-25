.. _module_hartree_fock: 
 
.. program:: hartree_fock 
 
.. default-role:: option 
 
============
hartree_fock
============


The Hartree-Fock module performs *Restricted* Hartree-Fock calculations (the
spatial part of the |MOs| is common for alpha and beta spinorbitals).

The Hartree-Fock in an SCF and therefore is based on the ``scf_utils`` structure.
It performs the following actions:

#. Compute/Read all the one- and two-electron integrals, and store them in memory

#. Check in the |EZFIO| database if there is a set of |MOs|. If there is, it
   will read them as initial guess. Otherwise, it will create a guess.
#. Perform the |SCF| iterations

The definition of the Fock matrix is in :file:`hartree_fock fock_matrix_hf.irp.f`
For the keywords related to the |SCF| procedure, see the ``scf_utils`` directory where you will find all options.
The main are:

# :option:`scf_utils thresh_scf`

# :option:`scf_utils level_shift`

At each iteration, the |MOs| are saved in the |EZFIO| database. Hence, if the calculation
crashes for any unexpected reason, the calculation can be restarted by running again
the |SCF| with the same |EZFIO| database.

The `DIIS`_ algorithm is implemented, as well as the `level-shifting`_ method.
If the |SCF| does not converge, try again with a higher value of :option:`level_shift`.

To start a calculation from scratch, the simplest way is to remove the
``mo_basis`` directory from the |EZFIO| database, and run the |SCF| again.




.. _DIIS: https://en.wikipedia.org/w/index.php?title=DIIS
.. _level-shifting: https://doi.org/10.1002/qua.560070407



 
 
 
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


    Alpha Fock matrix in AO basis set

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
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`hf_energy`

 
.. c:var:: ao_two_e_integral_beta


    File : :file:`hartree_fock/fock_matrix_hf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_integral_beta	(ao_num,ao_num)


    Alpha Fock matrix in AO basis set

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
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`hf_energy`

 
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

    Create a MO guess if no MOs are present in the EZFIO directory

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_coef`
       * :c:data:`mo_guess_type`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`ao_ortho_lowdin_coef`
       * :c:data:`mo_label`

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

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`

 
.. c:function:: run:


    File : :file:`hartree_fock/scf.irp.f`

    Run SCF calculation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`scf_energy`
       * :c:data:`mo_label`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2`
       * :c:func:`scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_hartree_fock_energy`
       * :c:func:`roothaan_hall_scf`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`level_shift`
       * :c:data:`mo_coef`

