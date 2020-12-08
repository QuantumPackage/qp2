.. _module_tools: 
 
.. program:: tools 
 
.. default-role:: option 
 
=====
tools
=====

Useful tools are grouped in this module.
 
 
 
Programs 
-------- 
 
 * :ref:`diagonalize_h` 
 * :ref:`fcidump` 
 * :ref:`four_idx_transform` 
 * :ref:`molden` 
 * :ref:`print_ci_vectors` 
 * :ref:`print_e_conv` 
 * :ref:`print_wf` 
 * :ref:`rotate_mos` 
 * :ref:`save_natorb` 
 * :ref:`save_one_e_dm` 
 * :ref:`save_ortho_mos` 
 * :ref:`sort_by_fock_energies` 
 * :ref:`swap_mos` 
 * :ref:`write_integrals_erf` 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: print_energy:


    File : :file:`print_energy.irp.f`

    .. code:: fortran

        subroutine print_energy


    Prints the energy of the wave function stored in the |EZFIO| directory.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`run`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`read_wf`

 
.. c:function:: print_hamiltonian:


    File : :file:`print_hamiltonian.irp.f`

    .. code:: fortran

        subroutine print_hamiltonian


    Prints the Hamiltonian matrix defined in the space of determinants
    present in the |EZFIO| directory.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`run`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`read_wf`

 
.. c:function:: routine:


    File : :file:`write_integrals_erf.irp.f`

    .. code:: fortran

        subroutine routine



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diagonalize_h`
       * :c:func:`print_ci_vectors`
       * :c:func:`print_wf`
       * :c:func:`write_integrals_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_erf_two_e_integrals_ao`
       * :c:func:`save_erf_two_e_integrals_mo`

 
.. c:function:: routine_e_conv:


    File : :file:`print_e_conv.irp.f`

    .. code:: fortran

        subroutine routine_e_conv


    routine called by :c:func:`print_e_conv`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_e_conv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_get_iterations_energy_iterations`
       * :c:func:`ezfio_get_iterations_n_det_iterations`
       * :c:func:`ezfio_get_iterations_n_iter`
       * :c:func:`ezfio_get_iterations_pt2_iterations`

 
.. c:function:: routine_save_one_e_dm:


    File : :file:`save_one_e_dm.irp.f`

    .. code:: fortran

        subroutine routine_save_one_e_dm


    routine called by :c:func:`save_one_e_dm`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_ao_alpha`
       * :c:data:`one_e_dm_mo_alpha`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`save_one_e_dm`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_aux_quantities_data_one_e_dm_alpha_ao`
       * :c:func:`ezfio_set_aux_quantities_data_one_e_dm_alpha_mo`
       * :c:func:`ezfio_set_aux_quantities_data_one_e_dm_beta_ao`
       * :c:func:`ezfio_set_aux_quantities_data_one_e_dm_beta_mo`

 
.. c:function:: run:


    File : :file:`print_hamiltonian.irp.f`

    .. code:: fortran

        subroutine run



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_energy`
       * :c:func:`print_hamiltonian`
       * :c:func:`pt2`
       * :c:func:`scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_h_j`

