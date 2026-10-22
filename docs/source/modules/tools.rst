.. _module_tools: 
 
.. program:: tools 
 
.. default-role:: option 
 
=====
tools
=====

Useful tools are grouped in this module.
 
 
 
Programs 
-------- 
 
 * :ref:`analyze_ccsd` 
 * :ref:`attachement_orb` 
 * :ref:`benchmark` 
 * :ref:`cas_complete` 
 * :ref:`check_orthonormality` 
 * :ref:`diagonalize_h` 
 * :ref:`fcidump` 
 * :ref:`fcidump_pyscf` 
 * :ref:`four_idx_transform` 
 * :ref:`molden` 
 * :ref:`pouet` 
 * :ref:`print_ci_vectors` 
 * :ref:`print_detweights` 
 * :ref:`print_dipole` 
 * :ref:`print_energy` 
 * :ref:`print_hamiltonian` 
 * :ref:`print_hmat` 
 * :ref:`print_sorted_wf_coef` 
 * :ref:`print_var_energy` 
 * :ref:`print_wf` 
 * :ref:`rotate_mos` 
 * :ref:`save_natorb` 
 * :ref:`save_natorb_no_ov_rot` 
 * :ref:`save_natorb_no_ref` 
 * :ref:`save_one_e_dm` 
 * :ref:`save_ortho_mos` 
 * :ref:`sort_by_fock_energies` 
 * :ref:`sort_wf` 
 * :ref:`swap_mos` 
 * :ref:`truncate_wf` 
 * :ref:`write_integrals_erf` 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: do_print:


    File : :file:`check_orthonormality.irp.f`

    .. code:: fortran

        subroutine do_print



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mo_overlap`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`check_orthonormality`

 
.. c:macro:: main:


    File : :file:`print_detweights.irp.f`

    .. code:: fortran

        subroutine main()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`
       * :c:func:`dsort`
       * :c:func:`get_excitation_degree`

 
.. c:macro:: molden_attachment:


    File : :file:`attachement_orb.irp.f`

    .. code:: fortran

        subroutine molden_attachment


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
       * :c:data:`attachment_numbers_sorted`
       * :c:data:`attachment_orbitals`
       * :c:data:`element_name`
       * :c:data:`ezfio_filename`
       * :c:data:`n_attachment`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_list_shell_aos`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`attachement_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`isort`

 
.. c:macro:: print_exc:


    File : :file:`print_detweights.irp.f`

    .. code:: fortran

        subroutine print_exc()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_detweights`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_degree`

 
.. c:macro:: routine:


    File : :file:`write_integrals_erf.irp.f`

    .. code:: fortran

        subroutine routine



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`analyze_ccsd`
       * :c:func:`diagonalize_h`
       * :c:func:`print_ci_vectors`
       * :c:func:`print_sorted_wf_coef`
       * :c:func:`print_var_energy`
       * :c:func:`print_wf`
       * :c:func:`sort_wf`
       * :c:func:`truncate_wf`
       * :c:func:`write_integrals_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_erf_two_e_integrals_ao`
       * :c:func:`save_erf_two_e_integrals_mo`

 
.. c:macro:: routine_s2:


    File : :file:`truncate_wf.irp.f`

    .. code:: fortran

        subroutine routine_s2



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_configuration`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`weight_configuration`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`truncate_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_wavefunction_general`

 
.. c:macro:: routine_save_one_e_dm:


    File : :file:`save_one_e_dm.irp.f`

    .. code:: fortran

        subroutine routine_save_one_e_dm


    routine called by :c:func:`save_one_e_dm`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_ao`
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
       * :c:func:`ezfio_set_aux_quantities_data_one_e_dm_tot_ao`

 
.. c:macro:: run:


    File : :file:`print_hmat.irp.f`

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

       * :c:func:`print`
       * :c:func:`print_energy`
       * :c:func:`print_hamiltonian`
       * :c:func:`print_hmat`
       * :c:func:`pt2`
       * :c:func:`scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_degree`
       * :c:func:`i_h_j`

