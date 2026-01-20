.. _module_tools: 
 
.. program:: tools 
 
.. default-role:: option 
 
=====
tools
=====

Useful tools are grouped in this module.
 
 
 
Programs 
-------- 
 
 * :ref:`attachement_orb` 
 * :ref:`cas_complete` 
 * :ref:`diagonalize_h` 
 * :ref:`fcidump` 
 * :ref:`fcidump_pyscf` 
 * :ref:`four_idx_transform` 
 * :ref:`guess_hcore` 
 * :ref:`guess_huckel` 
 * :ref:`molden` 
 * :ref:`print_ci_vectors` 
 * :ref:`print_detweights` 
 * :ref:`print_dipole` 
 * :ref:`print_energy` 
 * :ref:`print_hamiltonian` 
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
 
.. c:function:: molden_attachment:


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

 
.. c:function:: print_exc:


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

 
.. c:function:: routine_s2:


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

