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
 * :ref:`print_e_conv` 
 * :ref:`print_wf` 
 * :ref:`save_natorb` 
 * :ref:`save_one_e_dm` 
 * :ref:`save_ortho_mos` 
 * :ref:`write_integrals_erf` 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: routine:


    File : :file:`write_integrals_erf.irp.f`


    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diagonalize_h`
       * :c:func:`print_wf`
       * :c:func:`write_integrals_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_erf_two_e_integrals_ao`
       * :c:func:`save_erf_two_e_integrals_mo`

 
.. c:function:: routine_e_conv:


    File : :file:`print_e_conv.irp.f`

    routine called by :c:func:`print_e_conv`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`ezfio_filename`

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

    routine called by :c:func:`save_one_e_dm`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_mo_alpha`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`save_one_e_dm`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_aux_quantities_data_one_e_dm_alpha_mo`
       * :c:func:`ezfio_set_aux_quantities_data_one_e_dm_beta_mo`

 
.. c:function:: write_ao_basis:


    File : :file:`molden.irp.f`

    .. code:: fortran

        subroutine write_Ao_basis(i_unit_output)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_list_shell_aos`
       * :c:data:`ao_coef`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_charge`
       * :c:data:`ao_l`
       * :c:data:`ao_expo`
       * :c:data:`element_name`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`molden`

 
.. c:function:: write_geometry:


    File : :file:`molden.irp.f`

    .. code:: fortran

        subroutine write_geometry(i_unit_output)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_charge`
       * :c:data:`element_name`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`molden`

 
.. c:function:: write_intro_gamess:


    File : :file:`molden.irp.f`

    .. code:: fortran

        subroutine write_intro_gamess(i_unit_output)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`molden`

 
.. c:function:: write_mo_basis:


    File : :file:`molden.irp.f`

    .. code:: fortran

        subroutine write_Mo_basis(i_unit_output)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mo_coef`
       * :c:data:`ao_num`
       * :c:data:`ao_l_char_space`
       * :c:data:`nucl_charge`
       * :c:data:`ao_nucl`
       * :c:data:`element_name`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`molden`

