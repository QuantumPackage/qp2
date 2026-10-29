.. _module_mo_optimization: 
 
.. program:: mo_optimization 
 
.. default-role:: option 
 
mo_optimization
===============

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Programs 
-------- 
 
 * :ref:`cipsi_orb_opt` 
 * :ref:`debug_gradient_list_opt` 
 * :ref:`debug_gradient_opt` 
 * :ref:`debug_hessian_list_opt` 
 * :ref:`debug_hessian_opt` 
 * :ref:`orb_opt` 
 
Providers 
--------- 
 
.. c:var:: bool do_ddci


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
.. c:var:: bool do_only_1h1p


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
.. c:var:: bool do_only_cas


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void save_energy(void)


    File : :file:`save_energy.irp.f`

    .. code:: fortran

        subroutine save_energy(E,pt2)


    Saves the energy in |EZFIO|.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_fci_energy`
       * :c:func:`ezfio_set_fci_energy_pt2`

