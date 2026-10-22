.. _module_mrci: 
 
.. program:: mrci 
 
.. default-role:: option 
 
====
mrci
====


|CIPSI| algorithm in the multi-reference CI space (CAS + Singles and Doubles).


This module is the same as the :ref:`fci` module, except for the choice of the
generator and selector determinants.

The inactive, active and virtual |MOs| will need to be set with the
:ref:`qp_set_mo_class` program.

.. seealso::

    The documentation of the :ref:`fci` module.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: energy
 
    Calculated Selected CASSD energy
 
 
.. option:: energy_pt2
 
    Calculated CASSD energy + PT2
 
 
.. option:: do_ddci
 
    If true, remove purely inactive double excitations
 
    Default: False
 
.. option:: do_only_1h1p
 
    If true, do only one hole/one particle excitations
 
    Default: False
 
 
Programs 
-------- 
 
 * :ref:`mrci` 
 
Providers 
--------- 
 
.. c:var:: do_only_cas


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_cas	


    In the CAS+SD case, always false


 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: save_energy:


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

       * :c:func:`ezfio_set_mrci_energy`
       * :c:func:`ezfio_set_mrci_energy_pt2`

