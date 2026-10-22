.. _module_mrci: 
 
.. program:: mrci 
 
.. default-role:: option 
 
====
mrci
====

|CIPSI| algorithm in the multi-reference |CI| space (|CAS| + singles
and doubles).

The :ref:`mrci` program performs a |CIPSI|-like selected |CI| using
a stochastic scheme for both the selection of the important Slater
determinants and the computation of the |PT2| correction, for the
:option:`determinants n_states` lowest states of the variational
space. The inactive, active and virtual |MOs| are set with the
:ref:`qp_set_mo_class` program.

The boolean options `do_ddci` and `do_only_1h1p` restrict the space
of excitations: removing the purely inactive doubles, or keeping only
the one hole/one particle excitations.

.. seealso::

    The documentation of the :ref:`module_cipsi` module.
 
 
 
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

