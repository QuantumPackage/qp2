.. _module_fci: 
 
.. program:: fci 
 
.. default-role:: option 
 
===
fci
===


|CIPSI| algorithm in the full configuration interaction space.


The user point of view
----------------------

* :ref:`fci` performs |CIPSI| calculations using a stochastic scheme for both
  the selection and the |PT2| contribution,

* :ref:`pt2` computes the |PT2| contribution using the wave function stored in
  the |EZFIO| database.


The main keywords/options for this module are:

* :option:`determinants n_det_max` : maximum number of Slater determinants in
  the |CIPSI| wave function. The :ref:`fci` program will stop when the size of
  the |CIPSI| wave function will exceed :option:`determinants n_det_max`.

* :option:`perturbation pt2_max` : absolute value of the |PT2| to stop the
  |CIPSI| calculation. Once the abs(|PT2|) :math:`<` :option:`perturbation pt2_max`,
  the |CIPSI| calculation stops.

* :option:`determinants n_states` : number of states to consider in the |CIPSI|
  calculation.

* :option:`determinants read_wf` : if |false|, starts with a |ROHF|-like
  determinant, if |true|, starts with the current wave function(s) stored in
  the |EZFIO| directory.

.. note::
   For a multi-state calculation, it is recommended to start with :ref:`cis`
   or :ref:`cisd` wave functions as a guess.

* :option:`determinants expected_s2` : expected value of |S^2| for the
  desired spin multiplicity.

* :option:`determinants s2_eig` : if |true|, systematically add all the
  determinants needed to have a pure value of |S^2|. Also, if |true|, it
  tracks only the states having the good :option:`determinants expected_s2`.




The programmer's point of view
------------------------------

This module was created with the :ref:`module_cipsi` module.

.. seealso::

    The documentation of the :ref:`module_cipsi` module.


 
 
 
EZFIO parameters 
---------------- 
 
.. option:: energy
 
    Calculated Selected |FCI| energy
 
 
.. option:: energy_pt2
 
    Calculated |FCI| energy + |PT2|
 
 
 
Programs 
-------- 
 
 * :ref:`fci` 
 * :ref:`print` 
 * :ref:`pt2` 
 
Providers 
--------- 
 
.. c:var:: do_ddci


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
.. c:var:: do_only_1h1p


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
.. c:var:: do_only_cas


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the FCI case, all those are always false


 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: run:


    File : :file:`pt2.irp.f`

    .. code:: fortran

        subroutine run



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`do_pt2`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_configuration`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy_with_nucl_rep`
       * :c:data:`pt2_relative_error`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print`
       * :c:func:`pt2`
       * :c:func:`scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_summary`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`save_energy`
       * :c:func:`zmq_pt2`
       * :c:func:`zmq_selection`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`variance_match_weight`

 
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

       * :c:func:`run`
       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_fci_energy`
       * :c:func:`ezfio_set_fci_energy_pt2`

