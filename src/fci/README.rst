===
fci
===


|CIPSI| algorithm in the full configuration interaction space.


The user point of view
----------------------

* :c:func:`fci` performs |CIPSI| calculations using a stochastic scheme for both the selection and the |PT2| contribution,
* :c:func:`pt2` computes the |PT2| contribution using the wave function stored in the |EZFIO|
  database.


The main keywords/options for this module are:

* :option:`determinants n_det_max` : maximum number of Slater determinants in the CIPSI wave function. The :command:`fci` program will stop when the size of the CIPSI wave function will exceed :option:`determinants n_det_max`.

* :option:`perturbation pt2_max` : absolute value of the |PT2| to stop the CIPSI calculation. Once the |PT2| :math:`<` :option:`perturbation pt2_max`, the CIPSI calculation stops.

* :option:`determinants n_states` : number of states to consider in the CIPSI calculation.

* :option:`determinants read_wf` : if False, starts with a ROHF-like determinant, if True, starts with the current wave function(s) stored in the |EZFIO| folder.

.. note::
   For a multi-state calculation, it is recommended to start with :c:func:`cis` or :c:func:`cisd`
   wave functions as a guess.

* :option:`determinants s2_eig` : if True, systematically add all the determinants needed to have a pure value of :math:`S^2`. Also, if True, it tracks only the states having the good :option:`determinants expected_s2`.

.. note::
   For a multi-state calculation, it is recommended to start with :c:func:`cis` or :c:func:`cisd`
   wave functions as a guess.

* :option:`determinants expected_s2` : expected value of :math:`S^2` for the desired spin multiplicity.



The programmer point of view
----------------------------

This module have been created with the :ref:`cipsi` module.

.. seealso::

    The documentation of the :ref:`cipsi` module.


