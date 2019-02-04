===
cis
===

This module contains a |CIS| program.

The user point of view
----------------------

The :ref:`cis` program performs the CI to obtain the ROHF reference + all
single excitations on top of it. This program can be very useful to:

* **Ground state calculations**: generate a guess for the ground state wave
  function if one is not sure that the :ref:`scf` program gave the lowest |SCF|
  solution. In combination with :ref:`save_natorb` it can produce new |MOs| in
  order to reperform an :ref:`scf` optimization.

* **Excited states calculations**: generate guesses for all the
  :option:`determinants n_states` wave functions, that will be used by the
  :ref:`fci` program.


The main keywords/options to be used are:

* :option:`determinants n_states`: number of states to consider for the |CIS| calculation

* :option:`determinants s2_eig`: force all states to have the desired value of |S^2|

* :option:`determinants expected_s2`: desired value of |S^2|




The programmer's point of view
------------------------------

This module was built by setting the following rules:

* The only generator determinant is the Hartree-Fock (single-reference method)
* All generated singly excited determinants are included in the wave function (no perturbative
  selection)

These rules are set in the ``H_apply.irp.f`` file.


