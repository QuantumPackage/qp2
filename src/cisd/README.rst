====
cisd
====

Configuration interaction with singles and doubles.

The user point of view
----------------------

The :ref:`cisd` program performs the |CI| of all single and double
excitations on top of a reference Slater determinant with ROHF-like
occupancy, disregarding spatial symmetry, and computes the
:option:`determinants n_states` lowest eigenstates of that |CI|
matrix. This program can be very useful to:

* **Ground state calculations**: generate a guess for the ground
  state wave function if one is not sure that the :ref:`scf` program
  gave the lowest |SCF| solution. In combination with
  :ref:`save_natorb` it can produce new |MOs| in order to reperform a
  :ref:`scf` optimization.

* **Excited states calculations**: generate a guess for all the
  :option:`determinants n_states` wave functions, that will be used
  by the :ref:`fci` program.

The :ref:`lccsd` program performs the linearized |CCSD| computation
on the same reference.

The main keywords/options to be used are:

* :option:`determinants n_states` : number of states to consider for
  the |CISD| calculation
* :option:`determinants s2_eig` : force all states to have the
  desired value of :math:`S^2`
* :option:`determinants expected_s2` : desired value of :math:`S^2`

The programmer point of view
----------------------------

This module has been built by setting the following rules:

* The only generator determinant is the Hartree-Fock determinant
  (single-reference method)
* All generated determinants are included in the wave function (no
  perturbative selection)

These rules are set in the ``H_apply.irp.f`` file.
