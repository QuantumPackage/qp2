====
cisd
====

This module contains a CI of single and double excitations.

The user point of view
----------------------

The :command:`cisd` program performs the CI of the ROHF-like + all single and double excitations on top of it.
This program can be very useful to :

* **Ground state calculations**: generate a guess for the ground state wave function if one is not sure that the :c:func:`scf` program gave the lowest SCF solution. In combination with :c:func:`save_natorb` it can produce new |MOs| in order to reperform an :c:func:`scf` optimization.

* **Excited states calculations**: generate guess for all the :option:`determinants n_states` wave functions, that will be used by the :c:func:`fci` program.


The main keywords/options to be used are:

* :option:`determinants n_states` : number of states to consider for the |cisd| calculation

* :option:`determinants s2_eig` : force all states to have the desired value of :math:`S^2`

* :option:`determinants expected_s2` : desired value of :math:`S^2`

The programmer point of view
----------------------------

This module have been built by setting the following rules:


* The only generator determinant is the Hartree-Fock (single-reference method)
* All generated determinants are included in the wave function (no perturbative
  selection)

These rules are set in the ``H_apply.irp.f`` file.


