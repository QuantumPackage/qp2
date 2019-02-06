============
determinants
============

Contains everything for the computation of the Hamiltonian matrix elements in the basis of orthogonal Slater determinants built on a restricted spin-orbitals basis.

The main providers for this module are:

* :option:`determinants n_states`: number of states to be computed
* :c:data:`psi_det`: list of determinants in the wave function used in many routines/providers of the |QP|.
* :c:data:`psi_coef`: list of coefficients, for all :option:`determinants n_states` states, and all determinants.

The main routines for this module are:

* :c:func:`i_H_j`: computes the Hamiltonian matrix element between two arbitrary Slater determinants.
* :c:func:`i_H_j_s2`: computes the Hamiltonian and (|S^2|) matrix element between two arbitrary Slater determinants.
* :c:func:`i_H_j_verbose`: returns the decomposition in terms of one- and two-body components of the Hamiltonian matrix elements between two arbitrary Slater determinants. Also return the fermionic phase factor.
* :c:func:`i_H_psi`: computes the Hamiltonian matrix element between an arbitrary Slater determinant and a wave function composed of a sum of arbitrary Slater determinants.


For an example of how to use these routines and providers, take a look at :file:`example.irp.f`.
