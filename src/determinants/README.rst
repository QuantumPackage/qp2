============
determinants
============

Storage of the wave function in a basis of orthogonal Slater
determinants built on the restricted spin-orbital basis, and the basic
operations on it.

The main providers for this module are:

* :option:`determinants n_states`: number of states to be computed
* :c:data:`psi_det`: the list of determinants in the wave function,
  stored as spin bitmasks of the occupied spin orbitals, for all
  :option:`determinants n_states` states. It is used by many
  routines and providers of |QP|
* :c:data:`psi_coef`: the list of coefficients, for all states and
  all determinants
* the sorted variants :c:data:`psi_det_sorted` and
  :c:data:`psi_coef_sorted`, ordered by decreasing norm contribution,
  and the :c:data:`psi_bilinear_matrix` variants, which store the
  wave function in the bilinear form :math:`D_\alpha^\dagger C
  D_\beta` used by the distributed computations

The main routines for this module are:

* :c:func:`i_h_j`: computes the Hamiltonian matrix element between
  two arbitrary Slater determinants
* :c:func:`i_h_j_s2`: computes the Hamiltonian and |S^2| matrix
  element between two arbitrary Slater determinants
* :c:func:`i_h_j_verbose`: returns the decomposition in terms of
  one- and two-body components of the Hamiltonian matrix element
  between two arbitrary Slater determinants, and the fermionic phase
  factor
* :c:func:`i_h_psi`: computes the Hamiltonian matrix element between
  an arbitrary Slater determinant and a wave function composed of a
  sum of arbitrary Slater determinants

The other files of the module provide the creation of excitations
(`create_excitations.irp.f`), the application of the Hamiltonian to
the wave function (`h_apply.irp.f`), the |S^2| operator
(`s2.irp.f`), the one-, two- and three-electron density matrices
(`density_matrix.irp.f`, `two_e_density_matrix.irp.pouet`,
`tr_density_matrix.irp.f`), the dipole moments, the pruning of the
wave function (`prune_wf.irp.f`), the projection on a |CAS|
(`psi_cas.irp.f`), and the saving of the wave function in the |EZFIO|
database (`save_wavefunction`).

For an example of how to use these routines and providers, take a
look at :file:`example.irp.f`.
