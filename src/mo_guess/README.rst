========
mo_guess
========

Initial guess for the |MOs|, and nuclear attraction integrals in
orthonormal |AO| bases.

The subroutine :c:func:`hcore_guess` produces the initial set of
|MOs|: it diagonalizes the core Hamiltonian, the nuclear attraction
integrals in the |MO| basis, and stores the eigenvectors as the
initial |MO| coefficients.

The other providers build orthonormal |AO| bases and the nuclear
attraction integrals in them, as a starting point for the |SCF|:

* :c:data:`ao_ortho_lowdin_coef` and :c:data:`ao_ortho_lowdin_overlap`
  for the Lowdin :math:`S^{-1/2}` orthonormalization of the |AOs|, and
  :c:data:`ao_ortho_lowdin_nucl_elec_integrals` for the nuclear
  attraction in that basis
* :c:data:`ao_ortho_canonical_nucl_elec_integrals` for the nuclear
  attraction in the canonical orthonormalized |AO| basis
