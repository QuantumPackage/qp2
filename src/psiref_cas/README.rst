==========
psiref_cas
==========

Reference wave function, defined as a |CAS| wave function.

The reference is the |CAS| part of the wave function of the
:ref:`module_determinants` module: `psi_ref` and `psi_ref_coef` are
the determinants and the coefficients of the reference, `N_det_ref`
is its number of determinants and `idx_ref` gives the index of each
reference determinant in the full wave function.

This module is required for the methods that build on a reference
wave function: |CAS-SD|, |MRPT| and |MRCC|.
