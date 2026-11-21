.. _module_psiref_utils: 
 
.. program:: psiref_utils 
 
.. default-role:: option 
 
============
psiref_utils
============

Utilities related to the use of a reference wave function. This
module needs to be loaded with any `psi_ref_*` module.

The subroutine `extract_ref` replaces the total wave function by its
normalized projection on the reference. The complement of the
reference in the wave function, `psi_non_ref`, is provided with its
transposed coefficients, and the hole and particle operators that
connect the reference to its complement are
`holes_operators` and `particles_operators`.
 
