========
davidson
========

Abstract module for Davidson's diagonalization.
It contains everything required for the Davidson algorithm, dressed or
not. If a dressing is used, the dressing column should be defined and
the :ref:`module_davidson_dressed` module should be used. If no dressing
is required, the :ref:`module_davidson` module should be used, and it
has a default zero dressing vector.

The important providers for that module are:

#. :c:data:`psi_energy` which is the expectation value over the wave
   function (:c:data:`psi_det`, :c:data:`psi_coef`) of the Hamiltonian,
   dressed or not. It uses the general subroutine :c:func:`u_0_H_u_0`.

#. :c:data:`psi_energy_two_e` which is the expectation value over the
   wave function (:c:data:`psi_det`, :c:data:`psi_coef`) of the standard
   two-electron Coulomb operator. It uses the general routine
   :c:func:`u_0_H_u_0_two_e`.
