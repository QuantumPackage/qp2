========
davidson
========

Abstract module for Davidson's diagonalization.
It contains everything required for the Davidson algorithm, dressed or not. If
a dressing is used, the dressing column should be defined and the
:ref:`davidson_dressed` module should be used. If no dressing is required,
the :ref:`davidson` module should be used, and it has a default zero dressing vector.

The important providers for that module are:

# `psi_energy` which is the expectation value over the wave function (`psi_det`, `psi_coef`) of the Hamiltonian, dressed or not. It uses the general subroutine `u_0_H_u_0`.
# `psi_energy_two_e` which is the expectation value over the wave function (`psi_det`, `psi_coef`) of the standard two-electrons coulomb operator. It uses the general routine `u_0_H_u_0_two_e`.
