.. _module_davidson_dressed: 
 
.. program:: davidson_dressed 
 
.. default-role:: option 
 
================
davidson_dressed
================

Davidson diagonalization of the dressed Hamiltonian.

In the |CIPSI| iteration, the Hamiltonian is dressed (see the
:ref:`module_dressing` module) to bring its spectrum in the selected
space closer to the spectrum of the full Hamiltonian. This module
diagonalizes the dressed Hamiltonian with the Davidson algorithm and
provides the dressed CI energies `CI_energy_dressed`, the electronic
energies `CI_electronic_energy_dressed`, the eigenvectors
`CI_eigenvectors_dressed` and the diagonal matrix of the dressed
Hamiltonian `h_matrix_dressed`. The non symmetric variant, which does
not store the matrix, is in `nonsym_diagonalize_ci.irp.f`.
 
