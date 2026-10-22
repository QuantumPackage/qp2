===
mp2
===

Moller-Plesset perturbation theory to second order.

The :ref:`mp2` program computes the |MP2| correlation energy on top
of the reference determinant stored in the |EZFIO| database, and
prints the total energy :math:`E_\text{HF} + E_\text{MP2}` for each
of the :option:`determinants n_states` states. The
:c:func:`compute_mp2` subroutine is the entry point of the
computation, and the `h_apply_mp2` variants apply the Hamiltonian in
the |MP2| space, including the variants restricted to the selected
determinants.
