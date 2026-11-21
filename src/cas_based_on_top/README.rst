================
cas_based_on_top
================

On-top pair density of a |CAS|-based wave function.

For a wave function in which all the determinants belong to an active
space, this module computes the on-top pair density :math:`n_2(r, r)`
and the related quantities:

* :c:data:`total_cas_on_top_density`: :math:`n_2(r, r)` at each point
  of the Becke grid, with all the core, inactive and active
  contributions
* :c:func:`give_n2_cas`: :math:`\mu(r)`, :math:`f_\Psi` and
  :math:`n_2^\Psi` for a general |CAS| wave function
* the effective density matrices :c:data:`effective_alpha_dm`,
  :c:data:`effective_beta_dm` and :c:data:`effective_spin_dm`, with
  their gradients, used by the range separated |KS-DFT| calculation
  of the :ref:`module_kohn_sham_rs` module
* the one-electron densities of the core, inactive and active parts:
  :c:data:`core_density`, :c:data:`inact_density`,
  :c:data:`one_e_act_density_alpha` and the related providers

The programs :ref:`cas_based_on_top` and :ref:`cas_based_density`
are small examples of the use of these quantities.
