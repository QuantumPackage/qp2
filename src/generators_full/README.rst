===============
generators_full
===============

Generators for the perturbation energy, as all the determinants of the
variational space.

The generators are the determinants from which the external space is
generated. This module defines them as the determinants of the wave
function truncated to the largest norm contributions:
:c:data:`n_det_generators` is the smallest number of determinants, in
the sorted order, whose norm contribution reaches
:option:`determinants threshold_generators`, and it is at least 1.
:c:data:`psi_det_generators` and :c:data:`psi_coef_generators` are
the selected determinants and coefficients,
:c:data:`psi_det_sorted_gen` their sorted order, and
:c:data:`degree_max_generators` the largest excitation degree (with
respect to |HF|) reached by the selected generators.

This module is intended to be included in the :file:`NEED` file to
define a full set of generators.
