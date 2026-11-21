==============
selectors_full
==============

Selectors for the perturbation energy, as all the determinants of the
wave function.

The selectors are the determinants that enter the numerator of the
perturbation sum, in the matrix elements :math:`\langle i | H |
\psi \rangle`. This module defines them as the determinants of the
wave function truncated to the largest norm contributions:
:c:data:`n_det_selectors` is the smallest number of determinants, in
the sorted order, whose norm contribution reaches the threshold
:c:data:`threshold_selectors`, the square root of
:option:`determinants threshold_generators`, and it always includes
the generators. The selected determinants and their coefficients are
:c:data:`psi_selectors` and :c:data:`psi_selectors_coef`.
