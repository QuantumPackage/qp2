==========
iterations
==========

Iteration counters and energy history, for the extrapolation to the
|FCI| limit.

:c:data:`n_iter` and :c:data:`n_iter_max` count the iterations
(:c:func:`increment_n_iter`), and the energies of the iterations are
stored in :c:data:`energy_iterations` (the variational energies) and
:c:data:`pt2_iterations` (the |PT2| energies).
:c:data:`extrapolated_energy` is the value extrapolated to the |FCI|
limit, computed and printed by :c:func:`print_extrapolated_energy`,
and the summary of the calculation is printed by
:c:func:`print_summary` (or :c:func:`print_summary_tc` for the
transcorrelated variant). Both printers take the state irrep labels from
:ref:`module_symmetry_det` when they are usable (a ``# Irrep`` row, with a
``# Purity`` row when a state has no definite irrep, in the summary table;
an ``Irrep:`` field in the extrapolation headings) and are unchanged
otherwise; see the ``symmetry_det`` documentation.
