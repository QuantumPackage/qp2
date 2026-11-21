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
transcorrelated variant).
