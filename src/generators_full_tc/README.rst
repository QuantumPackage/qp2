==================
generators_full_tc
==================

Generators for the transcorrelated wave function, as all the
determinants of the variational space.

As in :ref:`module_generators_full`, the generators are the
determinants of the wave function truncated to the largest norm
contributions, with
:option:`determinants threshold_generators`. The difference is that
the norm contributions are those of the transcorrelated wave
function, so that the truncation applies to the transcorrelated
expansion.

This module is intended to be included in the :file:`NEED` file to
define a full set of generators for transcorrelated calculations.
