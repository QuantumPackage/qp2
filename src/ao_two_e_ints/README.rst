==================
ao_two_e_ints
==================

All the two-electron integrals in the |AO| basis.

The integrals :math:`(ij|kl)` are computed on the fly or read from
the |EZFIO| database. As they have four indices and many of them are
zero, they are stored in a map, as defined in
:file:`utils/map_module.f90`. To fetch an integral, use the
:c:func:`get_ao_two_e_integral` function with the map
:c:data:`ao_integrals_map`.

The computation of the integrals is controlled by the Schwarz
screening (:c:data:`ao_2e_cgtos_schwartz`) and by the Cholesky
decomposition of the Coulomb kernel (:c:data:`cholesky_ao`,
:c:data:`do_ao_cholesky`, with :c:data:`ao_cholesky_threshold`), and
the range separated (erf) integrals are treated by the
`*_erf` variants.

The conventions are:

* for |AO| integrals: :math:`(ij|kl) = (11|22) = \langle ik|jl
  \rangle = \langle 12|12 \rangle`
