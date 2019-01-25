==================
mo_two_e_ints
==================

Here, all two-electron integrals (:math:`1/r_{12}`) are computed.
As they have 4 indices and many are zero, they are stored in a map, as defined
in :file:`Utils/map_module.f90`.

To fetch an |AO| integral, use the
`get_ao_two_e_integral(i,j,k,l,ao_integrals_map)` function, and
to fetch an |MO| integral, use
`get_two_e_integral(i,j,k,l,mo_integrals_map)` or
`mo_two_e_integral(i,j,k,l)`.

The conventions are:

* For |AO| integrals : (ik|jl) = (11|22)
* For |MO| integrals : <ij|kl> = <12|12>



