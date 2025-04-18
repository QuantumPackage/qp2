==================
ao_cart_two_e_ints
==================

Here, all two-electron integrals (:math:`1/r_{12}`) are computed.
As they have 4 indices and many are zero, they are stored in a map, as defined
in :file:`utils/map_module.f90`.

To fetch an |AO| integral, use the
`get_ao_cart_two_e_integral(i,j,k,l,ao_cart_integrals_map)` function.


The conventions are:
* For |AO| integrals : (ij|kl) = (11|22) = <ik|jl> = <12|12>



