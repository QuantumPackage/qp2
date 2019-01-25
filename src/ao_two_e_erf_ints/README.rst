======================
ao_two_e_erf_ints
======================

Here, all two-electron integrals (:math:`erf(\mu r_{12})/r_{12}`) are computed.
As they have 4 indices and many are zero, they are stored in a map, as defined
in :file:`utils/map_module.f90`.

The main parameter of this module is :option:`ao_two_e_erf_ints mu_erf` which is the range-separation parameter.

To fetch an |AO| integral, use the
`get_ao_two_e_integral_erf(i,j,k,l,ao_integrals_erf_map)` function.


The conventions are:
* For |AO| integrals : (ij|kl) = (11|22) = <ik|jl> = <12|12>



