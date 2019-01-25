======================
mo_two_e_erf_ints
======================

Here, all two-electron integrals (:math:`erf({\mu}_{erf} * r_{12})/r_{12}`) are computed.
As they have 4 indices and many are zero, they are stored in a map, as defined
in :file:`Utils/map_module.f90`.

The range separation parameter :math:`{\mu}_{erf}` is the variable :option:`ao_two_e_erf_ints mu_erf`.

To fetch an |MO| integral, use
`get_mo_two_e_integral_erf(i,j,k,l,mo_integrals_map_erf)`

The conventions are:

* For |MO| integrals : <ij|kl> = <12|12>

Be aware that it might not be the same conventions for |MO| and |AO| integrals.


