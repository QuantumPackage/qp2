==================
ao_one_e_ints
==================

All the one-electron integrals in the |AO| basis.

The most important providers for usual quantum-chemistry calculations
are:

* :c:data:`ao_overlap`, the overlap matrix, with the
  :c:data:`s_inv` and :c:data:`s_half_inv` factors used to
  orthonormalize the |AOs|
* :c:data:`ao_kinetic_integrals`: the kinetic operator integrals
* :c:data:`ao_integrals_n_e`: the nuclear-electron integrals
* :c:data:`ao_one_e_integrals`: the core Hamiltonian integrals
* :c:data:`ao_pseudo_integrals`, with the local and non local parts
  :c:data:`ao_pseudo_integrals_local` and
  :c:data:`ao_pseudo_integrals_non_local`

The other files of the module provide the conversion between the
Cartesian and the spherical harmonic |AOs| (`ao_spherical.irp.f`),
the integrals of the position operator, dipole moments and spreads
(`spread_dipole_ao.irp.f`), the integrals with point charges
(`pot_pt_charges.irp.f`), the on-the-fly evaluation of the
integrals with the CGCOS scheme (`aos_cgtos.irp.f`,
`one_e_kin_integrals_cgtos.irp.f`,
`one_e_coul_integrals_cgtos.irp.f`), and the reading and writing of
the integrals in the |EZFIO| database
(`read_ao_integrals_*`, `write_ao_integrals_*`).
