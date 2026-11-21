==================
mo_one_e_ints
==================

All the one-electron integrals in the |MO| basis.

They are obtained from the |AO| integrals of the
:ref:`module_ao_one_e_ints` module with the |MO| transformation
(:c:func:`mo_to_ao` and its variants, and
:c:data:`ao_one_e_integrals_from_mo`):

* :c:data:`mo_kinetic_integrals`: the kinetic operator integrals
  (`kin_mo_ints.irp.f`)
* :c:data:`mo_integrals_n_e`: the nuclear-electron integrals
  (`pot_mo_ints.irp.f`)
* :c:data:`mo_one_e_integrals`: the core Hamiltonian integrals
  (`mo_one_e_ints.irp.f`)
* :c:data:`mo_pseudo_integrals`, with the local and non local parts
  (`pot_mo_pseudo_ints.irp.f`)

Other integrals related to the position operator, dipole moments and
spreads, are in `spread_dipole_mo.irp.f`, the orthonormalization of
the |MOs| in `orthonormalize.irp.f`, and the reading and writing of
the integrals in the |EZFIO| database are provided by the
`read_mo_integrals_*` and `write_mo_integrals_*` subroutines.
