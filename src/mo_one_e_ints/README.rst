==================
mo_one_e_integrals
==================

All the one-electron integrals in |MO| basis are defined here.

The most important providers for usual quantum-chemistry calculation are:

* `mo_kinetic_integrals` which are the kinetic operator integrals on the |AO| basis (see :file:`kin_mo_ints.irp.f`)
* `mo_integrals_n_e` which are the nuclear-elctron operator integrals on the |AO| basis (see :file:`pot_mo_ints.irp.f`)
* `mo_one_e_integrals` which are the the h_core operator integrals on the |AO| basis (see :file:`mo_mono_ints.irp.f`)

Note that you can find other interesting integrals related to the position operator in :file:`spread_dipole_mo.irp.f`.
