==================
ao_one_e_integrals
==================

All the one-electron integrals in the |AO| basis are here.

The most important providers for usual quantum-chemistry calculation are:

* `ao_kinetic_integral` which are the kinetic operator integrals on the |AO| basis (see :file:`kin_ao_ints.irp.f`)
* `ao_nucl_elec_integral` which are the nuclear-elctron operator integrals on the |AO| basis (see :file:`pot_ao_ints.irp.f`)
* `ao_one_e_integrals` which are the the h_core operator integrals on the |AO| basis (see :file:`ao_mono_ints.irp.f`)


Note that you can find other interesting integrals related to the position operator in :file:`spread_dipole_ao.irp.f`.
