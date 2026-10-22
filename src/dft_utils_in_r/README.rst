==============
dft_utils_in_r
==============

Quantities of the |DFT| calculation evaluated in real space, on the
Becke numerical grid.

The |AOs|, the |MOs|, the density and its derivatives are evaluated
at each point of the grid and stored, as these quantities are used
and re-used by the different |DFT| modules (see the
:ref:`module_becke_numerical_grid` module for the grid).

The main providers for this module are:

* `aos_in_r_array`: values of the |AO| basis on the grid points
  (`ao_in_r.irp.f`)
* `mos_in_r_array`: values of the |MO| basis on the grid points
  (`mo_in_r.irp.f`)
* `one_e_dm_and_grad_alpha_in_r` and
  `one_e_dm_and_grad_beta_in_r`: the spin densities and their
  gradients on the grid points (`dm_in_r.irp.f`)
* `alpha_dens_kin_in_r` and `beta_dens_kin_in_r`: the kinetic
  densities (`kin_dens.irp.f`)
* the dipole moments on the grid (`dipole_becke.irp.f`) and the
  gradient integrals (`ints_grad.irp.f`)
