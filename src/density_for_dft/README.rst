===============
density_for_dft
===============

One-electron density matrices used by the |DFT| calculations.

The density is selected by the :option:`density_for_dft density_for_dft`
parameter:

* `WFT`: the density of the wave function, which may be
  multi-determinant
* `KS`: the single-determinant density of the occupied |MOs|
* `input_density`: the one-body density matrix stored in the |EZFIO|
  database, already in the |MO| basis
* `input_density_ao`: as above, but stored in the |AO| basis and
  transformed to the |MO| basis here
* `damping_rs_dft`: the |WFT| density damped toward the input density
  with the factor :option:`density_for_dft damping_for_rs_dft`
* `state_average_dens`: the state-averaged density

The corresponding providers are :c:data:`one_e_dm_mo_alpha_for_dft`,
:c:data:`one_e_dm_mo_beta_for_dft` and :c:data:`one_e_dm_mo_for_dft`
in the |MO| basis, :c:data:`one_e_dm_alpha_ao_for_dft` and
:c:data:`one_e_dm_beta_ao_for_dft` in the |AO| basis, and
:c:data:`one_e_dm_average_mo_for_dft` (with the alpha and beta
variants) for the averages.

When :option:`density_for_dft no_core_density` is |true|, the elements
of the density matrix involving at least one core orbital are set to
zero; the providers with the `_no_core` suffix carry the result,
normalized to the valence electrons when `normalize_dm` is |true|.
