======
pseudo
======

Parameters of the effective core potentials.

When :option:`pseudo do_pseudo` is |true|, the core electrons are
replaced by pseudo-potentials, described by the |EZFIO| parameters of
this module:

* `nucl_charge_remove`: the nuclear charge removed per atom
* the local component, a sum of exponentials: the number of terms
  `pseudo_klocmax`, the powers of r-2 `pseudo_n_k`, the coefficients
  `pseudo_v_k` and the exponents `pseudo_dz_k`
* the non-local component, per angular momentum up to `pseudo_lmax`:
  up to `pseudo_kmax` functions with powers `pseudo_n_kl`,
  coefficients `pseudo_v_kl` and exponents `pseudo_dz_kl`
* the numerical grid on which the pseudo-potential is evaluated,
  `pseudo_grid_size` points up to `pseudo_grid_rmax`, and the grids
  stored in the |EZFIO| database for the |QMC| interface,
  `ao_pseudo_grid` and `mo_pseudo_grid`
