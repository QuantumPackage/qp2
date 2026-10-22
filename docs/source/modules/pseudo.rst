.. _module_pseudo: 
 
.. program:: pseudo 
 
.. default-role:: option 
 
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
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: nucl_charge_remove
 
    Nuclear charges removed per atom
 
 
.. option:: pseudo_klocmax
 
    Maximum value of k for the local component
 
 
.. option:: pseudo_n_k
 
    Powers of r - 2 in the local component
 
 
.. option:: pseudo_v_k
 
    Coefficients in the local component
 
 
.. option:: pseudo_dz_k
 
    Exponents in the local component
 
 
.. option:: pseudo_lmax
 
    Maximum angular momentum
 
 
.. option:: pseudo_kmax
 
    Maximum number of functions in the non-local component
 
 
.. option:: pseudo_n_kl
 
    Powers of r - 2 in the non-local component
 
 
.. option:: pseudo_v_kl
 
    Coefficients in the non-local component
 
 
.. option:: pseudo_dz_kl
 
    Exponents in the non-local component
 
 
.. option:: do_pseudo
 
    If `True`, pseudo-potentials are used.
 
    Default: False
 
.. option:: pseudo_grid_size
 
    Nb of points of the grid for the QMC interfaces
 
    Default: 1000
 
.. option:: pseudo_grid_rmax
 
    R_max of the QMC grid
 
    Default: 10.0
 
.. option:: ao_pseudo_grid
 
    Grid for the QMC interface
 
 
.. option:: mo_pseudo_grid
 
    Grid for the QMC interface
 
