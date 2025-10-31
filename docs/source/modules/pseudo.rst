.. _module_pseudo: 
 
.. program:: pseudo 
 
.. default-role:: option 
 
======
pseudo
======

This module defines the |EZFIO| parameters of the effective core potentials.
 
 
 
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
 
