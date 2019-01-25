==============
dft_utils_in_r
==============

This module contains most of the fundamental quantities (AOs, MOs or density derivatives) evaluated in real-space representation that are needed for the various DFT modules.

As these quantities might be used and re-used, the values at each point of the grid are stored (see ``becke_numerical_grid`` for more information on the grid).

The main providers for this module are:

* `aos_in_r_array`: values of the |AO| basis on the grid point.
* `mos_in_r_array`: values of the |MO| basis on the grid point.
* `one_e_dm_and_grad_alpha_in_r`: values of the density and its gradienst on the grid points.

