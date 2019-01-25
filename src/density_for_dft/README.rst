===============
density_for_dft
===============


This module defines the *provider* of the density used for the DFT related calculations.
This definition is done through the keyword :option:`density_for_dft density_for_dft`.
The density can be:

* WFT : the density is computed with a potentially multi determinant wave function (see variables `psi_det` and `psi_det`)# input_density : the density is set to a density previously stored in the |EZFIO| folder (see ``aux_quantities``)
* damping_rs_dft : the density is damped between the input_density and the WFT density, with a damping factor of :option:`density_for_dft damping_for_rs_dft`

