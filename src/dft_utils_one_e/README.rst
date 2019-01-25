===============
dft_utils_one_e
===============

This module contains all the one-body related quantities needed to perform DFT or RS-DFT calculations.
Therefore, it contains most of the properties which depends on the one-body density and density matrix.

The most important files and variables are:

* The general *providers* for the x/c energies in :file:`e_xc_general.irp.f`
* The general *providers* for the x/c potentials in :file:`pot_general.irp.f`
* The short-range hartree operator and all related quantities in :file:`sr_coulomb.irp.f`

These *providers* will be used in many DFT-related programs, such as :file:`ks_scf.irp.f` or :file:`rs_ks_scf.irp.f`.
It is also needed to compute the effective one-body operator needed in multi-determinant RS-DFT (see plugins by eginer).

Some other interesting quantities:

* The LDA and PBE *providers* for the x/c energies in :file:`e_xc.irp.f` and :file:`sr_exc.irp.f`
* The LDA and PBE *providers* for the x/c potentials on the AO basis in :file:`pot_ao.irp.f` and  :file:`sr_pot_ao.irp.f`
* The :math:`h_{core}` energy computed directly with the one-body density matrix in :file:`one_e_energy_dft.irp.f`
* LDA and PBE short-range functionals *subroutines* in :file:`exc_sr_lda.irp.f` and :file:`exc_sr_pbe.irp.f`


