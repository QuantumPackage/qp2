===============
dft_utils_one_e
===============

This module contains all the one-body related quantities needed to perform DFT or RS-DFT calculations with the LDA and PBE functionals.
Therefore, it contains most of the properties which depends on the one-body density and density matrix.

Some interesting quantities you might take a look at:

* The LDA and PBE *providers* for the x/c energies in :file:`e_xc.irp.f` and :file:`sr_exc.irp.f`
* The LDA and PBE *providers* for the x/c potentials on the AO basis in :file:`pot_ao.irp.f` and  :file:`sr_pot_ao.irp.f`
* The :math:`h_{core}` energy computed directly with the one-body density matrix in :file:`one_e_energy_dft.irp.f`
* LDA and PBE short-range functionals *subroutines* in :file:`exc_sr_lda.irp.f` and :file:`exc_sr_pbe.irp.f`


