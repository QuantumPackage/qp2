.. _module_aux_quantities: 
 
.. program:: aux_quantities 
 
.. default-role:: option 
 
==============
aux_quantities
==============


This module contains some global variables (such as densities and energies)
which are stored in the EZFIO folder in a different place than determinants.
This is used in practice to store density matrices which can be obtained from
any methods, as long as they are stored in the same MO basis which is used for
the calculations. In |RSDFT| calculations, this can be done to perform damping
on the density in order to speed up convergence.

The main providers of that module are:

* `data_one_e_dm_alpha_mo` and `data_one_e_dm_beta_mo` which are the
  one-body alpha and beta densities which are necessary read from the EZFIO
  folder.


Thanks to these providers you can use any density matrix that does not
necessary corresponds to that of the current wave function.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: data_energy_var
 
    Variational energy computed with the wave function
 
 
.. option:: data_energy_proj
 
    Projected energy computed with the wave function
 
 
.. option:: data_one_e_dm_alpha_mo
 
    Alpha one body density matrix on the |MO| basis computed with the wave function
 
 
.. option:: data_one_e_dm_beta_mo
 
    Beta one body density matrix on the |MO| basis computed with the wave function
 
