============================================
Developping new functionals for KS or RS-DFT
============================================

The very basics
===============

To develop new functionals for |DFT| (or |RSDFT|) to be used in the SCF programs (:ref:`ks_scf` or :ref:`rs_ks_scf` ) or in multi-configurational |RSDFT| (see `<https://gitlab.com/eginer/qp_plugins_eginer>`_), you need to specify, at the end of the day, two things:

* *exchange/correlation* **energy functionals** used to compute the **energy**

* *exchange/correlation* **potentials** for (alpha/beta electrons) used to **optimize the wave function**


The providers to define such quantities, and then used by all the all the |DFT| programs, are (see :ref:`module_dft_one_e`):  

* :c:data:`energy_x` and :c:data:`energy_c` : the *exchange* and *correlation* energy functionals (see :file:`e_xc_genera l.irp.f`)
 
* :c:data:`potential_x_alpha_ao` and :c:data:`potential_x_beta_ao` : the exchange potential for alpha/beta electrons on the |AO| basis (se e :file:`pot_general.irp.f`)
 
* :c:data:`potential_c_alpha_ao` and :c:data:`potential_c_beta_ao` : the correlation potential for alpha/beta electrons on the |AO| basis ( see :file:`pot_general.irp.f`)

From the |AO| basis, the providers for the exchange/correlation potentials of alpha/beta electrons on the |MO| basis are automatically obtained by a |AO| --> |MO| transformation. 

So, at the end of the day, adding a new functional consists only in **setting a new value to these providers**. 


The general philosphy
---------------------

The directory **functionals** contains only files ending with .irp.f whose name being the name of a specific functional. 
All files in *a_functional*.irp.f **must** contain **at least** the following providers

* :c:data:`energy_x_a_functional` and  :c:data:`energy_c_a_functional` which are of course the exchange and correlation energies

* :c:data:`potential_x_alpha_ao_a_functional` and :c:data:`potential_x_beta_ao_a_functional` which are the exchange alpha/beta potentials 

* :c:data:`potential_c_alpha_ao_a_functional` and :c:data:`potential_c_beta_ao_a_functional` which are the correlation alpha/beta potentials 

For instance, the file :file:`sr_lda.irp.f` contains the following providers 

* :c:data:`energy_x_sr_lda` and  :c:data:`energy_c_sr_lda` which are of course the exchange and correlation energies

* :c:data:`potential_x_alpha_ao_sr_lda` and :c:data:`potential_x_beta_ao_sr_lda` which are the exchange alpha/beta potentials 

* :c:data:`potential_c_alpha_ao_sr_lda` and :c:data:`potential_c_beta_ao_sr_lda` which are the correlation alpha/beta potentials 


Therefore, if you want to develop a new functional, just design a provider 

To use a functional 

Using the density for DFT calculations in the |QP|
==================================================

Different ways of defining the density for the DFT
--------------------------------------------------

There are many ways of defining a density, and the keyword to define it is :option:`density_for_dft density_for_dft`. 
Here are the following options for that keyword: 

* "KS" : density is obtained from **a single Slater determinant** 

* "WFT" : density is obtained from **the wave function** which is stored in the |EZFIO| data base

* "input_density" : a one-body density matrix on the |MO| basis is read from the |EZFIO| data base, and the density is built from there (see :c:data:`data_one_e_dm_alpha_mo`)  

* "damping_rs_dft" : damped density between "WFT" and "input_density" with the damping factor :option:`density_for_dft damping_for_rs_dft`. 

 .. note:: 
    If an |MO| basis is already defined in the  |EZFIO| data base, the one-body density matrices will be defined 
    according to this |MO| basis. For instance, if "KS", the density constructed will be density of a single Slater 
    determinant built with the current |MO| basis stored in the |EZFIO| data base. 

Once that you have defined how to define the density, you can easily access to the providers associated to it.  


Value of the density and its gradients in real space
----------------------------------------------------

The density and its gradients evaluated on all grid points are (see :ref:`module_dft_utils_in_r`): 

* :c:data:`one_e_dm_alpha_at_r` and :c:data:`one_e_dm_beta_at_r` : alpha/beta density at grid points 

* :c:data:`one_e_dm_and_grad_alpha_in_r`, :c:data:`one_e_dm_and_grad_beta_in_r`: alpha/beta gradients (and densities)

If you want to evaluate the density and its gradients at a given point in space, please refer to: 

* :c:func:`density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r`

If you use these *providers* and subroutines, the density computed will be coherent with the choice of density that you specified 
with :option:`density_for_dft density_for_dft`, and it will impact automatically the general providers of :ref:`module_dft_one_e`. 


