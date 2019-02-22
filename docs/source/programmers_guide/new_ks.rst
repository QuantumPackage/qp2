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

We have created a quite easy way to develop new functionals that is based on 

* the used of **your own external plugins** 

* the module :ref:`module_new_functionals` acting as **hub** 

A pictorial representation of the main dependencies can be seen here: 

 .. image:: /_static/dependencies_func.pdf
    :align: center
    :width: 200px
    :alt: Summary of dependencies


The main idea is the following:
 
 1. Develop *new providers* for the new functionals (energy and potential for exchange and correlation) in some **external plugin**.
 
   * Example:
 
     * In an **external plugin** named **fancy_functionals**, you create *e_c_new_fancy_func* for the energy and *pot_ao_alpha_new_func* for the alpha potential
 
     * If you want to be able to use the |DFT| programs already available in the |QP|, these *providers* must use the providers for the density defined in :ref:`module_density_for_dft` and :ref:`module_dft_utils_in_r` (see below).
 
 
 2. Add the name of your **external plugin** to the :file:`NEED` in order to link your new providers to **new_functionals**
 
   * Example:
 
     * add **fancy_functionals** to the NEED file of **new_functionals**
 
 3. Change the file :file:`e_xc_new_func.irp.f` and :file:`pot_xc_new_func.irp.f` to set the value of your new providers to the providers defined in **new_functionals**
 
   * Example:
 
     * for the exchange/correlation energy
 
 
 .. code:: fortran
  
        BEGIN_PROVIDER[double precision, energy_x_new_functional, (N_states) ]
       &BEGIN_PROVIDER[double precision, energy_c_new_functional, (N_states) ]
        implicit none
        BEGIN_DOC
       ! energy_x_new_functional = define here your functional 
       ! energy_c_new_functional = define here your functional 
        END_DOC
         energy_c_new_functional = e_c_new_fancy_func
         energy_x_new_functional = e_x_new_fancy_func
      
        END_PROVIDER 
 
 
4. Compile at the root of the |QP|
 
   * Example:
 
 
 .. code:: bash
 
     cd ${QP_ROOT}
     ninja 
 
 
5. When you want to execute a program with your new functional, just set the options :option:`dft_keywords exchange_functional`  and :option:`dft_keywords correlation_functional` to "my_functional".


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


