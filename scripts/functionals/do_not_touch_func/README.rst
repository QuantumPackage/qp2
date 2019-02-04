new_functionals
===============

This module allows to add new |DFT| or |RSDFT| functionals in the |QP|. 
It works as a **hub** between the *providers* that one can build in some external plugins and the main *providers* 
used by the various |DFT| or |RSDFT| programs. 


.. warning::
     This module is not tracked by GIT and therefore, the modifications performed on this module 
     can be easily lost. Keep in mind this important fact in order not to loose your own work. 


A pictorial representation of the main dependencies can be seen here: 

 .. image:: /_static/dependencies_func.pdf
    :align: center
    :width: 200px
    :alt: Summary of dependencies

The main idea is the following: 

1. Develop *new providers* for the new functionals (energy and potential for exchange and correlation) in some **external plugin**. 
   
  * Example: 

    * In an **external plugin** named **fancy_functionals**, you create *e_c_new_fancy_func* for the energy and *pot_ao_alpha_new_func* for the alpha potential

    * If you want to be able to use the |DFT| programs already available in the |QP|, these *providers* must use the providers for the density defined in :ref:`module_density_for_dft` and :ref:`module_dft_utils_in_r`. 


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

