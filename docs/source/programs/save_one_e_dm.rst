.. _save_one_e_dm: 
 
.. program:: save_one_e_dm 
 
============= 
save_one_e_dm 
============= 
 
 
 
 
 Program that computes the one body density on the |MO| and |AO| basis 
 for $\alpha$ and $\beta$ electrons from the wave function 
 stored in the |EZFIO| directory, and then saves it into the 
 :ref:`module_aux_quantities`. 
  
 Then, the global variable :option:`aux_quantities data_one_e_dm_alpha_mo` 
 and :option:`aux_quantities data_one_e_dm_beta_mo` (and the corresponding for |AO|) 
 will automatically ! read this density in the next calculation. 
 This can be used to perform damping on the density in |RSDFT| calculations (see 
 :ref:`module_density_for_dft`). 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`routine_save_one_e_dm` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
