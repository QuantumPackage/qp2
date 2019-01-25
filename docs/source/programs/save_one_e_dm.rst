.. _save_one_e_dm: 
 
.. program:: save_one_e_dm 
 
============= 
save_one_e_dm 
============= 
 
 
 Then, the global variable data_one_e_dm_alpha_mo and data_one_e_dm_beta_mo will automatically read this density in a further calculation. 
  
 This can be used to perform damping on the density in RS-DFT calculation (see the density_for_dft module). 
 
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
