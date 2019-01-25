.. _print_e_conv: 
 
.. program:: print_e_conv 
 
============ 
print_e_conv 
============ 
 
 
 
 
 program that prints in a human readable format the convergence of the CIPSI algorithm. 
  
 for all istate, this program produces 
  
 * a file "EZFIO.istate.conv" containing the variational and var+PT2 energies as a function of N_det 
  
 * for istate > 1, a file EZFIO.istate.delta_e.conv containing the energy difference (both var and var+PT2) with the ground state as a function of N_det 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`ezfio_filename` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`routine_e_conv` 
