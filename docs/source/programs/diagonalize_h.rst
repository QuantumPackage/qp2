.. _diagonalize_h: 
 
.. program:: diagonalize_h 
 
============= 
diagonalize_h 
============= 
 
 
 
 
 Program that extracts the :option:`determinants n_states` lowest 
 states of the Hamiltonian within the set of Slater determinants stored 
 in the |EZFIO| directory. 
  
 If :option:`determinants s2_eig` = |true|, it will retain only states 
 which correspond to the desired value of 
 :option:`determinants expected_s2`. 
  
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`routine` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
