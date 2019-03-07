.. _print_ci_vectors: 
 
.. program:: print_ci_vectors 
 
================ 
print_ci_vectors 
================ 
 
 
 
 
 Print the ground state wave function stored in the |EZFIO| directory 
 in the intermediate normalization. 
  
 It also prints a lot of information regarding the excitation 
 operators from the reference determinant ! and a first-order 
 perturbative analysis of the wave function. 
  
 If the wave function strongly deviates from the first-order analysis, 
 something funny is going on :) 
 
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
