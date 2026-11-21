.. _attachement_orb: 
 
.. program:: attachement_orb 
 
=============== 
attachement_orb 
=============== 
 
 
 
 
 Writes the dettachment and attachment orbitals of the second 
 state, eigenvectors of the difference between the spin-summed 
 one-electron density matrices of the first and second states, 
 to the Molden file <ezfio_filename>.attachement.mol. The wave 
 function of the states must be stored in the |EZFIO| directory. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`molden_attachment` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
