.. _symmetry_det: 
 
.. program:: symmetry_det 
 
============ 
symmetry_det 
============ 
 
 
 
 
 Determinant and state symmetry driver. 
  
 Reads the wave function stored in the EZFIO directory given on the 
 command line, labels every determinant with the irrep of the 
 MO projector, and projects every CI state on the irreps. 
  
 Prints the determinant table, the per-irrep counts, the per-state 
 weight vectors, and runs the gates D1-D6. Exits 1 if a gate fails, 0 
 otherwise, so that it can be chained in a test script. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`pg_det_report` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
