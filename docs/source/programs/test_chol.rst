.. _test_chol: 
 
.. program:: test_chol 
 
========= 
test_chol 
========= 
 
 
 
 
 Standalone test driver for the two-electron integral 
 functions of this file (the Cholesky-based bielec_*_no 
 functions vs. the stored arrays). Sets read_wf so that 
 the wave function is read from the EZFIO file; the test 
 subroutine calls are all commented out. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
