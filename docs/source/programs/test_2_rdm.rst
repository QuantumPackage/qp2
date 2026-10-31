.. _test_2_rdm: 
 
.. program:: test_2_rdm 
 
========== 
test_2_rdm 
========== 
 
 
 
 
 Test program for the two-electron density matrix (2-RDM) 
 providers of the two_body_rdm module: it sets read_wf so that 
 the wave function is read from the EZFIO file, then calls 
 routine_active_only, routine_full_mos and 
 routine_active_only_trans, which print the two-electron 
 repulsion energy obtained from the 2-RDM providers. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`routine_active_only` 
    * :c:func:`routine_active_only_trans` 
    * :c:func:`routine_full_mos` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
