.. _minimize_tc_angles: 
 
.. program:: minimize_tc_angles 
 
================== 
minimize_tc_angles 
================== 
 
 
 
 
 program that minimizes the angle between left- and right-orbitals when degeneracies are found in the TC-Fock matrix 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`my_grid_becke` 
    * :c:data:`my_n_pt_a_grid` 
    * :c:data:`my_n_pt_r_grid` 
    * :c:data:`tc_grid1_a` 
    * :c:data:`tc_grid1_r` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`minimize_tc_orb_angles` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`mo_l_coef` 
    * :c:data:`mo_r_coef` 
    * :c:data:`my_n_pt_a_grid` 
    * :c:data:`my_n_pt_r_grid` 
