.. _diagonalize_tc_h: 
 
.. program:: diagonalize_tc_h 
 
================ 
diagonalize_tc_h 
================ 
 
 
 
 
  
 TODO : Reads psi_det in the EZFIO folder and prints out the left- and right-eigenvectors together 
        with the energy. Saves the left-right wave functions at the end. 
  
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`my_extra_grid_becke` 
    * :c:data:`my_grid_becke` 
    * :c:data:`my_n_pt_a_extra_grid` 
    * :c:data:`my_n_pt_a_grid` 
    * :c:data:`my_n_pt_r_extra_grid` 
    * :c:data:`my_n_pt_r_grid` 
    * :c:data:`n_det` 
    * :c:data:`n_int` 
    * :c:data:`n_states` 
    * :c:data:`pruning` 
    * :c:data:`read_wf` 
    * :c:data:`tc_grid1_a` 
    * :c:data:`tc_grid1_r` 
    * :c:data:`tc_grid2_a` 
    * :c:data:`tc_grid2_r` 
    * :c:data:`tc_integ_type` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`routine_diag` 
    * :c:func:`write_int` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`my_extra_grid_becke` 
    * :c:data:`my_grid_becke` 
    * :c:data:`my_n_pt_a_extra_grid` 
    * :c:data:`my_n_pt_a_grid` 
    * :c:data:`my_n_pt_r_extra_grid` 
    * :c:data:`my_n_pt_r_grid` 
    * :c:data:`pruning` 
    * :c:data:`read_wf` 
