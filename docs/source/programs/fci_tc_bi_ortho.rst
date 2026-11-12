.. _fci_tc_bi_ortho: 
 
.. program:: fci_tc_bi_ortho 
 
=============== 
fci_tc_bi_ortho 
=============== 
 
 
 
 
 Selected Full Configuration Interaction with stochastic selection 
 and PT2. 
  
 This program performs a |CIPSI|-like selected |CI| using a 
 stochastic scheme for both the selection of the important Slater 
 determinants and the computation of the |PT2| correction. This 
 |CIPSI|-like algorithm will be performed for the lowest states of 
 the variational space (see :option:`determinants n_states`). The 
 |FCI| program will stop when reaching at least one the two following 
 conditions: 
  
 * number of Slater determinants > :option:`determinants n_det_max` 
 * abs(|PT2|) less than :option:`perturbation pt2_max` 
  
 The following other options can be of interest: 
  
 :option:`determinants read_wf` 
   When set to |false|, the program starts with a ROHF-like Slater 
   determinant as a guess wave function. When set to |true|, the 
   program starts with the wave function(s) stored in the |EZFIO| 
   directory as guess wave function(s). 
  
 :option:`determinants s2_eig` 
   When set to |true|, the selection will systematically add all the 
   necessary Slater determinants in order to have a pure spin wave 
   function with an |S^2| value corresponding to 
   :option:`determinants expected_s2`. 
  
 For excited states calculations, it is recommended to start with 
 :ref:`cis` or :ref:`cisd` guess wave functions, eventually in 
 a restricted set of |MOs|, and to set :option:`determinants s2_eig` 
 to |true|. 
  
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`my_extra_grid_becke` 
    * :c:data:`my_grid_becke` 
    * :c:data:`my_n_pt_a_extra_grid` 
    * :c:data:`my_n_pt_a_grid` 
    * :c:data:`my_n_pt_r_extra_grid` 
    * :c:data:`my_n_pt_r_grid` 
    * :c:data:`pruning` 
    * :c:data:`tc_grid1_a` 
    * :c:data:`tc_grid1_r` 
    * :c:data:`tc_grid2_a` 
    * :c:data:`tc_grid2_r` 
    * :c:data:`tc_integ_type` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`run_cipsi_tc` 
    * :c:func:`write_int` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`psi_configuration` 
    * :c:data:`n_det` 
    * :c:data:`c0_weight` 
    * :c:data:`eigval_right_tc_bi_orth` 
    * :c:data:`eigval_right_tc_bi_orth` 
    * :c:data:`eigval_right_tc_bi_orth` 
    * :c:data:`my_extra_grid_becke` 
    * :c:data:`my_grid_becke` 
    * :c:data:`my_n_pt_a_extra_grid` 
    * :c:data:`my_n_pt_a_grid` 
    * :c:data:`my_n_pt_r_extra_grid` 
    * :c:data:`my_n_pt_r_grid` 
    * :c:data:`eigval_right_tc_bi_orth` 
    * :c:data:`pruning` 
    * :c:data:`psi_coef` 
    * :c:data:`psi_det_sorted_bit` 
    * :c:data:`psi_det_sorted_bit` 
    * :c:data:`psi_configuration` 
    * :c:data:`psi_det` 
    * :c:data:`psi_det_size` 
    * :c:data:`psi_det_sorted_bit` 
    * :c:data:`psi_det_sorted_bit_tc` 
    * :c:data:`psi_energy` 
    * :c:data:`psi_l_coef_bi_ortho` 
    * :c:data:`psi_det_sorted_bit_tc` 
    * :c:data:`psi_r_coef_bi_ortho` 
    * :c:data:`psi_det_sorted_bit_tc` 
    * :c:data:`psi_energy` 
    * :c:data:`pt2_match_weight` 
    * :c:data:`pt2_overlap` 
    * :c:data:`pt2_stoch_istate` 
    * :c:data:`eigval_right_tc_bi_orth` 
    * :c:data:`selection_weight` 
    * :c:data:`state_average_weight` 
    * :c:data:`threshold_davidson_pt2` 
    * :c:data:`threshold_generators` 
    * :c:data:`variance_match_weight` 
