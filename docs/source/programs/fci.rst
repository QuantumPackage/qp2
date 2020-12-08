.. _fci: 
 
.. program:: fci 
 
=== 
fci 
=== 
 
 
 
 
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
 
    * :c:data:`do_pt2` 
    * :c:data:`is_zmq_slave` 
    * :c:data:`mo_two_e_integrals_in_map` 
    * :c:data:`psi_coef` 
    * :c:data:`psi_det` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`run_cipsi` 
    * :c:func:`run_slave_cipsi` 
    * :c:func:`run_stochastic_cipsi` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`ci_electronic_energy` 
    * :c:data:`ci_electronic_energy` 
    * :c:data:`ci_energy` 
    * :c:data:`ci_electronic_energy` 
    * :c:data:`n_det` 
    * :c:data:`n_iter` 
    * :c:data:`psi_occ_pattern` 
    * :c:data:`c0_weight` 
    * :c:data:`distributed_davidson` 
    * :c:data:`psi_coef` 
    * :c:data:`psi_det_sorted_bit` 
    * :c:data:`psi_det` 
    * :c:data:`psi_det_size` 
    * :c:data:`psi_det_sorted_bit` 
    * :c:data:`psi_energy` 
    * :c:data:`psi_occ_pattern` 
    * :c:data:`psi_energy` 
    * :c:data:`pt2_e0_denominator` 
    * :c:data:`pt2_match_weight` 
    * :c:data:`pt2_overlap` 
    * :c:data:`pt2_stoch_istate` 
    * :c:data:`read_wf` 
    * :c:data:`selection_weight` 
    * :c:data:`state_average_weight` 
    * :c:data:`threshold_davidson_pt2` 
    * :c:data:`threshold_generators` 
    * :c:data:`variance_match_weight` 
