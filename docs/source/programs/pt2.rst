.. _pt2: 
 
.. program:: pt2 
 
=== 
pt2 
=== 
 
 
 
 
 Second order perturbative correction to the wave function contained 
 in the |EZFIO| directory. 
  
 This programs runs the stochastic |PT2| correction on all 
 :option:`determinants n_states` wave functions stored in the |EZFIO| 
 directory. 
  
 The main option for the |PT2| correction is the 
 :option:`perturbation pt2_relative_error` which is the relative 
 stochastic error on the |PT2| to reach before stopping the 
 sampling. 
  
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`is_zmq_slave` 
    * :c:data:`mo_two_e_integrals_in_map` 
    * :c:data:`psi_energy` 
    * :c:data:`read_wf` 
    * :c:data:`threshold_generators` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`run` 
    * :c:func:`run_slave_cipsi` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`distributed_davidson` 
    * :c:data:`pt2_e0_denominator` 
    * :c:data:`pt2_stoch_istate` 
    * :c:data:`read_wf` 
    * :c:data:`selection_weight` 
    * :c:data:`state_average_weight` 
    * :c:data:`threshold_generators` 
