===========
cipsi_utils
===========

Shared routines of the |CIPSI| selection and of the stochastic |PT2|
computation.

* the selection buffer: :c:func:`create_selection_buffer`,
  :c:func:`add_to_selection_buffer`, :c:func:`merge_selection_buffers`
  and :c:func:`sort_selection_buffer` manage the determinants
  selected by the slaves, and :c:func:`push_selection_results` /
  :c:func:`pull_selection_results` exchange them with the master
* :c:func:`run_slave_cipsi` and :c:func:`run_slave_main` are the main
  programs of the selection slaves
* the |PT2| bookkeeping: :c:func:`pt2_alloc`, :c:func:`pt2_add`,
  :c:func:`pt2_add2` and :c:func:`pt2_dealloc` accumulate the energy,
  variance, |rPT2| and overlap of each state (:c:func:`pt2_serialize`
  and :c:func:`pt2_deserialize` carry them to the selection slaves),
  and :c:data:`pt2_mindetinfirstteeth`, :c:data:`pt2_n_tasks`,
  :c:data:`pt2_n_tasks_max` and :c:data:`pt2_r` describe the comb used
  to sample the generators
* the |PT2| driver itself is :c:func:`pt2_stochastic`, in this module:
  it samples the walk positions with OpenMP tasks in the master process,
  and no longer uses a |ZeroMQ| slave
