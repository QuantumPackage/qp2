===========
cipsi_utils
===========

Shared routines of the |CIPSI| selection and of the stochastic |PT2|
computation, for the master process and for the |ZeroMQ| slave
processes.

* the selection buffer: :c:func:`create_selection_buffer`,
  :c:func:`add_to_selection_buffer`, :c:func:`merge_selection_buffers`
  and :c:func:`sort_selection_buffer` manage the determinants
  selected by the slaves, and :c:func:`push_selection_results` /
  :c:func:`pull_selection_results` exchange them with the master
* the stochastic |PT2|: the `pt2_*` routines build the samples and
  accumulate the stochastic estimate of the |PT2| energy on the
  slaves (:c:func:`run_pt2_slave`), and
  :c:func:`push_pt2_results` / :c:func:`pull_pt2_results` exchange
  the partial results with the master
* :c:func:`run_slave_cipsi` and :c:func:`run_slave_main` are the main
  programs of the slave processes
