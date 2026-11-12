===========
cipsi_utils
===========

Shared routines of the |CIPSI| selection and of the stochastic |PT2|
computation.

* the selection buffer: :c:func:`create_selection_buffer`,
  :c:func:`add_to_selection_buffer`, :c:func:`merge_selection_buffers`
  and :c:func:`sort_selection_buffer` manage the determinants offered by
  the selection; :c:func:`omp_selection` is the in-process, fixed-order
  OpenMP driver (one task per generator/subset pair, a shared max-heap
  merged under a critical section), with :c:func:`omp_selection_task`
  as its per-task workhorse
* the |PT2| bookkeeping: :c:func:`pt2_alloc`, :c:func:`pt2_add`,
  :c:func:`pt2_add2` and :c:func:`pt2_dealloc` accumulate the energy,
  variance, |rPT2| and overlap of each state (:c:func:`pt2_serialize`
  and :c:func:`pt2_deserialize` move the bookkeeping to and from the
  transcorrelated plugins), and :c:data:`pt2_mindetinfirstteeth`,
  :c:data:`pt2_n_tasks`, :c:data:`pt2_n_tasks_max` and :c:data:`pt2_r`
  describe the comb used to sample the generators
* the |PT2| driver itself is :c:func:`pt2_stochastic`
  (:file:`pt2_omp.irp.f`): it walks the comb and evaluates each window
  with an OpenMP team, one (position, subset) per task, so the result is
  bit-identical for any thread count (see
  :file:`src/fci/BENCHMARK.md` for the contract)
* :c:func:`dress_find_sample` (:file:`sample_find.irp.f`) finds the
  sample of a dressed state; it moved here from the removed
  :file:`dressing` module
