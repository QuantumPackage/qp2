.. _module_cipsi_utils: 
 
.. program:: cipsi_utils 
 
.. default-role:: option 
 
===========
cipsi_utils
===========

Shared routines of the |CIPSI| selection and of the stochastic |PT2|
computation.

* the selection buffer: :c:func:`create_selection_buffer`,
  :c:func:`add_to_selection_buffer`, :c:func:`merge_selection_buffers`
  and :c:func:`sort_selection_buffer` manage the determinants offered by
  the selection; :c:func:`omp_selection` is the in-process OpenMP driver
  (one task per generator/subset pair, a shared max-heap whose merge
  under a critical section is order-independent), with
  :c:func:`omp_selection_task` as its per-task workhorse
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
 
 
 
Providers 
--------- 
 
.. c:var:: double pt2_cw


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pt2_w	(N_det_generators)
        double precision, allocatable	:: pt2_cw	(0:N_det_generators)
        double precision	:: pt2_w_t	
        double precision	:: pt2_u_0	
        integer, allocatable	:: pt2_n_0	(pt2_N_teeth+1)


    Weights of the Monte Carlo comb of the stochastic PT2.
    pt2_w(i): weight of norm-sorted generator i, rescaled so
    that each of the pt2_N_teeth teeth holds the total weight
    pt2_W_T. pt2_cW: cumulative sum of pt2_w. pt2_u_0: weight
    of the first pt2_n_0(1) generators, always included in the
    comb. pt2_n_0: index bounds of the teeth.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`psi_det_generators`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_f`
       * :c:data:`pt2_j`

 
.. c:var:: int pt2_f


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer, allocatable	:: pt2_f	(N_det_generators)
        integer	:: pt2_n_tasks_max	


    pt2_F(i) is the number of Monte Carlo tasks (combs)
    assigned to generator i (i indexes the norm-sorted
    generator list):
    pt2_n_tasks_max*pt2_min_parallel_tasks for the leading
    (largest-norm) generators, pt2_min_parallel_tasks for a
    middle range, and 1 for the rest.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mpi_master`
       * :c:data:`n_core_orb`
       * :c:data:`n_det_generators`
       * :c:data:`pt2_min_parallel_tasks`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`


 
.. c:var:: int pt2_j


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer, allocatable	:: pt2_j	(N_det_generators)
        integer, allocatable	:: pt2_r	(N_det_generators)


    pt2_J contains the list of generators after ordering them according to the
    Monte Carlo sampling.
    
    pt2_R(i) is the number of combs drawn when determinant i is computed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`pt2_n_tasks`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_u`
       * :c:data:`pt2_w`
       * :c:data:`qp_max_mem`


 
.. c:var:: double pt2_match_weight


    File : :file:`cipsi_utils/selection_weight.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pt2_match_weight	(N_states)


    Weights adjusted along the selection to make the PT2 contributions
    of each state coincide.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`selection_weight`

 
.. c:var:: int pt2_mindetinfirstteeth


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer	:: pt2_n_teeth	
        integer	:: pt2_mindetinfirstteeth	


    pt2_N_teeth: number of comb teeth of the stochastic PT2
    sampling.  1 means "no comb": the caller must then fall back
    to the deterministic selection (omp_selection), because a
    one-tooth comb samples one generator per comb and costs as
    much as the exact pass while adding noise.  This happens when
    N_det_generators < 1024, or when no comb can be built (a
    warning is printed in that case).
    Otherwise: the largest N from 100 down to 2 for which
    testTeethBuilding(pt2_minDetInFirstTeeth, N) is true.
    pt2_minDetInFirstTeeth: minimum number of determinants
    required in the first tooth (1 or min(5, N_det_generators)).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det_generators`
       * :c:data:`psi_det_generators`
       * :c:data:`pt2_stoch_istate`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_f`
       * :c:data:`pt2_j`
       * :c:data:`pt2_w`

 
.. c:var:: int pt2_n_0


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pt2_w	(N_det_generators)
        double precision, allocatable	:: pt2_cw	(0:N_det_generators)
        double precision	:: pt2_w_t	
        double precision	:: pt2_u_0	
        integer, allocatable	:: pt2_n_0	(pt2_N_teeth+1)


    Weights of the Monte Carlo comb of the stochastic PT2.
    pt2_w(i): weight of norm-sorted generator i, rescaled so
    that each of the pt2_N_teeth teeth holds the total weight
    pt2_W_T. pt2_cW: cumulative sum of pt2_w. pt2_u_0: weight
    of the first pt2_n_0(1) generators, always included in the
    comb. pt2_n_0: index bounds of the teeth.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`psi_det_generators`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_f`
       * :c:data:`pt2_j`

 
.. c:var:: int pt2_n_tasks


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer	:: pt2_n_tasks	


    Number of parallel tasks for the Monte Carlo

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_j`

 
.. c:var:: int pt2_n_tasks_max


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer, allocatable	:: pt2_f	(N_det_generators)
        integer	:: pt2_n_tasks_max	


    pt2_F(i) is the number of Monte Carlo tasks (combs)
    assigned to generator i (i indexes the norm-sorted
    generator list):
    pt2_n_tasks_max*pt2_min_parallel_tasks for the leading
    (largest-norm) generators, pt2_min_parallel_tasks for a
    middle range, and 1 for the rest.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mpi_master`
       * :c:data:`n_core_orb`
       * :c:data:`n_det_generators`
       * :c:data:`pt2_min_parallel_tasks`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`


 
.. c:var:: int pt2_n_teeth


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer	:: pt2_n_teeth	
        integer	:: pt2_mindetinfirstteeth	


    pt2_N_teeth: number of comb teeth of the stochastic PT2
    sampling.  1 means "no comb": the caller must then fall back
    to the deterministic selection (omp_selection), because a
    one-tooth comb samples one generator per comb and costs as
    much as the exact pass while adding noise.  This happens when
    N_det_generators < 1024, or when no comb can be built (a
    warning is printed in that case).
    Otherwise: the largest N from 100 down to 2 for which
    testTeethBuilding(pt2_minDetInFirstTeeth, N) is true.
    pt2_minDetInFirstTeeth: minimum number of determinants
    required in the first tooth (1 or min(5, N_det_generators)).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det_generators`
       * :c:data:`psi_det_generators`
       * :c:data:`pt2_stoch_istate`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_f`
       * :c:data:`pt2_j`
       * :c:data:`pt2_w`

 
.. c:var:: int pt2_omp_window


    File : :file:`cipsi_utils/pt2_omp.irp.f`

    .. code:: fortran

        integer	:: pt2_omp_window	


    Number of walk positions computed in one parallel region.  THIS IS
    PART OF THE DEFINITION OF THE RESULT: it must never depend on the
    thread count, otherwise the set of positions computed before a stop
    would change with the number of threads and the PT2 would stop being
    reproducible.  It is a provider only so that the memory budget in
    pt2_stochastic can charge the per-window scratch.


 
.. c:var:: double pt2_overlap


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pt2_overlap	(N_states,N_states)


    Overlap between the perturbed wave functions

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`


 
.. c:var:: int pt2_r


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer, allocatable	:: pt2_j	(N_det_generators)
        integer, allocatable	:: pt2_r	(N_det_generators)


    pt2_J contains the list of generators after ordering them according to the
    Monte Carlo sampling.
    
    pt2_R(i) is the number of combs drawn when determinant i is computed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`pt2_n_tasks`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_u`
       * :c:data:`pt2_w`
       * :c:data:`qp_max_mem`


 
.. c:var:: int pt2_stoch_istate


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer	:: pt2_stoch_istate	


    State for stochatsic PT2

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: double pt2_u


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pt2_u	(N_det_generators)


    Uniform random numbers in (0,1), one per Monte Carlo
    sample (comb) of the stochastic PT2. The random generator
    is seeded with (1, 2, ..., m) before being filled.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_j`

 
.. c:var:: double pt2_u_0


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pt2_w	(N_det_generators)
        double precision, allocatable	:: pt2_cw	(0:N_det_generators)
        double precision	:: pt2_w_t	
        double precision	:: pt2_u_0	
        integer, allocatable	:: pt2_n_0	(pt2_N_teeth+1)


    Weights of the Monte Carlo comb of the stochastic PT2.
    pt2_w(i): weight of norm-sorted generator i, rescaled so
    that each of the pt2_N_teeth teeth holds the total weight
    pt2_W_T. pt2_cW: cumulative sum of pt2_w. pt2_u_0: weight
    of the first pt2_n_0(1) generators, always included in the
    comb. pt2_n_0: index bounds of the teeth.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`psi_det_generators`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_f`
       * :c:data:`pt2_j`

 
.. c:var:: double pt2_w


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pt2_w	(N_det_generators)
        double precision, allocatable	:: pt2_cw	(0:N_det_generators)
        double precision	:: pt2_w_t	
        double precision	:: pt2_u_0	
        integer, allocatable	:: pt2_n_0	(pt2_N_teeth+1)


    Weights of the Monte Carlo comb of the stochastic PT2.
    pt2_w(i): weight of norm-sorted generator i, rescaled so
    that each of the pt2_N_teeth teeth holds the total weight
    pt2_W_T. pt2_cW: cumulative sum of pt2_w. pt2_u_0: weight
    of the first pt2_n_0(1) generators, always included in the
    comb. pt2_n_0: index bounds of the teeth.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`psi_det_generators`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_f`
       * :c:data:`pt2_j`

 
.. c:var:: double pt2_w_t


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pt2_w	(N_det_generators)
        double precision, allocatable	:: pt2_cw	(0:N_det_generators)
        double precision	:: pt2_w_t	
        double precision	:: pt2_u_0	
        integer, allocatable	:: pt2_n_0	(pt2_N_teeth+1)


    Weights of the Monte Carlo comb of the stochastic PT2.
    pt2_w(i): weight of norm-sorted generator i, rescaled so
    that each of the pt2_N_teeth teeth holds the total weight
    pt2_W_T. pt2_cW: cumulative sum of pt2_w. pt2_u_0: weight
    of the first pt2_n_0(1) generators, always included in the
    comb. pt2_n_0: index bounds of the teeth.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`psi_det_generators`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_f`
       * :c:data:`pt2_j`

 
.. c:var:: double selection_weight


    File : :file:`cipsi_utils/selection_weight.irp.f`

    .. code:: fortran

        double precision, allocatable	:: selection_weight	(N_states)


    Weights used in the selection criterion

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`c0_weight`
       * :c:data:`n_states`
       * :c:data:`pt2_match_weight`
       * :c:data:`state_average_weight`
       * :c:data:`variance_match_weight`
       * :c:data:`weight_selection`


 
.. c:var:: double variance_match_weight


    File : :file:`cipsi_utils/selection_weight.irp.f`

    .. code:: fortran

        double precision, allocatable	:: variance_match_weight	(N_states)


    Weights adjusted along the selection to make the variances
    of each state coincide.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`selection_weight`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void add_to_selection_buffer(void)


    File : :file:`cipsi_utils/selection_buffer.irp.f`

    .. code:: fortran

        subroutine add_to_selection_buffer(b, det, val)


    Adds the determinant det with weight val to the selection buffer
    b, but only if b%N > 0 and val <= b%mini, where b%mini is the
    current threshold of the buffer.
    If adding det fills the storage, the buffer is sorted by
    increasing val and the element with the largest val is dropped.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`fill_buffer_double`
       * :c:func:`fill_buffer_single`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`sort_selection_buffer`

 
.. c:function:: void create_selection_buffer(void)


    File : :file:`cipsi_utils/selection_buffer.irp.f`

    .. code:: fortran

        subroutine create_selection_buffer(N, size_in, res)


    Allocates the memory for a selection buffer.
    The arrays have dimension size_in and the maximum number of elements is N

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection`
       * :c:func:`pt2_stochastic`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`

 
.. c:function:: void delete_selection_buffer(void)


    File : :file:`cipsi_utils/selection_buffer.irp.f`

    .. code:: fortran

        subroutine delete_selection_buffer(b)


    Deallocates the storage of the selection buffer b (the det and
    val pointers), and resets b%cur, b%mini and b%N to zero.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection`
       * :c:func:`pt2_stochastic`

 
.. c:function:: void dress_find_sample(void)


    File : :file:`cipsi_utils/sample_find.irp.f`

    .. code:: fortran

        integer function dress_find_sample(v, w)


    Binary search on the cumulative weight array
    w(0:N_det_generators), assumed non-decreasing: returns the
    index i such that w(i-1) < v <= w(i), i.e. the generator
    selected by the sample value v. If several consecutive
    cumulative values are equal (zero-width bin), the last index
    of the run of equal values is returned.
    Probably to avoid a division by zero or a zero-width bin
    (author, 2026-10-31); the code shows the returned index is a
    lookup key into the generator arrays (edI(i), pt2_w(i)) and
    the selected weight pt2_w(i) is the divisor of the estimator
    x =! x +( edI(i) * pt2_W_T / pt2_w(i).)
    
    Moved verbatim from the removed dressing module (2026-11-12);
    the name is kept so that the external declaration in
    pt2_stoch_routines stays untouched.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`

 
.. c:function:: void make_selection_buffer_s2(void)


    File : :file:`cipsi_utils/selection_buffer.irp.f`

    .. code:: fortran

        subroutine make_selection_buffer_s2(b)


    Groups the determinants of the selection buffer b by
    configuration (the bitwise exclusive-or and and of the alpha and
    beta masks): duplicates are removed keeping the smallest val,
    and the result is sorted by increasing val.
    The configurations are then expanded back into determinants with
    configuration_to_dets, at most the original b%cur of them.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection`
       * :c:func:`pt2_stochastic`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`configuration_to_dets`
       * :c:func:`configuration_to_dets_size`
       * :c:func:`dsort`
       * :c:func:`i8sort`

 
.. c:function:: void merge_selection_buffers(void)


    File : :file:`cipsi_utils/selection_buffer.irp.f`

    .. code:: fortran

        subroutine merge_selection_buffers(b1, b2)


    Merges the selection buffers b1 and b2 into b2

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: void omp_selection(void)


    File : :file:`cipsi_utils/selection_omp.irp.f`

    .. code:: fortran

        subroutine omp_selection(N_in, pt2_data)


    Performs the deterministic determinant selection.
    
    In-process OpenMP port of the former distributed selection
    (the task server was removed, 2026-11-12); see the header of this
    file for the design and the determinism contract.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`do_pt2`
       * :c:data:`elec_alpha_num`
       * :c:data:`excitation_alpha_max`
       * :c:data:`excitation_beta_max`
       * :c:data:`excitation_max`
       * :c:data:`generators_bitmask`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_del_orb`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_inact_orb`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_virt_orb`
       * :c:data:`nproc`
       * :c:data:`nthreads_pt2`
       * :c:data:`pseudo_sym`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_coef_transp`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_f`
       * :c:data:`pt2_overlap`
       * :c:data:`qp_max_mem`
       * :c:data:`s2_eig`
       * :c:data:`selection_weight`
       * :c:data:`seniority_max`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_stochastic`
       * :c:func:`run_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_selection_buffer`
       * :c:func:`delete_selection_buffer`
       * :c:func:`fill_h_apply_buffer_no_selection`
       * :c:func:`make_selection_buffer_s2`
       * :c:func:`omp_selection_task`
       * :c:func:`pt2_add`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`pt2_omp_sort`
       * :c:func:`update_pt2_and_variance_weights`
       * :c:func:`write_double`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`variance_match_weight`

 
.. c:function:: void omp_selection_task(void)


    File : :file:`cipsi_utils/selection_omp.irp.f`

    .. code:: fortran

        subroutine omp_selection_task(g, j, b_workers, b, pt2_out)


    One task = one (generator g, subset j) pair: evaluate the subset
    into the caller-provided slot `pt2_out' (the window scratch), and
    merge the local selection buffer into the shared buffer `b' under
    the named critical region, exactly as pt2_omp_task merges into the
    shared heap.
    
    The local buffer is reset (cur = 0, mini = 0) at the start and at
    the end of the task: what a task offers must not depend on which
    thread ran it nor on the previous tasks it received, so the
    offered set is a pure function of (g, j) and the shared top-N
    under pt2_omp_lt is independent of the merge order.
    
    `b_workers(it)' is indexed by the OpenMP thread number clamped to
    0 (the serial warm-up runs outside a parallel region, where
    omp_get_thread_num() returns -1): it then uses slot 1, never
    b_workers(0), as in pt2_omp_task.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_f`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_merge`
       * :c:func:`select_connected`

 
.. c:function:: void pt2_add(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_add(p1, w, p2)


    p1 =! p1 +( w * p2)

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection`
       * :c:func:`pt2_omp_walk`

 
.. c:function:: void pt2_add2(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_add2(p1, w, p2)


    p1 =! p1 +( w * p2**2)

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_walk`

 
.. c:function:: void pt2_alloc(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_alloc(pt2_data,N)


    Allocates and zero-initializes the pt2, rpt2, variance and
    overlap arrays of a pt2_type for N states.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection`
       * :c:func:`pt2_omp_walk`
       * :c:func:`pt2_stochastic`
       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

 
.. c:function:: void pt2_dealloc(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_dealloc(pt2_data)


    Deallocates the pt2, rpt2, variance and overlap arrays of
    a pt2_type.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection`
       * :c:func:`pt2_omp_walk`
       * :c:func:`pt2_stochastic`
       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

 
.. c:function:: void pt2_deserialize(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_deserialize(pt2_data, n, x)


    Restores a pt2_type from the array x written by
    pt2_serialize: n entries each of pt2, rpt2 and variance,
    then the n x n overlap matrix.

 
.. c:function:: void pt2_find_sample(void)


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer function pt2_find_sample(v, w)


    Returns the index where the cumulative weight w crosses v,
    by calling pt2_find_sample_lr on the full range
    0..N_det_generators.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`

 
.. c:function:: void pt2_find_sample_lr(void)


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer function pt2_find_sample_lr(v, w, l_in, r_in)


    Binary search on the cumulative weight array w for the
    first index r such that w(r) >= v, within the range
    l_in..r_in. If w is flat there (determinants of zero
    weight), the last index of the flat run is returned.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`

 
.. c:function:: void pt2_omp_lt(void)


    File : :file:`cipsi_utils/pt2_omp.irp.f`

    .. code:: fortran

        logical function pt2_omp_lt(b1, i1, b2, i2)


    Total order used by the selection heap:
      (v1,d1) < (v2,d2)  iff  v1 < v2, or (v1 == v2 and d1 <lex d2),
    where lex compares the (N_int,2) bit patterns integer by integer,
    spin by spin (column-major order of b%det), the first difference
    deciding.  Signed comparison of the integer(8) patterns is a valid
    total order, and that is all the heap needs.  Ties are real
    (symmetry-equivalent determinants carry bit-identical weights), so
    this is what makes the retained set independent of the merge order.
    
    Entry i1 of buffer b1 is compared with entry i2 of buffer b2; the
    heap routines pass the same buffer twice.  (A single-buffer dummy
    cannot express the cross-buffer comparison that pt2_omp_merge needs.)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: void pt2_omp_merge(void)


    File : :file:`cipsi_utils/pt2_omp.irp.f`

    .. code:: fortran

        subroutine pt2_omp_merge(src, dst)


    Merge one local selection buffer into the shared max-heap, keeping
    the dst%N best candidates under the total order of pt2_omp_lt.
    No sort of src is needed and its order is irrelevant (lemma 3).
    
    Invariant on entry and exit: positions 1..dst%cur form a max-heap
    under pt2_omp_lt (the root is the worst kept candidate), and
    dst%mini = dst%val(1) when the buffer is full, 0 otherwise.
    
    (i) A rejected entry is not better than the current worst kept, so it
    cannot belong to the final dst%N smallest.  (ii) Replacing the root
    and sifting down restores the heap.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection_task`
       * :c:func:`pt2_omp_task`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_sift_down`
       * :c:func:`pt2_omp_sift_up`

 
.. c:function:: void pt2_omp_sift_down(void)


    File : :file:`cipsi_utils/pt2_omp.irp.f`

    .. code:: fortran

        subroutine pt2_omp_sift_down(b, i)


    Restore the max-heap invariant of b after the element at position i
    has been replaced (or decreased).  The heap is positions 1..b%cur;
    b%cur is therefore the heap size, which pt2_omp_sort shrinks as it
    grows the sorted suffix.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_merge`
       * :c:func:`pt2_omp_sort`

 
.. c:function:: void pt2_omp_sift_up(void)


    File : :file:`cipsi_utils/pt2_omp.irp.f`

    .. code:: fortran

        subroutine pt2_omp_sift_up(b, i)


    Restore the max-heap invariant of b (positions 1..b%cur) after the
    element at position i has been appended or increased.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_merge`

 
.. c:function:: void pt2_omp_sort(void)


    File : :file:`cipsi_utils/pt2_omp.irp.f`

    .. code:: fortran

        subroutine pt2_omp_sort(b)


    In-place heapsort of the b%cur entries of b by increasing order of
    pt2_omp_lt (ascending weight, ties ordered by determinant).  This is
    the replacement of sort_selection_buffer, which cannot be used on
    this buffer because it reads b%val(b%N) unguarded.
    
    b%cur is used as the live heap size: it is shrunk as the sorted
    suffix grows, and restored to its entry value before returning.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection`
       * :c:func:`pt2_stochastic`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_sift_down`

 
.. c:function:: void pt2_omp_task(void)


    File : :file:`cipsi_utils/pt2_omp.irp.f`

    .. code:: fortran

        subroutine pt2_omp_task(pos, j, b_workers, b, pt2_out)


    One task = one (walk position, subset) pair: evaluate subset `j' of
    the generator g = pt2_J(pos) into the caller-provided accumulator
    `pt2_out', and merge the local selection buffer into the shared one
    under the named critical region.
    
    `pt2_out' is either the generator's own slot (serial warm-up, running
    accumulation, exactly the sequential reference) or a private
    per-(position,subset) accumulator that the consume phase folds into
    the generator's slot in the fixed order 1..pt2_F(g).  The subset
    never sees the accumulator of another subset, and it starts from a
    zeroed accumulator, so its value depends on (g, j) only.
    
    The local buffer is reset at the start and at the end of every
    subset: what a subset offers must not depend on which thread ran it
    nor on the previous subsets (add_to_selection_buffer drops the entry
    with the largest kept weight when its storage fills, so a buffer that
    outlives a subset could discard a candidate the shared heap would
    have kept).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`nthreads_pt2`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_f`
       * :c:data:`pt2_j`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_walk`
       * :c:func:`pt2_stochastic`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_merge`
       * :c:func:`select_connected`

 
.. c:function:: void pt2_omp_walk(void)


    File : :file:`cipsi_utils/pt2_omp.irp.f`

    .. code:: fortran

        subroutine pt2_omp_walk(E, relative_error, pt2_data, pt2_data_err, b, b_workers, &
                                pt2_data_I, pt2_data_S, pt2_data_S2, pos_warm, nproc_target)


    The walk: the legacy collector's body minus every task-transport
    (arrival-order) branch.
    Runs on a single thread, computes a window of pt2_omp_window
    positions in one parallel-do region, then consumes the positions of
    the window in order.
    
    `E' is the total CI energy of the current state; `pos_warm' is the
    position already computed by the serial warm-up (skipped here);
    `nproc_target' is the size of the team, the walking thread included.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`n_states`
       * :c:data:`nthreads_pt2`
       * :c:data:`pt2_f`
       * :c:data:`pt2_j`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_omp_window`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`pt2_u`
       * :c:data:`pt2_w`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_stochastic`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`format_w_error`
       * :c:func:`pt2_add`
       * :c:func:`pt2_add2`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`pt2_omp_task`
       * :c:func:`wall_time`

 
.. c:function:: void pt2_serialize(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_serialize(pt2_data, n, x)


    Packs the first n entries of pt2, rpt2 and variance and
    the n x n overlap matrix of a pt2_type into x, in that
    order (3*n + n*n doubles).

 
.. c:function:: void pt2_stochastic(void)


    File : :file:`cipsi_utils/pt2_omp.irp.f`

    .. code:: fortran

        subroutine pt2_stochastic(E, pt2_data, pt2_data_err, relative_error, N_in)


    Stochastic |PT2| on OpenMP tasks: the only implementation (see the
    header of this file for the algorithm).  Loops over the states,
    allocates the per-state storage, runs the serial warm-up and the
    parallel region, then does the per-state epilogue and the final
    selection epilogue.
    
    `E' is the total CI energy of each state (psi_energy_with_nucl_rep);
    `relative_error' is the target relative stochastic error on PT2;
    `N_in' is the number of determinants to select (0 for qp_run pt2).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`banned_excitation`
       * :c:data:`det_to_configuration`
       * :c:data:`elec_alpha_num`
       * :c:data:`excitation_alpha_max`
       * :c:data:`excitation_beta_max`
       * :c:data:`excitation_max`
       * :c:data:`generators_bitmask`
       * :c:data:`h0_type`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_del_orb`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_inact_orb`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_virt_orb`
       * :c:data:`nproc`
       * :c:data:`nthreads_pt2`
       * :c:data:`pseudo_sym`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_configuration_hii`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_hii`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_coef_transp`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_f`
       * :c:data:`pt2_j`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_omp_window`
       * :c:data:`pt2_overlap`
       * :c:data:`pt2_j`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`pt2_u`
       * :c:data:`pt2_w`
       * :c:data:`qp_max_mem`
       * :c:data:`s2_eig`
       * :c:data:`selection_weight`
       * :c:data:`seniority_max`
       * :c:data:`state_average_weight`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`create_selection_buffer`
       * :c:func:`delete_selection_buffer`
       * :c:func:`fill_h_apply_buffer_no_selection`
       * :c:func:`make_selection_buffer_s2`
       * :c:func:`omp_selection`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`pt2_omp_sort`
       * :c:func:`pt2_omp_task`
       * :c:func:`pt2_omp_walk`
       * :c:func:`remove_duplicates_in_selection_buffer`
       * :c:func:`resident_memory`
       * :c:func:`update_pt2_and_variance_weights`
       * :c:func:`write_double`
       * :c:func:`write_int`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`variance_match_weight`

 
.. c:function:: void remove_duplicates_in_selection_buffer(void)


    File : :file:`cipsi_utils/selection_buffer.irp.f`

    .. code:: fortran

        subroutine remove_duplicates_in_selection_buffer(b)


    Removes duplicate determinants from the selection buffer b:
    determinants are sorted with det_search_key, and within each
    group of identical determinants all but one copy is dropped.
    If duplicates are found, b is compacted in place and b%cur and
    b%N are set to the resulting number of determinants.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_stochastic`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`i8sort`

 
.. c:function:: void sort_selection_buffer(void)


    File : :file:`cipsi_utils/selection_buffer.irp.f`

    .. code:: fortran

        subroutine sort_selection_buffer(b)


    Sorts the b%cur determinants of the selection buffer b by
    increasing val, and keeps the first min(b%N, b%cur) of them.
    b%mini is set to the smallest of its previous value and the
    largest kept val, b%val(b%N).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_to_selection_buffer`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsort`

 
.. c:function:: void testteethbuilding(void)


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        logical function testTeethBuilding(minF, N)


    Checks whether the comb can be built with N teeth so that
    the first tooth holds at least minF determinants. Tries
    leading offsets n0 of the norm-sorted generator list and
    returns .true. if the first minF determinants after n0 fit
    in one tooth of width Wt = (1 - u0)/N (u0: cumulative
    weight before n0); .false. if Wt <= 1e-3 or no offset works.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`psi_det_generators`
       * :c:data:`pt2_stoch_istate`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`

 
.. c:function:: void update_pt2_and_variance_weights(void)


    File : :file:`cipsi_utils/selection_weight.irp.f`

    .. code:: fortran

        subroutine update_pt2_and_variance_weights(pt2_data, N_st)


    Updates the PT2- and Variance- matching weights.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_relative_error`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`variance_match_weight`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`omp_selection`
       * :c:func:`pt2_stochastic`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_match_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`variance_match_weight`

