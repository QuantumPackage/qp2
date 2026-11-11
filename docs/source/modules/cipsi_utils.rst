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
       * :c:data:`psi_det_sorted_gen`
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
    to the deterministic selection (ZMQ_selection), because a
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
       * :c:data:`psi_det_sorted_gen`
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
       * :c:data:`psi_det_sorted_gen`
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
    to the deterministic selection (ZMQ_selection), because a
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
       * :c:data:`psi_det_sorted_gen`
       * :c:data:`pt2_stoch_istate`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_f`
       * :c:data:`pt2_j`
       * :c:data:`pt2_w`

 
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
       * :c:data:`psi_det_sorted_gen`
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
       * :c:data:`psi_det_sorted_gen`
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
       * :c:data:`psi_det_sorted_gen`
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
       * :c:func:`selection_collector`

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

       * :c:func:`pt2_stochastic`
       * :c:func:`run_selection_slave`
       * :c:func:`selection_collector`
       * :c:func:`zmq_selection`

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

       * :c:func:`pt2_stochastic`
       * :c:func:`run_selection_slave`
       * :c:func:`selection_collector`
       * :c:func:`zmq_selection`

 
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

       * :c:func:`pt2_stochastic`
       * :c:func:`zmq_selection`

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

 
.. c:function:: void provide_everything(void)


    File : :file:`cipsi_utils/slave_cipsi.irp.f`

    .. code:: fortran

        subroutine provide_everything


    Collects the IRPF90 entities required by the slave CIPSI
    process: the wave function (psi_det, psi_coef) and its
    generator and selector determinants, the selection and PT2
    parameters, and the ZeroMQ state. Called from run_slave_cipsi
    before the qp_run address is switched to the master.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`ci_energy`
       * :c:data:`generators_bitmask`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`pseudo_sym`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_det`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_selectors`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_generators`
       * :c:data:`zmq_context`
       * :c:data:`zmq_state`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_slave_cipsi`

 
.. c:function:: void pt2_add(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_add(p1, w, p2)


    p1 =! p1 +( w * p2)

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_walk`
       * :c:func:`selection_collector`

 
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

       * :c:func:`pt2_omp_walk`
       * :c:func:`pt2_stochastic`
       * :c:func:`run_cipsi`
       * :c:func:`run_selection_slave`
       * :c:func:`run_stochastic_cipsi`
       * :c:func:`selection_collector`

 
.. c:function:: void pt2_dealloc(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_dealloc(pt2_data)


    Deallocates the pt2, rpt2, variance and overlap arrays of
    a pt2_type.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_walk`
       * :c:func:`pt2_stochastic`
       * :c:func:`run_cipsi`
       * :c:func:`run_selection_slave`
       * :c:func:`run_stochastic_cipsi`
       * :c:func:`selection_collector`

 
.. c:function:: void pt2_deserialize(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_deserialize(pt2_data, n, x)


    Restores a pt2_type from the array x written by
    pt2_serialize: n entries each of pt2, rpt2 and variance,
    then the n x n overlap matrix.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pull_selection_results`

 
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

       * :c:func:`pt2_stochastic`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_omp_sift_down`

 
.. c:function:: void pt2_omp_task(void)


    File : :file:`cipsi_utils/pt2_omp.irp.f`

    .. code:: fortran

        subroutine pt2_omp_task(pos, b_workers, b, pt2_data_I)


    One task = one walk position: evaluate all pt2_F(g) subsets of the
    generator g = pt2_J(pos) into pt2_data_I(g), and merge the local
    selection buffer into the shared one after every subset, under the
    named critical region.  All subsets are accumulated in the same slot
    in the fixed order j = 1..pt2_F(g): the sums are bit-identical for any
    schedule (lemma 2).  The local buffer is empty at the start of every
    subset, so what a subset offers does not depend on which thread ran it
    nor on the previous subsets (lemma 3).

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
                                pt2_data_I, pt2_data_S, pt2_data_S2, pos_warm)


    The walk: the legacy collector's body minus every ZeroMQ/arrival
    branch.
    Runs on a single thread, spawns one task per walk position for a
    window of pt2_omp_window positions, waits for the window with a
    TASKGROUP, then consumes the positions of the window in order.
    
    `E' is the total CI energy of the current state; `pos_warm' is the
    position already computed by the serial warm-up (skipped here).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`n_states`
       * :c:data:`nthreads_pt2`
       * :c:data:`pt2_j`
       * :c:data:`pt2_n_teeth`
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

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`push_selection_results`

 
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
       * :c:data:`pt2_f`
       * :c:data:`pt2_j`
       * :c:data:`pt2_n_teeth`
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

       * :c:func:`run`
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
       * :c:func:`zmq_selection`

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

 
.. c:function:: void pull_selection_results(void)


    File : :file:`cipsi_utils/run_selection_slave.irp.f`

    .. code:: fortran

        subroutine pull_selection_results(zmq_socket_pull, pt2_data, val, det, N, task_id, ntasks)


    Master-side counterpart of push_selection_results: receives the
    partial PT2 data of a slave (deserialized into pt2_data) and
    the content of its selection buffer (N, and val and det when
    N > 0), followed by ntasks and the task_id list.
    When the socket is a REP, replies 'ok'.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`selection_collector`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_deserialize`

 
.. c:function:: void push_selection_results(void)


    File : :file:`cipsi_utils/run_selection_slave.irp.f`

    .. code:: fortran

        subroutine push_selection_results(zmq_socket_push, pt2_data, b, task_id, ntasks)


    Sends to the master, over a ZeroMQ socket, the partial PT2 data
    of the slave (the serialized pt2_data) and the content of the
    selection buffer b (b%cur, and b%val and b%det when b%cur > 0),
    followed by ntasks and the task_id list.
    When the socket is a REQ, waits for the 'ok' reply.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_selection_slave`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_serialize`

 
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

 
.. c:function:: void run_selection_slave(void)


    File : :file:`cipsi_utils/run_selection_slave.irp.f`

    .. code:: fortran

        subroutine run_selection_slave(thread,iproc,energy)


    Slave of the CIPSI determinant selection: for each task received
    from the task server, select_connected accumulates the PT2
    contribution of the excitations of one generator determinant and
    fills the local selection buffer, which is sent to the master
    with push_selection_results after each task.
    thread: 1 uses in-process (inproc) ZeroMQ sockets, 0 uses TCP.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pseudo_sym`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`pt2_f`
       * :c:data:`weight_selection`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_slave_main`
       * :c:func:`selection_slave_inproc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_selection_buffer`
       * :c:func:`delete_selection_buffer`
       * :c:func:`end_zmq_push_socket`
       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`provide_for_selection_slave`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`push_selection_results`
       * :c:func:`select_connected`
       * :c:func:`sort_selection_buffer`
       * :c:func:`sscanf_ddd`
       * :c:func:`usleep`

 
.. c:function:: void run_slave_cipsi(void)


    File : :file:`cipsi_utils/slave_cipsi.irp.f`

    .. code:: fortran

        subroutine run_slave_cipsi


    Helper program for distributed parallelism

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`distributed_davidson`
       * :c:data:`read_wf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`fci`
       * :c:func:`pt2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`provide_everything`
       * :c:func:`run_slave_main`
       * :c:func:`switch_qp_run_to_master`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`distributed_davidson`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`read_wf`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_generators`

 
.. c:function:: void run_slave_main(void)


    File : :file:`cipsi_utils/slave_cipsi.irp.f`

    .. code:: fortran

        subroutine run_slave_main


    Main loop of the slave process: the master rank waits for a
    state change from qp_run and broadcasts it to the MPI
    communicator. It then dispatches the matching step:
    'selection' runs run_selection_slave in every OpenMP thread,
    'davidson' runs davidson_slave_tcp.  The stochastic PT2 has no
    slave state any more: it runs in process on OpenMP tasks (see
    pt2_omp.irp.f).  Stops on 'Stopped'.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`mpi_rank`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_generators`
       * :c:data:`zmq_context`
       * :c:data:`zmq_state`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_slave_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_slave_tcp`
       * :c:func:`mpi_print`
       * :c:func:`run_selection_slave`
       * :c:func:`usleep`
       * :c:func:`wait_for_states`
       * :c:func:`wall_time`
       * :c:func:`write_double`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_generators`

 
.. c:function:: void selection_collector(void)


    File : :file:`cipsi_utils/zmq_selection.irp.f`

    .. code:: fortran

        subroutine selection_collector(zmq_socket_pull, b, N, pt2_data)


    Thread i == 0 of ZMQ_selection: pulls the results pushed by
    the worker threads, accumulates their PT2 data in pt2_data and
    adds their best determinants to b, deleting the done tasks
    from the qp_run task server until it replies that no task is
    left. Sorts b. zmq_socket_pull: pull socket of the workers.
    N: buffer size. pt2_data: accumulated PT2 data (inout).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`zmq_selection`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`add_to_selection_buffer`
       * :c:func:`check_mem`
       * :c:func:`create_selection_buffer`
       * :c:func:`delete_selection_buffer`
       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`pt2_add`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`pull_selection_results`
       * :c:func:`sort_selection_buffer`

 
.. c:function:: void selection_slave_inproc(void)


    File : :file:`cipsi_utils/zmq_selection.irp.f`

    .. code:: fortran

        subroutine selection_slave_inproc(i)


    In-process selection worker of ZMQ_selection: called by every
    OpenMP thread with i > 0, it runs one selection slave with
    in-process ZeroMQ sockets, passing i as the thread number to
    run_selection_slave. Thread i == 0 runs selection_collector.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_e0_denominator`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`zmq_selection`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`run_selection_slave`

 
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
       * :c:func:`run_selection_slave`
       * :c:func:`selection_collector`

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
       * :c:data:`psi_det_sorted_gen`
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

       * :c:func:`pt2_stochastic`
       * :c:func:`zmq_selection`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_match_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`variance_match_weight`

 
.. c:function:: void zmq_selection(void)


    File : :file:`cipsi_utils/zmq_selection.irp.f`

    .. code:: fortran

        subroutine ZMQ_selection(N_in, pt2_data)


    Performs the determinant selection using ZeroMQ

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`do_pt2`
       * :c:data:`elec_alpha_num`
       * :c:data:`excitation_alpha_max`
       * :c:data:`excitation_beta_max`
       * :c:data:`excitation_max`
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
       * :c:data:`psi_selectors`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_f`
       * :c:data:`pt2_overlap`
       * :c:data:`qp_max_mem`
       * :c:data:`s2_eig`
       * :c:data:`selection_weight`
       * :c:data:`seniority_max`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_generators`

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
       * :c:func:`end_parallel_job`
       * :c:func:`fill_h_apply_buffer_no_selection`
       * :c:func:`make_selection_buffer_s2`
       * :c:func:`new_parallel_job`
       * :c:func:`selection_collector`
       * :c:func:`selection_slave_inproc`
       * :c:func:`update_pt2_and_variance_weights`
       * :c:func:`write_double`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_match_weight`
       * :c:data:`pt2_overlap`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`variance_match_weight`

