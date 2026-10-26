.. _module_cipsi_utils: 
 
.. program:: cipsi_utils 
 
.. default-role:: option 
 
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
 
 
 
Providers 
--------- 
 
.. c:var:: int global_selection_buffer


    File : :file:`cipsi_utils/run_pt2_slave.irp.f`

    .. code:: fortran

        type(selection_buffer)	:: global_selection_buffer	


    Global buffer for the OpenMP selection

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`global_selection_buffer_lock`
       * :c:data:`n_det_generators`
       * :c:data:`n_int`


 
.. c:var:: int global_selection_buffer_lock


    File : :file:`cipsi_utils/run_pt2_slave.irp.f`

    .. code:: fortran

        integer(omp_lock_kind)	:: global_selection_buffer_lock	


    Global buffer for the OpenMP selection

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`global_selection_buffer`

 
.. c:var:: double pt2_cw


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pt2_w	(N_det_generators)
        double precision, allocatable	:: pt2_cw	(0:N_det_generators)
        double precision	:: pt2_w_t	
        double precision	:: pt2_u_0	
        integer, allocatable	:: pt2_n_0	(pt2_N_teeth+1)



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



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`fill_buffer_double`
       * :c:func:`fill_buffer_single`
       * :c:func:`pt2_collector`
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

       * :c:data:`global_selection_buffer`
       * :c:func:`pt2_collector`
       * :c:func:`run_pt2_slave_large`
       * :c:func:`run_pt2_slave_small`
       * :c:func:`run_selection_slave`
       * :c:func:`selection_collector`
       * :c:func:`zmq_pt2`
       * :c:func:`zmq_selection`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`

 
.. c:function:: void delete_selection_buffer(void)


    File : :file:`cipsi_utils/selection_buffer.irp.f`

    .. code:: fortran

        subroutine delete_selection_buffer(b)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`global_selection_buffer`
       * :c:func:`pt2_collector`
       * :c:func:`run_pt2_slave_large`
       * :c:func:`run_pt2_slave_small`
       * :c:func:`run_selection_slave`
       * :c:func:`selection_collector`
       * :c:func:`zmq_pt2`
       * :c:func:`zmq_selection`

 
.. c:function:: void make_selection_buffer_s2(void)


    File : :file:`cipsi_utils/selection_buffer.irp.f`

    .. code:: fortran

        subroutine make_selection_buffer_s2(b)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`zmq_pt2`
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

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_pt2_slave_large`

 
.. c:function:: void provide_everything(void)


    File : :file:`cipsi_utils/slave_cipsi.irp.f`

    .. code:: fortran

        subroutine provide_everything



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

       * :c:func:`pt2_collector`
       * :c:func:`selection_collector`

 
.. c:function:: void pt2_add2(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_add2(p1, w, p2)


    p1 =! p1 +( w * p2**2)

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_collector`

 
.. c:function:: void pt2_alloc(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_alloc(pt2_data,N)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_collector`
       * :c:func:`run_cipsi`
       * :c:func:`run_pt2_slave_large`
       * :c:func:`run_pt2_slave_small`
       * :c:func:`run_selection_slave`
       * :c:func:`run_stochastic_cipsi`
       * :c:func:`selection_collector`

 
.. c:function:: void pt2_collector(void)


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        subroutine pt2_collector(zmq_socket_pull, E, relative_error, pt2_data, pt2_data_err, b, N_)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`n_states`
       * :c:data:`pt2_f`
       * :c:data:`pt2_j`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`pt2_u`
       * :c:data:`pt2_w`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`zmq_pt2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`add_to_selection_buffer`
       * :c:func:`check_mem`
       * :c:func:`create_selection_buffer`
       * :c:func:`delete_selection_buffer`
       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`format_w_error`
       * :c:func:`pt2_add`
       * :c:func:`pt2_add2`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`pull_pt2_results`
       * :c:func:`sleep`
       * :c:func:`sort_selection_buffer`
       * :c:func:`wall_time`

 
.. c:function:: void pt2_dealloc(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_dealloc(pt2_data)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_collector`
       * :c:func:`run_cipsi`
       * :c:func:`run_pt2_slave_large`
       * :c:func:`run_pt2_slave_small`
       * :c:func:`run_selection_slave`
       * :c:func:`run_stochastic_cipsi`
       * :c:func:`selection_collector`

 
.. c:function:: void pt2_deserialize(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_deserialize(pt2_data, n, x)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pull_pt2_results`
       * :c:func:`pull_selection_results`

 
.. c:function:: void pt2_find_sample(void)


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer function pt2_find_sample(v, w)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`

 
.. c:function:: void pt2_find_sample_lr(void)


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        integer function pt2_find_sample_lr(v, w, l_in, r_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`

 
.. c:function:: void pt2_serialize(void)


    File : :file:`cipsi_utils/pt2_type.irp.f`

    .. code:: fortran

        subroutine pt2_serialize(pt2_data, n, x)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`push_pt2_results_async_send`
       * :c:func:`push_selection_results`

 
.. c:function:: void pt2_slave_inproc(void)


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        subroutine pt2_slave_inproc(i)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`global_selection_buffer`
       * :c:data:`pt2_e0_denominator`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`zmq_pt2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`run_pt2_slave`

 
.. c:function:: void pull_pt2_results(void)


    File : :file:`cipsi_utils/run_pt2_slave.irp.f`

    .. code:: fortran

        subroutine pull_pt2_results(zmq_socket_pull, index, pt2_data, task_id, n_tasks, b)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_collector`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_deserialize`

 
.. c:function:: void pull_selection_results(void)


    File : :file:`cipsi_utils/run_selection_slave.irp.f`

    .. code:: fortran

        subroutine pull_selection_results(zmq_socket_pull, pt2_data, val, det, N, task_id, ntasks)



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

 
.. c:function:: void push_pt2_results(void)


    File : :file:`cipsi_utils/run_pt2_slave.irp.f`

    .. code:: fortran

        subroutine push_pt2_results(zmq_socket_push, index, pt2_data, b, task_id, n_tasks)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_pt2_slave_small`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`push_pt2_results_async_recv`
       * :c:func:`push_pt2_results_async_send`

 
.. c:function:: void push_pt2_results_async_recv(void)


    File : :file:`cipsi_utils/run_pt2_slave.irp.f`

    .. code:: fortran

        subroutine push_pt2_results_async_recv(zmq_socket_push,mini,sending)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`push_pt2_results`
       * :c:func:`run_pt2_slave_large`

 
.. c:function:: void push_pt2_results_async_send(void)


    File : :file:`cipsi_utils/run_pt2_slave.irp.f`

    .. code:: fortran

        subroutine push_pt2_results_async_send(zmq_socket_push, index, pt2_data, b, task_id, n_tasks, sending)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`push_pt2_results`
       * :c:func:`run_pt2_slave_large`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_serialize`

 
.. c:function:: void push_selection_results(void)


    File : :file:`cipsi_utils/run_selection_slave.irp.f`

    .. code:: fortran

        subroutine push_selection_results(zmq_socket_push, pt2_data, b, task_id, ntasks)



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



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`zmq_pt2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`i8sort`

 
.. c:function:: void run_pt2_slave(void)


    File : :file:`cipsi_utils/run_pt2_slave.irp.f`

    .. code:: fortran

        subroutine run_pt2_slave(thread,iproc,energy)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states_diag`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_slave_inproc`
       * :c:func:`run_slave_main`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`run_pt2_slave_large`

 
.. c:function:: void run_pt2_slave_large(void)


    File : :file:`cipsi_utils/run_pt2_slave.irp.f`

    .. code:: fortran

        subroutine run_pt2_slave_large(thread,iproc,energy)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`global_selection_buffer`
       * :c:data:`global_selection_buffer_lock`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`pt2_f`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_pt2_slave`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_selection_buffer`
       * :c:func:`delete_selection_buffer`
       * :c:func:`end_zmq_push_socket`
       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`merge_selection_buffers`
       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`push_pt2_results_async_recv`
       * :c:func:`push_pt2_results_async_send`
       * :c:func:`select_connected`
       * :c:func:`sleep`
       * :c:func:`sort_selection_buffer`
       * :c:func:`sscanf_ddd`
       * :c:func:`wall_time`

 
.. c:function:: void run_pt2_slave_small(void)


    File : :file:`cipsi_utils/run_pt2_slave.irp.f`

    .. code:: fortran

        subroutine run_pt2_slave_small(thread,iproc,energy)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`pt2_f`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_selection_buffer`
       * :c:func:`delete_selection_buffer`
       * :c:func:`end_zmq_push_socket`
       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`pt2_alloc`
       * :c:func:`pt2_dealloc`
       * :c:func:`push_pt2_results`
       * :c:func:`select_connected`
       * :c:func:`sort_selection_buffer`
       * :c:func:`sscanf_ddd`
       * :c:func:`usleep`
       * :c:func:`wall_time`

 
.. c:function:: void run_selection_slave(void)


    File : :file:`cipsi_utils/run_selection_slave.irp.f`

    .. code:: fortran

        subroutine run_selection_slave(thread,iproc,energy)



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

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`provide_everything`
       * :c:func:`run_slave_main`
       * :c:func:`switch_qp_run_to_master`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`distributed_davidson`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`read_wf`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_generators`

 
.. c:function:: void run_slave_main(void)


    File : :file:`cipsi_utils/slave_cipsi.irp.f`

    .. code:: fortran

        subroutine run_slave_main



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_configuration`
       * :c:data:`elec_alpha_num`
       * :c:data:`global_selection_buffer`
       * :c:data:`h0_type`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`mpi_rank`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nthreads_pt2`
       * :c:data:`pseudo_sym`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_selectors_coef_transp`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_f`
       * :c:data:`pt2_min_parallel_tasks`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`qp_max_mem`
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

       * :c:func:`check_mem`
       * :c:func:`davidson_slave_tcp`
       * :c:func:`mpi_print`
       * :c:func:`provide_for_zmq_pt2`
       * :c:func:`resident_memory`
       * :c:func:`run_pt2_slave`
       * :c:func:`run_selection_slave`
       * :c:func:`usleep`
       * :c:func:`wait_for_states`
       * :c:func:`wall_time`
       * :c:func:`write_double`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`pt2_e0_denominator`
       * :c:data:`pt2_stoch_istate`
       * :c:data:`selection_weight`
       * :c:data:`state_average_weight`
       * :c:data:`threshold_generators`

 
.. c:function:: void selection_collector(void)


    File : :file:`cipsi_utils/zmq_selection.irp.f`

    .. code:: fortran

        subroutine selection_collector(zmq_socket_pull, b, N, pt2_data)



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



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_to_selection_buffer`
       * :c:func:`pt2_collector`
       * :c:func:`run_pt2_slave_large`
       * :c:func:`run_pt2_slave_small`
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

       * :c:func:`zmq_pt2`
       * :c:func:`zmq_selection`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_match_weight`
       * :c:data:`threshold_davidson_pt2`
       * :c:data:`variance_match_weight`

 
.. c:function:: void zmq_pt2(void)


    File : :file:`cipsi_utils/pt2_stoch_routines.irp.f`

    .. code:: fortran

        subroutine ZMQ_pt2(E, pt2_data, pt2_data_err, relative_error, N_in)


    Computes the PT2 energy using ZMQ

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`det_to_configuration`
       * :c:data:`elec_alpha_num`
       * :c:data:`excitation_alpha_max`
       * :c:data:`excitation_beta_max`
       * :c:data:`excitation_max`
       * :c:data:`global_selection_buffer`
       * :c:data:`h0_type`
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
       * :c:data:`psi_selectors`
       * :c:data:`pt2_e0_denominator`
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
       * :c:data:`threshold_generators`

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
       * :c:func:`end_parallel_job`
       * :c:func:`fill_h_apply_buffer_no_selection`
       * :c:func:`make_selection_buffer_s2`
       * :c:func:`new_parallel_job`
       * :c:func:`provide_for_zmq_pt2`
       * :c:func:`pt2_collector`
       * :c:func:`pt2_slave_inproc`
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

       * :c:func:`run_cipsi`
       * :c:func:`zmq_pt2`

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

