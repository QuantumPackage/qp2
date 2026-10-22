.. _module_zmq: 
 
.. program:: zmq 
 
.. default-role:: option 
 
===
zmq
===

|ZeroMQ| sockets and messages, used to parallelize a |QP| calculation
over several processes.

The master process runs the calculation. The slave processes, started
with `qp_run -slave`, connect to the master through the address
:c:data:`qp_run_address`. The context of the |ZeroMQ| library is
:c:data:`zmq_context`, and the socket addresses and ports are provided
by :c:data:`zmq_socket_pull_tcp_address` and the related providers.
Sockets of the pair, pull, push and sub types are created by
:c:func:`new_zmq_pair_socket` and the related subroutines, and closed
by the `end_zmq_*` subroutines.

Two levels of parallelism are provided:

* parallel jobs, started with :c:func:`new_parallel_job` and ended
  with :c:func:`end_parallel_job`: during the job, the slave tasks
  execute a designated subroutine, and
  :c:func:`zmq_set_running` marks the job as running. The current job
  is identified by the state :c:data:`zmq_state`
* a task server, for a finer granularity of the work: the workers
  connect with :c:func:`connect_to_taskserver`, which returns their
  worker id, and then fetch and acknowledge tasks through the
  `get_task_from_taskserver` and `task_done_to_taskserver` functions

Data is exchanged between the master and the workers with the
`zmq_put_*` and `zmq_get_*` functions, which move vectors, matrices
and scalars through the qp_run socket: for example
:c:func:`zmq_put_dvector` and :c:func:`zmq_get_dvector` for double
precision vectors.
 
 
 
Providers 
--------- 
 
.. c:var:: is_zmq_slave


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        logical	:: is_zmq_slave	


    If |true|, the current process is a |ZeroMQ| slave.


 
.. c:var:: qp_run_address


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        character*(128)	:: qp_run_address	
        integer	:: zmq_port_start	


    Address of the qp_run socket
    Example : tcp://130.120.229.139:12345

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_socket_pull_tcp_address`

 
.. c:var:: zmq_context


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer(ZMQ_PTR)	:: zmq_context	
        integer(omp_lock_kind)	:: zmq_lock	


    Context for the ZeroMQ library


 
.. c:var:: zmq_lock


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer(ZMQ_PTR)	:: zmq_context	
        integer(omp_lock_kind)	:: zmq_lock	


    Context for the ZeroMQ library


 
.. c:var:: zmq_port_start


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        character*(128)	:: qp_run_address	
        integer	:: zmq_port_start	


    Address of the qp_run socket
    Example : tcp://130.120.229.139:12345

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_socket_pull_tcp_address`

 
.. c:var:: zmq_socket_pair_inproc_address


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        character*(128)	:: zmq_socket_pull_tcp_address	
        character*(128)	:: zmq_socket_pair_inproc_address	
        character*(128)	:: zmq_socket_push_tcp_address	
        character*(128)	:: zmq_socket_pull_inproc_address	
        character*(128)	:: zmq_socket_push_inproc_address	
        character*(128)	:: zmq_socket_sub_tcp_address	


    Socket which pulls the results (2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`


 
.. c:var:: zmq_socket_pull_inproc_address


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        character*(128)	:: zmq_socket_pull_tcp_address	
        character*(128)	:: zmq_socket_pair_inproc_address	
        character*(128)	:: zmq_socket_push_tcp_address	
        character*(128)	:: zmq_socket_pull_inproc_address	
        character*(128)	:: zmq_socket_push_inproc_address	
        character*(128)	:: zmq_socket_sub_tcp_address	


    Socket which pulls the results (2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`


 
.. c:var:: zmq_socket_pull_tcp_address


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        character*(128)	:: zmq_socket_pull_tcp_address	
        character*(128)	:: zmq_socket_pair_inproc_address	
        character*(128)	:: zmq_socket_push_tcp_address	
        character*(128)	:: zmq_socket_pull_inproc_address	
        character*(128)	:: zmq_socket_push_inproc_address	
        character*(128)	:: zmq_socket_sub_tcp_address	


    Socket which pulls the results (2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`


 
.. c:var:: zmq_socket_push_inproc_address


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        character*(128)	:: zmq_socket_pull_tcp_address	
        character*(128)	:: zmq_socket_pair_inproc_address	
        character*(128)	:: zmq_socket_push_tcp_address	
        character*(128)	:: zmq_socket_pull_inproc_address	
        character*(128)	:: zmq_socket_push_inproc_address	
        character*(128)	:: zmq_socket_sub_tcp_address	


    Socket which pulls the results (2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`


 
.. c:var:: zmq_socket_push_tcp_address


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        character*(128)	:: zmq_socket_pull_tcp_address	
        character*(128)	:: zmq_socket_pair_inproc_address	
        character*(128)	:: zmq_socket_push_tcp_address	
        character*(128)	:: zmq_socket_pull_inproc_address	
        character*(128)	:: zmq_socket_push_inproc_address	
        character*(128)	:: zmq_socket_sub_tcp_address	


    Socket which pulls the results (2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`


 
.. c:var:: zmq_socket_sub_tcp_address


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        character*(128)	:: zmq_socket_pull_tcp_address	
        character*(128)	:: zmq_socket_pair_inproc_address	
        character*(128)	:: zmq_socket_push_tcp_address	
        character*(128)	:: zmq_socket_pull_inproc_address	
        character*(128)	:: zmq_socket_push_inproc_address	
        character*(128)	:: zmq_socket_sub_tcp_address	


    Socket which pulls the results (2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`


 
.. c:var:: zmq_state


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        character*(128)	:: zmq_state	


    Threads executing work through the ZeroMQ interface


 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: add_task_to_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function add_task_to_taskserver(zmq_to_qp_run_socket,task)


    Get a task from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: connect_to_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread)


    Connect to the task server and obtain the worker ID

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`sscanf_ssds`

 
.. c:macro:: disconnect_from_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function disconnect_from_taskserver(zmq_to_qp_run_socket, worker_id)


    Disconnect from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: disconnect_from_taskserver_state:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function disconnect_from_taskserver_state(zmq_to_qp_run_socket, worker_id, state)


    Disconnect from the task server

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`sscanf_ss`

 
.. c:macro:: end_parallel_job:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine end_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,name_in)


    End a new parallel job with name 'name'. The slave tasks execute subroutine 'slave'

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_context`
       * :c:data:`zmq_state`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_zmq`
       * :c:func:`h_u_0_nstates_zmq`
       * :c:func:`zmq_pt2`
       * :c:func:`zmq_selection`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`end_zmq_pull_socket`
       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`lowercase`
       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`
       * :c:func:`sleep`

 
.. c:macro:: end_zmq_pair_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine end_zmq_pair_socket(zmq_socket_pair)


    Terminate socket on which the results are sent.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_context`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:macro:: end_zmq_pull_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine end_zmq_pull_socket(zmq_socket_pull)


    Terminate socket on which the results are sent.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_context`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`end_parallel_job`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:macro:: end_zmq_push_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine end_zmq_push_socket(zmq_socket_push,thread)


    Terminate socket on which the results are sent.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_context`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_csf_run_slave`
       * :c:func:`davidson_nos2_run_slave`
       * :c:func:`davidson_run_slave`
       * :c:func:`run_pt2_slave_large`
       * :c:func:`run_pt2_slave_small`
       * :c:func:`run_selection_slave`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:macro:: end_zmq_sub_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine end_zmq_sub_socket(zmq_socket_sub)


    Terminate socket on which the results are sent.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_context`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`wait_for_next_state`
       * :c:func:`wait_for_state`
       * :c:func:`wait_for_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:macro:: end_zmq_to_qp_run_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)


    Terminate the socket from the application to qp_run

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_csf_run_slave`
       * :c:func:`davidson_nos2_run_slave`
       * :c:func:`davidson_run_slave`
       * :c:func:`end_parallel_job`
       * :c:func:`pt2_collector`
       * :c:func:`run_pt2_slave_large`
       * :c:func:`run_pt2_slave_small`
       * :c:func:`run_selection_slave`
       * :c:func:`selection_collector`

 
.. c:macro:: get_task_from_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function get_task_from_taskserver(zmq_to_qp_run_socket,worker_id,task_id,task)


    Get a task from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`sscanf_sd`

 
.. c:macro:: get_tasks_from_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function get_tasks_from_taskserver(zmq_to_qp_run_socket,worker_id,task_id,task,n_tasks)


    Get multiple tasks from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: new_parallel_job:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,name_in)


    Start a new parallel job with name 'name'. The slave tasks execute subroutine 'slave'

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_context`
       * :c:data:`zmq_socket_pull_tcp_address`
       * :c:data:`zmq_state`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_zmq`
       * :c:func:`h_u_0_nstates_zmq`
       * :c:func:`zmq_pt2`
       * :c:func:`zmq_selection`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lowercase`
       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:macro:: new_zmq_pair_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function new_zmq_pair_socket(bind)


    Socket on which the collector and the main communicate

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_context`
       * :c:data:`zmq_socket_pull_tcp_address`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:macro:: new_zmq_pull_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function new_zmq_pull_socket()


    Socket on which the results are sent. If thread is 1, use inproc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`
       * :c:data:`zmq_context`
       * :c:data:`zmq_socket_pull_tcp_address`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`
       * :c:func:`sleep`

 
.. c:macro:: new_zmq_push_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function new_zmq_push_socket(thread)


    Socket on which the results are sent. If thread is 1, use inproc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_context`
       * :c:data:`zmq_socket_pull_tcp_address`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:macro:: new_zmq_sub_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function new_zmq_sub_socket()


    Socket to read the state published by the Task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_context`
       * :c:data:`zmq_socket_pull_tcp_address`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:macro:: new_zmq_to_qp_run_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function new_zmq_to_qp_run_socket()


    Socket on which the qp_run process replies

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`
       * :c:data:`zmq_context`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:macro:: reset_zmq_addresses:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine reset_zmq_addresses


    Socket which pulls the results (2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`
       * :c:data:`zmq_socket_pull_tcp_address`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`switch_qp_run_to_master`

 
.. c:macro:: switch_qp_run_to_master:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine switch_qp_run_to_master


    Address of the master qp_run socket
    Example : tcp://130.120.229.139:12345

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`is_zmq_slave`
       * :c:data:`qp_run_address`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_slave_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`getenv`
       * :c:func:`reset_zmq_addresses`

 
.. c:macro:: task_done_to_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function task_done_to_taskserver(zmq_to_qp_run_socket, worker_id, task_id)


    Get a task from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: tasks_done_to_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function tasks_done_to_taskserver(zmq_to_qp_run_socket, worker_id, task_id, n_tasks)


    Get a task from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: wait_for_next_state:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine wait_for_next_state(state)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`end_zmq_sub_socket`

 
.. c:macro:: wait_for_state:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine wait_for_state(state_wait,state)


    Wait for the ZMQ state to be ready

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`end_zmq_sub_socket`

 
.. c:macro:: wait_for_states:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine wait_for_states(state_wait,state,n)


    Wait for the ZMQ state to be ready

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_slave_main`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`end_zmq_sub_socket`

 
.. c:macro:: zmq_abort:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_abort(zmq_to_qp_run_socket)


    Aborts a running parallel computation

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`sleep`

 
.. c:macro:: zmq_delete_task:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_delete_task(zmq_to_qp_run_socket,zmq_socket_pull,task_id,more)


    When a task is done, it has to be removed from the list of tasks on the qp_run
    queue. This guarantees that the results have been received in the pull.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_delete_task_async_recv:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_delete_task_async_recv(zmq_to_qp_run_socket,more,sending)


    When a task is done, it has to be removed from the list of tasks on the qp_run
    queue. This guarantees that the results have been received in the pull.

 
.. c:macro:: zmq_delete_task_async_send:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_delete_task_async_send(zmq_to_qp_run_socket,task_id,sending)


    When a task is done, it has to be removed from the list of tasks on the qp_run
    queue. This guarantees that the results have been received in the pull.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_delete_tasks:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_delete_tasks(zmq_to_qp_run_socket,zmq_socket_pull,task_id,n_tasks,more)


    When a task is done, it has to be removed from the list of tasks on the qp_run
    queue. This guarantees that the results have been received in the pull.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_delete_tasks_async_recv:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_delete_tasks_async_recv(zmq_to_qp_run_socket,more,sending)


    When a task is done, it has to be removed from the list of tasks on the qp_run
    queue. This guarantees that the results have been received in the pull.

 
.. c:macro:: zmq_delete_tasks_async_send:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_delete_tasks_async_send(zmq_to_qp_run_socket,task_id,n_tasks,sending)


    When a task is done, it has to be removed from the list of tasks on the qp_run
    queue. This guarantees that the results have been received in the pull.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_get8_dvector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get8_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Get a float vector from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`zmq_state`

 
.. c:macro:: zmq_get8_ivector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get8_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Get a vector of integers from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`zmq_state`

 
.. c:macro:: zmq_get_dmatrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_dmatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Get a float vector from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`zmq_state`

 
.. c:macro:: zmq_get_dvector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Get a float vector from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`zmq_state`

 
.. c:macro:: zmq_get_i8matrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_i8matrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Get a float vector from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`zmq_state`

 
.. c:macro:: zmq_get_imatrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_imatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Get a float vector from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`zmq_state`

 
.. c:macro:: zmq_get_int:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_int(zmq_to_qp_run_socket, worker_id, name, x)


    Get a vector of integers from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`zmq_state`

 
.. c:macro:: zmq_get_int_nompi:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_int_nompi(zmq_to_qp_run_socket, worker_id, name, x)


    Get a vector of integers from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_get_ivector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Get a vector of integers from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`zmq_state`

 
.. c:macro:: zmq_port:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function zmq_port(ishift)


    Return the value of the ZMQ port from the corresponding integer

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`

 
.. c:macro:: zmq_put8_dvector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put8_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Put a float vector on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_put8_ivector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put8_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Put a vector of integers on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_put_dmatrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_dmatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Put a float vector on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_put_dvector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Put a float vector on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_put_i8matrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_i8matrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Put a float vector on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_put_imatrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_imatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Put a float vector on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_put_int:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_int(zmq_to_qp_run_socket, worker_id, name, x)


    Put a vector of integers on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_put_ivector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Put a vector of integers on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:macro:: zmq_set_running:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_set_running(zmq_to_qp_run_socket)


    Set the job to Running in QP-run

