.. _module_zmq: 
 
.. program:: zmq 
 
.. default-role:: option 
 
===
zmq
===

Definition of |ZeroMQ| sockets and messages.


 
 
 
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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`

 
.. c:var:: zmq_lock


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer(ZMQ_PTR)	:: zmq_context	
        integer(omp_lock_kind)	:: zmq_lock	


    Context for the ZeroMQ library

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`

 
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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`

 
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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`

 
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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`

 
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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`

 
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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`

 
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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`

 
.. c:var:: zmq_state


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        character*(128)	:: zmq_state	


    Threads executing work through the ZeroMQ interface

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: add_task_to_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function add_task_to_taskserver(zmq_to_qp_run_socket,task)


    Get a task from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: connect_to_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread)


    Connect to the task server and obtain the worker ID

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: disconnect_from_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function disconnect_from_taskserver(zmq_to_qp_run_socket, worker_id)


    Disconnect from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: disconnect_from_taskserver_state:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function disconnect_from_taskserver_state(zmq_to_qp_run_socket, worker_id, state)


    Disconnect from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: end_parallel_job:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine end_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,name_in)


    End a new parallel job with name 'name'. The slave tasks execute subroutine 'slave'

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`
       * :c:data:`zmq_context`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:func:`h_s2_u_0_nstates_zmq`
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

 
.. c:function:: end_zmq_pair_socket:


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

 
.. c:function:: end_zmq_pull_socket:


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

 
.. c:function:: end_zmq_push_socket:


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

       * :c:func:`ao_two_e_integrals_erf_in_map_slave`
       * :c:func:`ao_two_e_integrals_in_map_slave`
       * :c:func:`davidson_run_slave`
       * :c:func:`run_pt2_slave`
       * :c:func:`run_selection_slave`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:function:: end_zmq_sub_socket:


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

 
.. c:function:: end_zmq_to_qp_run_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)


    Terminate the socket from the application to qp_run

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_erf_in_map_collector`
       * :c:func:`ao_two_e_integrals_erf_in_map_slave`
       * :c:func:`ao_two_e_integrals_in_map_collector`
       * :c:func:`ao_two_e_integrals_in_map_slave`
       * :c:func:`davidson_run_slave`
       * :c:func:`end_parallel_job`
       * :c:func:`pt2_collector`
       * :c:func:`run_pt2_slave`
       * :c:func:`run_selection_slave`
       * :c:func:`selection_collector`

 
.. c:function:: get_task_from_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function get_task_from_taskserver(zmq_to_qp_run_socket,worker_id,task_id,task)


    Get a task from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: get_tasks_from_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function get_tasks_from_taskserver(zmq_to_qp_run_socket,worker_id,task_id,task,n_tasks)


    Get multiple tasks from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: new_parallel_job:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,name_in)


    Start a new parallel job with name 'name'. The slave tasks execute subroutine 'slave'

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`
       * :c:data:`zmq_socket_pull_tcp_address`
       * :c:data:`zmq_context`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:func:`h_s2_u_0_nstates_zmq`
       * :c:func:`zmq_pt2`
       * :c:func:`zmq_selection`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lowercase`
       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:function:: new_zmq_pair_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function new_zmq_pair_socket(bind)


    Socket on which the collector and the main communicate

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_socket_pull_tcp_address`
       * :c:data:`zmq_context`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:function:: new_zmq_pull_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function new_zmq_pull_socket()


    Socket on which the results are sent. If thread is 1, use inproc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`
       * :c:data:`zmq_socket_pull_tcp_address`
       * :c:data:`zmq_context`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`
       * :c:func:`sleep`

 
.. c:function:: new_zmq_push_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function new_zmq_push_socket(thread)


    Socket on which the results are sent. If thread is 1, use inproc

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_socket_pull_tcp_address`
       * :c:data:`zmq_context`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:function:: new_zmq_sub_socket:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function new_zmq_sub_socket()


    Socket to read the state published by the Task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_socket_pull_tcp_address`
       * :c:data:`zmq_context`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`

 
.. c:function:: new_zmq_to_qp_run_socket:


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

 
.. c:function:: reset_zmq_addresses:


    File : :file:`zmq/utils.irp.f`

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

 
.. c:function:: switch_qp_run_to_master:


    File : :file:`zmq/utils.irp.f`

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

 
.. c:function:: task_done_to_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function task_done_to_taskserver(zmq_to_qp_run_socket, worker_id, task_id)


    Get a task from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: tasks_done_to_taskserver:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function tasks_done_to_taskserver(zmq_to_qp_run_socket, worker_id, task_id, n_tasks)


    Get a task from the task server

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: wait_for_next_state:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine wait_for_next_state(state)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`end_zmq_sub_socket`

 
.. c:function:: wait_for_state:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        subroutine wait_for_state(state_wait,state)


    Wait for the ZMQ state to be ready

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`end_zmq_sub_socket`

 
.. c:function:: wait_for_states:


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

 
.. c:function:: zmq_abort:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_abort(zmq_to_qp_run_socket)


    Aborts a running parallel computation

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`sleep`

 
.. c:function:: zmq_delete_task:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_delete_task(zmq_to_qp_run_socket,zmq_socket_pull,task_id,more)


    When a task is done, it has to be removed from the list of tasks on the qp_run
    queue. This guarantees that the results have been received in the pull.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_delete_tasks:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_delete_tasks(zmq_to_qp_run_socket,zmq_socket_pull,task_id,n_tasks,more)


    When a task is done, it has to be removed from the list of tasks on the qp_run
    queue. This guarantees that the results have been received in the pull.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_delete_tasks_async_recv:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_delete_tasks_async_recv(zmq_to_qp_run_socket,zmq_socket_pull,task_id,n_tasks,more)


    When a task is done, it has to be removed from the list of tasks on the qp_run
    queue. This guarantees that the results have been received in the pull.

 
.. c:function:: zmq_delete_tasks_async_send:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_delete_tasks_async_send(zmq_to_qp_run_socket,zmq_socket_pull,task_id,n_tasks,more)


    When a task is done, it has to be removed from the list of tasks on the qp_run
    queue. This guarantees that the results have been received in the pull.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_get8_dvector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get8_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Get a float vector from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`
       * :c:data:`mpi_master`

 
.. c:function:: zmq_get8_ivector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get8_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Get a vector of integers from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`
       * :c:data:`mpi_master`

 
.. c:function:: zmq_get_dmatrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_dmatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Get a float vector from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`
       * :c:data:`mpi_master`

 
.. c:function:: zmq_get_dvector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Get a float vector from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`
       * :c:data:`mpi_master`

 
.. c:function:: zmq_get_i8matrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_i8matrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Get a float vector from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`
       * :c:data:`mpi_master`

 
.. c:function:: zmq_get_imatrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_imatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Get a float vector from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`
       * :c:data:`mpi_master`

 
.. c:function:: zmq_get_int:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_int(zmq_to_qp_run_socket, worker_id, name, x)


    Get a vector of integers from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`
       * :c:data:`mpi_master`

 
.. c:function:: zmq_get_int_nompi:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_int_nompi(zmq_to_qp_run_socket, worker_id, name, x)


    Get a vector of integers from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_get_ivector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_get_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Get a vector of integers from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`
       * :c:data:`mpi_master`

 
.. c:function:: zmq_port:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        function zmq_port(ishift)


    Return the value of the ZMQ port from the corresponding integer

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_run_address`

 
.. c:function:: zmq_put8_dvector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put8_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Put a float vector on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_put8_ivector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put8_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Put a vector of integers on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_dmatrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_dmatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Put a float vector on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_dvector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Put a float vector on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_i8matrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_i8matrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Put a float vector on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_imatrix:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_imatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)


    Put a float vector on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_int:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_int(zmq_to_qp_run_socket, worker_id, name, x)


    Put a vector of integers on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_put_ivector:


    File : :file:`zmq/put_get.irp.f`

    .. code:: fortran

        integer function zmq_put_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)


    Put a vector of integers on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`zmq_state`

 
.. c:function:: zmq_set_running:


    File : :file:`zmq/utils.irp.f`

    .. code:: fortran

        integer function zmq_set_running(zmq_to_qp_run_socket)


    Set the job to Running in QP-run

