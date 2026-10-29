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
