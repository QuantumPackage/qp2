===========
Parallelism
===========

|qp| uses two independent levels of parallelism. Inside one process,
openMP threads accelerate the linear algebra and the stochastic sampling.
Between processes, a task server distributes the elementary tasks over
one or several machines. Both levels can be combined, and the same
calculation runs unchanged on a laptop, on one node or on a cluster.


Threads
=======

The programs are compiled with openMP and the number of threads is fixed
by the usual environment variables before the calculation starts:

.. code:: bash

    export OMP_NUM_THREADS=8
    export MKL_NUM_THREADS=8

The first one sets the number of openMP threads, and the second one the
number of threads used by the |BLAS| library. Setting them to the number
of physical cores of the machine is usually a good choice, and the
comments in :file:`${QP_ROOT}/etc/local.rc` note that hyperthreading is
profitable for some steps.

Two variables tune the two most expensive parts of a |CIPSI| calculation
independently of the number of openMP threads: :envvar:`QP_NTHREADS_PT2`
for the stochastic evaluation of the |PT2| correction,
:envvar:`QP_NTHREADS_DAVIDSON` for the matrix-vector products of the
Davidson diagonalization. Both default to the number of openMP threads,
and the values in use are printed in the output of the programs, for
example ``Number of threads for PT2``.

Set these variables in your job script or in
:file:`${QP_ROOT}/etc/local.rc`, which is read by
:file:`quantum_package.rc` together with the other files of the
:file:`etc` directory.


Memory
======

Before a large allocation, the programs check that the memory needed by
the new arrays, added to the memory already in use, fits in a limit given
in gigabytes by :envvar:`QP_MAXMEM`. When the limit is exceeded the
program prints ``Not enough memory: aborting in <routine>`` and stops.
The default limit is the total memory of the machine minus one gigabyte,
with a minimum of three gigabytes.

The molecular integrals are stored in the :file:`work` directory of
the |EZFIO| database as memory-mapped files. This directory is
written during the calculation, it can become large, and it must be
visible from all the processes that take part in the calculation. On
a cluster, keep the |EZFIO| directory on a scratch file system for
this reason, and remove :file:`work` before archiving a database.


Distributed calculations
========================

A distributed calculation involves one master process, started with
:ref:`qp_run` as usual, and any number of slave processes started with
the ``--slave`` option. The master checks the input, opens a |ZeroMQ|
socket on a free TCP port and writes its address, for example
``tcp://130.120.229.139:42379``, in the file
:file:`$EZFIO/work/qp_run_address`. The slaves read that file, connect to
the master and print the address they use::

   Switched to qp_run master : tcp://130.120.229.139:42379

The master then sends the elementary tasks of the calculation, which are
the selection of determinants, part of the Davidson diagonalization and
the sampling of the |PT2| correction, and collects the results. Only the
master updates the |EZFIO| database, so the slaves only need to be able to
read it, including its :file:`work` directory. The address file is
removed when the master terminates.

To start a calculation by hand on two machines that share the |EZFIO|
directory, run the master first and then the slave:

.. code:: bash

    qp run fci &
    qp_run --slave fci hcn.ezfio

The master spends its first seconds computing the integrals and the
initial wave function, so starting the slaves a little later is harmless:
they wait for their first task after connecting.

The ``distributed_davidson`` option must be set before a calculation that
runs on several nodes, because the distributed algorithm is only used
when the option is |true|. It takes over when the Davidson space contains
more than 100 000 determinants, which is the typical regime of a large
|CIPSI| calculation:

.. code:: bash

    qp set davidson_keywords distributed_davidson true

Slaves can be started at any moment during a calculation: they connect to
the running master and receive tasks as soon as they are ready.


SLURM
=====

On a cluster that uses |SLURM|, the :command:`qp_srun` script does the
bookkeeping. It reads the list of nodes from ``srun hostname``, starts the
master on one node and one slave on each of the other nodes. The job
script must request one task per node, so that the script sees exactly
one host name per node:

.. code:: bash

    #!/bin/bash
    #SBATCH --nodes=4
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=32
    #SBATCH --mem=120G
    #SBATCH --signal=B:USR1@120

    source ${QP_ROOT}/quantum_package.rc
    export OMP_NUM_THREADS=32
    export QP_NTHREADS_PT2=32
    export QP_NTHREADS_DAVIDSON=32
    export QP_MAXMEM=110

    qp_srun fci hcn.ezfio

The script checks that there is a single task per host, and it warns when
the two-electron integrals are not stored on disk, because each process
would then compute them separately. Running

.. code:: bash

    qp_run four_idx_transform hcn.ezfio

once as a single-node job stores the |MO| integrals in
:file:`$EZFIO/work` and subsequent runs read them directly.

When the scheduler sends the ``USR1`` signal to :ref:`qp_run`, the program
calls :ref:`qp_stop`, which asks the calculation to terminate cleanly. The
option ``--signal=B:USR1@120`` of the example sends this signal two
minutes before the time limit. The :command:`qp_stop` command can also be
run by hand from another shell: it leaves a marker in the :file:`work`
directory, and the running calculation checks it at the end of each
iteration before stopping with the wave function saved on disk. The
marker stays in the database, so it must be cancelled with
``qp_stop --cancel`` before the next run; the ``--soft`` variant instead
removes the marker as soon as the running program reads it.


Other launchers
===============

:command:`qp_mpirun` does the same job as :command:`qp_srun` with an MPI
launcher instead of |SLURM|. It is written for Intel MPI and calls
``mpiexec.hydra``, and it also expects one rank per node.


Separate networks
=================

When the machines are behind different firewalls, the slaves cannot reach
the master directly. The :ref:`qp_tunnel` command opens an SSH tunnel
between the two front-end nodes, so that the slaves of the second cluster
can connect to the master of the first one. Its manual page describes the
two commands to run and the case where the |EZFIO| directory is not on a
shared file system.


Troubleshooting
===============

If a slave stops with ``No master is not running``, the master has not
written the address file yet or the slave does not see the same |EZFIO|
directory. If the connection to the master fails, the address printed by
the slave is a good indication: the master chooses the network interface
by resolving its host name, unless :envvar:`QP_NIC` names one, as in
``export QP_NIC=ib0`` for an InfiniBand network. The messages exchanged
between the two sides are printed when :envvar:`QP_TASK_DEBUG` is set to
``1``.

.. seealso::

    The manual pages of :ref:`qp_run`, :ref:`qp_stop` and :ref:`qp_tunnel`,
    and the ``etc/local.rc`` file of the installation for the environment
    variables of the local machine.
