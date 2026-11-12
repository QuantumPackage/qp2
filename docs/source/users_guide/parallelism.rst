===========
Parallelism
===========

|qp| runs on a single machine and parallelizes with openMP threads inside
each executable: the threads accelerate the linear algebra, the stochastic
sampling and the determinant selection. The same calculation runs
unchanged on a laptop or on one node of a cluster.


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
written during the calculation, it can become large, and it should sit
on a fast local file system. Remove :file:`work` before archiving a
database.


SLURM
=====

On a cluster that uses |SLURM|, a single-node job is submitted as usual.
The job script requests one node and the thread count of the machine:

.. code:: bash

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=32
    #SBATCH --mem=120G
    #SBATCH --signal=B:USR1@120

    source ${QP_ROOT}/quantum_package.rc
    export OMP_NUM_THREADS=32
    export QP_NTHREADS_PT2=32
    export QP_NTHREADS_DAVIDSON=32
    export QP_MAXMEM=110

    qp run fci hcn.ezfio

Running

.. code:: bash

    qp run four_idx_transform hcn.ezfio

once stores the |MO| integrals in :file:`$EZFIO/work` and subsequent runs
read them directly.

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

.. seealso::

    The manual pages of :ref:`qp_run` and :ref:`qp_stop`, and the
    ``etc/local.rc`` file of the installation for the environment
    variables of the local machine.
