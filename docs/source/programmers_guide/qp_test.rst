.. _qp_test:

=======
qp_test
=======

.. program:: qp_test

This command runs the consistency tests of |qp|.
The tests are run with the |Bats| shell testing environment.
If the name of a test or its number is specified on the command line, only this
test will be run.

Usage
-----

.. code:: bash

    qp_test [FLAGS] [TEST]
    Flags :
      [-a]          Run all the tests
      [-v]          Verbose mode: shows the output of the runs


.. option:: -a

   Runs the tests of all the installed modules. By default, run only the
   tests of the current directory, and the directories below.

.. option:: -v

   Verbose mode. Print the output of the running executions of |qp|.


Testing the command-line tools
-------------------------------

:command:`qp_test` exercises the tools end to end: the |Bats| tests of
the modules call :command:`qp create_ezfio`, :command:`qp edit`,
:command:`qp run` and the other tools on real |EZFIO| directories, so the
tools must have been built first (see :ref:`rust_tooling_build`).

The tool layer has its own test suite, which runs without a Fortran
build:

.. code:: bash

    source quantum_package.rc
    cd $QP_ROOT/rust
    cargo test --workspace -j4



