.. _qp_test:

=======
qp_test
=======

.. program:: qp_test

This command runs the consistency test of |qp|.
The tests are run with the |Bats| shell testing environment.
If the name of a test of its number is specified on the command line, only this
test will be run.

Usage
-----

.. code:: bash

    qp_test [FLAGS] [TEST]
    Flags :
      [-a]          Run all the tests
      [-v]          Verbose mode: shows the output of the runs


.. option:: -a

   Runs all the tests. By default, run only the tests of the current
   directory, and the directories below.

.. option:: -v

   Verbose mode. Print the output of the running executions of |qp|.



