.. _qp_stop:

=======
qp_stop
=======

.. program:: qp_stop

Requests for a clean termination of the program.

This will have the effect to exit the Davidson diagonalization, the
|SCF| procedure or the determinant selection to save the current wave
function and exit the program.

Usage
-----

.. code:: bash

  qp_stop [-chq] EZFIO_DIR

.. option:: -c, --cancel

  Cancel the qp_stop order.

.. option:: -q, --query

  Ask if :file:`EZFIO_DIR` was requested to stop.


