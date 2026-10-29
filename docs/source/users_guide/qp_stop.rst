.. _qp_stop:

=======
qp_stop
=======

.. program:: qp_stop

Requests for a clean termination of the program.

This will have the effect to exit the Davidson diagonalization, the
|SCF| procedure or the determinant selection to save the current wave
function and exit the program.

The program checks two marker files in ``$EZFIO_DIR/work``: ``qpkill``,
written by ``qp_stop`` without option, and ``qpstop``, written by
``qp_stop -s``. The marker is read at the end of the current iteration,
so the wave function and the |EZFIO| directory are left in a consistent
state.

Usage
-----

.. code:: bash

  qp_stop [-chqs] EZFIO_DIR

.. option:: -c, --cancel

  Cancel the qp_stop order by removing the marker file. Without ``-s``,
  ``qpkill`` is removed; with ``-s``, ``qpstop`` is removed.

.. option:: -h, --help

  Print the help message.

.. option:: -q, --query

  Ask if :file:`EZFIO_DIR` was requested to stop. The exit status is 0
  when a stop was requested and 1 otherwise. As for ``-c``, the file
  checked depends on the presence of ``-s``.

.. option:: -s, --soft

  Terminate the current step, but don't kill the program. The ``qpstop``
  marker is deleted as soon as the program reads it, so the next
  iterations continue normally. In contrast, the marker written by
  ``qp_stop`` without ``-s`` stays in place and every subsequent run on
  this |EZFIO| directory stops immediately, until it is removed with
  ``qp_stop --cancel``.
