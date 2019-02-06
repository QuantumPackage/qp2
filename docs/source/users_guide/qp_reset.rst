.. _qp_reset:

========
qp_reset
========

.. program:: qp_reset

This command resets parts of the |EZFIO| directory.


Usage
-----

.. code:: bash

  qp_reset [-adhm] EZFIO_DIR

.. option:: -a, --all

  Reset to the state in which the directory is after after running :ref:`qp_create_ezfio`.

.. option:: -d, --dets

  Deletes the determinants and CI coefficients.

.. option:: -m, --mos

  Deletes the |MOs|, and consequently the determinants and CI coefficients.



