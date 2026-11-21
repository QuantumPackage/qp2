.. _qp_reset:

========
qp_reset
========

.. program:: qp_reset

This command resets parts of the |EZFIO| directory. Whatever the
options, :ref:`qp_edit` is run with ``--check`` at the end to restore a
consistent directory.

Usage
-----

.. code:: bash

  qp_reset [-adhm] EZFIO_DIR

.. option:: -a, --all

  Reset to the state in which the directory is after running :ref:`qp_create_ezfio`.
  This is the same as ``--dets --mos``.

.. option:: -d, --dets

  Deletes the determinants and CI coefficients: ``determinants/n_det``,
  ``determinants/psi_det.gz``, ``determinants/psi_coef.gz`` and the
  ``qp_edit`` copies, as well as the bi-orthogonal coefficients.

.. option:: -h, --help

  Print the help message.

.. option:: -m, --mos

  Deletes the |MOs| (the ``mo_basis`` and ``bi_ortho_mos`` directories)
  and the memory-mapped |MO| integrals in ``work/mo_ints_*``, and resets
  the ``io_*`` flags of the one- and two-electron integrals. If the
  active space was changed, a warning reminds you that the |MO| classes
  have to be defined again.

  The determinants are *not* deleted by this option alone: pass
  ``--dets`` as well, or use ``--all``, to start again from scratch.
