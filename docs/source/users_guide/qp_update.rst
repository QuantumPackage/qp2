.. _qp_update:

=========
qp_update
=========

.. program:: qp_update

This command makes an update of the |QP| to the latest version of the
``master`` branch.

Usage
-----

.. code:: bash

    qp_update [-h]

.. option:: -h, --help

    Print the help message and exit.

After asking for confirmation, the command runs the following steps in
:file:`$QP_ROOT`:

.. code:: bash

    ninja clean
    git pull origin master
    ninja

The update therefore requires :file:`$QP_ROOT` to be a git clone with a
``master`` branch, and it rebuilds the whole code from scratch. Local
changes that conflict with the upstream version make ``git pull`` fail;
the rebuild still runs, on the unchanged tree. As with any update, the
|EZFIO| files
created by a previous version should be checked with :ref:`qp_edit`
with ``--check``, or simply re-created.
