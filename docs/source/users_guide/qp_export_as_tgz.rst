================
qp_export_as_tgz
================

.. program:: qp_export_as_tgz

In some HPC facilities, the access to the internet is limited for
security reasons. In such an environment, the installation of |QP| is
sometimes very painful because the OCaml compiler and the libraries
can't be installed by a non-root user.

This command creates a self-contained binary distribution in the form of
a `tar.gz` file that can be copied on another machine.

Usage
-----

.. code:: bash

    qp_export_as_tgz [-h|--help]

.. option:: -h, --help

    Prints the help message

.. note::
   There can be conflicts due to the version of glibc. The machine on which |QP| is
   compiled should be the oldest one.
   

