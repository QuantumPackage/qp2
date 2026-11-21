.. _qpsh:

====
qpsh
====

.. program:: qpsh


:command:`qpsh` is the |qp| shell. It is a Bash shell with all the
required environment variables loaded, a modified prompt, and the
:ref:`qp` command.


.. _qp:

.. program:: qp

qp
==

This command is a hub to the most used command within |qp|. The power
of the :ref:`qpsh` shell is the auto-completion that comes when the
:kbd:`<Tab>` key is pressed with the :ref:`qp` command.

:ref:`qp` is a shell function defined in :file:`$QP_ROOT/etc/qp.rc`, so
it is only available after sourcing :file:`quantum_package.rc`, either
directly or through :command:`qpsh`.


EZFIO access
------------

.. option:: set_file

   .. code:: bash

     qp set_file EZFIO_DIR

   Sets the current |EZFIO| directory. All the following instructions will be
   relative to this directory.

   It behaves like :command:`ezfio set_file EZFIO_DIR`, except that the
   directory must already exist, while :command:`ezfio set_file` creates
   it.


.. option:: unset_file

   .. code:: bash

     qp unset_file

   Unsets the current |EZFIO| directory.

   This command is equivalent to :command:`ezfio unset_file`.


.. option:: has

   .. code:: bash

     qp has <module> <parameter>

   Returns 0 if `<module>/<parameter>` is set in the |EZFIO| directory,
   and 1 otherwise. This makes it usable directly in a shell conditional:

   .. code:: bash

     if qp has davidson n_states_diag ; then
       echo "n_states_diag is set"
     fi

   This command is equivalent to :command:`ezfio has <module> <parameter>`.


.. option:: get

   .. code:: bash

     qp get <module> <parameter>

   Returns the value of `<module>/<parameter>`.

   This command is equivalent to :command:`ezfio get <module> <parameter>`.


.. option:: set

   .. code:: bash

     qp set <module> <parameter> [<value>]

   Sets the value of `<module>/<parameter>`. If the value is not given in
   the command line it is read from the standard input.

   This command is equivalent to
   :command:`ezfio set <module> <parameter> [<value>]`.


Running programs
----------------

.. option:: run

   .. code:: bash

     qp (run|srun|mpirun) [options] <program>

   Runs :ref:`qp_run`, :command:`qp_srun`, or :command:`qp_mpirun` using the
   current |EZFIO| directory.

.. option:: stop

   :command:`qp stop` : runs :ref:`qp_stop`

Getting help
------------

.. option:: man

   .. code:: bash

     qp man (<program>|<qp_command>)

   Displays a man page for a |qp| program or a |qp| command.


Running quantum package commands
--------------------------------

The ``qp_`` commands can be run without specifying the |EZFIO| directory:

.. option:: convert_output_to_ezfio

   :command:`qp convert_output_to_ezfio` : runs :ref:`qp_convert_output_to_ezfio`

.. option:: create_ezfio

   :command:`qp create_ezfio` : runs :ref:`qp_create_ezfio`

.. option:: plugins

   :command:`qp plugins` : runs :ref:`qp_plugins`

.. option:: update

   :command:`qp update` : runs :ref:`qp_update`

Any other ``qp_`` command can be called directly from :ref:`qpsh`,
since :file:`quantum_package.rc` puts ``$QP_ROOT/bin`` in the
:envvar:`PATH`.
