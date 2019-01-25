.. _qpsh:

====
qpsh
====

.. program:: qpsh


:command:`qpsh` is the |qp| shell. It is a Bash shell with all the
required evironment variables loaded, a modified prompt, and the
:ref:`qp` command.


.. _qp:

.. program:: qp

qp
==

This command is a hub to the most used command within |qp|. The power
of the :ref:`qpsh` shell is the auto-completion that comes when the
:kbd:`<Tab>` key is pressed with the :ref:`qp` command.


EZFIO access
------------

.. option:: set_file 

   .. code:: bash

     qp set_file EZFIO_DIR 

   Sets the current |EZFIO| directory. All the following instruction will be
   relative to this directory.

   This command is equivalent to :command:`ezfio set_file EZFIO_DIR`.


.. option:: unset_file

   .. code:: bash

     qp unset_file

   Unsets the current |EZFIO| directory.

   This command is equivalent to :command:`ezfio unset_file`.


.. option:: has

   .. code:: bash

     qp has <module> <parameter> 

   If the `<module>/<parameter>` is set in the |EZFIO| directory, returns 1.
   Otherwise returns 0.

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

  Runs :ref:`qp_run`, :ref:`qp_srun`, or :ref:`qp_mpirun` using the current 
  |EZFIO| directory.

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

.. option:: reset

  :command:`qp reset` : runs :ref:`qp_reset`

.. option:: set_frozen_core

  :command:`qp set_frozen_core` : runs :ref:`qp_set_frozen_core`

.. option:: set_mo_class

  :command:`qp set_mo_class` : runs :ref:`qp_set_mo_class`

.. option:: update

  :command:`qp update` : runs :ref:`qp_update`





