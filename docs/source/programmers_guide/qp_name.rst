.. _qp_name:

qp_name
=======

.. program:: qp_name

Displays the names of all the files in which the provider/subroutine/function
given as argument is used. With the `-r` flag, the name can be changed in the
whole quantum package.

Usage
-----

.. code:: bash

   qp_name <name> [-r <new_name> | --rename=<new_name>]


.. option:: -h

   Prints the help message.


.. option:: -r <new_name> --rename=<new_name>

   Renames the provider/subroutine/function and all its occurences.


.. note::

    It is safe to create a commit before renaming a provider, and then to
    check what has changed using git diff.


