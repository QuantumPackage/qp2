.. _qp_edit:

=======
qp_edit
=======

.. program:: qp_edit


This command reads the content of the |EZFIO| directory and creates
a temporary file containing the data. The data is presented as a
*ReStructured Text* (rst) document, where each section corresponds to
the corresponding |qp| module. The content of the file can be modified
to change the input parameters. When the text editor is closed, the
updated data is saved into the |EZFIO| directory.

.. note::
   The text editor which will be opened is defined by the :envvar:`EDITOR`
   environment variable. If this variable is not set, the :command:`vi`
   text editor will be used by default.
   
.. warning::
   When the wave function is too large (more than 10 000 determinants), the
   determinants are not displayed.

.. note::
   On some machines the terminal will be stuck in inverted colors after using
   qp_edit. To Avoid this problem, put in your :file:`$HOME/.vimrc`::

        set t_ti=
        set t_te=


Usage
-----

.. code:: bash

  qp_edit [-c] [-h] [-n <int>] [-s <range>] [--] EZFIO_DIR          

.. option:: -c, --check

   Checks the input data

.. option:: -h, --help

   Print the help text and exits

.. option:: -n, --ndet=<int>

   Truncates the wavefunction to the target number of determinants

.. option:: -s, --state=<range>

   Select the states to extract from the |EZFIO| directory, using the same conventions
   as :ref:`qp_set_mo_class`. See example below.


Example
-------

.. code:: bash

   qp_edit --state="[1,3-5]" test.ezfio

Removes all states except states 1,3,4 and 5 from :file:`test.ezfio`.
The resulting |EZFIO| directory has 4 states.

