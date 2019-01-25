.. _qp_set_frozen_core:

==================
qp_set_frozen_core
==================

.. program:: qp_set_frozen_core

Automatically finds *n*, the number of core electrons. Calls
:ref:`qp_set_mo_class` setting all |MOs| as ``Active``, except the
:math:`n/2` first ones which are set as ``Core``. If pseudo-potentials
are used, all the |MOs| are set as ``Active``.

For elements on the right of the periodic table, `qp_set_frozen_core`
will work as expected. But for elements on the left, a small core will
be chosen. For example, a Carbon atom will have 2 core electrons, but a
Lithium atom will have zero.

Usage 
-----

.. code:: bash

      qp_set_frozen_core [-q]  EZFIO_DIR


.. option:: -q 

    Prints in the standard output the number of core electrons.


