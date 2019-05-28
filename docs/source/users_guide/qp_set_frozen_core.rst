.. _qp_set_frozen_core:

==================
qp_set_frozen_core
==================

.. program:: qp_set_frozen_core

Automatically finds *n*, the number of core electrons. Calls
:ref:`qp_set_mo_class` setting all |MOs| as ``Active``, except the
:math:`n/2` first ones which are set as ``Core``. If pseudo-potentials
are used, all the |MOs| are set as ``Active``.


========== ========= ======= ======= 
  Range     Default   Small   Large  
========== ========= ======= ======= 
 H  -> He      0        0        0   
 Li -> Be      0        0        2   
 B  -> Ne      2        2        2   
 Na -> Mg      2        2       10   
 Al -> Ar     10        2       10   
 K  -> Ca     10       10       18   
 Sc -> Zn     10       10       18   
 Ga -> Kr     18       10       18   
 Rb -> Sr     18       18       36   
 Y  -> Cd     18       18       36   
 In -> Xe     36       18       36   
 Cs -> Ba     36       36       54   
 La -> Hg     36       36       54   
 Tl -> Rn     54       36       54   
 Fr -> Ra     54       54       86   
 Ac -> Cn     54       54       86   
 Nh -> Og     86       54       86   
========== ========= ======= ======= 

For elements on the right of the periodic table, `qp_set_frozen_core`
will work as expected. But for elements on the left, a small core will
be chosen. For example, a Carbon atom will have 2 core electrons, but a
Lithium atom will have zero.

Usage 
-----

.. code:: bash

      qp_set_frozen_core [-q|--query] [(-l|-s|--large|--small)  EZFIO_DIR


.. option:: -q, --query 

    Prints in the standard output the number of core electrons.

.. option:: -s, --small

    Use a small core.

.. option:: -l, --large

    Use a large core.


