.. _qp_create_ezfio:

qp_create_ezfio
===============

.. program:: qp_create_ezfio

This command creates an |EZFIO| directory from a standard `xyz` file or
from a `z-matrix` file in Gaussian format.

Usage 
-----

.. code:: bash

   qp_create_ezfio [-a] -b <string> [-c <int>] [-d <float>] 
      [-h] [-m <int>] [-o EZFIO_DIR] [-p <string>] [-x] [--] FILE


.. option:: -a, --au

   If present, input geometry is in atomic units.


.. option:: -b, --basis=<string>

   Name of basis set. The basis set is defined as a single string if
   all the atoms are taken from the same basis set, otherwise specific
   elements can be defined as follows::

      -b "cc-pcvdz | H:cc-pvdz | C:6-31g"
      -b "cc-pvtz | 1,H:sto-3g | 3,H:6-31g"

   By default, the basis set is obtained from the local database of the.
   |qp| This option is mandatory                                       .


.. option:: -c, --charge=<int>

   Total charge of the molecule. Default is 0.


.. option:: -d, --dummy=<float>

   Add dummy atoms (X) between atoms when the distance between two atoms
   is less than :math:`x \times \sum R_\mathrm{cov}`, the covalent radii
   of the atoms. The default is x=0, so no dummy atom is added.


.. option:: -h, --help

   Print the help text and exit


.. option:: -m, --multiplicity=<int>

   Spin multiplicity :math:`2S+1` of the molecule. Default is 1.


.. option:: -o, --output=EZFIO_DIR

   Name of the created |EZFIO| directory.

.. option:: -p <string>, --pseudo=<string>

   Name of the pseudo-potential. Follows the same conventions as the basis set.

.. option:: -x, --cart

   Compute |AOs| in the Cartesian basis set (6d, 10f, ...)


Using custom atomic basis sets
------------------------------

If a file with the same name as the basis set exists, this file will
be read. For example, if the file containing the basis set is named
``custom.basis``, and the *xyz* geometry is in ``molecule.xyz``, the
following should be used::

    qp_create_ezfio -b custom.basis molecule.xyz

Basis set files should be given in |GAMESS| format, without combined sp, spd, ...
contractions. The full names of the atoms are given, and the basis sets for each element are
separated by a blank line. Here is an example ::

      HYDROGEN
      S   3
      1     13.0100000              0.0196850
      2      1.9620000              0.1379770
      3      0.4446000              0.4781480
      S   1
      1      0.1220000              1.0000000
      P   1
      1      0.7270000              1.0000000

      BORON
      S   8
      1   4570.0000000              0.0006960
      2    685.9000000              0.0053530
      3    156.5000000              0.0271340
      4     44.4700000              0.1013800
      5     14.4800000              0.2720550
      6      5.1310000              0.4484030
      7      1.8980000              0.2901230
      8      0.3329000              0.0143220
      S   8
      1   4570.0000000             -0.0001390
      2    685.9000000             -0.0010970
      3    156.5000000             -0.0054440
      4     44.4700000             -0.0219160
      5     14.4800000             -0.0597510
      6      5.1310000             -0.1387320
      7      1.8980000             -0.1314820
      8      0.3329000              0.5395260
      S   1
      1      0.1043000              1.0000000
      P   3
      1      6.0010000              0.0354810
      2      1.2410000              0.1980720
      3      0.3364000              0.5052300
      P   1
      1      0.0953800              1.0000000
      D   1
      1      0.3430000              1.0000000


Files can be extracted from the Basis Set Exchange database 
https://www.basissetexchange.org , with the ``qp_basis`` tool.




Using custom pseudo-potentials
------------------------------

As for the basis set, if a file with the same name as the
pseudo-potential exists, this file will be read. For example, if the
file containing the custom pseudo-potential is named ``custom.pseudo``,
the basis set is named ``custom.basis``, and the *xyz* geometry is in
``molecule.xyz``, the following command should be used

.. code:: bash

    qp_create_ezfio -b custom.basis -p custom.pseudo molecule.xyz

Pseudo-potential files should be given in a format very close to
|GAMESS| format. The first line should be formatted as ``%s GEN %d %d``
where the first string is the chemical symbol, the first integer is
the number of core electrons to be removed and the second integer is
LMAX+1 as in |GAMESS| format. The pseudo-potential for each element are
separated by a blank line. Here is an example ::

      Ne GEN 2 1
      3
      8.00000000 1 10.74945199
      85.99561593 3 10.19801460
      -56.79004456 2 10.18694048
      1
      55.11144535 2 12.85042963

      F GEN 2 1
      3
      7.00000000 1 11.39210685
      79.74474797 3 10.74911370
      -49.45159098 2 10.45120693
      1
      50.25646328 2 11.30345826




