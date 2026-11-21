=================
Quick-start guide
=================

This tutorial walks through a complete frozen-core |CIPSI| calculation on
the HCN molecule in the 6-31G basis set, and it ends with an estimate of
the |FCI| energy in this basis. It introduces the commands that the other
pages of this guide rely on, and it assumes that |qp| is already
installed and compiled (see :doc:`../intro/install`).


Demo video
==========

The same tutorial is available as a video:
`<https://www.youtube.com/watch?v=4nmdCAPkZlc>`_


Loading the |QP| environment
============================

The commands described below are shell functions and programs that need
the environment variables defined in :file:`quantum_package.rc`. Load it
once in the shell where you work:

.. code:: bash

    source ${QP_ROOT}/quantum_package.rc

where :file:`${QP_ROOT}` is the installation directory of |qp|. The
script adds the |qp| programs to :envvar:`PATH` and defines the ``qp``
shell function described below.

Another way to obtain the same environment is to start a |qpsh|:

.. code:: bash

    ${QP_ROOT}/bin/qpsh

This opens a new interactive shell that has already loaded
:file:`quantum_package.rc`, and whose prompt shows the |EZFIO| directory
currently in use.

Each command of |qp| exists as an executable named ``qp_<command>`` that
takes the |EZFIO| directory as its last argument. The ``qp`` shell
function calls that executable with the current database, so both lines
below run the |SCF| program on the directory ``hcn`` when
:envvar:`EZFIO_FILE` is set to ``hcn``:

.. code:: bash

    qp_run scf hcn
    qp run scf

The second form is shorter, it completes the command and option names
with the :kbd:`Tab` key, and it gives access to the manual pages:

.. code:: bash

    qp --help
    qp man scf


Create the |EZFIO| database
===========================

All the data of a calculation are stored in an |EZFIO| database. |EZFIO|
is a hierarchical data format that uses the file system hierarchy to
organize the data: each quantity is a small file in a directory named
after the module that defines it, plain text for the scalars and
gzipped for the arrays, so the content of a database can be inspected
with any text editor. Programs should access it
through the APIs provided by |EZFIO|, in Fortran, |Python|, |OCaml| or
Bash, or through tools built on these APIs such as :ref:`qp_edit`.

First write an ``xyz`` file with the geometry of the molecule in
:file:`hcn.xyz`::

    3
    HCN molecule
    C    0.0    0.0    0.0
    H    0.0    0.0    1.064
    N    0.0    0.0    -1.156

Then create the database with the 6-31G basis set:

.. code:: bash

    qp create_ezfio -b "6-31g" hcn.xyz -o hcn

The command creates the directory :file:`hcn` and makes it the current
|EZFIO| database. The other options of :ref:`qp_create_ezfio` set the
charge, the multiplicity, a pseudopotential or ghost atoms.

The database now contains the geometry, the basis set and the number of
electrons, which can be queried from the shell:

.. code:: bash

    $ qp get nuclei nucl_num
    3
    $ qp get electrons elec_alpha_num
    7
    $ qp get ao_basis ao_num
    20

In a new shell, the same database is selected with

.. code:: bash

    qp set_file hcn


Run a Hartree-Fock calculation
==============================

The program :ref:`qp_run` is the driver that runs all the other programs.
A restricted Hartree-Fock calculation is obtained with

.. code:: bash

    qp run scf

The calculation converges in nine iterations to ``-92.82785666`` au. The
:ref:`scf` program stores the |MOs| in the |EZFIO| database, where all
the post-Hartree-Fock methods read them, and prints their Fock
eigenvalues. The |MOs| of this calculation are already ordered by
increasing Fock energy; if they are not,
:ref:`sort_by_fock_energies` reorders them.

.. seealso::

    The documentation of the :ref:`module_hartree_fock` module and of the
    :ref:`scf` program.


Choose the target |MO| space
============================

The next step restricts the |CIPSI| calculation to the valence |MOs| and
keeps the core |MOs| frozen. The command :ref:`qp_set_frozen_core` finds
the number of core electrons of each atom and sets the orbital classes
accordingly:

.. code:: bash

    qp set_frozen_core

For HCN in 6-31G the 20 |MOs| are split into 2 core orbitals and 18
active orbitals. The result is stored in the |EZFIO| database and can be
checked with the :ref:`qp_set_mo_class` command:

.. code:: bash

    $ qp set_mo_class -q
    Electrons: 7 7
    MO  : 20
    n_int: 1
    Core     : (1, 2)
    Inactive
    Active   : (3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
    Virtual
    Deleted

The same classes are obtained explicitly with the ranges

.. code:: bash

    qp set_mo_class --core="[1-2]" --act="[3-20]"

where ``--core`` orbitals are always doubly occupied, ``--inact``
orbitals can only receive holes, ``--virt`` orbitals can only receive
particles, and ``--del`` orbitals are never occupied. Ranges are written
as ``[36-53,72-107]``.


Run the |CIPSI| calculation
===========================

The |CIPSI| algorithm estimates the |FCI| energy by iteratively enlarging
a variational space. The calculation starts from the Hartree-Fock
determinant and, at each iteration, selects the determinants of the
external space that contribute most to the |PT2| correction and adds them
to the internal space. All the determinants needed to make the wave
function an eigenfunction of |S^2| are added as well, because
:option:`determinants s2_eig` is |true| by default. The Hamiltonian is
then diagonalized in the enlarged space, the |PT2| correction is
estimated by stochastic sampling, and the variational energy is
extrapolated to the |FCI| limit by fitting
:math:`E = E_\text{FCI} - \alpha E_\text{PT2}`.

.. code:: bash

    qp run fci | tee hcn.fci.out

The trajectory of the calculation is printed in the output file, with one
summary block per iteration:

.. code:: text

    -----------------------------------
    Summary at N_det =        20154
    -----------------------------------
    ...
    # E             -93.04419463
    # PT2            -0.00567047
    # E+PT2         -93.04986510

The loop stops when the number of determinants reaches
:option:`determinants n_det_max` (one million by default) or when the
|PT2| energy becomes smaller than :option:`perturbation pt2_max` (10\ :sup:`-4`
au by default). The stochastic evaluation of |PT2| itself stops when its
relative error reaches :option:`perturbation pt2_relative_error`.

The run of this tutorial ends at 1 290 340 determinants, where both
criteria are met. Because the |PT2| estimate is sampled stochastically,
the digits printed by a new run will differ slightly; the |PT2|-corrected
energy and the fit are stable to about five decimal places. A few of the
summary blocks show how the energy converges:

.. list-table::
   :header-rows: 1

   * - Number of determinants
     - Variational energy
     - |PT2|
     - E + |PT2|
   * - 1 256
     - -93.02808189
     - -0.02229222
     - -93.05037410
   * - 20 154
     - -93.04419463
     - -0.00567047
     - -93.04986510
   * - 161 277
     - -93.04909906
     - -0.00093748
     - -93.05003653
   * - 645 164
     - -93.04987142
     - -0.00017428
     - -93.05004571
   * - 1 290 340
     - -93.04997937
     - -0.00006896
     - -93.05004833

The last iteration saves the final wave function and prints its
properties, including the spin expectation value that confirms the
singlet character of the state:

.. code:: text

    * State            1
    < S^2 >         =    2.4553002300743978E-011
    E               =   -93.049979367027902
    PT2             =   -6.8962142235121408E-005  +/-  1.3802497317171859E-007
    E+PT2           =   -93.050048329170139       +/-  1.3802497317171859E-007

The fit of the last iterations gives ``-93.05005`` au for the smallest
|PT2| cutoffs, in agreement with the |PT2| corrected energy. This is the
estimated |FCI| energy of HCN in 6-31G with a frozen core, to be compared
with the |HF| energy of ``-92.82785666`` au.

Because the wave function is stored in the |EZFIO| database after each
iteration, the calculation can be continued from where it stopped by
setting :option:`determinants read_wf`, and the :ref:`pt2` program
recomputes its |PT2| correction without performing any selection:

.. code:: bash

    qp set determinants read_wf true
    qp run fci
    qp run pt2

.. seealso::

    The documentation of the :ref:`module_fci` module and of the
    :ref:`fci` program. :doc:`cipsi` describes the selection algorithm
    and its parameters in more detail.


What next
=========

The same workflow applies to any molecule and basis set, and the
parameters described in :doc:`cipsi` control the accuracy and the cost of
the calculation. :doc:`excited_states` explains how to compute several
states, :doc:`parallelism` how to run the calculation on several cores or
nodes, and :doc:`interfaces` how to export the wave function to other
quantum-chemistry codes.
