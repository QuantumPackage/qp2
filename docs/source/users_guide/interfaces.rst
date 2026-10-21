Interfaces
==========

|qp| reads geometries and integral files produced by other programs and
writes files that other programs can read. The programs described here
are part of the standard distribution unless their name says otherwise;
the ones that come from plugins are installed by :ref:`qp_plugins`.


Reading data from other codes
-----------------------------

A geometry in ``xyz`` format or as a z-matrix is enough to create a
database with :ref:`qp_create_ezfio`, which fills the nuclei, the
electron numbers, the basis set and, when requested, the
pseudopotential. The nuclear repulsion and the integrals are computed
on demand by the programs that need them. Basis sets are read from the
files of :file:`${QP_ROOT}/data/basis`, and a pseudopotential is
selected with the ``-p`` option.

The output of Gaussian can be converted into an |EZFIO| database with

.. code:: bash

    qp convert_output_to_ezfio gaussian.out -o h2o.ezfio

This command reads the geometry, the basis set, the |MOs| and the
integrals through the TREXIO Python packages, so it needs the ``trexio``
and ``trexio-tools`` modules of Python, which are not part of the |qp|
requirements.

Two programs import the content of a TREXIO file: ``import_trexio_integrals``
reads the nuclear repulsion and the |AO| integrals, and
``import_trexio_determinants`` reads a wave function. Both read the file
named by the ``trexio_file`` option, so the database must already contain
the geometry and the basis set:

.. code:: bash

    qp set trexio trexio_file h2o.h5
    qp run import_trexio_integrals
    qp run import_trexio_determinants

The ``import_integrals_ao`` and ``import_integrals_mo`` programs, from a
development plugin, read plain-text integral files instead.


Writing files for other codes
-----------------------------

The :ref:`molden` program writes the current |MOs| in the Molden format,
in the file :file:`$EZFIO.mol`, for visualization with `Molden
<http://cheminf.cmbi.ru.nl/molden>`_ or any program that reads this
format.

The :ref:`fcidump` program writes the Hamiltonian in the
:file:`$EZFIO.FCIDUMP` file, the format used by `NECI`_, `Dice`_ and
other |FCI| - |QMC| codes. Only the active |MOs| are written, the core
orbitals contribute to the core energy stored in the last line, and the
orbital classes set with :ref:`qp_set_mo_class` define the active space.
The ``fcidump_pyscf`` program writes the same information with the
``&end`` terminator that the PySCF reader expects, in the same file name,
so the two programs are alternative choices rather than successive steps.

The ``export_trexio`` program writes a TREXIO file with the geometry, the
basis set and the |MOs|, and optionally the integrals and the two-body
density matrix:

.. code:: bash

    qp set trexio trexio_file h2o.h5
    qp run export_trexio

The format of the file is HDF5 by default and plain text when the
``backend`` option is set to 1, and the ``export_*`` options of the same
module select the quantities to write. For a multi-state calculation,
one file is produced per state: the first state goes in the named file,
and the following states in files suffixed ``.001``, ``.002``, and so
on, counting the first state as zero.

The ``export_integrals_ao`` and ``export_integrals_mo`` programs, from the
same development plugin as their import counterparts, write the integrals
in a plain-text format.

For a human-readable wave function, :ref:`print_wf` prints the
determinants in intermediate normalization, with the excitation
operators relative to the reference determinant and a first-order
perturbative estimate of each coefficient; :ref:`print_ci_vectors`
prints the same determinants and coefficients without the perturbative
analysis.
:ref:`save_one_e_dm` stores the one-body density matrix of the alpha and
beta electrons, in the |AO| and |MO| bases, in the ``aux_quantities``
section of the database.


Quantum Monte Carlo
-------------------

Trial wave functions for |QMC| are produced by three plugin programs.
``save_for_qmcchem`` writes the spin determinants and the
pseudopotential in the layout expected by `QMC=Chem`_, and stores the
reference energy taken from the |CIPSI| calculation.
``save_for_qmcpack`` prepares the wave function, the |MOs| and the
pseudopotential for `QMCPack`_, and ``save_for_champ`` does the same for
`CHAMP <https://www.utwente.nl/en/tnw/ccp/research/CHAMP.html>`_.
The three programs are installed with the |qp| plugins of the
corresponding codes, described in the `external plugins`_ repository, and
they must be run after the |CIPSI| calculation has converged:

.. code:: bash

    qp run save_for_qmcchem

A few auxiliary programs complete these interfaces.
``truncate_wf_qmcchem`` removes the alpha and beta strings whose weight is
below the ``qmcchem/ci_threshold`` option, 10\ :sup:`-8` by default, and
saves the truncated wave function. ``expand_spindets_qmcchem`` generates
all the alpha-beta determinant products of the stored wave function and
diagonalizes it again, which is the dense representation expected by some
|QMC| codes, and ``qmc_create_wf`` does the same for `QMCPack`_.
``qmc_e_curve`` prints the energy of the wave function truncated at a
series of thresholds, which helps choosing the cutoff of the |QMC|
calculation.

.. seealso::

    The :ref:`qp_plugins` command, which installs and removes the plugins,
    and the documentation of the :ref:`molden`, :ref:`fcidump`,
    :ref:`print_wf` and :ref:`print_ci_vectors` programs.
