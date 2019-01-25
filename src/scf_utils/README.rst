=========
scf_utils
=========



The scf_utils module is an abstract module which contains the basics to perform *Restricted* SCF calculations (the
spatial part of the |MOs| is common for alpha and beta spinorbitals) based on a single-determinant wave function.

This module does not produce any executable *and must not do*, but instead it contains everything one needs to perform an orbital optimization based on an Fock matrix.
The ``scf_utils`` module is meant to be included in the :file:`NEED` of the various single determinant SCF procedures, such as ``hartree_fock`` or ``kohn_sham``, where a specific definition of the Fock matrix is given (see :file:`hartree_fock fock_matrix_hf.irp.f` for an example).

All SCF programs perform the following actions:


#. Compute/Read all the one- and two-electron integrals, and store them in memory

#. Check in the |EZFIO| database if there is a set of |MOs|. If there is, it
   will read them as initial guess. Otherwise, it will create a guess.
#. Perform the |SCF| iterations based on the definition of the Fock matrix


The main keywords/options are:

* :option:`scf_utils thresh_scf`
* :option:`scf_utils level_shift`

At each iteration, the |MOs| are saved in the |EZFIO| database. Hence, if the calculation
crashes for any unexpected reason, the calculation can be restarted by running again
the |SCF| with the same |EZFIO| database.

The `DIIS`_ algorithm is implemented, as well as the `level-shifting`_ method.
If the |SCF| does not converge, try again with a higher value of :option:`level_shift`.

To start a calculation from scratch, the simplest way is to remove the
``mo_basis`` directory from the |EZFIO| database, and run the |SCF| again.

.. _DIIS: https://en.wikipedia.org/w/index.php?title=DIIS
.. _level-shifting: https://doi.org/10.1002/qua.560070407

