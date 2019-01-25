============
kohn_sham_rs
============


The Range-separated Kohn-Sham module performs *Restricted* Kohn-Sham calculations (the
spatial part of the |MOs| is common for alpha and beta spinorbitals) where the coulomb interaction is partially treated using exact exchange.
The splitting of the interaction between long- and short-range is determined by the range-separation parameter :option:`ao_two_e_erf_ints mu_erf`. The long-range part of the interaction is explicitly treated with exact exchange, and the short-range part of the interaction is treated with appropriate DFT functionals.

The Range-separated Kohn-Sham in an SCF and therefore is based on the ``scf_utils`` structure.
It performs the following actions:

#. Compute/Read all the one- and two-electron integrals, and store them in memory
#. Check in the |EZFIO| database if there is a set of |MOs|. If there is, it
   will read them as initial guess. Otherwise, it will create a guess.
#. Perform the |SCF| iterations

The definition of the Fock matrix is in :file:`kohn_sham_rs fock_matrix_rs_ks.irp.f`
For the keywords related to the |SCF| procedure, see the ``scf_utils`` directory where you will find all options.
The main are:
# :option:`scf_utils thresh_scf`
# :option:`scf_utils level_shift`


At each iteration, the |MOs| are saved in the |EZFIO| database. Hence, if the calculation
crashes for any unexpected reason, the calculation can be restarted by running again
the |SCF| with the same |EZFIO| database.

The `DIIS`_ algorithm is implemented, as well as the `level-shifting`_ method.
If the |SCF| does not converge, try again with a higher value of :option:`level_shift`.

To start a calculation from scratch, the simplest way is to remove the
``mo_basis`` directory from the |EZFIO| database, and run the |SCF| again.


.. _DIIS: https://en.wikipedia.org/w/index.php?title=DIIS
.. _level-shifting: https://doi.org/10.1002/qua.560070407



