.. _scf: 
 
.. program:: scf 
 
=== 
scf 
=== 
 
 
 
 
  
 The :ref:`scf` program performs *Restricted* Hartree-Fock 
 calculations (the spatial part of the |MOs| is common for alpha and beta 
 spinorbitals). 
  
 It performs the following actions: 
  
 #. Compute/Read all the one- and two-electron integrals, and store them 
    in memory 
 #. Check in the |EZFIO| database if there is a set of |MOs|. 
    If there is, it will read them as initial guess. Otherwise, it will 
    create a guess. 
 #. Perform the |SCF| iterations 
  
 For the keywords related to the |SCF| procedure, see the ``scf_utils`` 
 directory where you will find all options. 
  
 At each iteration, the |MOs| are saved in the |EZFIO| database. Hence, 
 if the calculation crashes for any unexpected reason, the calculation 
 can be restarted by running again the |SCF| with the same |EZFIO| 
 database. 
  
 To start again a fresh |SCF| calculation, the |MOs| can be reset by 
 running the :ref:`qp_reset` command. 
  
 The `DIIS`_ algorithm is implemented, as well as the `level-shifting`_ 
 method. If the |SCF| does not converge, try again with a higher value of 
 :option:`level_shift`. 
  
 .. _DIIS: https://en.wikipedia.org/w/index.php?title=DIIS 
 .. _level-shifting: https://doi.org/10.1002/qua.560070407 
  
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`create_guess` 
    * :c:func:`orthonormalize_mos` 
    * :c:func:`run` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`fock_matrix_ao_alpha` 
    * :c:data:`fock_matrix_ao_alpha` 
    * :c:data:`mo_coef` 
    * :c:data:`mo_label` 
