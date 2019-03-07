.. _cis: 
 
.. program:: cis 
 
=== 
cis 
=== 
 
 
 
 
  
 Configuration Interaction with Single excitations. 
  
 This program takes a reference Slater determinant of ROHF-like 
 occupancy, and performs all single excitations on top of it. 
 Disregarding spatial symmetry, it computes the `n_states` lowest 
 eigenstates of that CI matrix. (see :option:`determinants n_states`) 
  
 This program can be useful in many cases: 
  
  
 1. Ground state calculation 
  
    To be sure to have the lowest |SCF| solution, perform an :ref:`scf` 
    (see the :ref:`module_hartree_fock` module), then a :ref:`cis`, save the 
    natural orbitals (see :ref:`save_natorb`) and re-run an :ref:`scf` 
    optimization from this |MO| guess. 
  
  
 2. Excited states calculations 
  
    The lowest excited states are much likely to be dominated by 
    single-excitations. Therefore, running a :ref:`cis` will save the 
    `n_states` lowest states within the |CIS| space in the |EZFIO| 
    directory, which can afterwards be used as guess wave functions for 
    a further multi-state |FCI| calculation if :option:`determinants 
    read_wf` is set to |true| before running the :ref:`fci` executable. 
  
  
 If :option:`determinants s2_eig` is set to |true|, the |CIS| 
 will only retain states having the expected |S^2| value (see 
 :option:`determinants expected_s2`). Otherwise, the |CIS| will take 
 the lowest :option:`determinants n_states`, whatever multiplicity 
 they are. 
  
 .. note:: 
  
   To discard some orbitals, use the :ref:`qp_set_mo_class` 
   command to specify: 
  
   * *core* orbitals which will be always doubly occupied 
  
   * *act* orbitals where an electron can be either excited from or to 
  
   * *del* orbitals which will be never occupied 
  
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`run` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`fock_matrix_ao_alpha` 
    * :c:data:`fock_matrix_ao_alpha` 
    * :c:data:`mo_coef` 
    * :c:data:`level_shift` 
    * :c:data:`mo_coef` 
    * :c:data:`read_wf` 
