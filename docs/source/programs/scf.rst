.. _scf: 
 
.. program:: scf 
 
=== 
scf 
=== 
 
 
  
 output: hartree_fock.energy 
  
 optional: mo_basis.mo_coef 
 
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
    * :c:data:`fock_matrix_ao_alpha` 
    * :c:data:`fock_matrix_ao_alpha` 
    * :c:data:`mo_coef` 
    * :c:data:`level_shift` 
    * :c:data:`mo_coef` 
    * :c:data:`mo_label` 
