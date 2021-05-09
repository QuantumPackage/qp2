.. _sort_by_fock_energies: 
 
.. program:: sort_by_fock_energies 
 
===================== 
sort_by_fock_energies 
===================== 
 
 
 
 
 Program that saves the current |MOs| ordered by diagonal element of the Fock operator. 
  
 Warning : the Fock operator, and therefore its matrix elements, depends on the occupancy. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`ao_num` 
    * :c:data:`fock_matrix_mo` 
    * :c:data:`mo_coef` 
    * :c:data:`mo_num` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`dsort` 
    * :c:func:`save_mos` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`mo_coef` 
