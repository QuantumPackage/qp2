.. _save_natorb: 
 
.. program:: save_natorb 
 
=========== 
save_natorb 
=========== 
 
 
 
 
 Save natural |MOs| into the |EZFIO|. 
  
 This program reads the wave function stored in the |EZFIO| directory, 
 extracts the corresponding natural orbitals and setd them as the new 
 |MOs|. 
  
 If this is a multi-state calculation, the density matrix that produces 
 the natural orbitals is obtained from an average of the density 
 matrices of each state with the corresponding 
 :option:`determinants state_average_weight` 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`ezfio_set_mo_one_e_ints_io_mo_integrals_kinetic` 
    * :c:func:`ezfio_set_mo_one_e_ints_io_mo_integrals_n_e` 
    * :c:func:`ezfio_set_mo_one_e_ints_io_mo_integrals_pseudo` 
    * :c:func:`ezfio_set_mo_one_e_ints_io_mo_one_e_integrals` 
    * :c:func:`ezfio_set_mo_two_e_ints_io_mo_two_e_integrals` 
    * :c:func:`save_natural_mos` 
    * :c:func:`save_ref_determinant` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`mo_coef` 
    * :c:data:`mo_occ` 
    * :c:data:`read_wf` 
