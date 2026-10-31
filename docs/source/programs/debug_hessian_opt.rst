.. _debug_hessian_opt: 
 
.. program:: debug_hessian_opt 
 
================= 
debug_hessian_opt 
================= 
 
 
 
 
 Debug program for the orbital hessian, for all mo_num orbitals 
 (n = mo_num*(mo_num-1)/2 rotation pairs), in the 2D and 4D 
 forms. For optimization_method 'full' or 'diag', compares the 
 current implementation (hessian_opt or diag_hessian_opt) with the 
 original one (first_hessian_opt or first_diag_hessian_opt) and 
 lists the elements differing by more than 1d-12. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`all_mo_integrals` 
    * :c:data:`mo_num` 
    * :c:data:`optimization_method` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`abort` 
    * :c:func:`diag_hessian_opt` 
    * :c:func:`first_diag_hessian_opt` 
    * :c:func:`first_hessian_opt` 
    * :c:func:`hessian_opt` 
