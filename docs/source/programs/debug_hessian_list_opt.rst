.. _debug_hessian_list_opt: 
 
.. program:: debug_hessian_list_opt 
 
====================== 
debug_hessian_list_opt 
====================== 
 
 
 
 
 Debug program for the orbital hessian, restricted to the active 
 orbitals list_act (m = dim_list_act_orb orbitals, n = 
 m*(m-1)/2 rotation pairs). For optimization_method 'full' or 
 'diag', compares the hessian computed by the current code 
 (hessian_list_opt or diag_hessian_list_opt) with the original 
 implementation (first_hessian_list_opt or first_diag_hessian_list_opt). 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`all_mo_integrals` 
    * :c:data:`dim_list_act_orb` 
    * :c:data:`list_act` 
    * :c:data:`optimization_method` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`diag_hessian_list_opt` 
    * :c:func:`first_diag_hessian_list_opt` 
    * :c:func:`first_hessian_list_opt` 
    * :c:func:`hessian_list_opt` 
