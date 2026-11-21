.. _debug_gradient_list_opt: 
 
.. program:: debug_gradient_list_opt 
 
======================= 
debug_gradient_list_opt 
======================= 
 
 
 
 
 Debug program for the orbital gradient, same check as 
 debug_gradient_opt but restricted to the active orbitals list_act 
 (m = dim_list_act_orb orbitals, n = m*(m-1)/2 rotation pairs), 
 using gradient_list_opt and first_gradient_list_opt. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`all_mo_integrals` 
    * :c:data:`dim_list_act_orb` 
    * :c:data:`list_act` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`diagonalize_ci` 
    * :c:func:`first_gradient_list_opt` 
    * :c:func:`gradient_list_opt` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`ci_electronic_energy` 
    * :c:data:`ci_electronic_energy` 
    * :c:data:`ci_energy` 
    * :c:data:`ci_electronic_energy` 
    * :c:data:`psi_coef` 
    * :c:data:`psi_energy` 
    * :c:data:`psi_energy` 
