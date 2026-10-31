.. _debug_gradient_opt: 
 
.. program:: debug_gradient_opt 
 
================== 
debug_gradient_opt 
================== 
 
 
 
 
 Debug program for the orbital gradient, for all mo_num orbitals 
 (n = mo_num*(mo_num-1)/2 rotation pairs). Compares the gradient 
 of the energy with respect to orbital rotations computed by the 
 current implementation (gradient_opt) with the original one 
 (first_gradient_opt), and lists the elements whose difference 
 exceeds the 1d-12 threshold. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`all_mo_integrals` 
    * :c:data:`mo_num` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`diagonalize_ci` 
    * :c:func:`first_gradient_opt` 
    * :c:func:`gradient_opt` 
 
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
