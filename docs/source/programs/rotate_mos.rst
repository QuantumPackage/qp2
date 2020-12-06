.. _rotate_mos: 
 
.. program:: rotate_mos 
 
========== 
rotate_mos 
========== 
 
 
 
 
 Rotates molecular orbitals i and j by combining them as 
 $1/\sqrt{2} ( \phi_i + \phi_j )$ and 
 $1/\sqrt{2} ( \phi_i - \phi_j )$. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`ao_num` 
    * :c:data:`mo_coef` 
    * :c:data:`mo_num` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`save_mos` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`mo_coef` 
