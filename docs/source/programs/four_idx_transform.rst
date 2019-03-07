.. _four_idx_transform: 
 
.. program:: four_idx_transform 
 
================== 
four_idx_transform 
================== 
 
 
 
 
 4-index transformation of two-electron integrals from |AO| to |MO| 
 integrals. 
  
 This program will compute the two-electron integrals on the |MO| basis 
 and store it into the |EZFIO| directory. 
  
 This program can be useful if the AO --> MO transformation is an 
 expensive step by itself. 
  
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`io_mo_two_e_integrals` 
    * :c:data:`mo_two_e_integrals_in_map` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`io_mo_two_e_integrals` 
