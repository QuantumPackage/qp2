.. _fcidump: 
 
.. program:: fcidump 
 
======= 
fcidump 
======= 
 
 
  
 The fcidump program supports 3 types of MO_class : 
  
 * the "core" orbitals which are always doubly occupied in the calculation 
  
 * the "del" orbitals that are never occupied in the calculation 
  
 * the "act" orbitals that will be occupied by a varying number of electrons 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`elec_beta_num` 
    * :c:data:`list_act` 
    * :c:data:`ezfio_filename` 
    * :c:data:`core_fock_operator` 
    * :c:data:`core_bitmask` 
    * :c:data:`elec_num` 
    * :c:data:`mo_two_e_integrals_in_map` 
    * :c:data:`elec_alpha_num` 
    * :c:data:`mo_one_e_integrals` 
    * :c:data:`n_act_orb` 
    * :c:data:`mo_integrals_threshold` 
    * :c:data:`mo_integrals_map` 
    * :c:data:`core_energy` 
