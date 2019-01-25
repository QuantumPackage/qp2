.. _fcidump: 
 
.. program:: fcidump 
 
======= 
fcidump 
======= 
 
 
 
 
 Produce a regular FCIDUMP file from the |MOs| stored in the |EZFIO| folder. 
  
 To specify an active space, the class of the mos have to set in the |EZFIO| folder (see :ref:`qp_set_mo_class`). 
  
 The fcidump program supports 3 types of MO_class : 
  
 * the "core" orbitals which are always doubly occupied in the calculation 
  
 * the "del" orbitals that are never occupied in the calculation 
  
 * the "act" orbitals that will be occupied by a varying number of electrons 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`elec_beta_num` 
    * :c:data:`ezfio_filename` 
    * :c:data:`core_fock_operator` 
    * :c:data:`elec_num` 
    * :c:data:`mo_two_e_integrals_in_map` 
    * :c:data:`elec_alpha_num` 
    * :c:data:`mo_one_e_integrals` 
    * :c:data:`n_core_orb` 
    * :c:data:`mo_integrals_threshold` 
    * :c:data:`list_inact` 
    * :c:data:`mo_integrals_map` 
    * :c:data:`core_energy` 
