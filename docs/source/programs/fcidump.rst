.. _fcidump: 
 
.. program:: fcidump 
 
======= 
fcidump 
======= 
 
 
 
 
 Produce a regular `FCIDUMP` file from the |MOs| stored in the |EZFIO| 
 directory. 
  
 To specify an active space, the class of the |MOs| have to set in the 
 |EZFIO| directory (see :ref:`qp_set_mo_class`). 
  
 The :ref:`fcidump` program supports 3 types of |MO| classes : 
  
 * the *core* orbitals which are always doubly occupied in the 
   calculation 
  
 * the *deleted* orbitals that are never occupied in the calculation 
  
 * the *active* orbitals that are occupied with a varying number of 
   electrons 
  
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`core_energy` 
    * :c:data:`core_fock_operator` 
    * :c:data:`elec_alpha_num` 
    * :c:data:`elec_beta_num` 
    * :c:data:`elec_num` 
    * :c:data:`ezfio_filename` 
    * :c:data:`list_act` 
    * :c:data:`mo_integrals_map` 
    * :c:data:`mo_integrals_threshold` 
    * :c:data:`mo_one_e_integrals` 
    * :c:data:`mo_two_e_integrals_in_map` 
    * :c:data:`n_act_orb` 
    * :c:data:`n_core_orb` 
