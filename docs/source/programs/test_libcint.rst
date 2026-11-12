.. _test_libcint: 
 
.. program:: test_libcint 
 
============ 
test_libcint 
============ 
 
 
 
 
 Print every |AO| integral produced by the active engine, with tags, so 
 that two runs differing only by `ao_basis.use_libcint` can be diffed. 
  
 All integral lines are of the form 
  
     TAG i j value        (one-electron) 
     TAG i j k l value    (two-electron) 
  
 and the accompanying `INFO` lines are not compared. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`ao_integrals_n_e` 
    * :c:data:`ao_kinetic_integrals` 
    * :c:data:`ao_l` 
    * :c:data:`ao_l` 
    * :c:data:`ao_num` 
    * :c:data:`ao_overlap` 
    * :c:data:`ao_overlap_qp2` 
    * :c:data:`libcint_available` 
    * :c:data:`libcint_nuc_available` 
    * :c:data:`mu_erf` 
    * :c:data:`nucl_num` 
    * :c:data:`shell_num` 
    * :c:data:`use_libcint` 
    * :c:data:`use_libcint_effective` 
