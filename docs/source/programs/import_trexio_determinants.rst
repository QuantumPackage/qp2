.. _import_trexio_determinants: 
 
.. program:: import_trexio_determinants 
 
========================== 
import_trexio_determinants 
========================== 
 
 
 
 
 Reads the determinants and coefficients of a wave 
 function from the TREXIO file named by the 
 trexio_file option and saves them in the |EZFIO| 
 database, which must already contain the geometry 
 and basis set. Stops if the file declares a number 
 of states. 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`run` 
