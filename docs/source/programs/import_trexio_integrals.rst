.. _import_trexio_integrals: 
 
.. program:: import_trexio_integrals 
 
======================= 
import_trexio_integrals 
======================= 
 
 
 
 
 Reads the nuclear repulsion and the one- and 
 two-electron integrals (AO and MO) from the TREXIO 
 file named by the trexio_file option and stores 
 them in the |EZFIO| database. The two-electron 
 integrals may come as Cholesky vectors or as the 
 list of nonzero integrals. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`trexio_filename` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`run` 
    * :c:func:`trexio_assert` 
