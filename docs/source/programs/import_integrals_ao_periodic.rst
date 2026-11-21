.. _import_integrals_ao_periodic: 
 
.. program:: import_integrals_ao_periodic 
 
============================ 
import_integrals_ao_periodic 
============================ 
 
 
 
 
 Asks for the number of AOs, then imports the one- 
 electron integrals of a periodic calculation into the 
 |EZFIO| database from text files in the current 
 directory: kinetic (T.qp), overlap (S.qp), 
 pseudopotential (P.qp), nucleus-electron (V.qp), plus 
 the nuclear repulsion from the first value of E.qp. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`ao_num` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`run` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`ao_num` 
