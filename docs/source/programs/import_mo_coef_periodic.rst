.. _import_mo_coef_periodic: 
 
.. program:: import_mo_coef_periodic 
 
======================= 
import_mo_coef_periodic 
======================= 
 
 
 
 
 Imports the MO coefficients (real and imaginary parts) 
 from the text file C.qp in the current directory into 
 the |EZFIO| database named by ezfio_filename. 
 TODO(doc): mo_num, mo_occ, mo_class and ao_md5 are not 
 set in this program; are they read from the database by 
 save_mos? 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`ezfio_filename` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`run` 
