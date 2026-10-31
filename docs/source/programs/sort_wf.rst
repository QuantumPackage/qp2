.. _sort_wf: 
 
.. program:: sort_wf 
 
======= 
sort_wf 
======= 
 
 
 
 
 Writes to the file <ezfio_filename>.wf_sorted the amplitudes 
 |c_i/c_1| of the wave function sorted by contribution to the 
 norm, read from the |EZFIO| directory, one line per determinant. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`routine` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`read_wf` 
