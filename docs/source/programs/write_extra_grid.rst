.. _write_extra_grid: 
 
.. program:: write_extra_grid 
 
================ 
write_extra_grid 
================ 
 
 
 
 
 Program that writes the extra grid to the EZFIO file: the 
 number of points, the points and the weights. It sets 
 becke_numerical_grid/io_extra_grid to "Write", forces the 
 computation of final_grid_points_extra, and resets 
 io_extra_grid to "Read" so that later calculations read 
 the grid back from disk. 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`io_extra_grid` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`routine` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`io_extra_grid` 
