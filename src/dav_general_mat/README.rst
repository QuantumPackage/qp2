===============
dav_general_mat
===============

This modules allows to use the Davidson Algorithm for general squared symmetric matrices 
You have two options : 
 a)  the routine "davidson_general" to whom you pass the matrix you want to diagonalize
 b)  the routine "davidson_general_ext_rout" to whom you pass the subroutine that realizes v = H u 

See the routines in "test_dav.irp.f" for a clear example. 



