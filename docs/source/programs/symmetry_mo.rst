.. _symmetry_mo: 
 
.. program:: symmetry_mo 
 
=========== 
symmetry_mo 
=========== 
 
 
 
 
 Step-2 driver of the symmetry module. 
  
 Reads the nuclei from the EZFIO directory given on the command line 
 (through the standard :c:data:`ezfio_filename` provider), builds the canonical 
 frame, enumerates the point group, names it, prints the classes, and 
 prints the character table shipped in ``QP_ROOT/data/Symmetry``. 
  
 The EZFIO directory is the first command-line argument, as for every 
 other QP2 program. 
  
 An optional environment variable ``PG_WRITE_MO_SYMMETRY=1`` writes the 
 labels into ``mo_basis/mo_symmetry``. Nothing is 
 written when a gate fails, so a written field is always a checked one. 
 An environment variable is used because :c:data:`ezfio_filename` stops unless the 
 command line has exactly one argument. 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`getenv` 
    * :c:func:`pg_report` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`mo_symmetry` 
