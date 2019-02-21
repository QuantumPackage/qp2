============
hartree_fock
============


Quick description
-----------------

The :ref:`scf` program performs *Restricted* Hartree-Fock
calculations (the spatial part of the |MOs| is common for alpha and beta
spinorbitals).

.. seealso:: 
   To see the keywords/options associated to the |SCF| algorithm itself,  
   see the documentation of the :ref:`module_scf_utils` module. 


More advanced description
-------------------------

The Hartree-Fock algorithm is a |SCF| and therefore is based on the
:ref:`module_scf_utils` module. 

The Fock matrix is defined in :file:`fock_matrix_hf.irp.f`.


.. seealso:: 
   For a more detailed description of the |SCF| structure, 
   see the documentation of the :ref:`module_scf_utils` module. 


