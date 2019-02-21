=========
kohn_sham
=========

Quick description
-----------------

The Kohn-Sham module performs *Restricted* Kohn-Sham calculations (the
spatial part of the |MOs| is common for alpha and beta spinorbitals). 
The program associated to it is the :ref:`ks_scf` executable. 

.. seealso:: 
 
     The documentation of the :ref:`module_dft_keywords` module for the various keywords 
     such as the exchange/correlation functionals or the amount of |HF| exchange. 


.. seealso:: 
   To see the keywords/options associated to the |SCF| algorithm itself,  
   see the documentation of the :ref:`module_scf_utils` module. 


More advanced description
-------------------------

The Kohn-Sham in an SCF and therefore is based on the :ref:`module_scf_utils` structure.

The definition of the Fock matrix is in :file:`kohn_sham fock_matrix_ks.irp.f`


.. seealso:: 
   For a more detailed description of the |SCF| structure, 
   see the documentation of the :ref:`module_scf_utils` module. 


