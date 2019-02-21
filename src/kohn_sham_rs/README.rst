============
kohn_sham_rs
============


Quick description
-----------------

The Range-separated Kohn-Sham module performs *Restricted* range-separated Hybrid calculation, 
which means that only the long-range part of the *exact* exchange is taken into account. 

The program associated to it is the :ref:`rs_ks_scf` executable. 

.. seealso:: 
 
     The documentation of the :ref:`module_dft_keywords` module for the various keywords 
     such as the exchange/correlation functionals or the range-separation parameter. 


.. seealso:: 
   To see the keywords/options associated to the |SCF| algorithm itself,  
   see the documentation of the :ref:`module_scf_utils` module. 


More advanced description
-------------------------

The splitting of the interaction between long- and short-range is determined by the range-separation parameter :option:`ao_two_e_erf_ints mu_erf`. The long-range part of the interaction is explicitly treated with exact exchange, and the short-range part of the interaction is treated with appropriate DFT functionals.

The Range-separated Kohn-Sham in an SCF and therefore is based on the :ref:`module_scf_utils` structure.

The definition of the Fock matrix is in :file:`kohn_sham_rs fock_matrix_rs_ks.irp.f`


.. seealso:: 
   For a more detailed description of the |SCF| structure, 
   see the documentation of the :ref:`module_scf_utils` module. 


