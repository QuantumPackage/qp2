.. _module_kohn_sham: 
 
.. program:: kohn_sham 
 
.. default-role:: option 
 
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


 
 
 
Programs 
-------- 
 
 * :ref:`ks_scf` 
 
Providers 
--------- 
 
.. c:var:: ks_energy


    File : :file:`ks_enery.irp.f`

    .. code:: fortran

        double precision	:: ks_energy	
        double precision	:: two_e_energy	
        double precision	:: one_e_energy	
        double precision	:: fock_matrix_energy	
        double precision	:: trace_potential_xc	


    Kohn-Sham energy containing the nuclear repulsion energy, and the various components of this quantity.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`e_correlation_dft`
       * :c:data:`e_exchange_dft`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`nuclear_repulsion`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`extra_e_contrib_density`

 
 
Subroutines / functions 
----------------------- 
