.. _module_kohn_sham_rs: 
 
.. program:: kohn_sham_rs 
 
.. default-role:: option 
 
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


 
 
 
EZFIO parameters 
---------------- 
 
.. option:: energy
 
    Energy range separated hybrid
 
 
 
Programs 
-------- 
 
 * :ref:`rs_ks_scf` 
 
Providers 
--------- 
 
.. c:var:: rs_ks_energy


    File : :file:`rs_ks_energy.irp.f`

    .. code:: fortran

        double precision	:: rs_ks_energy	
        double precision	:: two_e_energy	
        double precision	:: one_e_energy	
        double precision	:: fock_matrix_energy	
        double precision	:: trace_potential_xc	


    Range-separated Kohn-Sham energy containing the nuclear repulsion energy, and the various components of this quantity.

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
