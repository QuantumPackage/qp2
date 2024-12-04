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
 * :ref:`print_mos` 
 
Providers 
--------- 
 
.. c:var:: ao_potential_alpha_xc


    File : :file:`pot_functionals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_potential_alpha_xc	(ao_num,ao_num)
        double precision, allocatable	:: ao_potential_beta_xc	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`potential_c_alpha_ao`
       * :c:data:`potential_x_alpha_ao`
       * :c:data:`potential_xc_alpha_ao`
       * :c:data:`same_xc_func`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`ks_energy`

 
.. c:var:: ao_potential_beta_xc


    File : :file:`pot_functionals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_potential_alpha_xc	(ao_num,ao_num)
        double precision, allocatable	:: ao_potential_beta_xc	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`potential_c_alpha_ao`
       * :c:data:`potential_x_alpha_ao`
       * :c:data:`potential_xc_alpha_ao`
       * :c:data:`same_xc_func`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`ks_energy`

 
.. c:var:: e_correlation_dft


    File : :file:`pot_functionals.irp.f`

    .. code:: fortran

        double precision	:: e_correlation_dft	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`energy_c`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`extra_e_contrib_density`
       * :c:data:`ks_energy`

 
.. c:var:: e_exchange_dft


    File : :file:`pot_functionals.irp.f`

    .. code:: fortran

        double precision	:: e_exchange_dft	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`extra_e_contrib_density`
       * :c:data:`ks_energy`

 
.. c:var:: fock_matrix_alpha_no_xc_ao


    File : :file:`fock_matrix_ks.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_alpha_no_xc_ao	(ao_num,ao_num)
        double precision, allocatable	:: fock_matrix_beta_no_xc_ao	(ao_num,ao_num)


    Mono electronic an Coulomb matrix in ao basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_two_e_integral_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`

 
.. c:var:: fock_matrix_beta_no_xc_ao


    File : :file:`fock_matrix_ks.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_alpha_no_xc_ao	(ao_num,ao_num)
        double precision, allocatable	:: fock_matrix_beta_no_xc_ao	(ao_num,ao_num)


    Mono electronic an Coulomb matrix in ao basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_two_e_integral_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`

 
.. c:var:: fock_matrix_energy


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

 
.. c:var:: one_e_energy


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

 
.. c:var:: trace_potential_xc


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

 
.. c:var:: two_e_energy


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
 
.. c:function:: check_coherence_functional:


    File : :file:`ks_scf.irp.f`

    .. code:: fortran

        subroutine check_coherence_functional



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ks_scf`

