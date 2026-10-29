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

The splitting of the interaction between long- and short-range is determined by the range-separation parameter :option:`mu_erf`. The long-range part of the interaction is explicitly treated with exact exchange, and the short-range part of the interaction is treated with appropriate DFT functionals.

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
 
.. c:var:: double ao_potential_alpha_xc


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
       * :c:data:`rs_ks_energy`

 
.. c:var:: double ao_potential_beta_xc


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
       * :c:data:`rs_ks_energy`

 
.. c:var:: double ao_two_e_integral_alpha


    File : :file:`fock_matrix_rs_ks.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_integral_beta	(ao_num,ao_num)


    Alpha Fock matrix in ao basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_alpha_no_xc_ao`
       * :c:data:`rs_ks_energy`

 
.. c:var:: double ao_two_e_integral_beta


    File : :file:`fock_matrix_rs_ks.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_integral_alpha	(ao_num,ao_num)
        double precision, allocatable	:: ao_two_e_integral_beta	(ao_num,ao_num)


    Alpha Fock matrix in ao basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_alpha_no_xc_ao`
       * :c:data:`rs_ks_energy`

 
.. c:var:: double e_correlation_dft


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
       * :c:data:`rs_ks_energy`

 
.. c:var:: double e_exchange_dft


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
       * :c:data:`rs_ks_energy`

 
.. c:var:: double extra_e_contrib_density


    File : :file:`rs_ks_energy.irp.f`

    .. code:: fortran

        double precision	:: extra_e_contrib_density	


    Extra contribution to the SCF energy coming from the density.
    
    For a Hartree-Fock calculation: extra_e_contrib_density = 0
    
    For a Kohn-Sham or Range-separated Kohn-Sham: the exchange/correlation - 1/2 trace of the V_xc potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_correlation_dft`
       * :c:data:`e_exchange_dft`
       * :c:data:`rs_ks_energy`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`scf_energy`

 
.. c:var:: double fock_matrix_alpha_no_xc_ao


    File : :file:`fock_matrix_rs_ks.irp.f`

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

 
.. c:var:: double fock_matrix_ao_alpha


    File : :file:`fock_matrix_rs_ks.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_ao_alpha	(ao_num,ao_num)
        double precision, allocatable	:: fock_matrix_ao_beta	(ao_num,ao_num)


    Alpha Fock matrix in ao basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`fock_matrix_alpha_no_xc_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao`
       * :c:data:`fock_matrix_mo_alpha`
       * :c:data:`fock_matrix_mo_beta`
       * :c:data:`rs_ks_energy`
       * :c:data:`scf_energy`

 
.. c:var:: double fock_matrix_ao_beta


    File : :file:`fock_matrix_rs_ks.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_ao_alpha	(ao_num,ao_num)
        double precision, allocatable	:: fock_matrix_ao_beta	(ao_num,ao_num)


    Alpha Fock matrix in ao basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`fock_matrix_alpha_no_xc_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao`
       * :c:data:`fock_matrix_mo_alpha`
       * :c:data:`fock_matrix_mo_beta`
       * :c:data:`rs_ks_energy`
       * :c:data:`scf_energy`

 
.. c:var:: double fock_matrix_beta_no_xc_ao


    File : :file:`fock_matrix_rs_ks.irp.f`

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

 
.. c:var:: double fock_matrix_energy


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

 
.. c:var:: double one_e_energy


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

 
.. c:var:: double rs_ks_energy


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

 
.. c:var:: double trace_potential_xc


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

 
.. c:var:: double two_e_energy


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
 
.. c:function:: void check_coherence_functional(void)


    File : :file:`rs_ks_scf.irp.f`

    .. code:: fortran

        subroutine check_coherence_functional



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`

 
.. c:function:: void create_guess(void)


    File : :file:`rs_ks_scf.irp.f`

    .. code:: fortran

        subroutine create_guess


    Create a MO guess if no MOs are present in the EZFIO directory

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_lowdin_coef`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_coef`
       * :c:data:`mo_guess_type`
       * :c:data:`mo_label`
       * :c:data:`mo_one_e_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`rs_ks_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_has_mo_basis_mo_coef`
       * :c:func:`huckel_guess`
       * :c:func:`mo_as_eigvectors_of_mo_matrix`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`

 
.. c:function:: void run(void)


    File : :file:`rs_ks_scf.irp.f`

    .. code:: fortran

        subroutine run


    Run SCF calculation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`e_correlation_dft`
       * :c:data:`e_exchange_dft`
       * :c:data:`json_unit`
       * :c:data:`level_shift`
       * :c:data:`mo_label`
       * :c:data:`rs_ks_energy`
       * :c:data:`scf_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`rs_ks_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_kohn_sham_rs_energy`
       * :c:func:`json_close`
       * :c:func:`roothaan_hall_scf`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`level_shift`
       * :c:data:`mo_coef`

