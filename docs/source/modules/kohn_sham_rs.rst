.. _module_kohn_sham_rs: 
 
.. program:: kohn_sham_rs 
 
.. default-role:: option 
 
============
kohn_sham_rs
============


The Range-separated Kohn-Sham module performs *Restricted* Kohn-Sham calculations (the
spatial part of the |MOs| is common for alpha and beta spinorbitals) where the coulomb interaction is partially treated using exact exchange.
The splitting of the interaction between long- and short-range is determined by the range-separation parameter :option:`ao_two_e_erf_ints mu_erf`. The long-range part of the interaction is explicitly treated with exact exchange, and the short-range part of the interaction is treated with appropriate DFT functionals.

The Range-separated Kohn-Sham in an SCF and therefore is based on the ``scf_utils`` structure.
It performs the following actions:

#. Compute/Read all the one- and two-electron integrals, and store them in memory
#. Check in the |EZFIO| database if there is a set of |MOs|. If there is, it
   will read them as initial guess. Otherwise, it will create a guess.
#. Perform the |SCF| iterations

The definition of the Fock matrix is in :file:`kohn_sham_rs fock_matrix_rs_ks.irp.f`
For the keywords related to the |SCF| procedure, see the ``scf_utils`` directory where you will find all options.
The main are:
# :option:`scf_utils thresh_scf`
# :option:`scf_utils level_shift`


At each iteration, the |MOs| are saved in the |EZFIO| database. Hence, if the calculation
crashes for any unexpected reason, the calculation can be restarted by running again
the |SCF| with the same |EZFIO| database.

The `DIIS`_ algorithm is implemented, as well as the `level-shifting`_ method.
If the |SCF| does not converge, try again with a higher value of :option:`level_shift`.

To start a calculation from scratch, the simplest way is to remove the
``mo_basis`` directory from the |EZFIO| database, and run the |SCF| again.


.. _DIIS: https://en.wikipedia.org/w/index.php?title=DIIS
.. _level-shifting: https://doi.org/10.1002/qua.560070407



 
 
 
EZFIO parameters 
---------------- 
 
.. option:: energy
 
    Energy range separated hybrid
 
 
 
Programs 
-------- 
 
 * :ref:`rs_ks_scf` 
 
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
       * :c:data:`potential_x_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`rs_ks_energy`

 
.. c:var:: ao_potential_beta_xc


    File : :file:`pot_functionals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_potential_alpha_xc	(ao_num,ao_num)
        double precision, allocatable	:: ao_potential_beta_xc	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`potential_x_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`rs_ks_energy`

 
.. c:var:: e_correlation_dft


    File : :file:`pot_functionals.irp.f`

    .. code:: fortran

        double precision	:: e_correlation_dft	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`extra_e_contrib_density`
       * :c:data:`rs_ks_energy`

 
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
       * :c:data:`rs_ks_energy`

 
.. c:var:: fock_matrix_alpha_no_xc_ao


    File : :file:`fock_matrix_rs_ks.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_alpha_no_xc_ao	(ao_num,ao_num)
        double precision, allocatable	:: fock_matrix_beta_no_xc_ao	(ao_num,ao_num)


    Mono electronic an Coulomb matrix in AO basis set

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


    File : :file:`fock_matrix_rs_ks.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_alpha_no_xc_ao	(ao_num,ao_num)
        double precision, allocatable	:: fock_matrix_beta_no_xc_ao	(ao_num,ao_num)


    Mono electronic an Coulomb matrix in AO basis set

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

 
.. c:var:: one_e_energy


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

 
.. c:var:: trace_potential_xc


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

 
.. c:var:: two_e_energy


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
 
.. c:function:: check_coherence_functional:


    File : :file:`rs_ks_scf.irp.f`

    .. code:: fortran

        subroutine check_coherence_functional



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`exchange_functional`
       * :c:data:`correlation_functional`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`rs_ks_scf`

