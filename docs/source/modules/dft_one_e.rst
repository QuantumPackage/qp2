.. _module_dft_one_e: 
 
.. program:: dft_one_e 
 
.. default-role:: option 
 
dft_one_e
=========

This module defines the most important providers needed for the |DFT| and |RSDFT| calculations: 

* :c:data:`energy_x` and :c:data:`energy_c` : the *exchange* and *correlation* energy functionals (see :file:`e_xc_general.irp.f`)

* :c:data:`potential_x_alpha_ao` and :c:data:`potential_x_beta_ao` : the exchange potential for alpha/beta electrons  (see :file:`pot_general.irp.f`)

* :c:data:`potential_c_alpha_ao` and :c:data:`potential_c_beta_ao` : the correlation potential for alpha/beta electrons (see :file:`pot_general.irp.f`)  


These providers are then used in the :ref:`ks_scf` and :ref:`rs_ks_scf` programs, together within some |RSDFT| external 
plugins (see `<https://gitlab.com/eginer/qp_plugins_eginer>`_). 

The flexibility of the functionals is handle by the two following keywords (see :ref:`module_dft_keywords`): 

* :option:`dft_keywords exchange_functional` : defines which *exchange* functionals will be set 

* :option:`dft_keywords correlation_functional` : defines which *correlation* functionals will be set 


In the core modules of the |QP|, two functionals are implemented: 

 * "LDA" or "short_range_LDA" for, respectively the |LDA| and its short-range version

 * "PBE" or "short_range_PBE" for, respectively the |PBE| and its short-range version


 
 
 
Providers 
--------- 
 
.. c:var:: ao_effective_one_e_potential


    File : :file:`dft_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_effective_one_e_potential	(ao_num,ao_num,N_states)
        double precision, allocatable	:: ao_effective_one_e_potential_without_kin	(ao_num,ao_num,N_states)


    Effective_one_e_potential(i,j) = :math:`\rangle i_{AO}| v_{H}^{sr} |j_{AO}\rangle  + \rangle i_{AO}| h_{core} |j_{AO}\rangle  + \rangle i_{AO}|v_{xc} |j_{AO}\rangle` 
    
    on the |MO| basis
    
    Taking the expectation value does not provide any energy, but
    
    ao_effective_one_e_potential(i,j) is the potential coupling DFT and WFT parts
    
    and it is used in any RS-DFT based calculations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`effective_one_e_potential`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`s_mo_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_effective_one_e_potential_sa`

 
.. c:var:: ao_effective_one_e_potential_sa


    File : :file:`dft_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_effective_one_e_potential_sa	(ao_num,ao_num)
        double precision, allocatable	:: ao_effective_one_e_potential_without_kin_sa	(ao_num,ao_num)


    State-averaged potential in AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_effective_one_e_potential`
       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`state_average_weight`


 
.. c:var:: ao_effective_one_e_potential_without_kin


    File : :file:`dft_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_effective_one_e_potential	(ao_num,ao_num,N_states)
        double precision, allocatable	:: ao_effective_one_e_potential_without_kin	(ao_num,ao_num,N_states)


    Effective_one_e_potential(i,j) = :math:`\rangle i_{AO}| v_{H}^{sr} |j_{AO}\rangle  + \rangle i_{AO}| h_{core} |j_{AO}\rangle  + \rangle i_{AO}|v_{xc} |j_{AO}\rangle` 
    
    on the |MO| basis
    
    Taking the expectation value does not provide any energy, but
    
    ao_effective_one_e_potential(i,j) is the potential coupling DFT and WFT parts
    
    and it is used in any RS-DFT based calculations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`effective_one_e_potential`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`s_mo_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_effective_one_e_potential_sa`

 
.. c:var:: ao_effective_one_e_potential_without_kin_sa


    File : :file:`dft_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_effective_one_e_potential_sa	(ao_num,ao_num)
        double precision, allocatable	:: ao_effective_one_e_potential_without_kin_sa	(ao_num,ao_num)


    State-averaged potential in AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_effective_one_e_potential`
       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`state_average_weight`


 
.. c:var:: effective_one_e_potential


    File : :file:`dft_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_one_e_potential	(mo_num,mo_num,N_states)
        double precision, allocatable	:: effective_one_e_potential_without_kin	(mo_num,mo_num,N_states)


    Effective_one_e_potential(i,j) = :math:`\rangle i_{MO}| v_{H}^{sr} |j_{MO}\rangle  + \rangle i_{MO}| h_{core} |j_{MO}\rangle  + \rangle i_{MO}|v_{xc} |j_{MO}\rangle` 
    
    on the |MO| basis
    
    Taking the expectation value does not provide any energy, but
    
    effective_one_e_potential(i,j) is the potential coupling DFT and WFT parts
    
    and it is used in any RS-DFT based calculations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_c_alpha_mo`
       * :c:data:`potential_x_alpha_mo`
       * :c:data:`short_range_hartree_operator`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_effective_one_e_potential`
       * :c:data:`effective_one_e_potential_sa`

 
.. c:var:: effective_one_e_potential_sa


    File : :file:`dft_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_one_e_potential_sa	(mo_num,mo_num)
        double precision, allocatable	:: effective_one_e_potential_without_kin_sa	(mo_num,mo_num)


    State-averaged potential in MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`state_average_weight`


 
.. c:var:: effective_one_e_potential_without_kin


    File : :file:`dft_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_one_e_potential	(mo_num,mo_num,N_states)
        double precision, allocatable	:: effective_one_e_potential_without_kin	(mo_num,mo_num,N_states)


    Effective_one_e_potential(i,j) = :math:`\rangle i_{MO}| v_{H}^{sr} |j_{MO}\rangle  + \rangle i_{MO}| h_{core} |j_{MO}\rangle  + \rangle i_{MO}|v_{xc} |j_{MO}\rangle` 
    
    on the |MO| basis
    
    Taking the expectation value does not provide any energy, but
    
    effective_one_e_potential(i,j) is the potential coupling DFT and WFT parts
    
    and it is used in any RS-DFT based calculations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_c_alpha_mo`
       * :c:data:`potential_x_alpha_mo`
       * :c:data:`short_range_hartree_operator`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_effective_one_e_potential`
       * :c:data:`effective_one_e_potential_sa`

 
.. c:var:: effective_one_e_potential_without_kin_sa


    File : :file:`dft_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_one_e_potential_sa	(mo_num,mo_num)
        double precision, allocatable	:: effective_one_e_potential_without_kin_sa	(mo_num,mo_num)


    State-averaged potential in MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`state_average_weight`


 
.. c:var:: energy_c


    File : :file:`dft_one_e/e_xc_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_c	(N_states)


    correlation and exchange energies general providers.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`correlation_functional`
       * :c:data:`energy_c_lda`
       * :c:data:`energy_c_none`
       * :c:data:`energy_c_sr_lda`
       * :c:data:`energy_x_pbe`
       * :c:data:`energy_x_sr_pbe`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_correlation_dft`

 
.. c:var:: energy_x


    File : :file:`dft_one_e/e_xc_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x	(N_states)


    correlation energies general providers.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x_lda`
       * :c:data:`energy_x_none`
       * :c:data:`energy_x_pbe`
       * :c:data:`energy_x_sr_lda`
       * :c:data:`energy_x_sr_pbe`
       * :c:data:`exchange_functional`
       * :c:data:`hf_exchange`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_exchange_dft`

 
.. c:var:: mu_erf_dft


    File : :file:`dft_one_e/mu_erf_dft.irp.f`

    .. code:: fortran

        double precision	:: mu_erf_dft	


    range separation parameter used in RS-DFT.
    
    It is set to mu_erf in order to be consistent with the module "hamiltonian"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_dft`

 
.. c:var:: mu_grad_rho


    File : :file:`dft_one_e/mu_erf_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_grad_rho	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_dft`

 
.. c:var:: mu_of_r_dft


    File : :file:`dft_one_e/mu_erf_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_of_r_dft	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_dft_type`
       * :c:data:`mu_erf_dft`
       * :c:data:`mu_grad_rho`
       * :c:data:`mu_of_r_hf`
       * :c:data:`mu_rsc_of_r`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`energy_c_sr_lda`
       * :c:data:`energy_x_sr_lda`
       * :c:data:`energy_x_sr_pbe`
       * :c:data:`mu_of_r_dft_average`

 
.. c:var:: mu_of_r_dft_average


    File : :file:`dft_one_e/mu_erf_dft.irp.f`

    .. code:: fortran

        double precision	:: mu_of_r_dft_average	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`


 
.. c:var:: mu_rsc_of_r


    File : :file:`dft_one_e/mu_erf_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_rsc_of_r	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_dft`

 
.. c:var:: potential_c_alpha_ao


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_c_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao	(ao_num,ao_num,N_states)


    general providers for the alpha/beta correlation potentials on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`correlation_functional`
       * :c:data:`n_states`
       * :c:data:`potential_c_alpha_ao_lda`
       * :c:data:`potential_c_alpha_ao_none`
       * :c:data:`potential_c_alpha_ao_sr_lda`
       * :c:data:`potential_c_beta_ao_none`
       * :c:data:`potential_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_sr_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`potential_c_alpha_mo`

 
.. c:var:: potential_c_alpha_mo


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_c_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_c_beta_mo	(mo_num,mo_num,N_states)


    general providers for the alpha/beta correlation potentials on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_c_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`trace_v_xc`

 
.. c:var:: potential_c_beta_ao


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_c_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao	(ao_num,ao_num,N_states)


    general providers for the alpha/beta correlation potentials on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`correlation_functional`
       * :c:data:`n_states`
       * :c:data:`potential_c_alpha_ao_lda`
       * :c:data:`potential_c_alpha_ao_none`
       * :c:data:`potential_c_alpha_ao_sr_lda`
       * :c:data:`potential_c_beta_ao_none`
       * :c:data:`potential_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_sr_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`potential_c_alpha_mo`

 
.. c:var:: potential_c_beta_mo


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_c_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_c_beta_mo	(mo_num,mo_num,N_states)


    general providers for the alpha/beta correlation potentials on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_c_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`trace_v_xc`

 
.. c:var:: potential_x_alpha_ao


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao	(ao_num,ao_num,N_states)


    general providers for the alpha/beta exchange potentials on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`exchange_functional`
       * :c:data:`hf_exchange`
       * :c:data:`n_states`
       * :c:data:`potential_x_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_none`
       * :c:data:`potential_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_sr_lda`
       * :c:data:`potential_x_alpha_ao_sr_pbe`
       * :c:data:`potential_x_beta_ao_none`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`potential_x_alpha_mo`

 
.. c:var:: potential_x_alpha_mo


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_x_beta_mo	(mo_num,mo_num,N_states)


    general providers for the alpha/beta exchange potentials on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_x_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`trace_v_xc`

 
.. c:var:: potential_x_beta_ao


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao	(ao_num,ao_num,N_states)


    general providers for the alpha/beta exchange potentials on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`exchange_functional`
       * :c:data:`hf_exchange`
       * :c:data:`n_states`
       * :c:data:`potential_x_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_none`
       * :c:data:`potential_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_sr_lda`
       * :c:data:`potential_x_alpha_ao_sr_pbe`
       * :c:data:`potential_x_beta_ao_none`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`potential_x_alpha_mo`

 
.. c:var:: potential_x_beta_mo


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_x_beta_mo	(mo_num,mo_num,N_states)


    general providers for the alpha/beta exchange potentials on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_x_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`trace_v_xc`

 
.. c:var:: potential_xc_alpha_ao


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_xc_beta_ao	(ao_num,ao_num,N_states)


    general providers for the alpha/beta exchange/correlation potentials on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`exchange_functional`
       * :c:data:`n_states`
       * :c:data:`potential_xc_alpha_ao_lda`
       * :c:data:`potential_xc_alpha_ao_none`
       * :c:data:`potential_xc_alpha_ao_pbe`
       * :c:data:`potential_xc_alpha_ao_sr_lda`
       * :c:data:`potential_xc_alpha_ao_sr_pbe`
       * :c:data:`potential_xc_beta_ao_none`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`potential_xc_alpha_mo`

 
.. c:var:: potential_xc_alpha_mo


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_xc_beta_mo	(mo_num,mo_num,N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_xc_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`trace_v_xc_new`

 
.. c:var:: potential_xc_beta_ao


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_xc_beta_ao	(ao_num,ao_num,N_states)


    general providers for the alpha/beta exchange/correlation potentials on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`exchange_functional`
       * :c:data:`n_states`
       * :c:data:`potential_xc_alpha_ao_lda`
       * :c:data:`potential_xc_alpha_ao_none`
       * :c:data:`potential_xc_alpha_ao_pbe`
       * :c:data:`potential_xc_alpha_ao_sr_lda`
       * :c:data:`potential_xc_alpha_ao_sr_pbe`
       * :c:data:`potential_xc_beta_ao_none`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`potential_xc_alpha_mo`

 
.. c:var:: potential_xc_beta_mo


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_xc_beta_mo	(mo_num,mo_num,N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_xc_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`trace_v_xc_new`

 
.. c:var:: psi_dft_energy_h_core


    File : :file:`dft_one_e/one_e_energy_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_dft_energy_kinetic	(N_states)
        double precision, allocatable	:: psi_dft_energy_nuclear_elec	(N_states)
        double precision, allocatable	:: psi_dft_energy_h_core	(N_states)


    kinetic, electron-nuclear and total h_core energy computed with the density matrix one_e_dm_mo_beta_for_dft+one_e_dm_mo_alpha_for_dft

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`


 
.. c:var:: psi_dft_energy_kinetic


    File : :file:`dft_one_e/one_e_energy_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_dft_energy_kinetic	(N_states)
        double precision, allocatable	:: psi_dft_energy_nuclear_elec	(N_states)
        double precision, allocatable	:: psi_dft_energy_h_core	(N_states)


    kinetic, electron-nuclear and total h_core energy computed with the density matrix one_e_dm_mo_beta_for_dft+one_e_dm_mo_alpha_for_dft

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`


 
.. c:var:: psi_dft_energy_nuclear_elec


    File : :file:`dft_one_e/one_e_energy_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_dft_energy_kinetic	(N_states)
        double precision, allocatable	:: psi_dft_energy_nuclear_elec	(N_states)
        double precision, allocatable	:: psi_dft_energy_h_core	(N_states)


    kinetic, electron-nuclear and total h_core energy computed with the density matrix one_e_dm_mo_beta_for_dft+one_e_dm_mo_alpha_for_dft

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`


 
.. c:var:: regular_range_hartree


    File : :file:`dft_one_e/sr_coulomb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: regular_range_hartree_operator	(mo_num,mo_num,N_states)
        double precision, allocatable	:: regular_range_hartree	(N_states)


    regular_range_Hartree_operator(i,j) = :math:`\int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}` 
    
    regular_range_Hartree = :math:`1/2  \sum_{i,j} \rho_{ij} \mathtt{regular_range_Hartree_operator}(i,j)` 
    
                        = :math:`1/2  \int dr \int r' \rho(r) \rho(r') W_{ee}^{sr}` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_average_mo_for_dft`
       * :c:data:`one_e_dm_mo_for_dft`


 
.. c:var:: regular_range_hartree_operator


    File : :file:`dft_one_e/sr_coulomb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: regular_range_hartree_operator	(mo_num,mo_num,N_states)
        double precision, allocatable	:: regular_range_hartree	(N_states)


    regular_range_Hartree_operator(i,j) = :math:`\int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}` 
    
    regular_range_Hartree = :math:`1/2  \sum_{i,j} \rho_{ij} \mathtt{regular_range_Hartree_operator}(i,j)` 
    
                        = :math:`1/2  \int dr \int r' \rho(r) \rho(r') W_{ee}^{sr}` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_average_mo_for_dft`
       * :c:data:`one_e_dm_mo_for_dft`


 
.. c:var:: short_range_hartree


    File : :file:`dft_one_e/sr_coulomb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: short_range_hartree_operator	(mo_num,mo_num,N_states)
        double precision, allocatable	:: short_range_hartree	(N_states)


    short_range_Hartree_operator(i,j) = :math:`\int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}` 
    
    short_range_Hartree = :math:`1/2  \sum_{i,j} \rho_{ij} \mathtt{short_range_Hartree_operator}(i,j)` 
    
                        = :math:`1/2  \int dr \int r' \rho(r) \rho(r') W_{ee}^{sr}` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_average_mo_for_dft`
       * :c:data:`one_e_dm_mo_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`trace_v_xc`

 
.. c:var:: short_range_hartree_operator


    File : :file:`dft_one_e/sr_coulomb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: short_range_hartree_operator	(mo_num,mo_num,N_states)
        double precision, allocatable	:: short_range_hartree	(N_states)


    short_range_Hartree_operator(i,j) = :math:`\int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}` 
    
    short_range_Hartree = :math:`1/2  \sum_{i,j} \rho_{ij} \mathtt{short_range_Hartree_operator}(i,j)` 
    
                        = :math:`1/2  \int dr \int r' \rho(r) \rho(r') W_{ee}^{sr}` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_average_mo_for_dft`
       * :c:data:`one_e_dm_mo_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`trace_v_xc`

 
.. c:var:: trace_v_h


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: trace_v_xc	(N_states)
        double precision, allocatable	:: trace_v_h	(N_states)
        double precision, allocatable	:: trace_v_hxc	(N_states)


    Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
    Trace_v_Hxc = \sum_{i,j} v^{H}_{ij} (rho_{ij}_\alpha + rho_{ij}_\beta)
    Trace_v_Hxc = \sum_{i,j} rho_{ij} v^{Hxc}_{ij}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`potential_c_alpha_mo`
       * :c:data:`potential_x_alpha_mo`
       * :c:data:`short_range_hartree_operator`


 
.. c:var:: trace_v_hxc


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: trace_v_xc	(N_states)
        double precision, allocatable	:: trace_v_h	(N_states)
        double precision, allocatable	:: trace_v_hxc	(N_states)


    Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
    Trace_v_Hxc = \sum_{i,j} v^{H}_{ij} (rho_{ij}_\alpha + rho_{ij}_\beta)
    Trace_v_Hxc = \sum_{i,j} rho_{ij} v^{Hxc}_{ij}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`potential_c_alpha_mo`
       * :c:data:`potential_x_alpha_mo`
       * :c:data:`short_range_hartree_operator`


 
.. c:var:: trace_v_xc


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: trace_v_xc	(N_states)
        double precision, allocatable	:: trace_v_h	(N_states)
        double precision, allocatable	:: trace_v_hxc	(N_states)


    Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
    Trace_v_Hxc = \sum_{i,j} v^{H}_{ij} (rho_{ij}_\alpha + rho_{ij}_\beta)
    Trace_v_Hxc = \sum_{i,j} rho_{ij} v^{Hxc}_{ij}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`potential_c_alpha_mo`
       * :c:data:`potential_x_alpha_mo`
       * :c:data:`short_range_hartree_operator`


 
.. c:var:: trace_v_xc_new


    File : :file:`dft_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: trace_v_xc_new	(N_states)


    Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`potential_xc_alpha_mo`


