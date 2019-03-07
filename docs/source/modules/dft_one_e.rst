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
       * :c:data:`energy_c_pbe`
       * :c:data:`energy_c_sr_lda`
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
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_exchange_dft`

 
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


