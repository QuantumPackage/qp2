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


If you have designed your own exchange/correlation functionals (see the documentation of the :ref:`module_new_functionals`), 
you can use them in all |DFT|-based programs by  setting the :option:`dft_keywords exchange_functional` and :option:`dft_keywords correlation_functional` keywords to "my_functional". 
