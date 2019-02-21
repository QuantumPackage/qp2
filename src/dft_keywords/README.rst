============
dft_keywords
============

This module contains the main keywords related to a DFT calculation or RS-DFT calculation. 
These keywords are related to the following programs of the |QP| core modules:

* :ref:`ks_scf` : Kohn-Sham |DFT|
* :ref:`rs_ks_scf` : Range separated Hybrids |DFT|


Modifying the exchange/correlation functionals
----------------------------------------------
To modify the exchange/correlation functionals, see the following keywords: 

* :option:`dft_keywords exchange_functional`: type of exchange functionals
* :option:`dft_keywords correlation_functional`: type of correlation functionals 

Each of these keywords can have the following value: 
* "LDA" or "short_range_LDA" for, respectively the |LDA| and its short-range version 
* "PBE" or "short_range_PBE" for, respectively the |PBE| and its short-range version 


Modifying the amount of |HF| exchange
-------------------------------------
* :option:`dft_keywords HF_exchange`  : only relevent for the :ref:`ks_scf` program


Other related keywords not defined in :ref:`module_dft_keywords`
----------------------------------------------------------------
The keyword for the **range separation parameter**  :math:`\mu` is the :option:`ao_two_e_erf_ints mu_erf` keyword.

The keyword for the **type of density used in RS-DFT** calculation with a **multi-configurational wave function** is the :option:`density_for_dft density_for_dft` keyword.
