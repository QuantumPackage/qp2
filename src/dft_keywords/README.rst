============
dft_keywords
============

This module contains the main keywords related to a DFT calculation or RS-DFT calculation, such as:

* :option:`dft_keywords exchange_functional`
* :option:`dft_keywords correlation_functional`
* :option:`dft_keywords HF_exchange`  : only relevent for the :c:func:`rs_ks_scf` program

The keyword for the **range separation parameter**  :math:`\mu` is the :option:`ao_two_e_erf_ints mu_erf` keyword.

The keyword for the type of density used in RS-DFT calculation with a multi-configurational wave function is the :option:`density_for_dft density_for_dft` keyword.
