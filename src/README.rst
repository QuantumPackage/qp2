==========================
The core modules of the QP
==========================

*** How are handled the DFT functionals in QP2 ?
================================================
    The Exchange and Correlation energies/potentials can be accessed by the following providers 
    energy_x
    energy_c
    potential_x_alpha_ao
    potential_c_alpha_ao
    potential_x_beta_ao
    potential_c_beta_ao

    These providers are automatically linked to the providers of the actual exchange/correlation energies of a given functional 
    through the character keywords 
    "exchange_functional"
    "correlation_functional"

    All the providers for the available functionals are in the folder "functionals", with one file "my_functional.irp.f" per functional.  

    Ex : if "exchange_functional" == "sr_pbe", then energy_x will contain the exchange correlation functional defined in "functiona/sr_pbe.irp.f", which corresponds to the short-range PBE functional (at the value mu_erf for the range separation parameter) 


*** How are handled the DFT functionals in QP2 ?
================================================

    Creating a new functional and propagating it through the whole QP2 programs is easy as all dependencies are handled by a script. 

    To do so, let us assume that the name of your functional is "my_func". 
    Then you just have to create the file "my_func.irp.f" in the folder "functional" which shoud contain 

    +) if you're adding an exchange functional, then create the provider "energy_x_my_func" 

    +) if you're adding a correlation functional, create the provider "energy_c_my_func" 
    
    +) if you want to add the echange potentials, create the providers "potential_x_alpha_ao_my_func", "potential_x_beta_ao_my_func"  which are the exchange potentials on the AO basis for the alpha/beta electrons 

    +) if you want to add the correlation potentials, create the providers "potential_c_alpha_ao_my_func", "potential_c_beta_ao_my_func"  which are the correlation potentials on the AO basis for the alpha/beta electrons 
    
    That's all :) 
   
    Then, when running whatever DFT calculation or accessing/using the providers: 
    energy_x
    energy_c
    potential_x_alpha_ao
    potential_c_alpha_ao
    potential_x_beta_ao
    potential_c_beta_ao
 
    if exchange_functional = mu_func, then you will automatically have access to what you need, such as kohn sham orbital optimization and so on ... 
