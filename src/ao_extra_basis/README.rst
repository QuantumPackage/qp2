===========
extra_basis
===========

Plugin to handle an extra basis, which is attached to the extra_nuclei.
It is essentially a duplication of all important quantities (coefficients, exponents and so on) of the usual |AO| basis. 

An interesting feature is the possibility to fit any basis made at most with "p" functions onto a purely "s" basis. 
This is done with the various scripts here: 
 
 - qp_fit_1s_basis : script that creates an |EZFIO| folder corresponding to an .xyz file and a basis fitted with only "s" functions
 - qp_add_extra_fit_system : script that takes as input an |EZFIO| folder and an .xyz file and add an extra_basis and extra_nuclei with a 1s fitted basis
 
Ex: 
qp_add_extra_fit_system LiH.ezfio/ h2o.xyz # takes the EZFIO folder "LiH.ezfio" and creates all necessary additional basis and nuclei based on h2o.xyz, but only with 1s functions. 
