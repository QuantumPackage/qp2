==================
ao_many_one_e_ints
==================

This module contains A LOT of one-electron integrals of the type 
A_ij( r ) = \int dr' phi_i(r') w(r,r') phi_j(r') 
where r is a point in real space. 

+) ao_gaus_gauss.irp.f: w(r,r') is a exp(-(r-r')^2)  , and can be multiplied by x/y/z
+) ao_erf_gauss.irp.f : w(r,r') is a exp(-(r-r')^2) erf(mu * |r-r'|)/|r-r'| , and can be multiplied by x/y/z
+) ao_erf_gauss_grad.irp.f: w(r,r') is a exp(-(r-r')^2) erf(mu * |r-r'|)/|r-r'| , and can be multiplied by x/y/z, but evaluated with also one gradient of an AO function. 

Fit of a Slater function and corresponding integrals
----------------------------------------------------
The file fit_slat_gauss.irp.f contains many useful providers/routines to fit a Slater function with 20 gaussian. 
+) coef_fit_slat_gauss : coefficients of the gaussians to fit e^(-x)
+) expo_fit_slat_gauss : exponents of the gaussians to fit e^(-x)

Integrals involving Slater functions : stg_gauss_int.irp.f

Taylor expansion of full correlation factor
-------------------------------------------
In taylor_exp.irp.f you might find interesting integrals of the type 
\int dr' exp( e^{-alpha |r-r|' - beta |r-r'|^2}) phi_i(r') phi_j(r') 
evaluated as a Taylor expansion of the exponential. 
