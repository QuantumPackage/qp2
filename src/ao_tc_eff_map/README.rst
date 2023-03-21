ao_tc_eff_map
=============

This is a module to obtain the integrals on the AO basis of the SCALAR HERMITIAN 
effective potential defined in Eq. 32 of JCP 154, 084119 (2021)
It also contains the modification by a one-body Jastrow factor.  

The main routine/providers are

+) ao_tc_sym_two_e_pot_map : map of the SCALAR PART of total effective two-electron on the AO basis in PHYSICIST notations. It might contain the two-electron term coming from the one-e correlation factor. 
+) get_ao_tc_sym_two_e_pot(i,j,k,l,ao_tc_sym_two_e_pot_map) : routine to get the integrals from ao_tc_sym_two_e_pot_map. 
+) ao_tc_sym_two_e_pot(i,j,k,l) : FUNCTION that returns the scalar part of TC-potential EXCLUDING the erf(mu r12)/r12. See two_e_ints_gauss.irp.f for more details. 
