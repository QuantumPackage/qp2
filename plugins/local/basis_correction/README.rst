================
basis_correction
================

This module proposes the various flavours of the DFT-based basis set correction originally proposed in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714. 

This basis set correction relies mainy on : 

  +) The definition of a range-separation function \mu(r) varying in space to mimic the incompleteness of the basis set used to represent the coulomb interaction. This procedure needs a two-body rdm representing qualitatively the spacial distribution of the opposite spin electron pairs. 
    Two types of \mu(r) are proposed, according to the strength of correlation, through the keyword "mu_of_r_potential" in the module "mu_of_r": 
    a) "mu_of_r_potential = hf" uses the two-body rdm of a HF-like wave function (i.e. a single Slater determinant developped with the MOs stored in the EZFIO folder). 
       When HF is a qualitative representation of the electron pairs (i.e. weakly correlated systems), such an approach for \mu(r) is OK. 
       See for instance JPCL, 10, 2931-2937 (2019) for typical flavours of the results. 
       Thanks to the trivial nature of such a two-body rdm, the equation (22) of J. Chem. Phys. 149, 194301 (2018) can be rewritten in a very efficient way, and therefore the limiting factor of such an approach is the AO->MO four-index transformation of the two-electron integrals.  
    b) "mu_of_r_potential = cas_ful" uses the two-body rdm of CAS-like wave function (i.e. linear combination of Slater determinants developped in an active space with the MOs stored in the EZFIO folder). 
      If the CAS is properly chosen (i.e. the CAS-like wave function qualitatively represents the wave function of the systems), then such an approach is OK for \mu(r) even in the case of strong correlation. 

  +) The use of DFT correlation functionals with multi-determinant reference (Ecmd). These functionals are originally defined in the RS-DFT framework (see for instance Theor.  Chem.  Acc.114,  305(2005)) and design to capture short-range correlation effects. A important quantity arising in the Ecmd is the exact on-top pair density of the system, and the main differences of approximated Ecmd relies on different approximations for the exact on-top pair density.  

     The two main flavours of Ecmd depends on the strength of correlation in the system: 

     a) for weakly correlated systems, the ECMD PBE-UEG functional based on the seminal work of in RSDFT (see JCP, 150, 084103 1-10 (2019)) and adapted for the basis set correction in JPCL, 10, 2931-2937 (2019) uses the exact on-top pair density of the UEG at large mu and the PBE correlation functional at mu = 0. As shown in JPCL, 10, 2931-2937 (2019), such a functional is more accurate than the ECMD LDA for weakly correlated systems. 

     b) for strongly correlated systems, the ECMD PBE-OT, which uses the extrapolated on-top pair density of the CAS wave function thanks to the large \mu behaviour of the on-top pair density, is accurate, but suffers from S_z dependence (i.e. is not invariant with respect to S_z) because of the spin-polarization dependence of the PBE correlation functional entering at mu=0.   

        An alternative is ECMD SU-PBE-OT which uses the same on-top pair density that ECMD PBE-OT but a ZERO spin-polarization to remove the S_z dependence. As shown in ???????????, this strategy is one of the more accurate and respects S_z invariance and size consistency if the CAS wave function is correctly chosen. 

