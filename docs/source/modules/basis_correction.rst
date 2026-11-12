.. _module_basis_correction: 
 
.. program:: basis_correction 
 
.. default-role:: option 
 
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
    b) "mu_of_r_potential = cas_full" uses the two-body rdm of CAS-like wave function (i.e. linear combination of Slater determinants developped in an active space with the MOs stored in the EZFIO folder). 
      If the CAS is properly chosen (i.e. the CAS-like wave function qualitatively represents the wave function of the systems), then such an approach is OK for \mu(r) even in the case of strong correlation. 

  +) The use of DFT correlation functionals with multi-determinant reference (Ecmd). These functionals are originally defined in the RS-DFT framework (see for instance Theor.  Chem.  Acc.114,  305(2005)) and design to capture short-range correlation effects. A important quantity arising in the Ecmd is the exact on-top pair density of the system, and the main differences of approximated Ecmd relies on different approximations for the exact on-top pair density.  

     The two main flavours of Ecmd depends on the strength of correlation in the system: 

     a) for weakly correlated systems, the ECMD PBE-UEG functional based on the seminal work of in RSDFT (see JCP, 150, 084103 1-10 (2019)) and adapted for the basis set correction in JPCL, 10, 2931-2937 (2019) uses the exact on-top pair density of the UEG at large mu and the PBE correlation functional at mu = 0. As shown in JPCL, 10, 2931-2937 (2019), such a functional is more accurate than the ECMD LDA for weakly correlated systems. 

     b) for strongly correlated systems, the ECMD PBE-OT, which uses the extrapolated on-top pair density of the CAS wave function thanks to the large \mu behaviour of the on-top pair density, is accurate, but suffers from S_z dependence (i.e. is not invariant with respect to S_z) because of the spin-polarization dependence of the PBE correlation functional entering at mu=0.   

        An alternative is ECMD SU-PBE-OT which uses the same on-top pair density that ECMD PBE-OT but a ZERO spin-polarization to remove the S_z dependence. As shown in ???????????, this strategy is one of the more accurate and respects S_z invariance and size consistency if the CAS wave function is correctly chosen. 

 
 
 
Programs 
-------- 
 
 * :ref:`basis_correction` 
 * :ref:`print_n2_related_stuff` 
 * :ref:`print_su_pbe_ot` 
 * :ref:`test_chol_bas` 
 
Providers 
--------- 
 
.. c:var:: double average_extrapolated_on_top


    File : :file:`pbe_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: extrapolated_on_top	(n_points_final_grid,N_states)
        double precision, allocatable	:: average_extrapolated_on_top	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_potential`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`total_cas_on_top_density`


 
.. c:var:: double ecmd_lda_eff_xi_mu_of_r


    File : :file:`eff_xi_based_func.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ecmd_lda_eff_xi_mu_of_r	(N_states)


    ecmd_lda_eff_xi_mu_of_r = multi-determinantal Ecmd within the LDA and effective spin polarization approximation with mu(r),
    
    corresponds to equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
    
    !!!! BUT !!!! with an EFFECTIVE SPIN POLARIZATION DEPENDING ON THE ON-TOP PAIR DENSITY
    
    See P. Perdew, A. Savin, and K. Burke, Phys. Rev. A 51, 4531 (1995). for original Ref., and Eq. 29 in ???????????

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`effective_spin_dm`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double ecmd_lda_mu_of_r


    File : :file:`weak_corr_func.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ecmd_lda_mu_of_r	(N_states)


    ecmd_lda_mu_of_r = multi-determinantal Ecmd within the LDA approximation with mu(r) ,
    
    see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double ecmd_pbe_on_top_mu_of_r


    File : :file:`pbe_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ecmd_pbe_on_top_mu_of_r	(N_states)


    
    Ecmd functional evaluated with mu(r) and depending on
       +) the on-top pair density
    
       +) the total density, density gradients
    
       +) the spin density
    
    Defined originally in Eq. (25) of JCP, 150, 084103 1-10 (2019) for RS-DFT calculations, but evaluated with mu(r).
    
    Such a functional is built by interpolating between two regimes :
    
       +) the large mu behaviour in cst/(\mu^3) \int dr on-top(r) where on-top(r) is supposed to be the exact on-top of the system
    
       +) mu= 0 with the usal ec_pbe(rho_a,rho_b,grad_rho_a,grad_rho_b)
    
    Here the approximation to the exact on-top is done through the assymptotic expansion (in \mu) of the exact on-top pair density (see Eq. 29) but with a mu(r) instead of a constant mu
    
    Such an asymptotic expansion was introduced in P. Gori-Giorgi and A. Savin, Phys. Rev. A73, 032506 (2006)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_potential`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`
       * :c:data:`total_cas_on_top_density`


 
.. c:var:: double ecmd_pbe_on_top_no_extrap_su_mu_of_r


    File : :file:`pbe_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ecmd_pbe_on_top_no_extrap_su_mu_of_r	(N_states)


    
    Ecmd functional evaluated with mu(r) and depending on
       +) the on-top pair density
    
       +) the total density, density gradients
    
       +) !!!!! NO SPIN POLAIRIZATION !!!!!
    
    Defined originally in Eq. (25) of JCP, 150, 084103 1-10 (2019) for RS-DFT calculations, but evaluated with mu(r).
    
    Such a functional is built by interpolating between two regimes :
    
       +) the large mu behaviour in cst/(\mu^3) \int dr on-top(r) where on-top(r) is supposed to be the exact on-top of the system
    
       +) mu= 0 with the usal ec_pbe(rho_a,rho_b,grad_rho_a,grad_rho_b)
    
    Here the approximation to the exact on-top is done through the assymptotic expansion (in \mu) of the exact on-top pair density (see Eq. 29) but with a mu(r) instead of a constant mu
    
    Such an asymptotic expansion was introduced in P. Gori-Giorgi and A. Savin, Phys. Rev. A73, 032506 (2006)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_potential`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`
       * :c:data:`total_cas_on_top_density`


 
.. c:var:: double ecmd_pbe_on_top_su_mu_of_r


    File : :file:`pbe_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ecmd_pbe_on_top_su_mu_of_r	(N_states)


    
    Ecmd functional evaluated with mu(r) and depending on
       +) the on-top pair density
    
       +) the total density, density gradients
    
       +) !!!!! NO SPIN POLAIRIZATION !!!!!
    
    Defined originally in Eq. (25) of JCP, 150, 084103 1-10 (2019) for RS-DFT calculations, but evaluated with mu(r).
    
    Such a functional is built by interpolating between two regimes :
    
       +) the large mu behaviour in cst/(\mu^3) \int dr on-top(r) where on-top(r) is supposed to be the exact on-top of the system
    
       +) mu= 0 with the usal ec_pbe(rho_a,rho_b,grad_rho_a,grad_rho_b)
    
    Here the approximation to the exact on-top is done through the assymptotic expansion (in \mu) of the exact on-top pair density (see Eq. 29) but with a mu(r) instead of a constant mu
    
    Such an asymptotic expansion was introduced in P. Gori-Giorgi and A. Savin, Phys. Rev. A73, 032506 (2006)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_potential`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`
       * :c:data:`total_cas_on_top_density`


 
.. c:var:: double ecmd_pbe_ueg_eff_xi_mu_of_r


    File : :file:`eff_xi_based_func.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ecmd_pbe_ueg_eff_xi_mu_of_r	(N_states)


    ecmd_pbe_ueg_eff_xi_mu_of_r = multi-determinantal Ecmd within the PBE-UEG and effective spin polarization approximation with mu(r),
    
    see Eqs. 30 in ???????????
    
    Based on the PBE-on-top functional (see Eqs. 26, 27 of J. Chem. Phys.150, 084103 (2019); doi: 10.1063/1.5082638)
    
    and replaces the approximation of the exact on-top pair density by the exact on-top of the UEG
    
    !!!! BUT !!!! with an EFFECTIVE SPIN POLARIZATION DEPENDING ON THE ON-TOP PAIR DENSITY
    
    See P. Perdew, A. Savin, and K. Burke, Phys. Rev. A 51, 4531 (1995). for original Ref., and Eq. 29 in ???????????

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`effective_spin_dm`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double ecmd_pbe_ueg_mu_of_r


    File : :file:`weak_corr_func.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ecmd_pbe_ueg_mu_of_r	(N_states)


    ecmd_pbe_ueg_mu_of_r = multi-determinantal Ecmd within the PBE-UEG approximation with mu(r) ,
    
    see Eqs. 13-14b in Phys.Chem.Lett.2019, 10, 2931   2937;  https://pubs.acs.org/doi/10.1021/acs.jpclett.9b01176
    
    Based on the PBE-on-top functional (see Eqs. 26, 27 of J. Chem. Phys.150, 084103 (2019); doi: 10.1063/1.5082638)
    
    but it the on-top pair density of the UEG as an approximation of the exact on-top pair density

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double extrapolated_on_top


    File : :file:`pbe_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: extrapolated_on_top	(n_points_final_grid,N_states)
        double precision, allocatable	:: average_extrapolated_on_top	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_potential`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`total_cas_on_top_density`


 
.. c:var:: double sr_correction_on_top3_mu_of_r


    File : :file:`sr_correction_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: sr_correction_on_top_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_on_top3_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu3_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho3_of_r	(N_states)


    :math:`\frac{1}{2} \int P_2(R) \cdot a_3 \cdot (4\epsilon(R))^{3/2} \frac{1}{1 + \frac{4}{\sqrt{\pi}}\sqrt{\epsilon(R)}} \, dR` 
    
    :math:`\frac{1}{2} \int P_2(R) \cdot \left(1 + \frac{2}{\sqrt{\pi}\mu}\right)^{-1}\cdot a_3 \mu^{-3}  \, dR` 
    
    :math:`a_3 = -\frac43 \sqrt{\pi} (\sqrt{2}-1)` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double sr_correction_on_top_mu_of_r


    File : :file:`sr_correction_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: sr_correction_on_top_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_on_top3_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu3_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho3_of_r	(N_states)


    :math:`\frac{1}{2} \int P_2(R) \cdot a_3 \cdot (4\epsilon(R))^{3/2} \frac{1}{1 + \frac{4}{\sqrt{\pi}}\sqrt{\epsilon(R)}} \, dR` 
    
    :math:`\frac{1}{2} \int P_2(R) \cdot \left(1 + \frac{2}{\sqrt{\pi}\mu}\right)^{-1}\cdot a_3 \mu^{-3}  \, dR` 
    
    :math:`a_3 = -\frac43 \sqrt{\pi} (\sqrt{2}-1)` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double sr_correction_rho3_of_r


    File : :file:`sr_correction_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: sr_correction_on_top_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_on_top3_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu3_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho3_of_r	(N_states)


    :math:`\frac{1}{2} \int P_2(R) \cdot a_3 \cdot (4\epsilon(R))^{3/2} \frac{1}{1 + \frac{4}{\sqrt{\pi}}\sqrt{\epsilon(R)}} \, dR` 
    
    :math:`\frac{1}{2} \int P_2(R) \cdot \left(1 + \frac{2}{\sqrt{\pi}\mu}\right)^{-1}\cdot a_3 \mu^{-3}  \, dR` 
    
    :math:`a_3 = -\frac43 \sqrt{\pi} (\sqrt{2}-1)` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double sr_correction_rho_mu3_of_r


    File : :file:`sr_correction_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: sr_correction_on_top_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_on_top3_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu3_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho3_of_r	(N_states)


    :math:`\frac{1}{2} \int P_2(R) \cdot a_3 \cdot (4\epsilon(R))^{3/2} \frac{1}{1 + \frac{4}{\sqrt{\pi}}\sqrt{\epsilon(R)}} \, dR` 
    
    :math:`\frac{1}{2} \int P_2(R) \cdot \left(1 + \frac{2}{\sqrt{\pi}\mu}\right)^{-1}\cdot a_3 \mu^{-3}  \, dR` 
    
    :math:`a_3 = -\frac43 \sqrt{\pi} (\sqrt{2}-1)` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double sr_correction_rho_mu_of_r


    File : :file:`sr_correction_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: sr_correction_on_top_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_on_top3_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu3_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho3_of_r	(N_states)


    :math:`\frac{1}{2} \int P_2(R) \cdot a_3 \cdot (4\epsilon(R))^{3/2} \frac{1}{1 + \frac{4}{\sqrt{\pi}}\sqrt{\epsilon(R)}} \, dR` 
    
    :math:`\frac{1}{2} \int P_2(R) \cdot \left(1 + \frac{2}{\sqrt{\pi}\mu}\right)^{-1}\cdot a_3 \mu^{-3}  \, dR` 
    
    :math:`a_3 = -\frac43 \sqrt{\pi} (\sqrt{2}-1)` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double sr_correction_rho_of_r


    File : :file:`sr_correction_on_top.irp.f`

    .. code:: fortran

        double precision, allocatable	:: sr_correction_on_top_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_on_top3_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho_mu3_of_r	(N_states)
        double precision, allocatable	:: sr_correction_rho3_of_r	(N_states)


    :math:`\frac{1}{2} \int P_2(R) \cdot a_3 \cdot (4\epsilon(R))^{3/2} \frac{1}{1 + \frac{4}{\sqrt{\pi}}\sqrt{\epsilon(R)}} \, dR` 
    
    :math:`\frac{1}{2} \int P_2(R) \cdot \left(1 + \frac{2}{\sqrt{\pi}\mu}\right)^{-1}\cdot a_3 \mu^{-3}  \, dR` 
    
    :math:`a_3 = -\frac43 \sqrt{\pi} (\sqrt{2}-1)` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void g0_gori(void)


    File : :file:`sr_correction_on_top.irp.f`

    .. code:: fortran

        double precision function g0_gori(rs, rho_a, rho_b)


    Spin-polarized approximate on-top pair density of the UEG.
    
    Equation (47) of DOI: 10.1103/PhysRevB.64.155102
    P. Gori-Giorgi and J. P. Perdew, PHYSICAL REVIEW B, VOLUME 64, 155102

 
.. c:function:: void g0_gori_mu(void)


    File : :file:`sr_correction_on_top.irp.f`

    .. code:: fortran

        double precision function g0_gori_mu(rs, rho_a, rho_b,mu)


    Spin-polarized approximate on-top pair density of the UEG.
    
    Equation (47) of DOI: 10.1103/PhysRevB.64.155102
    P. Gori-Giorgi and J. P. Perdew, PHYSICAL REVIEW B, VOLUME 64, 155102

 
.. c:function:: void g0_gori_su(void)


    File : :file:`sr_correction_on_top.irp.f`

    .. code:: fortran

        double precision function g0_gori_su(rs)


    Spin-unpolarized approximate on-top pair density of the UEG.
    
    Equation (30) of DOI: 10.1103/PhysRevB.64.155102
    P. Gori-Giorgi and J. P. Perdew, PHYSICAL REVIEW B, VOLUME 64, 155102

 
.. c:function:: void g0_gori_su_mu(void)


    File : :file:`sr_correction_on_top.irp.f`

    .. code:: fortran

        double precision function g0_gori_su_mu(rs, mu)


    Eq (47) in http://dx.doi.org/10.1103/PhysRevA.73.032506
    Without the -1/2 at the end

 
.. c:function:: void print_basis_correction(void)


    File : :file:`print_routine.irp.f`

    .. code:: fortran

        subroutine print_basis_correction



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ecmd_lda_mu_of_r`
       * :c:data:`ecmd_pbe_on_top_mu_of_r`
       * :c:data:`ecmd_pbe_on_top_su_mu_of_r`
       * :c:data:`ecmd_pbe_ueg_mu_of_r`
       * :c:data:`mu_average_prov`
       * :c:data:`mu_average_prov2`
       * :c:data:`mu_average_prov3`
       * :c:data:`mu_average_prov_ot`
       * :c:data:`mu_of_r_potential`
       * :c:data:`n_states`
       * :c:data:`sr_correction_on_top_mu_of_r`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`basis_correction`

 
.. c:function:: void print_su_pbe_ot_subroutine(void)


    File : :file:`print_su_pbe_ot.irp.f`

    .. code:: fortran

        subroutine print_su_pbe_ot_subroutine



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ecmd_pbe_on_top_su_mu_of_r`
       * :c:data:`ecmd_pbe_ueg_mu_of_r`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_su_pbe_ot`

 
.. c:function:: void routine(void)


    File : :file:`print_n2_related_stuff.irp.f`

    .. code:: fortran

        subroutine routine



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`average_on_top`
       * :c:data:`extrapolated_on_top`
       * :c:data:`mu_average_prov`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`cas_based_on_top`
       * :c:func:`print_n2_related_stuff`

 
.. c:function:: void test(void)


    File : :file:`test_chol_bas.irp.f`

    .. code:: fortran

        subroutine test



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_hf_cholesky`
       * :c:data:`f_hf_cholesky_sparse`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_chol_bas`
       * :c:func:`test_proj_op`

