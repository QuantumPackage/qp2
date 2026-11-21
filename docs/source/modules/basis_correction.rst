.. _module_basis_correction: 
 
.. program:: basis_correction 
 
.. default-role:: option 
 
basis_correction
================

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Providers 
--------- 
 
.. c:var:: ecmd_lda_eff_xi_mu_of_r


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


 
.. c:var:: ecmd_lda_mu_of_r


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


 
.. c:var:: ecmd_pbe_on_top_mu_of_r


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


 
.. c:var:: ecmd_pbe_on_top_no_extrap_su_mu_of_r


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


 
.. c:var:: ecmd_pbe_on_top_su_mu_of_r


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


 
.. c:var:: ecmd_pbe_ueg_eff_xi_mu_of_r


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


 
.. c:var:: ecmd_pbe_ueg_mu_of_r


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


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: basis_corr_su_pbe_ot:


    File : :file:`print_su_pbe_ot.irp.f`

    .. code:: fortran

        subroutine basis_corr_su_pbe_ot


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`io_mo_two_e_integrals`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`no_core_density`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_su_pbe_ot`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`no_core_density`
       * :c:data:`read_wf`

 
.. c:function:: basis_correction:


    File : :file:`basis_correction.irp.f`

    .. code:: fortran

        subroutine basis_correction


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`io_mo_two_e_integrals`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`no_core_density`
       * :c:data:`read_wf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_basis_correction`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`no_core_density`
       * :c:data:`read_wf`

 
.. c:function:: print_basis_correction:


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
       * :c:data:`mu_of_r_potential`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`basis_correction`

 
.. c:function:: print_su_pbe_ot:


    File : :file:`print_su_pbe_ot.irp.f`

    .. code:: fortran

        subroutine print_su_pbe_ot



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ecmd_pbe_on_top_su_mu_of_r`
       * :c:data:`ecmd_pbe_ueg_mu_of_r`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`basis_corr_su_pbe_ot`

