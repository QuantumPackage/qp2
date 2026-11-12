.. _module_dft_utils_func: 
 
.. program:: dft_utils_func 
 
.. default-role:: option 
 
===============
dft_utils_func
===============

Exchange and correlation functionals of |DFT| and |RSDFT|.

* the LDA exchange and correlation, and their short range
  variants: `ecmd_lda.irp.f` and `routines_exc_sr_lda.irp.f`
  (:c:data:`ex_lda`, :c:data:`ec_lda`, :c:data:`ex_lda_sr`,
  :c:data:`ec_lda_sr`, and the range separated correlation
  :c:data:`ecorrlr` / :c:data:`vcorrlr`)
* the PBE and the SR-PBE functionals: `ecmd_pbe_general.irp.f` and
  `routines_exc_sr_pbe.irp.f` (:c:data:`ex_pbe_sr`,
  :c:data:`ec_pbe_sr`)
* the meta-GGA |PBE| on-top correlation, with its implementation in
  terms of the uniform electron gas: `ecmd_pbe_on_top.irp.f`,
  `ecmd_pbe_ueg.irp.f` and `on_top_from_ueg.irp.f`
  (:c:data:`ecmd_pbe_ueg_at_r`, :c:data:`g0_ueg_mu`,
  :c:data:`mu_correction_of_on_top`)
* the range separation functions of the interaction:
  `mu_of_r_dft.irp.f` (:c:data:`mu_rs_c`,
  :c:data:`mu_grad_rho_func`) and the error function and screened
  Coulomb functions in `utils.irp.f` (:c:data:`g0`, :c:data:`g0f`,
  :c:data:`g0s`, :c:data:`g0d` and their derivatives)
* the transformations between the spin resolved and the total
  density, for the screened Coulomb interaction:
  `rho_ab_to_rho_tot.irp.f` (:c:data:`rho_oc_to_rho_ab`,
  :c:data:`v_rho_oc_to_v_rho_ab` and their variants)
 
 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void berf(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        function berf(a)


    Second-order exchange gradient expansion coefficient for the erf
    interaction, as a function of a = mu / (2 kF). The value is
    multiplied by the damping factor 2.540118935556 * exp(-19 a**2);
    series for a < 0.075 and for a > 50, exact expression (with
    derf) in between.

 
.. c:function:: void beta_rs(void)


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        double precision function beta_rs(rs)



 
.. c:function:: void correction_to_on_top_from_ueg(void)


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        double precision function correction_to_on_top_from_UEG(mu,r,istate)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dm_dft_alpha_beta_at_r`

 
.. c:function:: void dberfda(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        function dberfda(a)


    Derivative with respect to a = mu / (2 kF) of the coefficient
    berf for the erf interaction. Same branches and damping factor
    as berf.

 
.. c:function:: void delta_lrsr_ldaerf(void)


    File : :file:`dft_utils_func/ecmd_lda.irp.f`

    .. code:: fortran

        subroutine DELTA_LRSR_LDAERF (rhot,rhos,mu,dospin,e)


    LDA approximation of the Delta_(LR-SR) term of the
    multideterminant short-range LDA correlation. rhot is the total
    density, rhos the spin density, mu the range-separation
    parameter. If dospin is false, z = 0. e is the contribution to
    the energy density (delta * rhot). Not tested for z != 0.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`esrc_md_ldaerf`

 
.. c:function:: void dpol(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function dpol(rs)



 
.. c:function:: void dpold(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function dpold(rs)



 
.. c:function:: void dpoldd(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function dpoldd(rs)



 
.. c:function:: void ec_lda(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ec_lda(rho_a,rho_b,ec,vc_a,vc_b)


    LDA correlation: energy density ec and potentials vc_a, vc_b
    computed from the alpha and beta spin densities. The correlation
    per particle is taken from ecPW (uniform gas, Coulomb
    interaction); ec is that value times rho = rho_a + rho_b. For
    total densities below 1.d-10 the outputs are set to 1.d-15.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_pbe_only`
       * :c:func:`ec_pbe_sr`
       * :c:data:`energy_c_lda`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ecpw`

 
.. c:function:: void ec_lda_sr(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ec_lda_sr(mu,rho_a,rho_b,ec,vc_a,vc_b)


    Short-range LDA correlation: energy density ec and potentials
    vc_a, vc_b for the spin densities, with range-separation
    parameter mu. Each quantity is the full LDA value minus the
    long-range part: ec from ecPW and ecorrlr, vc_a and vc_b from
    the LDA correlation potential and vcorrlr. ec is multiplied by
    rho = rho_a + rho_b.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:func:`ec_pbe_only`
       * :c:func:`ec_pbe_sr`
       * :c:data:`energy_c_sr_lda`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ecorrlr`
       * :c:func:`ecpw`
       * :c:func:`vcorrlr`

 
.. c:function:: void ec_lyp2(void)


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        double precision function ec_lyp2(RhoA,RhoB,GA,GB,GAB)



 
.. c:function:: void ec_lyp_88(void)


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        double precision function ec_lyp_88(rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)


    LYP functional of the Lee, Yan, Parr, Phys. Rev B 1988, Vol 37, page 785.
    The expression used is the one by Miehlich, Savin, Stoll, Preuss, CPL, 1989 which gets rid of the laplacian of the density

 
.. c:function:: void ec_md_on_top_pbe_mu_corrected(void)


    File : :file:`dft_utils_func/ecmd_pbe_on_top.irp.f`

    .. code:: fortran

        subroutine ec_md_on_top_PBE_mu_corrected(mu,r,two_dm,eps_c_md_on_top_PBE)


    enter with "r(3)", and "two_dm(N_states)" which is the on-top pair density at "r" for each states
    
    you get out with the energy density defined in J. Chem. Phys.150, 084103 (2019); doi: 10.1063/1.508263
    
    by Eq. (26), which includes the correction of the on-top pair density of Eq. (29).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r`
       * :c:func:`ec_pbe_only`
       * :c:func:`grad_rho_ab_to_grad_rho_oc`
       * :c:func:`rho_ab_to_rho_oc`

 
.. c:function:: void ec_md_pbe_on_top_general(void)


    File : :file:`dft_utils_func/ecmd_pbe_general.irp.f`

    .. code:: fortran

        subroutine ec_md_pbe_on_top_general(mu,rho_a,rho_b,grad_rho_a,grad_rho_b,on_top,eps_c_md_on_top_PBE)


    
    General e_cmd functional interpolating between :
    
                                                     +) the large mu behaviour in cst/(\mu^3) on-top
    
                                                     +) mu= 0 with the usal ec_pbe at
    
    Depends on : mu, the density (rho_a,rho_b), the square of the density gradient (grad_rho_a,grad_rho_b)
    
                 the flavour of on-top densiyt (on_top) you fill in: in principle it should be the exact on-top
    
    The form of the functional was originally introduced in JCP, 150, 084103 1-10 (2019)
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_pbe_only`
       * :c:func:`grad_rho_ab_to_grad_rho_oc`
       * :c:func:`rho_ab_to_rho_oc`

 
.. c:function:: void ec_only_lda_sr(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ec_only_lda_sr(mu,rho_a,rho_b,ec)


    Short-range LDA correlation energy density only (no potentials),
    for the spin densities rho_a and rho_b with range-separation
    parameter mu: the full LDA correlation (ecPW) minus the
    long-range part (ecorrlr), times rho = rho_a + rho_b. ec is 0
    for total densities below 1.d-12.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_scan`
       * :c:func:`esrc_md_ldaerf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ecorrlr`
       * :c:func:`ecpw`

 
.. c:function:: void ec_pbe_only(void)


    File : :file:`dft_utils_func/routines_exc_sr_pbe.irp.f`

    .. code:: fortran

        subroutine ec_pbe_only(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec)


    Short-range pbe correlation energy functional for erf interaction
    
    input : ==========
    
    mu = range separated parameter
    
    rhoc, rhoo = total density and spin density
    
    sigmacc    = square of the gradient of the total density
    
    sigmaco    = square of the gradient of the spin density
    
    sigmaoo    = scalar product between the gradient of the total density and the one of the spin density
    
    output: ==========
    
    ec         = correlation energy
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_md_on_top_pbe_mu_corrected`
       * :c:func:`ec_md_pbe_on_top_general`
       * :c:func:`ecmd_pbe_ueg_at_r`
       * :c:func:`eps_c_md_pbe_at_grid_pt`
       * :c:func:`eps_c_md_pbe_from_density`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_lda`
       * :c:func:`ec_lda_sr`

 
.. c:function:: void ec_pbe_sr(void)


    File : :file:`dft_utils_func/routines_exc_sr_pbe.irp.f`

    .. code:: fortran

        subroutine ec_pbe_sr(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)


    Short-range pbe correlation energy functional for erf interaction
    
    input : ==========
    
    mu = range separated parameter
    
    rhoc, rhoo = total density and spin density
    
    sigmacc    = square of the gradient of the total density
    
    sigmaco    = square of the gradient of the spin density
    
    sigmaoo    = scalar product between the gradient of the total density and the one of the spin density
    
    output: ==========
    
    ec         = correlation energy
    
    all variables v** are energy derivatives with respect to components of the density
    
    vrhoc      = derivative with respect to the total density
    
    vrhoo      = derivative with respect to spin density
    
    vsigmacc   = derivative with respect to the square of the gradient of the total density
    
    vsigmaco   = derivative with respect to scalar product between the gradients of total and spin densities
    
    vsigmaoo   = derivative with respect to the square of the gradient of the psin density

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gga_sr_type_functionals`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_lda`
       * :c:func:`ec_lda_sr`

 
.. c:function:: void ec_scan(void)


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        double precision function ec_scan(rho_a,rho_b,tau,grad_rho_2)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_only_lda_sr`

 
.. c:function:: void ecmd_pbe_ueg_at_r(void)


    File : :file:`dft_utils_func/ecmd_pbe_ueg.irp.f`

    .. code:: fortran

        subroutine ecmd_pbe_ueg_at_r(mu,r,eps_c_md_PBE)


    provides the integrand of Eq. (13) of Phys.Chem.Lett.2019, 10, 2931   2937
    
    !!! WARNING !!! This is the total integrand of Eq. (13), which is e_cmd * n
    
    such a function is based on the exact behaviour of the Ecmd at large mu
    
    but with the exact on-top estimated with that of the UEG
    
    You enter with r(3), you get out with eps_c_md_PBE(1:N_states)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r`
       * :c:func:`ec_pbe_only`
       * :c:func:`grad_rho_ab_to_grad_rho_oc`
       * :c:func:`rho_ab_to_rho_oc`

 
.. c:function:: void ecorrlr(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ecorrlr(rs,z,mu,eclr)


    Long-range LDA correlation energy per particle of the uniform
    gas for the screened interaction, in Hartree, given rs (density
    parameter), z (spin polarization) and mu (range-separation
    parameter). The Coulomb part comes from ecPW, the screened part
    from Qrpa evaluated at x = mu * sqrt(rs) / phi. eclr is the
    result.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_lda_sr`
       * :c:func:`ec_only_lda_sr`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ecpw`

 
.. c:function:: void ecpw(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ecPW(x,y,ec,ecd,ecz,ecdd,eczd)


    Correlation energy per particle of the uniform gas in the limit
    mu = infinity, in Hartree: ec = ec(rs, zeta) with x -> rs and
    y -> zeta. ecd and ecz are the first derivatives with respect
    to rs and zeta, ecdd the second derivative with respect to rs,
    eczd is set to 0.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_lda`
       * :c:func:`ec_lda_sr`
       * :c:func:`ec_only_lda_sr`
       * :c:func:`ecorrlr`
       * :c:func:`vcorrlr`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpw`

 
.. c:function:: void eps_c_md_pbe_at_grid_pt(void)


    File : :file:`dft_utils_func/ecmd_pbe_ueg.irp.f`

    .. code:: fortran

        subroutine eps_c_md_PBE_at_grid_pt(mu,i_point,eps_c_md_PBE)


    provides the integrand of Eq. (13) of Phys.Chem.Lett.2019, 10, 2931   2937
    
    !!! WARNING !!! This is the total integrand of Eq. (13), which is e_cmd * n
    
    such a function is based on the exact behaviour of the Ecmd at large mu
    
    but with the exact on-top estimated with that of the UEG
    
    You enter with the alpha/beta density and density gradients
    
    You get out with eps_c_md_PBE(1:N_states)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_pbe_only`

 
.. c:function:: void eps_c_md_pbe_from_density(void)


    File : :file:`dft_utils_func/ecmd_pbe_ueg.irp.f`

    .. code:: fortran

        subroutine eps_c_md_PBE_from_density(mu,rho_a,rho_b, grad_rho_a, grad_rho_b,eps_c_md_PBE)


    provides the integrand of Eq. (13) of Phys.Chem.Lett.2019, 10, 2931   2937
    
    !!! WARNING !!! This is the total integrand of Eq. (13), which is e_cmd * n
    
    such a function is based on the exact behaviour of the Ecmd at large mu
    
    but with the exact on-top estimated with that of the UEG
    
    You enter with the alpha/beta density and density gradients
    
    You get out with eps_c_md_PBE(1:N_states)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_pbe_only`
       * :c:func:`grad_rho_ab_to_grad_rho_oc`
       * :c:func:`rho_ab_to_rho_oc`

 
.. c:function:: void esrc_md_ldaerf(void)


    File : :file:`dft_utils_func/ecmd_lda.irp.f`

    .. code:: fortran

        subroutine ESRC_MD_LDAERF (mu,rho_a,rho_b,dospin,e)


    Short-range spin-dependent LDA correlation energy density with a
    multideterminant reference, for OEP calculations. e is the sum of
    ec_only_lda_sr and the multideterminant correction
    DELTA_LRSR_LDAERF. mu is the range-separation parameter; rho_a and
    rho_b are the spin densities. If dospin is false, the Delta term is
    evaluated at zero spin polarization.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`delta_lrsr_ldaerf`
       * :c:func:`ec_only_lda_sr`

 
.. c:function:: void ex_lda(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ex_lda(rho_a,rho_b,ex,vx_a,vx_b)


    LDA exchange: energy density ex and potentials vx_a, vx_b
    computed from the alpha and beta spin densities rho_a and rho_b.
    ex is the sum of the cst_lda * rho**(4/3) terms for the two
    spins, and vx_a and vx_b are its derivatives with respect to
    rho_a and rho_b.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x_lda`

 
.. c:function:: void ex_lda_sr(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ex_lda_sr(mu,rho_a,rho_b,ex,vx_a,vx_b)


    Short-range LDA exchange: energy density ex and potentials
    vx_a, vx_b for the spin densities, with range-separation
    parameter mu. For each spin the expression branches on
    a = mu / (2 kF), kF = (6 pi**2 rho_s)**(1/3) from the spin
    density rho_s: small-a limit, exact expression, large-a
    expansion. ex is half of the sum of the two spin contributions.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:data:`energy_x_sr_lda`
       * :c:func:`ex_pbe_sr`
       * :c:func:`ex_pbe_sr_only`

 
.. c:function:: void ex_pbe_sr(void)


    File : :file:`dft_utils_func/routines_exc_sr_pbe.irp.f`

    .. code:: fortran

        subroutine ex_pbe_sr(mu,rho_a,rho_b,grd_rho_a_2,grd_rho_b_2,grd_rho_a_b,ex,vx_rho_a,vx_rho_b,vx_grd_rho_a_2,vx_grd_rho_b_2,vx_grd_rho_a_b)


    mu    = range separation parameter
    
    rho_a = density alpha
    
    rho_b = density beta
    
    grd_rho_a_2 = (gradient rho_a)^2
    
    grd_rho_b_2 = (gradient rho_b)^2
    
    grd_rho_a_b = (gradient rho_a).(gradient rho_b)
    
    ex = exchange energy density at the density and corresponding gradients of the density
    
    vx_rho_a = d ex / d rho_a
    
    vx_rho_b = d ex / d rho_b
    
    vx_grd_rho_a_2 = d ex / d grd_rho_a_2
    
    vx_grd_rho_b_2 = d ex / d grd_rho_b_2
    
    vx_grd_rho_a_b = d ex / d grd_rho_a_b

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gga_sr_type_functionals`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ex_lda_sr`

 
.. c:function:: void ex_pbe_sr_only(void)


    File : :file:`dft_utils_func/routines_exc_sr_pbe.irp.f`

    .. code:: fortran

        subroutine ex_pbe_sr_only(mu,rho_a,rho_b,grd_rho_a_2,grd_rho_b_2,grd_rho_a_b,ex)


    rho_a = density alpha
    
    rho_b = density beta
    
    grd_rho_a_2 = (gradient rho_a)^2
    
    grd_rho_b_2 = (gradient rho_b)^2
    
    grd_rho_a_b = (gradient rho_a).(gradient rho_b)
    
    ex = exchange energy density at point r

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ex_lda_sr`

 
.. c:function:: void g0_dg0(void)


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        subroutine g0_dg0(rho, rho_a, rho_b, g0, dg0drho)


    Give the on-top pair distribution function g0 and its derivative according to rho dg0drho

 
.. c:function:: void g0_dg0_d2g0(void)


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        subroutine g0_dg0_d2g0(rho, rho_a, rho_b, g0, dg0drho, d2g0drho2)


    Give the on-top pair distribution function g0 second derivative according to rho d2g0drho2

 
.. c:function:: void g0_ueg_mu(void)


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        double precision function g0_UEG_mu(mu,rho_a,rho_b)


    Pair distribution function g0(n_alpha,n_beta) of the UEG interacting with the long range interaction erf(mu r12)/r12
    
    Taken from P. Gori-Giorgi and A. Savin, Phys. Rev. A 73, 032506 (2006).

 
.. c:function:: void g0_ueg_mu_inf(void)


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        double precision function g0_UEG_mu_inf(rho_a,rho_b)


    Pair distribution function g0(n_alpha,n_beta) of the Colombic UEG
    
    Taken from Eq. (46)  P. Gori-Giorgi and A. Savin, Phys. Rev. A 73, 032506 (2006).

 
.. c:function:: void g0d(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function g0d(rs)



 
.. c:function:: void g0dd(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function g0dd(rs)



 
.. c:function:: void g0f(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function g0f(x)



 
.. c:function:: void g0s(void)


    File : :file:`dft_utils_func/ecmd_lda.irp.f`

    .. code:: fortran

        double precision function g0s(rs)



 
.. c:function:: void gga_sr_type_functionals(void)


    File : :file:`dft_utils_func/utils.irp.f`

    .. code:: fortran

        subroutine GGA_sr_type_functionals(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
                                        ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &
                                        ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b )


    routine that helps in building the x/c potentials on the AO basis for a GGA functional with a short-range interaction

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`energy_x_pbe`
       * :c:data:`energy_x_sr_pbe`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_pbe_sr`
       * :c:func:`ex_pbe_sr`
       * :c:func:`grad_rho_ab_to_grad_rho_oc`
       * :c:func:`rho_ab_to_rho_oc`
       * :c:func:`v_grad_rho_oc_to_v_grad_rho_ab`
       * :c:func:`v_rho_oc_to_v_rho_ab`

 
.. c:function:: void give_all_stuffs_in_r_for_lyp_88(void)


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        subroutine give_all_stuffs_in_r_for_lyp_88(r,rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)


    At the point r, gives the density and gradient quantities needed
    by the LYP functional ec_lyp_88: the spin densities rho_a and
    rho_b, the total density rho = rho_a + rho_b, the squared norms
    of the gradients of the spin densities, and grad_rho_2, the
    squared norm of the gradient of the total density. The arrays are
    indexed by the state istate = 1..N_states.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r`

 
.. c:function:: void gpw(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G,Gd,Gdd)


    Fit form G(x) used by ecPW, with x -> rs:
    G = -2 Ac (1 + alfa1 x) log(1 + 1/(2 Ac A(x))), where
    A(x) = beta1 x**0.5 + beta2 x + beta3 x**1.5 + beta4 x**2.
    Gd is dG/dx, Gdd is d2G/dx**2.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ecpw`

 
.. c:function:: void grad_rho_ab_to_grad_rho_oc(void)


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_o_2,grad_rho_c_2,grad_rho_o_c)


    convert (grad_rho_a_2, grad_rho_b_2, grad_rho_a.grad_rho_b, )
    
    to      (grad_rho_c_2, grad_rho_o_2, grad_rho_o.grad_rho_c)
    
    rho_c = total density, rho_o spin density

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_md_on_top_pbe_mu_corrected`
       * :c:func:`ec_md_pbe_on_top_general`
       * :c:func:`ecmd_pbe_ueg_at_r`
       * :c:func:`eps_c_md_pbe_from_density`
       * :c:func:`gga_sr_type_functionals`

 
.. c:function:: void h_func(void)


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        double precision function h_func(zeta)



 
.. c:function:: void mu_correction_of_on_top(void)


    File : :file:`dft_utils_func/ecmd_pbe_on_top.irp.f`

    .. code:: fortran

        double precision function mu_correction_of_on_top(mu,on_top)


    mu-based correction to the on-top pair density provided by the assymptotic expansion of
    
    P. Gori-Giorgi and A. Savin, Phys. Rev. A73, 032506 (2006)
    
    This is used in J. Chem. Phys.150, 084103 (2019); Eq. (29).

 
.. c:function:: void mu_grad_rho_func(void)


    File : :file:`dft_utils_func/mu_of_r_dft.irp.f`

    .. code:: fortran

        double precision function mu_grad_rho_func(r)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`density_and_grad_alpha_beta`

 
.. c:function:: void mu_rs_c(void)


    File : :file:`dft_utils_func/mu_of_r_dft.irp.f`

    .. code:: fortran

        double precision function mu_rs_c(rho)



 
.. c:function:: void qrpa(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function Qrpa(x)



 
.. c:function:: void qrpad(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function Qrpad(x)



 
.. c:function:: void qrpadd(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function Qrpadd(x)



 
.. c:function:: void rho_ab_to_rho_oc(void)


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)


    convert rho_alpha, rho_beta to rho_c, rho_o
    
    rho_c = total density, rho_o spin density

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_md_on_top_pbe_mu_corrected`
       * :c:func:`ec_md_pbe_on_top_general`
       * :c:func:`ecmd_pbe_ueg_at_r`
       * :c:func:`eps_c_md_pbe_from_density`
       * :c:func:`gga_sr_type_functionals`

 
.. c:function:: void rho_oc_to_rho_ab(void)


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine rho_oc_to_rho_ab(rho_o,rho_c,rho_a,rho_b)


    convert rho_c, rho_o to rho_alpha, rho_beta
    
    rho_c = total density, rho_o spin density

 
.. c:function:: void step_f(void)


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        double precision function step_f(x)



 
.. c:function:: void v_grad_rho_ab_to_v_grad_rho_oc(void)


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine v_grad_rho_ab_to_v_grad_rho_oc(v_grad_rho_a_2,v_grad_rho_b_2,v_grad_rho_a_b,v_grad_rho_o_2,v_grad_rho_c_2,v_grad_rho_o_c)


    convert  (v_grad_rho_a_2, v_grad_rho_b_2, v_grad_rho_a.grad_rho_b)
    
    to       (v_grad_rho_c_2, v_grad_rho_o_2, v_grad_rho_o.grad_rho_c)
    
    rho_c = total density, rho_o spin density

 
.. c:function:: void v_grad_rho_oc_to_v_grad_rho_ab(void)


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine v_grad_rho_oc_to_v_grad_rho_ab(v_grad_rho_o_2,v_grad_rho_c_2,v_grad_rho_o_c,v_grad_rho_a_2,v_grad_rho_b_2,v_grad_rho_a_b)


    convert  (v_grad_rho_c_2, v_grad_rho_o_2, v_grad_rho_o.grad_rho_c)
    
    to       (v_grad_rho_a_2, v_grad_rho_b_2, v_grad_rho_a.grad_rho_b)
    
    rho_c = total density, rho_o spin density

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gga_sr_type_functionals`

 
.. c:function:: void v_rho_ab_to_v_rho_oc(void)


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine v_rho_ab_to_v_rho_oc(v_rho_a,v_rho_b,v_rho_o,v_rho_c)


    convert v_rho_alpha, v_rho_beta to v_rho_c, v_rho_o
    
    rho_c = total density, rho_o spin density

 
.. c:function:: void v_rho_oc_to_v_rho_ab(void)


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine v_rho_oc_to_v_rho_ab(v_rho_o,v_rho_c,v_rho_a,v_rho_b)


    convert v_rho_alpha, v_rho_beta  to v_rho_c, v_rho_o
    
    rho_c = total density, rho_o spin density

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gga_sr_type_functionals`

 
.. c:function:: void vcorrlr(void)


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine vcorrlr(rs,z,mu,vclrup,vclrdown,vclrupd,vclrdownd)


    Long-range LDA correlation potentials for the screened
    interaction, in Hartree: vclrup (spin-up) and vclrdown
    (spin-down) for given rs, z and mu. vclrupd and vclrdownd are
    built from d(eclr)/drs and d2(eclr)/dr2.
    vclrupd and vclrdownd have no use (author, 2026-10-31); they
    are kept because external plugins may use them (author
    policy, 2026-10-31).

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_lda_sr`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ecpw`

