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
 
.. c:macro:: berf:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        function berf(a)



 
.. c:macro:: beta_rs:


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        double precision function beta_rs(rs)



 
.. c:macro:: correction_to_on_top_from_ueg:


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

 
.. c:macro:: dberfda:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        function dberfda(a)



 
.. c:macro:: delta_lrsr_ldaerf:


    File : :file:`dft_utils_func/ecmd_lda.irp.f`

    .. code:: fortran

        subroutine DELTA_LRSR_LDAERF (rhot,rhos,mu,dospin,e)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`esrc_md_ldaerf`

 
.. c:macro:: dpol:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function dpol(rs)



 
.. c:macro:: dpold:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function dpold(rs)



 
.. c:macro:: dpoldd:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function dpoldd(rs)



 
.. c:macro:: ec_lda:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ec_lda(rho_a,rho_b,ec,vc_a,vc_b)



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

 
.. c:macro:: ec_lda_sr:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ec_lda_sr(mu,rho_a,rho_b,ec,vc_a,vc_b)



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

 
.. c:macro:: ec_lyp2:


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        double precision function ec_lyp2(RhoA,RhoB,GA,GB,GAB)



 
.. c:macro:: ec_lyp_88:


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        double precision function ec_lyp_88(rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)


    LYP functional of the Lee, Yan, Parr, Phys. Rev B 1988, Vol 37, page 785.
    The expression used is the one by Miehlich, Savin, Stoll, Preuss, CPL, 1989 which gets rid of the laplacian of the density

 
.. c:macro:: ec_md_on_top_pbe_mu_corrected:


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

 
.. c:macro:: ec_md_pbe_on_top_general:


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

 
.. c:macro:: ec_only_lda_sr:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ec_only_lda_sr(mu,rho_a,rho_b,ec)



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

 
.. c:macro:: ec_pbe_only:


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

 
.. c:macro:: ec_pbe_sr:


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

 
.. c:macro:: ec_scan:


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        double precision function ec_scan(rho_a,rho_b,tau,grad_rho_2)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_only_lda_sr`

 
.. c:macro:: ecmd_pbe_ueg_at_r:


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

 
.. c:macro:: ecorrlr:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ecorrlr(rs,z,mu,eclr)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_lda_sr`
       * :c:func:`ec_only_lda_sr`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ecpw`

 
.. c:macro:: ecpw:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ecPW(x,y,ec,ecd,ecz,ecdd,eczd)



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

 
.. c:macro:: eps_c_md_pbe_at_grid_pt:


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

 
.. c:macro:: eps_c_md_pbe_from_density:


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

 
.. c:macro:: esrc_md_ldaerf:


    File : :file:`dft_utils_func/ecmd_lda.irp.f`

    .. code:: fortran

        subroutine ESRC_MD_LDAERF (mu,rho_a,rho_b,dospin,e)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`delta_lrsr_ldaerf`
       * :c:func:`ec_only_lda_sr`

 
.. c:macro:: ex_lda:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ex_lda(rho_a,rho_b,ex,vx_a,vx_b)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x_lda`

 
.. c:macro:: ex_lda_sr:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ex_lda_sr(mu,rho_a,rho_b,ex,vx_a,vx_b)



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

 
.. c:macro:: ex_pbe_sr:


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

 
.. c:macro:: ex_pbe_sr_only:


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

 
.. c:macro:: g0_dg0:


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        subroutine g0_dg0(rho, rho_a, rho_b, g0, dg0drho)


    Give the on-top pair distribution function g0 and its derivative according to rho dg0drho

 
.. c:macro:: g0_dg0_d2g0:


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        subroutine g0_dg0_d2g0(rho, rho_a, rho_b, g0, dg0drho, d2g0drho2)


    Give the on-top pair distribution function g0 second derivative according to rho d2g0drho2

 
.. c:macro:: g0_ueg_mu:


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        double precision function g0_UEG_mu(mu,rho_a,rho_b)


    Pair distribution function g0(n_alpha,n_beta) of the UEG interacting with the long range interaction erf(mu r12)/r12
    
    Taken from P. Gori-Giorgi and A. Savin, Phys. Rev. A 73, 032506 (2006).

 
.. c:macro:: g0_ueg_mu_inf:


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        double precision function g0_UEG_mu_inf(rho_a,rho_b)


    Pair distribution function g0(n_alpha,n_beta) of the Colombic UEG
    
    Taken from Eq. (46)  P. Gori-Giorgi and A. Savin, Phys. Rev. A 73, 032506 (2006).

 
.. c:macro:: g0d:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function g0d(rs)



 
.. c:macro:: g0dd:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function g0dd(rs)



 
.. c:macro:: g0f:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function g0f(x)



 
.. c:macro:: g0s:


    File : :file:`dft_utils_func/ecmd_lda.irp.f`

    .. code:: fortran

        double precision function g0s(rs)



 
.. c:macro:: gga_sr_type_functionals:


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

 
.. c:macro:: give_all_stuffs_in_r_for_lyp_88:


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        subroutine give_all_stuffs_in_r_for_lyp_88(r,rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r`

 
.. c:macro:: gpw:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G,Gd,Gdd)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ecpw`

 
.. c:macro:: grad_rho_ab_to_grad_rho_oc:


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

 
.. c:macro:: h_func:


    File : :file:`dft_utils_func/on_top_from_ueg.irp.f`

    .. code:: fortran

        double precision function h_func(zeta)



 
.. c:macro:: mu_correction_of_on_top:


    File : :file:`dft_utils_func/ecmd_pbe_on_top.irp.f`

    .. code:: fortran

        double precision function mu_correction_of_on_top(mu,on_top)


    mu-based correction to the on-top pair density provided by the assymptotic expansion of
    
    P. Gori-Giorgi and A. Savin, Phys. Rev. A73, 032506 (2006)
    
    This is used in J. Chem. Phys.150, 084103 (2019); Eq. (29).

 
.. c:macro:: mu_grad_rho_func:


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

 
.. c:macro:: mu_rs_c:


    File : :file:`dft_utils_func/mu_of_r_dft.irp.f`

    .. code:: fortran

        double precision function mu_rs_c(rho)



 
.. c:macro:: qrpa:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function Qrpa(x)



 
.. c:macro:: qrpad:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function Qrpad(x)



 
.. c:macro:: qrpadd:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function Qrpadd(x)



 
.. c:macro:: rho_ab_to_rho_oc:


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

 
.. c:macro:: rho_oc_to_rho_ab:


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine rho_oc_to_rho_ab(rho_o,rho_c,rho_a,rho_b)


    convert rho_c, rho_o to rho_alpha, rho_beta
    
    rho_c = total density, rho_o spin density

 
.. c:macro:: step_f:


    File : :file:`dft_utils_func/garbage_func.irp.f`

    .. code:: fortran

        double precision function step_f(x)



 
.. c:macro:: v_grad_rho_ab_to_v_grad_rho_oc:


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine v_grad_rho_ab_to_v_grad_rho_oc(v_grad_rho_a_2,v_grad_rho_b_2,v_grad_rho_a_b,v_grad_rho_o_2,v_grad_rho_c_2,v_grad_rho_o_c)


    convert  (v_grad_rho_a_2, v_grad_rho_b_2, v_grad_rho_a.grad_rho_b)
    
    to       (v_grad_rho_c_2, v_grad_rho_o_2, v_grad_rho_o.grad_rho_c)
    
    rho_c = total density, rho_o spin density

 
.. c:macro:: v_grad_rho_oc_to_v_grad_rho_ab:


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

 
.. c:macro:: v_rho_ab_to_v_rho_oc:


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine v_rho_ab_to_v_rho_oc(v_rho_a,v_rho_b,v_rho_o,v_rho_c)


    convert v_rho_alpha, v_rho_beta to v_rho_c, v_rho_o
    
    rho_c = total density, rho_o spin density

 
.. c:macro:: v_rho_oc_to_v_rho_ab:


    File : :file:`dft_utils_func/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine v_rho_oc_to_v_rho_ab(v_rho_o,v_rho_c,v_rho_a,v_rho_b)


    convert v_rho_alpha, v_rho_beta  to v_rho_c, v_rho_o
    
    rho_c = total density, rho_o spin density

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gga_sr_type_functionals`

 
.. c:macro:: vcorrlr:


    File : :file:`dft_utils_func/routines_exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine vcorrlr(rs,z,mu,vclrup,vclrdown,vclrupd,vclrdownd)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_lda_sr`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ecpw`

