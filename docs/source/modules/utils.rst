.. _module_utils: 
 
.. program:: utils 
 
.. default-role:: option 
 
=====
utils
=====

Contains general purpose utilities (sorting, maps, etc).

 
 
 
EZFIO parameters 
---------------- 
 
.. option:: restore_symm
 
    If true, try to find symmetry in the MO coefficient matrices
 
    Default: False
 
 
Providers 
--------- 
 
.. c:var:: au_to_d


    File : :file:`utils/units.irp.f`

    .. code:: fortran

        double precision	:: ha_to_ev	
        double precision	:: au_to_d	
        double precision	:: planck_cte	
        double precision	:: light_speed	
        double precision	:: ha_to_j	
        double precision	:: ha_to_nm	


    Some conversion between different units


 
.. c:var:: binom


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision, allocatable	:: binom	(0:40,0:40)
        double precision, allocatable	:: binom_transp	(0:40,0:40)


    Binomial coefficients

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`
       * :c:data:`dettocsftransformationmatrix`
       * :c:data:`nsomomax`
       * :c:data:`psi_csf_coef`

 
.. c:var:: binom_int


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        integer*8, allocatable	:: binom_int	(0:40,0:40)
        integer*8, allocatable	:: binom_int_transp	(0:40,0:40)


    Binomial coefficients, as integers*8

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`dominant_dets_of_cfgs`
       * :c:data:`n_dominant_dets_of_cfgs`
       * :c:data:`psi_configuration_to_psi_det`

 
.. c:var:: binom_int_transp


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        integer*8, allocatable	:: binom_int	(0:40,0:40)
        integer*8, allocatable	:: binom_int_transp	(0:40,0:40)


    Binomial coefficients, as integers*8

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`dominant_dets_of_cfgs`
       * :c:data:`n_dominant_dets_of_cfgs`
       * :c:data:`psi_configuration_to_psi_det`

 
.. c:var:: binom_transp


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision, allocatable	:: binom	(0:40,0:40)
        double precision, allocatable	:: binom_transp	(0:40,0:40)


    Binomial coefficients

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`
       * :c:data:`dettocsftransformationmatrix`
       * :c:data:`nsomomax`
       * :c:data:`psi_csf_coef`

 
.. c:var:: degree_max_integration_lebedev


    File : :file:`utils/angular_integration.irp.f`

    .. code:: fortran

        integer	:: degree_max_integration_lebedev	


    integrate correctly a polynom of order "degree_max_integration_lebedev"
    needed for the angular integration according to LEBEDEV formulae

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_integration_angular_lebedev`
       * :c:data:`theta_angular_integration_lebedev`

 
.. c:function:: dtranspose:


    File : :file:`utils/transpose.irp.f`

    .. code:: fortran

        recursive subroutine dtranspose(A,LDA,B,LDB,d1,d2)


    Transpose input matrix A into output matrix B

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`dtranspose`
       * :c:func:`h_s2_u_0_nstates_openmp`
       * :c:func:`h_s2_u_0_nstates_zmq`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp`
       * :c:func:`h_u_0_nstates_openmp`
       * :c:func:`h_u_0_nstates_zmq`
       * :c:func:`orb_range_2_rdm_openmp`
       * :c:func:`orb_range_2_rdm_state_av_openmp`
       * :c:func:`orb_range_2_trans_rdm_openmp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dtranspose`

 
.. c:var:: fact_inv


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fact_inv	(128)


    1/n!


 
.. c:function:: give_explicit_cpoly_and_cgaussian:


    File : :file:`utils/cgtos_utils.irp.f`

    .. code:: fortran

        subroutine give_explicit_cpoly_and_cgaussian(P_new, P_center, p, fact_k, iorder, &
               alpha, beta, a, b, Ae_center, Be_center, Ap_center, Bp_center, dim)


    
    Transforms the product of
    
      (x - x_Ap)^a(1) (x - x_Bp)^b(1) exp(-alpha (x - x_Ae)^2) exp(-beta (x - x_Be)^2) x
      (y - y_Ap)^a(2) (y - y_Bp)^b(2) exp(-alpha (y - y_Ae)^2) exp(-beta (y - y_Be)^2) x
      (z - z_Ap)^a(3) (z - z_Bp)^b(3) exp(-alpha (z - z_Ae)^2) exp(-beta (z - z_Be)^2)
    
    into
      fact_k * [sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x] exp (-p (x-P_center(1))^2)
             * [sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y] exp (-p (y-P_center(2))^2)
             * [sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z] exp (-p (z-P_center(3))^2)
    
    WARNING ::: IF fact_k is too smal then:
    returns a "s" function centered in zero
    with an inifinite exponent and a zero polynom coef
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_2e_cgtos_schwartz_accel`
       * :c:func:`ao_two_e_integral_cgtos`
       * :c:func:`overlap_cgaussian_xyz`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`cgaussian_product`
       * :c:func:`multiply_cpoly`
       * :c:func:`recentered_cpoly2`

 
.. c:function:: give_explicit_cpoly_and_cgaussian_x:


    File : :file:`utils/cgtos_utils.irp.f`

    .. code:: fortran

        subroutine give_explicit_cpoly_and_cgaussian_x(P_new, P_center, p, fact_k, iorder, &
                                               alpha, beta, a, b, Ae_center, Be_center, Ap_center, Bp_center, dim)


    
    Transform the product of
    
        (x - x_Ap)^a (x - x_Bp)^b exp(-alpha (r - Ae)^2) exp(-beta (r - Be)^2)
    
    into
    
        fact_k \sum_{i=0}^{iorder} (x - x_P)^i exp(-p (r - P)^2)
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_cgaussian_x`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`multiply_cpoly`
       * :c:func:`recentered_cpoly2`

 
.. c:var:: ha_to_ev


    File : :file:`utils/units.irp.f`

    .. code:: fortran

        double precision	:: ha_to_ev	
        double precision	:: au_to_d	
        double precision	:: planck_cte	
        double precision	:: light_speed	
        double precision	:: ha_to_j	
        double precision	:: ha_to_nm	


    Some conversion between different units


 
.. c:var:: ha_to_j


    File : :file:`utils/units.irp.f`

    .. code:: fortran

        double precision	:: ha_to_ev	
        double precision	:: au_to_d	
        double precision	:: planck_cte	
        double precision	:: light_speed	
        double precision	:: ha_to_j	
        double precision	:: ha_to_nm	


    Some conversion between different units


 
.. c:var:: ha_to_nm


    File : :file:`utils/units.irp.f`

    .. code:: fortran

        double precision	:: ha_to_ev	
        double precision	:: au_to_d	
        double precision	:: planck_cte	
        double precision	:: light_speed	
        double precision	:: ha_to_j	
        double precision	:: ha_to_nm	


    Some conversion between different units


 
.. c:var:: inv_int


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision, allocatable	:: inv_int	(128)


    1/i


 
.. c:var:: light_speed


    File : :file:`utils/units.irp.f`

    .. code:: fortran

        double precision	:: ha_to_ev	
        double precision	:: au_to_d	
        double precision	:: planck_cte	
        double precision	:: light_speed	
        double precision	:: ha_to_j	
        double precision	:: ha_to_nm	


    Some conversion between different units


 
.. c:var:: n_points_integration_angular_lebedev


    File : :file:`utils/angular_integration.irp.f`

    .. code:: fortran

        integer	:: n_points_integration_angular_lebedev	


    Number of points needed for the angular integral

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_integration_lebedev`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`theta_angular_integration_lebedev`

 
.. c:var:: nproc


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        integer	:: nproc	


    Number of current OpenMP threads

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`cholesky_ao_num`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`n_det`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`

 
.. c:function:: overlap_cgaussian_xyz:


    File : :file:`utils/cgtos_one_e.irp.f`

    .. code:: fortran

        subroutine overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                 Ap_center, Bp_center, overlap_x, overlap_y, overlap_z, overlap, dim)


    
    S_x = \int (x - Ap_x)^{a_x} exp(-\alpha (x - Ae_x)^2)
               (x - Bp_x)^{b_x} exp(-\beta  (x - Be_x)^2) dx
    
    S = S_x S_y S_z
    
    for  complex arguments
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos`
       * :c:data:`ao_deriv2_cgtos_x`
       * :c:data:`ao_overlap_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_cpoly_and_cgaussian`

 
.. c:var:: phi_angular_integration_lebedev


    File : :file:`utils/angular_integration.irp.f`

    .. code:: fortran

        double precision, allocatable	:: theta_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: phi_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: weights_angular_integration_lebedev	(n_points_integration_angular_lebedev)


    Theta phi values together with the weights values for the angular integration :
    integral [dphi,dtheta] f(x,y,z) = 4 * pi * sum (1<i<n_points_integration_angular_lebedev) f(xi,yi,zi)
    Note that theta and phi are in DEGREES !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_integration_lebedev`
       * :c:data:`n_points_integration_angular_lebedev`


 
.. c:var:: planck_cte


    File : :file:`utils/units.irp.f`

    .. code:: fortran

        double precision	:: ha_to_ev	
        double precision	:: au_to_d	
        double precision	:: planck_cte	
        double precision	:: light_speed	
        double precision	:: ha_to_j	
        double precision	:: ha_to_nm	


    Some conversion between different units


 
.. c:var:: qp_max_mem


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        integer	:: qp_max_mem	


    Maximum memory in Gb

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`file_lock`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha_chol`
       * :c:data:`cholesky_ao_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`pt2_j`
       * :c:data:`pt2_w`

 
.. c:var:: shiftfact_op5_inv


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision, allocatable	:: shiftfact_op5_inv	(128)


    
    1 / Gamma(n + 0.5)
    


 
.. c:var:: theta_angular_integration_lebedev


    File : :file:`utils/angular_integration.irp.f`

    .. code:: fortran

        double precision, allocatable	:: theta_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: phi_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: weights_angular_integration_lebedev	(n_points_integration_angular_lebedev)


    Theta phi values together with the weights values for the angular integration :
    integral [dphi,dtheta] f(x,y,z) = 4 * pi * sum (1<i<n_points_integration_angular_lebedev) f(xi,yi,zi)
    Note that theta and phi are in DEGREES !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_integration_lebedev`
       * :c:data:`n_points_integration_angular_lebedev`


 
.. c:function:: transpose:


    File : :file:`utils/transpose.irp.f`

    .. code:: fortran

        recursive subroutine transpose(A,LDA,B,LDB,d1,d2)


    Transpose input matrix A into output matrix B

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`transpose`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`transpose`

 
.. c:var:: weights_angular_integration_lebedev


    File : :file:`utils/angular_integration.irp.f`

    .. code:: fortran

        double precision, allocatable	:: theta_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: phi_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: weights_angular_integration_lebedev	(n_points_integration_angular_lebedev)


    Theta phi values together with the weights values for the angular integration :
    integral [dphi,dtheta] f(x,y,z) = 4 * pi * sum (1<i<n_points_integration_angular_lebedev) f(xi,yi,zi)
    Note that theta and phi are in DEGREES !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_integration_lebedev`
       * :c:data:`n_points_integration_angular_lebedev`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: add_cpoly:


    File : :file:`utils/cgtos_utils.irp.f`

    .. code:: fortran

        subroutine add_cpoly(b, nb, c, nc, d, nd)


    Add two complex polynomials
    D(t) =! D(t) +( B(t) + C(t))

 
.. c:function:: add_cpoly_multiply:


    File : :file:`utils/cgtos_utils.irp.f`

    .. code:: fortran

        subroutine add_cpoly_multiply(b, nb, cst, d, nd)


    Add a complex polynomial multiplied by a complex constant
    D(t) =! D(t) +( cst * B(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`general_primitive_integral_cgtos`

 
.. c:function:: add_poly:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine add_poly(b,nb,c,nc,d,nd)


    Add two polynomials
    D(t) =! D(t) +( B(t)+C(t))

 
.. c:function:: add_poly_multiply:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine add_poly_multiply(b,nb,cst,d,nd)


    Add a polynomial multiplied by a constant
    D(t) =! D(t) +( cst * B(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`general_primitive_integral`
       * :c:func:`general_primitive_integral_erf`

 
.. c:function:: apply_rotation:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine apply_rotation(A,LDA,R,LDR,B,LDB,m,n)


    Apply the rotation found by find_rotation

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: approx_dble:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function approx_dble(a,n)



 
.. c:function:: binom_func:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function binom_func(i,j)


    .. math                       ::
    
      \frac{i!}{j!(i-j)!}
    

 
.. c:function:: cgaussian_product:


    File : :file:`utils/cgtos_utils.irp.f`

    .. code:: fortran

        subroutine cgaussian_product(a, xa, b, xb, k, p, xp)


    complex Gaussian product
    e^{-a (r-r_A)^2} e^{-b (r-r_B)^2} = k e^{-p (r-r_P)^2}

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_cpoly_and_cgaussian`

 
.. c:function:: cgaussian_product_x:


    File : :file:`utils/cgtos_utils.irp.f`

    .. code:: fortran

        subroutine cgaussian_product_x(a, xa, b, xb, k, p, xp)


    complex Gaussian product in 1D.
    e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K e^{-p (x-x_P)^2}

 
.. c:function:: check_mem:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        subroutine check_mem(rss_in,routine)


    Checks if n gigabytes can be allocated. If not, exit the run.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_max_mem`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha_chol`
       * :c:func:`create_selection_buffer`
       * :c:func:`dav_double_dressed`
       * :c:func:`davidson_diag_csf_hjj`
       * :c:func:`davidson_diag_hjj`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`davidson_general`
       * :c:func:`davidson_general_diag_dressed_ext_rout_nonsym_b1space`
       * :c:func:`davidson_general_ext_rout`
       * :c:func:`davidson_general_ext_rout_diag_dressed`
       * :c:func:`davidson_general_ext_rout_dressed`
       * :c:func:`davidson_general_ext_rout_nonsym_b1space`
       * :c:func:`make_selection_buffer_s2`
       * :c:func:`pt2_collector`
       * :c:data:`pt2_j`
       * :c:data:`pt2_w`
       * :c:func:`remove_duplicates_in_selection_buffer`
       * :c:func:`run_cipsi`
       * :c:func:`run_slave_main`
       * :c:func:`run_stochastic_cipsi`
       * :c:func:`selection_collector`
       * :c:func:`testteethbuilding`
       * :c:func:`zmq_pt2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_memory_usage`
       * :c:func:`resident_memory`

 
.. c:function:: check_sym:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine check_sym(A, n)



 
.. c:function:: cpx_erf:


    File : :file:`utils/cpx_erf.irp.f`

    .. code:: fortran

        complex*16 function cpx_erf(x, y)


    
    compute erf(z) for z = x + i y
    
    REF: Abramowitz and Stegun
    

 
.. c:function:: cpx_erf_1:


    File : :file:`utils/cpx_erf.irp.f`

    .. code:: fortran

        complex*16 function cpx_erf_1(x, y)


    
    compute erf(z) for z = x + i y
    
    REF: Abramowitz and Stegun
    

 
.. c:function:: crint:


    File : :file:`utils/cgtos_utils.irp.f`

    .. code:: fortran

        complex*16 function crint(n, rho)



 
.. c:function:: crint_1:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        complex*16 function crint_1(n, rho)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zboysfun00_1`

 
.. c:function:: crint_1_vec:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        subroutine crint_1_vec(n_max, rho, vals)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`crint_sum`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`crint_smallz_vec`
       * :c:func:`zboysfun00_1`

 
.. c:function:: crint_2:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        complex*16 function crint_2(n, rho)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zboysfun`
       * :c:func:`zboysfunnrp`

 
.. c:function:: crint_2_vec:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        subroutine crint_2_vec(n_max, rho, vals)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`crint_smallz_vec`
       * :c:func:`zboysfun`
       * :c:func:`zboysfunnrp`

 
.. c:function:: crint_quad_1:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        subroutine crint_quad_1(n, rho, n_quad, crint_quad)



 
.. c:function:: crint_quad_12:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        subroutine crint_quad_12(n, rho, n_quad, crint_quad)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`crint_quad_12_vec`

 
.. c:function:: crint_quad_12_vec:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        subroutine crint_quad_12_vec(n_max, rho, vals)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`crint_quad_12`

 
.. c:function:: crint_quad_2:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        subroutine crint_quad_2(n, rho, n_quad, crint_quad)



 
.. c:function:: crint_smallz:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        complex*16 function crint_smallz(n, rho)


    Standard version of rint

 
.. c:function:: crint_smallz_vec:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        subroutine crint_smallz_vec(n_max, rho, vals)


    Standard version of rint

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`crint_1_vec`
       * :c:func:`crint_2_vec`

 
.. c:function:: crint_sum:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        complex*16 function crint_sum(n_pt_out, rho, d1)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`crint_1_vec`

 
.. c:function:: dble_fact:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function dble_fact(n)



 
.. c:function:: dble_fact_even:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function dble_fact_even(n) result(fact2)


    n!!

 
.. c:function:: dble_fact_odd:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function dble_fact_odd(n) result(fact2)


    n!!

 
.. c:function:: dble_logfact:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function dble_logfact(n) result(logfact2)


    n!!

 
.. c:function:: derf_mu_x:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function derf_mu_x(mu,x)



 
.. c:function:: diag_mat_per_fock_degen:


    File : :file:`utils/block_diag_degen.irp.f`

    .. code:: fortran

        subroutine diag_mat_per_fock_degen(fock_diag, mat_ref, n, thr_d, thr_nd, thr_deg, leigvec, reigvec, eigval)


    
    subroutine that diagonalizes a matrix mat_ref BY BLOCK
    
    the blocks are defined by the elements having the SAME DEGENERACIES in the entries "fock_diag"
    
    examples : all elements having degeneracy 1 in fock_diag (i.e. not being degenerated) will be treated together
    
             : all elements having degeneracy 2 in fock_diag (i.e. two elements are equal) will be treated together
    
             : all elements having degeneracy 3 in fock_diag (i.e. two elements are equal) will be treated together
    
    etc... the advantage is to guarentee no spurious mixing because of numerical problems.
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsort`
       * :c:func:`give_degen_full_list`
       * :c:func:`isort`
       * :c:func:`non_hrmt_bieig`

 
.. c:function:: diag_mat_per_fock_degen_core:


    File : :file:`utils/block_diag_degen_core.irp.f`

    .. code:: fortran

        subroutine diag_mat_per_fock_degen_core(fock_diag, mat_ref, listcore,ncore, n, thr_d, thr_nd, thr_deg, leigvec, reigvec, eigval)


    
    subroutine that diagonalizes a matrix mat_ref BY BLOCK
    
    the blocks are defined by the elements having the SAME DEGENERACIES in the entries "fock_diag"
    
    the elements of listcore are untouched
    
    examples : all elements having degeneracy 1 in fock_diag (i.e. not being degenerated) will be treated together
    
             : all elements having degeneracy 2 in fock_diag (i.e. two elements are equal) will be treated together
    
             : all elements having degeneracy 3 in fock_diag (i.e. two elements are equal) will be treated together
    
    etc... the advantage is to guarentee no spurious mixing because of numerical problems.
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsort`
       * :c:func:`give_degen_full_listcore`
       * :c:func:`isort`
       * :c:func:`non_hrmt_bieig`

 
.. c:function:: diag_nonsym_right:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine diag_nonsym_right(n, A, A_ldim, V, V_ldim, energy, E_ldim)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`davidson_general_diag_dressed_ext_rout_nonsym_b1space`
       * :c:func:`davidson_general_ext_rout_nonsym_b1space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgeevx`
       * :c:func:`dsort`

 
.. c:function:: diagonalize_sym_matrix:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine diagonalize_sym_matrix(N, A, e)


    
    Diagonalize a symmetric matrix
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsyev`

 
.. c:function:: dset_order:


    File : :file:`utils/sort.irp.f_template_90`

    .. code:: fortran

        subroutine dset_order(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos_ord`
       * :c:data:`ao_coef_normalized_ordered`
       * :c:func:`h_s2_u_0_nstates_openmp`
       * :c:func:`h_s2_u_0_nstates_zmq`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp`
       * :c:func:`h_u_0_nstates_openmp`
       * :c:func:`h_u_0_nstates_zmq`
       * :c:func:`orb_range_2_rdm_openmp`
       * :c:func:`orb_range_2_rdm_state_av_openmp`
       * :c:func:`orb_range_2_trans_rdm_openmp`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:func:`restore_symmetry`

 
.. c:function:: dset_order_big:


    File : :file:`utils/sort.irp.f_template_90`

    .. code:: fortran

        subroutine dset_order_big(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: eigsvd:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine eigSVD(A,LDA,U,LDU,D,Vt,LDVt,m,n)


    Algorithm 3 of https://arxiv.org/pdf/1810.06860.pdf
    
    A(m,n) = U(m,n) D(n) Vt(n,n) with m>n

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dscal`
       * :c:func:`lapack_diagd`
       * :c:func:`svd`

 
.. c:function:: erf_e:


    File : :file:`utils/cpx_erf.irp.f`

    .. code:: fortran

        complex*16 function erf_E(x, yabs)



 
.. c:function:: erf_f:


    File : :file:`utils/cpx_erf.irp.f`

    .. code:: fortran

        double precision function erf_F(x, yabs)



 
.. c:function:: erf_g:


    File : :file:`utils/cpx_erf.irp.f`

    .. code:: fortran

        complex*16 function erf_G(x, yabs)



 
.. c:function:: erf_h:


    File : :file:`utils/cpx_erf.irp.f`

    .. code:: fortran

        complex*16 function erf_H(x, yabs)



 
.. c:function:: exp_matrix:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine exp_matrix(X,n,exp_X)


    exponential of the matrix X: X has to be ANTI HERMITIAN !!
    
    taken from Hellgaker, jorgensen, Olsen book
    
    section evaluation of matrix exponential (Eqs. 3.1.29 to 3.1.31)

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_a_squared`
       * :c:func:`lapack_diagd`

 
.. c:function:: exp_matrix_taylor:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine exp_matrix_taylor(X,n,exp_X,converged)


    exponential of a general real matrix X using the Taylor expansion of exp(X)
    
    returns the logical converged which checks the convergence
    exponential of X using Taylor expansion

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`umat`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: extrapolate_data:


    File : :file:`utils/extrapolation.irp.f`

    .. code:: fortran

        subroutine extrapolate_data(N_data, data, pt2, output)


    Extrapolate the data to the FCI limit

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`increment_n_iter`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_pseudo_inverse`

 
.. c:function:: f_integral:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function F_integral(n,p)


    function that calculates the following integral
    \int_{\-infty}^{+\infty} x^n \exp(-p x^2) dx

 
.. c:function:: fact:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function fact(n)


    n!

 
.. c:function:: fc_integral:


    File : :file:`utils/cgtos_utils.irp.f`

    .. code:: fortran

        complex*16 function Fc_integral(n, inv_sq_p)


    function that calculates the following integral
    \int_{\-infty}^{+\infty} x^n \exp(-p x^2) dx
    for complex valued p

 
.. c:function:: find_rotation:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine find_rotation(A,LDA,B,m,C,n)


    Find A.C = B

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_pseudo_inverse`

 
.. c:function:: format_w_error:


    File : :file:`utils/format_w_error.irp.f`

    .. code:: fortran

        subroutine format_w_error(value,error,size_nb,max_nb_digits,format_value,str_error)


    Format for double precision, value(error)

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pt2_collector`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lock_io`
       * :c:func:`unlock_io`

 
.. c:function:: gaussian_product:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine gaussian_product(a,xa,b,xb,k,p,xp)


    Gaussian product in 1D.
    e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian`
       * :c:func:`give_explicit_poly_and_gaussian_double`

 
.. c:function:: gaussian_product_v:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine gaussian_product_v(a, xa, LD_xa, b, xb, k, p, xp, n_points)


    
    Gaussian product in 1D.
    e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
    
    Using multiple A centers
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian_v`

 
.. c:function:: gaussian_product_x:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine gaussian_product_x(a,xa,b,xb,k,p,xp)


    Gaussian product in 1D.
    e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_gaussian_xyz`

 
.. c:function:: gaussian_product_x_v:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine gaussian_product_x_v(a,xa,b,xb,k,p,xp,n_points)


    Gaussian product in 1D with multiple xa
    e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}

 
.. c:function:: get_a_squared:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine get_A_squared(A,n,A2)


    A2 = A A where A is n x n matrix. Use the dgemm routine

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`exp_matrix`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: get_ab_prod:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine get_AB_prod(A,n,m,B,l,AB)


    AB = A B where A is n x m, B is m x l. Use the dgemm routine

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: get_inverse:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine get_inverse(A,LDA,m,C,LDC)


    Returns the inverse of the square matrix A

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_canonical_coef_inv`
       * :c:data:`overlap_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgetrf`
       * :c:func:`dgetri`

 
.. c:function:: get_inverse_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine get_inverse_complex(A,LDA,m,C,LDC)


    Returns the inverse of the square matrix A

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zgetrf`
       * :c:func:`zgetri`

 
.. c:function:: get_pseudo_inverse:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine get_pseudo_inverse(A, LDA, m, n, C, LDC, cutoff)


    Find C = A^-1

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_inv`
       * :c:func:`extrapolate_data`
       * :c:func:`find_rotation`
       * :c:data:`s_inv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dgesvd`

 
.. c:function:: get_pseudo_inverse_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine get_pseudo_inverse_complex(A,LDA,m,n,C,LDC,cutoff)


    Find C = A^-1

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`s_inv_complex`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zgesvd`

 
.. c:function:: get_total_available_memory:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        integer function get_total_available_memory() result(res)


    Returns the total available memory on the current machine

 
.. c:function:: give_degen:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine give_degen(A, n, shift, list_degen, n_degen_list)


    returns n_degen_list :: the number of degenerated SET of elements (i.e. with |A(i)-A(i+1)| below shift)
    
    for each of these sets, list_degen(1,i) = first degenerate element of the set i,
    
                            list_degen(2,i) = last degenerate element of the set i.

 
.. c:function:: give_degen_full_list:


    File : :file:`utils/block_diag_degen.irp.f`

    .. code:: fortran

        subroutine give_degen_full_list(A, n, thr, list_degen, n_degen_list)


    you enter with an array A(n) and spits out all the elements degenerated up to thr
    
    the elements of A(n) DON'T HAVE TO BE SORTED IN THE ENTRANCE: TOTALLY GENERAL
    
    list_degen(i,0) = number of degenerate entries
    
    list_degen(i,1) = index of the first degenerate entry
    
    list_degen(i,2:list_degen(i,0)) = list of all other dengenerate entries
    
    if list_degen(i,0) == 1 it means that there is no degeneracy for that element

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_mat_per_fock_degen`

 
.. c:function:: give_degen_full_listcore:


    File : :file:`utils/block_diag_degen_core.irp.f`

    .. code:: fortran

        subroutine give_degen_full_listcore(A, n, listcore, ncore, thr, list_degen, n_degen_list)


    you enter with an array A(n) and spits out all the elements degenerated up to thr
    
    the elements of A(n) DON'T HAVE TO BE SORTED IN THE ENTRANCE: TOTALLY GENERAL
    
    list_degen(i,0) = number of degenerate entries
    
    list_degen(i,1) = index of the first degenerate entry
    
    list_degen(i,2:list_degen(i,0)) = list of all other dengenerate entries
    
    if list_degen(i,0) == 1 it means that there is no degeneracy for that element
    
    if list_degen(i,0) >= 1000 it means that it is core orbitals

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_mat_per_fock_degen_core`

 
.. c:function:: give_explicit_poly_and_gaussian:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine give_explicit_poly_and_gaussian(P_new,P_center,p,fact_k,iorder,alpha,beta,a,b,A_center,B_center,dim)


    Transforms the product of
             (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
    into
           fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
                  * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
                  * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
    
    WARNING ::: IF fact_k is too smal then:
    returns a "s" function centered in zero
    with an inifinite exponent and a zero polynom coef

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integral`
       * :c:func:`ao_two_e_integral_erf`
       * :c:func:`ao_two_e_integral_schwartz_accel`
       * :c:func:`ao_two_e_integral_schwartz_accel_erf`
       * :c:func:`give_explicit_poly_and_gaussian_double`
       * :c:func:`overlap_gaussian_xyz`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gaussian_product`
       * :c:func:`multiply_poly`
       * :c:func:`recentered_poly2`

 
.. c:function:: give_explicit_poly_and_gaussian_double:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine give_explicit_poly_and_gaussian_double(P_new,P_center,p,fact_k,iorder,alpha,beta,gama,a,b,A_center,B_center,Nucl_center,dim)


    Transforms the product of
             (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3)
             exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta) exp(-(r-Nucl_center)^2 gama
    
    into
           fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
                  * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
                  * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gaussian_product`
       * :c:func:`give_explicit_poly_and_gaussian`

 
.. c:function:: give_explicit_poly_and_gaussian_v:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine give_explicit_poly_and_gaussian_v(P_new, ldp, P_center, p, fact_k, iorder, alpha, beta, a, b, A_center, LD_A, B_center, n_points)


    Transforms the product of
             (x-x_A)^a(1) (x-x_B)^b(1) (y-y_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
    into
           fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
                  * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
                  * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
    
    WARNING                      :: : IF fact_k is too smal then:
    returns a "s" function centered in zero
    with an inifinite exponent and a zero polynom coef

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_gaussian_xyz_v`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gaussian_product_v`
       * :c:func:`multiply_poly_v`
       * :c:func:`recentered_poly2_v`
       * :c:func:`recentered_poly2_v0`

 
.. c:function:: give_explicit_poly_and_gaussian_x:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine give_explicit_poly_and_gaussian_x(P_new,P_center,p,fact_k,iorder,alpha,beta,a,b,A_center,B_center,dim)


    Transform the product of
             (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
    into
           fact_k  (x-x_P)^iorder(1)  (y-y_P)^iorder(2)  (z-z_P)^iorder(3) exp(-p(r-P)^2)

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_gaussian_x`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`multiply_poly`
       * :c:func:`recentered_poly2`

 
.. c:function:: give_pol_in_r:


    File : :file:`utils/prim_in_r.irp.f`

    .. code:: fortran

        double precision function give_pol_in_r(r,pol,center, alpha,iorder, max_dim)



 
.. c:function:: hermite:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function hermite(n,x)


    Hermite polynomial

 
.. c:function:: i2set_order:


    File : :file:`utils/sort.irp.f_template_90`

    .. code:: fortran

        subroutine i2set_order(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.

 
.. c:function:: i2set_order_big:


    File : :file:`utils/sort.irp.f_template_90`

    .. code:: fortran

        subroutine i2set_order_big(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: i8set_order:


    File : :file:`utils/sort.irp.f_template_90`

    .. code:: fortran

        subroutine i8set_order(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.

 
.. c:function:: i8set_order_big:


    File : :file:`utils/sort.irp.f_template_90`

    .. code:: fortran

        subroutine i8set_order_big(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: is_same_spin:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        logical function is_same_spin(sigma_1, sigma_2)


    
    true if sgn(sigma_1) = sgn(sigma_2)
    

 
.. c:function:: iset_order:


    File : :file:`utils/sort.irp.f_template_90`

    .. code:: fortran

        subroutine iset_order(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:func:`restore_symmetry`

 
.. c:function:: iset_order_big:


    File : :file:`utils/sort.irp.f_template_90`

    .. code:: fortran

        subroutine iset_order_big(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: kronecker_delta:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        function Kronecker_delta(i, j) result(delta)


    Kronecker Delta

 
.. c:function:: lapack_diag:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine lapack_diag(eigvalues,eigvectors,H,nmax,n)


    Diagonalize matrix H
    
    H is untouched between input and ouptut
    
    eigevalues(i) = ith lowest eigenvalue of the H matrix
    
    eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:func:`davidson_general`
       * :c:func:`davidson_general_ext_rout`
       * :c:func:`davidson_general_ext_rout_diag_dressed`
       * :c:data:`difference_dm_eigvect`
       * :c:func:`mo_as_eigvectors_of_mo_matrix`
       * :c:data:`multi_s_x_dipole_moment_eigenvec`
       * :c:data:`psi_coef_cas_diagonalized`
       * :c:data:`sxeigenvec`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsyev`

 
.. c:function:: lapack_diag_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine lapack_diag_complex(eigvalues,eigvectors,H,nmax,n)


    Diagonalize matrix H (complex)
    
    H is untouched between input and ouptut
    
    eigevalues(i) = ith lowest eigenvalue of the H matrix
    
    eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zheev`

 
.. c:function:: lapack_diagd:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine lapack_diagd(eigvalues,eigvectors,H,nmax,n)


    Diagonalize matrix H
    
    H is untouched between input and ouptut
    
    eigevalues(i) = ith lowest eigenvalue of the H matrix
    
    eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`eigsvd`
       * :c:func:`exp_matrix`
       * :c:data:`inertia_tensor_eigenvectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsyevd`

 
.. c:function:: lapack_diagd_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine lapack_diagd_complex(eigvalues,eigvectors,H,nmax,n)


    Diagonalize matrix H(complex)
    
    H is untouched between input and ouptut
    
    eigevalues(i) = ith lowest eigenvalue of the H matrix
    
    eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zheevd`

 
.. c:function:: lapack_diagd_diag_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine lapack_diagd_diag_complex(eigvalues,eigvectors,H,nmax,n)


    Diagonalize matrix H(complex)
    
    H is untouched between input and ouptut
    
    eigevalues(i) = ith lowest eigenvalue of the H matrix
    
    eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zheev`
       * :c:func:`zheevd`

 
.. c:function:: lapack_diagd_diag_in_place_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine lapack_diagd_diag_in_place_complex(eigvalues,eigvectors,nmax,n)


    Diagonalize matrix H(complex)
    
    H is untouched between input and ouptut
    
    eigevalues(i) = ith lowest eigenvalue of the H matrix
    
    eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zheev`
       * :c:func:`zheevd`

 
.. c:function:: logfact:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function logfact(n)


    n!

 
.. c:function:: lowercase:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine lowercase(txt,n)


    Transform to lower case

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`end_parallel_job`
       * :c:func:`new_parallel_job`

 
.. c:function:: map_load_from_disk:


    File : :file:`utils/map_functions.irp.f`

    .. code:: fortran

        subroutine map_load_from_disk(filename,map)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`c_f_pointer`
       * :c:func:`mmap`

 
.. c:function:: map_save_to_disk:


    File : :file:`utils/map_functions.irp.f`

    .. code:: fortran

        subroutine map_save_to_disk(filename,map)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:func:`save_erf_two_e_integrals_ao`
       * :c:func:`save_erf_two_e_integrals_mo`
       * :c:func:`save_erf_two_e_ints_ao_into_ints_ao`
       * :c:func:`save_erf_two_e_ints_mo_into_ints_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`c_f_pointer`
       * :c:func:`map_sort`
       * :c:func:`mmap`
       * :c:func:`msync`

 
.. c:function:: matrix_vector_product_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine matrix_vector_product_complex(u0,u1,matrix,sze,lda)


    performs u1 =! performs u1 +( u0 * matrix)

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zhemv`

 
.. c:function:: memory_of_double:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        double precision function memory_of_double(n)


    Computes the memory required for n double precision elements in gigabytes.

 
.. c:function:: memory_of_double8:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        double precision function memory_of_double8(n)


    Computes the memory required for n double precision elements in gigabytes.

 
.. c:function:: memory_of_int:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        double precision function memory_of_int(n)


    Computes the memory required for n double precision elements in gigabytes.

 
.. c:function:: memory_of_int8:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        double precision function memory_of_int8(n)


    Computes the memory required for n double precision elements in gigabytes.

 
.. c:function:: multiply_cpoly:


    File : :file:`utils/cgtos_utils.irp.f`

    .. code:: fortran

        subroutine multiply_cpoly(b, nb, c, nc, d, nd)


    Multiply two complex polynomials
    D(t) =! D(t) +( B(t) * C(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`general_primitive_integral_cgtos`
       * :c:func:`give_cpolynomial_mult_center_one_e`
       * :c:func:`give_explicit_cpoly_and_cgaussian`
       * :c:func:`give_explicit_cpoly_and_cgaussian_x`
       * :c:func:`i_x1_pol_mult_a1_cgtos`
       * :c:func:`i_x1_pol_mult_a2_cgtos`
       * :c:func:`i_x1_pol_mult_one_e_cgtos`
       * :c:func:`i_x1_pol_mult_recurs_cgtos`
       * :c:func:`i_x2_pol_mult_cgtos`
       * :c:func:`i_x2_pol_mult_one_e_cgtos`

 
.. c:function:: multiply_poly:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine multiply_poly(b,nb,c,nc,d,nd)


    Multiply two polynomials
    D(t) =! D(t) +( B(t)*C(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`general_primitive_integral_erf`
       * :c:func:`give_explicit_poly_and_gaussian`
       * :c:func:`give_explicit_poly_and_gaussian_x`
       * :c:func:`give_polynomial_mult_center_one_e`
       * :c:func:`give_polynomial_mult_center_one_e_erf`
       * :c:func:`give_polynomial_mult_center_one_e_erf_opt`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`multiply_poly_b0`
       * :c:func:`multiply_poly_b1`
       * :c:func:`multiply_poly_b2`
       * :c:func:`multiply_poly_c0`
       * :c:func:`multiply_poly_c1`
       * :c:func:`multiply_poly_c2`

 
.. c:function:: multiply_poly_b0:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine multiply_poly_b0(b,c,nc,d,nd)


    Multiply two polynomials
    D(t) =! D(t) +( B(t)*C(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`multiply_poly`

 
.. c:function:: multiply_poly_b1:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine multiply_poly_b1(b,c,nc,d,nd)


    Multiply two polynomials
    D(t) =! D(t) +( B(t)*C(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`multiply_poly`

 
.. c:function:: multiply_poly_b2:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine multiply_poly_b2(b,c,nc,d,nd)


    Multiply two polynomials
    D(t) =! D(t) +( B(t)*C(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`multiply_poly`

 
.. c:function:: multiply_poly_c0:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine multiply_poly_c0(b,nb,c,d,nd)


    Multiply two polynomials
    D(t) =! D(t) +( B(t)*C(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`multiply_poly`

 
.. c:function:: multiply_poly_c1:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine multiply_poly_c1(b,nb,c,d,nd)


    Multiply two polynomials
    D(t) =! D(t) +( B(t)*C(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`multiply_poly`

 
.. c:function:: multiply_poly_c2:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine multiply_poly_c2(b,nb,c,d,nd)


    Multiply two polynomials
    D(t) =! D(t) +( B(t)*C(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e`
       * :c:func:`i_x2_pol_mult_one_e`
       * :c:func:`multiply_poly`

 
.. c:function:: multiply_poly_v:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine multiply_poly_v(b,nb,c,nc,d,nd,n_points)


    Multiply pairs of polynomials
    D(t) =! D(t) +( B(t)*C(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian_v`

 
.. c:function:: normalize:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine normalize(u,sze)


    Normalizes vector u

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`copy_h_apply_buffer_to_wf`
       * :c:func:`dav_double_dressed`
       * :c:func:`davidson_diag_csf_hjj`
       * :c:func:`davidson_diag_hjj`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`davidson_general`
       * :c:func:`davidson_general_diag_dressed_ext_rout_nonsym_b1space`
       * :c:func:`davidson_general_ext_rout`
       * :c:func:`davidson_general_ext_rout_diag_dressed`
       * :c:func:`davidson_general_ext_rout_dressed`
       * :c:func:`davidson_general_ext_rout_nonsym_b1space`
       * :c:func:`save_wavefunction_general`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dscal`

 
.. c:function:: nullify_small_elements:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine nullify_small_elements(m,n,A,LDA,thresh)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:func:`dav_double_dressed`
       * :c:func:`davidson_diag_csf_hjj`
       * :c:func:`davidson_diag_hjj`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`davidson_general_ext_rout_dressed`
       * :c:func:`hcore_guess`
       * :c:data:`mo_one_e_integrals`
       * :c:func:`save_natural_mos`
       * :c:func:`save_natural_mos_canon_label`
       * :c:func:`save_natural_mos_no_ov_rot`
       * :c:func:`save_wavefunction_truncated`

 
.. c:function:: ortho_canonical:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_canonical(overlap,LDA,N,C,LDC,m,cutoff)


    Compute C_new=C_old.U.s^-1/2 canonical orthogonalization.
    
    overlap : overlap matrix
    
    LDA : leftmost dimension of overlap array
    
    N : Overlap matrix is NxN (array is (LDA,N) )
    
    C : Coefficients of the vectors to orthogonalize. On exit,
        orthogonal vectors
    
    LDC : leftmost dimension of C
    
    m : Coefficients matrix is MxN, ( array is (LDC,N) )
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_canonical_coef`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`svd`

 
.. c:function:: ortho_canonical_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_canonical_complex(overlap,LDA,N,C,LDC,m,cutoff)


    Compute C_new=C_old.U.s^-1/2 canonical orthogonalization.
    
    overlap : overlap matrix
    
    LDA : leftmost dimension of overlap array
    
    N : Overlap matrix is NxN (array is (LDA,N) )
    
    C : Coefficients of the vectors to orthogonalize. On exit,
        orthogonal vectors
    
    LDC : leftmost dimension of C
    
    m : Coefficients matrix is MxN, ( array is (LDC,N) )
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`svd_complex`
       * :c:func:`zgemm`

 
.. c:function:: ortho_lowdin:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_lowdin(overlap,LDA,N,C,LDC,m,cutoff)


    Compute C_new=C_old.S^-1/2 orthogonalization.
    
    overlap : overlap matrix
    
    LDA : leftmost dimension of overlap array
    
    N : Overlap matrix is NxN (array is (LDA,N) )
    
    C : Coefficients of the vectors to orthogonalize. On exit,
        orthogonal vectors
    
    LDC : leftmost dimension of C
    
    M : Coefficients matrix is MxN, ( array is (LDC,N) )
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_lowdin_coef`
       * :c:func:`orthonormalize_mos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`svd`

 
.. c:function:: ortho_lowdin_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_lowdin_complex(overlap,LDA,N,C,LDC,m,cutoff)


    Compute C_new=C_old.S^-1/2 orthogonalization.
    
    overlap : overlap matrix
    
    LDA : leftmost dimension of overlap array
    
    N : Overlap matrix is NxN (array is (LDA,N) )
    
    C : Coefficients of the vectors to orthogonalize. On exit,
        orthogonal vectors
    
    LDC : leftmost dimension of C
    
    M : Coefficients matrix is MxN, ( array is (LDC,N) )
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`svd_complex`
       * :c:func:`zgemm`

 
.. c:function:: ortho_qr:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_qr(A,LDA,m,n)


    Orthogonalization using Q.R factorization
    
    A : matrix to orthogonalize
    
    LDA : leftmost dimension of A
    
    m : Number of rows of A
    
    n : Number of columns of A
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`davidson_general`
       * :c:func:`davidson_general_diag_dressed_ext_rout_nonsym_b1space`
       * :c:func:`davidson_general_ext_rout`
       * :c:func:`davidson_general_ext_rout_diag_dressed`
       * :c:func:`davidson_general_ext_rout_dressed`
       * :c:func:`davidson_general_ext_rout_nonsym_b1space`
       * :c:func:`ortho_svd`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgeqrf`
       * :c:func:`dorgqr`

 
.. c:function:: ortho_qr_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_qr_complex(A,LDA,m,n)


    Orthogonalization using Q.R factorization
    
    A : matrix to orthogonalize
    
    LDA : leftmost dimension of A
    
    n : Number of rows of A
    
    m : Number of columns of A
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zgeqrf`
       * :c:func:`zungqr`

 
.. c:function:: ortho_qr_unblocked:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_qr_unblocked(A,LDA,m,n)


    Orthogonalization using Q.R factorization
    
    A : matrix to orthogonalize
    
    LDA : leftmost dimension of A
    
    n : Number of rows of A
    
    m : Number of columns of A
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgeqr2`
       * :c:func:`dorg2r`

 
.. c:function:: ortho_qr_unblocked_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_qr_unblocked_complex(A,LDA,m,n)


    Orthogonalization using Q.R factorization
    
    A : matrix to orthogonalize
    
    LDA : leftmost dimension of A
    
    n : Number of rows of A
    
    m : Number of columns of A
    

 
.. c:function:: ortho_svd:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_svd(A,LDA,m,n)


    Orthogonalization via fast SVD
    
    A : matrix to orthogonalize
    
    LDA : leftmost dimension of A
    
    m : Number of rows of A
    
    n : Number of columns of A
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`randomized_svd`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ortho_qr`
       * :c:func:`svd`

 
.. c:function:: overlap_cgaussian_x:


    File : :file:`utils/cgtos_one_e.irp.f`

    .. code:: fortran

        complex*16 function overlap_cgaussian_x(Ae_center, Be_center, alpha, beta, power_A, power_B, Ap_center, Bp_center, dim)


    
    \int_{-infty}^{+infty} (x - Ap_x)^ax (x - Bp_x)^bx exp(-alpha (x - Ae_x)^2) exp(-beta (x - Be_X)^2) dx
    
    with complex arguments
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_cpoly_and_cgaussian_x`

 
.. c:function:: overlap_gaussian_x:


    File : :file:`utils/one_e_integration.irp.f`

    .. code:: fortran

        double precision function overlap_gaussian_x(A_center,B_center,alpha,beta,power_A,power_B,dim)


    .. math::
    
     \sum_{-infty}^{+infty} (x-A_x)^ax (x-B_x)^bx exp(-alpha(x-A_x)^2) exp(-beta(x-B_X)^2) dx
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian_x`

 
.. c:function:: overlap_gaussian_xyz:


    File : :file:`utils/one_e_integration.irp.f`

    .. code:: fortran

        subroutine overlap_gaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, overlap_y, overlap_z, overlap, dim)


    .. math::
    
       S_x = \int (x-A_x)^{a_x} exp(-\alpha(x-A_x)^2)  (x-B_x)^{b_x} exp(-beta(x-B_x)^2) dx \\
       S = S_x S_y S_z
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_deriv_1_x`
       * :c:data:`ao_dipole_x`
       * :c:data:`ao_overlap`
       * :c:data:`ao_spread_x`
       * :c:data:`prim_normalization_factor`
       * :c:data:`shell_normalization_factor`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gaussian_product_x`
       * :c:func:`give_explicit_poly_and_gaussian`

 
.. c:function:: overlap_gaussian_xyz_v:


    File : :file:`utils/one_e_integration.irp.f`

    .. code:: fortran

        subroutine overlap_gaussian_xyz_v(A_center, B_center, alpha, beta, power_A, power_B, overlap, n_points)


    .. math::
    
       S_x = \int (x-A_x)^{a_x} exp(-\alpha(x-A_x)^2) (x-B_x)^{b_x} exp(-beta(x-B_x)^2) dx \\
       S = S_x S_y S_z
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian_v`

 
.. c:function:: overlap_x_abs:


    File : :file:`utils/one_e_integration.irp.f`

    .. code:: fortran

        subroutine overlap_x_abs(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, lower_exp_val, dx, nx)


    .. math                      ::
    
     \int_{-infty}^{+infty} (x-A_center)^(power_A) * (x-B_center)^power_B * exp(-alpha(x-A_center)^2) * exp(-beta(x-B_center)^2) dx
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap_abs`

 
.. c:function:: pivoted_cholesky:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine pivoted_cholesky( A, rank, tol, ndim, U)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`roothaan_hall_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dpstrf`

 
.. c:function:: pol_modif_center:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine pol_modif_center(A_center, B_center, iorder, A_pol, B_pol)


    
    Transform the pol centerd on A:
          [ \sum_i ax_i (x-x_A)^i ] [ \sum_j ay_j (y-y_A)^j ] [ \sum_k az_k (z-z_A)^k ]
    to a pol centered on B
          [ \sum_i bx_i (x-x_B)^i ] [ \sum_j by_j (y-y_B)^j ] [ \sum_k bz_k (z-z_B)^k ]
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pol_modif_center_x`

 
.. c:function:: pol_modif_center_x:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine pol_modif_center_x(A_center, B_center, iorder, A_pol, B_pol)


    
    Transform the pol centerd on A:
          [ \sum_i ax_i (x-x_A)^i ]
    to a pol centered on B
          [ \sum_i bx_i (x-x_B)^i ]
    
    bx_i = \sum_{j=i}^{iorder} ax_j (x_B - x_A)^(j-i) j! / [ i! (j-i)! ]
         = \sum_{j=i}^{iorder} ax_j (x_B - x_A)^(j-i) binom_func(j,i)
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pol_modif_center`

 
.. c:function:: primitive_value_explicit:


    File : :file:`utils/prim_in_r.irp.f`

    .. code:: fortran

        double precision function primitive_value_explicit(power_prim,center_prim,alpha,r)


    Evaluates at "r" a primitive of type :
    (x - center_prim(1))**power_prim(1)  (y - center_prim(2))**power_prim(2) * (z - center_prim(3))**power_prim(3)
    
    exp(-alpha * [(x - center_prim(1))**2 + (y - center_prim(2))**2 + (z - center_prim(3))**2] )

 
.. c:function:: print_memory_usage:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        subroutine print_memory_usage()


    Prints the memory usage in the output

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:data:`cholesky_ao_num`
       * :c:func:`write_time`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`resident_memory`
       * :c:func:`total_memory`

 
.. c:function:: randomized_svd:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine randomized_svd(A,LDA,U,LDU,D,Vt,LDVt,m,n,q,r)


    Randomized SVD: rank r, q power iterations
    
    1. Sample column space of A with P: Z = A.P where P is random with r+p columns.
    
    2. Power iterations : Z <- X . (Xt.Z)
    
    3. Z = Q.R
    
    4. Compute SVD on projected Qt.X = U' . S. Vt
    
    5. U = Q U'

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`ortho_svd`
       * :c:func:`random_number`
       * :c:func:`svd`

 
.. c:function:: recentered_cpoly2:


    File : :file:`utils/cgtos_utils.irp.f`

    .. code:: fortran

        subroutine recentered_cpoly2(P_A, x_A, x_P, a, P_B, x_B, x_Q, b)


    
    write two complex polynomials (x-x_A)^a (x-x_B)^b
    as P_A(x-x_P) and P_B(x-x_Q)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_cpoly_and_cgaussian`
       * :c:func:`give_explicit_cpoly_and_cgaussian_x`

 
.. c:function:: recentered_poly2:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine recentered_poly2(P_new, x_A, x_P, a, Q_new, x_B, x_Q, b)


    
    Recenter two polynomials:
    
      (x - x_A)^a -> \sum_{i=0}^{a} \binom{a}{i} (x_A - x_P)^{a-i} (x - x_P)^i
      (x - x_B)^b -> \sum_{i=0}^{b} \binom{b}{i} (x_B - x_Q)^{b-i} (x - x_Q)^i
    
    where:
      \binom{a}{i} = \binom{a}{a-i} = a! / [i! (a-i)!]
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian`
       * :c:func:`give_explicit_poly_and_gaussian_x`

 
.. c:function:: recentered_poly2_v:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine recentered_poly2_v(P_new, lda, x_A, LD_xA, x_P, a, P_new2, ldb, x_B, x_Q, b, n_points)


    Recenter two polynomials

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian_v`

 
.. c:function:: recentered_poly2_v0:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine recentered_poly2_v0(P_new, lda, x_A, LD_xA, x_P, a, n_points)


    
    Recenter two polynomials. Special case for b=(0,0,0)
    
    (x - A)^a (x - B)^0 = (x - P + P - A)^a  (x - Q + Q - B)^0
                        = (x - P + P - A)^a
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian_v`

 
.. c:function:: resident_memory:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        subroutine resident_memory(value)


    Returns the current used memory in gigabytes used by the current process.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha_chol`
       * :c:func:`check_mem`
       * :c:data:`cholesky_ao_num`
       * :c:func:`dav_double_dressed`
       * :c:func:`davidson_diag_csf_hjj`
       * :c:func:`davidson_diag_hjj`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`davidson_diag_nonsym_hjj`
       * :c:func:`davidson_general`
       * :c:func:`davidson_general_diag_dressed_ext_rout_nonsym_b1space`
       * :c:func:`davidson_general_ext_rout`
       * :c:func:`davidson_general_ext_rout_diag_dressed`
       * :c:func:`davidson_general_ext_rout_dressed`
       * :c:func:`davidson_general_ext_rout_nonsym_b1space`
       * :c:func:`print_memory_usage`
       * :c:func:`run_slave_main`
       * :c:func:`zmq_pt2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lock_io`
       * :c:func:`unlock_io`
       * :c:func:`usleep`

 
.. c:function:: restore_symmetry:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine restore_symmetry(m,n,A,LDA,thresh)


    Tries to find the matrix elements that are the same, and sets them
    to the average value.
    If restore_symm is False, only nullify small elements

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_to_mo`
       * :c:func:`create_guess`
       * :c:func:`huckel_guess`
       * :c:func:`roothaan_hall_scf`
       * :c:func:`svd_symm`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dset_order`
       * :c:func:`dsort`
       * :c:func:`iset_order`

 
.. c:function:: rint:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function rint(n,rho)


    .. math::
    
      \int_0^1 dx \exp(-p x^2) x^n
    

 
.. c:function:: rint1:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function rint1(n,rho)


    Standard version of rint

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fact_inv`
       * :c:data:`inv_int`

 
.. c:function:: rint_large_n:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function rint_large_n(n,rho)


    Version of rint for large values of n

 
.. c:function:: rint_sum:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function rint_sum(n_pt_out,rho,d1)


    Needed for the calculation of two-electron integrals.

 
.. c:function:: set_multiple_levels_omp:


    File : :file:`utils/set_multiple_levels_omp.irp.f`

    .. code:: fortran

        subroutine set_multiple_levels_omp(activate)


    If true, activate OpenMP nested parallelism. If false, deactivate.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map_cholesky`
       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_mo`
       * :c:func:`h_s2_u_0_nstates_zmq`
       * :c:func:`h_u_0_nstates_zmq`
       * :c:data:`mo_integrals_cache`
       * :c:func:`run_slave_cipsi`
       * :c:func:`run_slave_main`
       * :c:func:`zmq_pt2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_max_active_levels`
       * :c:func:`omp_set_nested`

 
.. c:function:: set_order:


    File : :file:`utils/sort.irp.f_template_90`

    .. code:: fortran

        subroutine set_order(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.

 
.. c:function:: set_order_big:


    File : :file:`utils/sort.irp.f_template_90`

    .. code:: fortran

        subroutine set_order_big(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: shank:


    File : :file:`utils/shank.irp.f`

    .. code:: fortran

        subroutine shank(array,n,nmax,shank1)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`shank_general`

 
.. c:function:: shank_function:


    File : :file:`utils/shank.irp.f`

    .. code:: fortran

        double precision function shank_function(array,i,n)



 
.. c:function:: shank_general:


    File : :file:`utils/shank.irp.f`

    .. code:: fortran

        double precision function shank_general(array,n,nmax)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`shank`

 
.. c:function:: sorted_dnumber:


    File : :file:`utils/sort.irp.f_template_32`

    .. code:: fortran

        subroutine sorted_dnumber(x,isize,n)


    Returns the number of sorted elements

 
.. c:function:: sorted_i2number:


    File : :file:`utils/sort.irp.f_template_32`

    .. code:: fortran

        subroutine sorted_i2number(x,isize,n)


    Returns the number of sorted elements

 
.. c:function:: sorted_i8number:


    File : :file:`utils/sort.irp.f_template_32`

    .. code:: fortran

        subroutine sorted_i8number(x,isize,n)


    Returns the number of sorted elements

 
.. c:function:: sorted_inumber:


    File : :file:`utils/sort.irp.f_template_32`

    .. code:: fortran

        subroutine sorted_inumber(x,isize,n)


    Returns the number of sorted elements

 
.. c:function:: sorted_number:


    File : :file:`utils/sort.irp.f_template_32`

    .. code:: fortran

        subroutine sorted_number(x,isize,n)


    Returns the number of sorted elements

 
.. c:function:: sub_a_at:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine sub_A_At(A, N)



 
.. c:function:: sum_a_at:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine sum_A_At(A, N)



 
.. c:function:: svd:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine svd(A,LDA,U,LDU,D,Vt,LDVt,m,n)


    Compute A = U.D.Vt
    
    LDx : leftmost dimension of x
    
    Dimension of A is m x n
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`eigsvd`
       * :c:func:`mo_as_svd_vectors_of_mo_matrix`
       * :c:func:`mo_as_svd_vectors_of_mo_matrix_eig`
       * :c:func:`mo_coef_new_as_svd_vectors_of_mo_matrix_eig`
       * :c:data:`natorbsci`
       * :c:func:`ortho_canonical`
       * :c:func:`ortho_lowdin`
       * :c:func:`ortho_svd`
       * :c:func:`randomized_svd`
       * :c:data:`s_half`
       * :c:data:`s_half_inv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgesvd`

 
.. c:function:: svd_complex:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine svd_complex(A,LDA,U,LDU,D,Vt,LDVt,m,n)


    Compute A = U.D.Vt
    
    LDx : leftmost dimension of x
    
    Dimension of A is m x n
    A,U,Vt are complex*16
    D is double precision

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ortho_canonical_complex`
       * :c:func:`ortho_lowdin_complex`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zgesvd`

 
.. c:function:: svd_symm:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine svd_symm(A,LDA,U,LDU,D,Vt,LDVt,m,n)


    Compute A = U.D.Vt
    
    LDx : leftmost dimension of x
    
    Dimension of A is m x n
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dgesvd`
       * :c:func:`restore_symmetry`

 
.. c:function:: total_memory:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        subroutine total_memory(value)


    Returns the current used memory in gigabytes used by the current process.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_memory_usage`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lock_io`
       * :c:func:`unlock_io`

 
.. c:function:: u_dot_u:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function u_dot_u(u,sze)


    Compute <u|u>

 
.. c:function:: u_dot_v:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function u_dot_v(u,v,sze)


    Compute <u|v>

 
.. c:function:: v2_over_x:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine v2_over_x(v,x,res)



 
.. c:function:: v_phi:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function V_phi(n, m)


    Computes the angular $\phi$ part of the nuclear attraction integral:
    
    $\int_{0}^{2 \pi} \cos(\phi)^n \sin(\phi)^m d\phi$.

 
.. c:function:: v_theta:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function V_theta(n, m)


    Computes the angular $\theta$ part of the nuclear attraction integral:
    
    $\int_{0}^{\pi} \cos(\theta)^n \sin(\theta)^m d\theta$

 
.. c:function:: wall_time:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine wall_time(t)


    The equivalent of cpu_time, but for the wall time.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`act_2_rdm_trans_spin_trace_mo`
       * :c:func:`add_integrals_to_map`
       * :c:func:`add_integrals_to_map_erf`
       * :c:data:`ao_pseudo_integrals_local`
       * :c:data:`ao_pseudo_integrals_non_local`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:func:`apply_mo_rotation`
       * :c:data:`bielecci`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:func:`diag_hessian_list_opt`
       * :c:func:`diag_hessian_opt`
       * :c:func:`diagonalization_hessian`
       * :c:func:`first_diag_hessian_list_opt`
       * :c:func:`first_diag_hessian_opt`
       * :c:func:`first_hessian_list_opt`
       * :c:func:`first_hessian_opt`
       * :c:func:`gradient_list_opt`
       * :c:func:`gradient_opt`
       * :c:func:`hessian_list_opt`
       * :c:func:`hessian_opt`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`output_wall_time_0`
       * :c:func:`pt2_collector`
       * :c:func:`rotation_matrix`
       * :c:func:`rotation_matrix_iterative`
       * :c:func:`run_pt2_slave_large`
       * :c:func:`run_pt2_slave_small`
       * :c:func:`run_slave_main`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`
       * :c:func:`trust_region_expected_e`
       * :c:func:`trust_region_optimal_lambda`
       * :c:func:`trust_region_step`
       * :c:func:`write_time`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`system_clock`

 
.. c:function:: wallis:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function Wallis(n)


    Wallis integral:
    
    $\int_{0}^{\pi} \cos(\theta)^n d\theta$.

 
.. c:function:: write_git_log:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine write_git_log(iunit)


    Write the last git commit in file iunit.

 
.. c:function:: zboysfun:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        subroutine zboysfun(n_max, x, vals)


    
    Computes values of the Boys function for n = 0, 1, ..., n_max
    for a complex valued argument
    
    Input: x --- argument, complex*16, Re(x) >= 0
    Output: vals  --- values of the Boys function, n = 0, 1, ..., n_max
    
    Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
    https://doi.org/10.1063/5.0062444
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`crint_2`
       * :c:func:`crint_2_vec`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zboysfun00_2`

 
.. c:function:: zboysfun00_1:


    File : :file:`utils/cpx_erf.irp.f`

    .. code:: fortran

        subroutine zboysfun00_1(rho, F0)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`crint_1`
       * :c:func:`crint_1_vec`

 
.. c:function:: zboysfun00_2:


    File : :file:`utils/cpx_erf.irp.f`

    .. code:: fortran

        subroutine zboysfun00_2(z, val)


    
    Computes values of the Boys function for n=0
    for a complex valued argument
    
    Input: z  --- argument, complex*16, Real(z) >= 0
    Output: val  --- value of the Boys function n=0
    
    Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
    https://doi.org/10.1063/5.0062444
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`zboysfun`

 
.. c:function:: zboysfun00nrp:


    File : :file:`utils/cpx_erf.irp.f`

    .. code:: fortran

        subroutine zboysfun00nrp(z, val)


    
    Computes values of the exp(z) F(0,z)
    (where F(0,z) is the Boys function)
    for a complex valued argument with Real(z)<=0
    
    Input: z  --- argument, complex*16, !!! Real(z)<=0 !!!
    Output: val  --- value of the function !!! exp(z) F(0,z) !!!, where F(0,z) is the Boys function
    
    Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
    https://doi.org/10.1063/5.0062444
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`zboysfunnrp`

 
.. c:function:: zboysfunnrp:


    File : :file:`utils/cpx_boys.irp.f`

    .. code:: fortran

        subroutine zboysfunnrp(n_max, x, vals)


    
    Computes values of e^z F(n,z) for n = 0, 1, ..., n_max
    (where F(n,z) are the Boys functions)
    for a complex valued argument WITH NEGATIVE REAL PART
    
    Input: x  --- argument, complex *16 Re(x)<=0
    Output: vals  --- values of e^z F(n,z), n = 0, 1, ..., n_max
    
    Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
    https://doi.org/10.1063/5.0062444
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`crint_2`
       * :c:func:`crint_2_vec`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`zboysfun00nrp`

