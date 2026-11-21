.. _module_non_h_ints_mu: 
 
.. program:: non_h_ints_mu 
 
.. default-role:: option 
 
non_h_ints_mu
=============

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Providers 
--------- 
 
.. c:var:: ao_non_hermit_term_chemist


    File : :file:`grad_tc_int.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_non_hermit_term_chemist	(ao_num,ao_num,ao_num,ao_num)


                               1 1 2 2      1 2                                1 2
    
    ao_non_hermit_term_chemist(k,i,l,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the AO basis
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_erf_rk_cst_mu`
       * :c:data:`x_v_ij_erf_rk_cst_mu`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_non_hermit_term_chemist`

 
.. c:var:: ao_tc_int_chemist


    File : :file:`total_tc_int.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_tc_int_chemist	(ao_num,ao_num,ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_tc_int_chemist_test`
       * :c:data:`ao_two_e_coul`
       * :c:data:`tc_grad_and_lapl_ao`
       * :c:data:`tc_grad_square_ao`
       * :c:data:`test_cycle_tc`


 
.. c:var:: ao_tc_int_chemist_test


    File : :file:`total_tc_int.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_tc_int_chemist_test	(ao_num,ao_num,ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_coul`
       * :c:data:`tc_grad_and_lapl_ao_test`
       * :c:data:`tc_grad_square_ao_test`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_tc_int_chemist`

 
.. c:var:: ao_two_e_coul


    File : :file:`total_tc_int.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_coul	(ao_num,ao_num,ao_num,ao_num)


    
    ao_two_e_coul(k,i,l,j) = ( k i | 1/r12 | l j ) = < l k | 1/r12 | j i >
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache`
       * :c:data:`ao_integrals_cache_min`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_tc_int_chemist`
       * :c:data:`ao_tc_int_chemist_test`
       * :c:data:`ao_vartc_int_chemist`

 
.. c:var:: ao_vartc_int_chemist


    File : :file:`total_tc_int.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_vartc_int_chemist	(ao_num,ao_num,ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_coul`
       * :c:data:`tc_grad_square_ao`
       * :c:data:`tc_grad_square_ao_test`
       * :c:data:`test_cycle_tc`


 
.. c:var:: grad12_j12


    File : :file:`grad_squared.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grad12_j12	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`expo_gauss_1_erf_x_2`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1u2_grad2u2_j1b2`
       * :c:data:`j1b_type`
       * :c:data:`n_max_fit_slat`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`v_1b`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_square_ao`
       * :c:data:`tc_grad_square_ao_loop`

 
.. c:var:: grad12_j12_test


    File : :file:`grad_squared_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grad12_j12_test	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`expo_gauss_1_erf_x_2`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1u2_grad2u2_j1b2_test`
       * :c:data:`j1b_type`
       * :c:data:`n_max_fit_slat`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`v_1b`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_square_ao_test`
       * :c:data:`tc_grad_square_ao_test_ref`

 
.. c:var:: gradu_squared_u_ij_mu


    File : :file:`grad_squared.irp.f`

    .. code:: fortran

        double precision, allocatable	:: gradu_squared_u_ij_mu	(ao_num,ao_num,n_points_final_grid)


    
    if J(r1,r2) = u12:
    
    gradu_squared_u_ij_mu = -0.50 x \int r2 [ (grad_1 u12)^2 + (grad_2 u12^2)] \phi_i(2) \phi_j(2)
                          = -0.25 x \int r2       (1 - erf(mu*r12))^2          \phi_i(2) \phi_j(2)
    and
      (1 - erf(mu*r12))^2 = \sum_i coef_gauss_1_erf_x_2(i) * exp(-expo_gauss_1_erf_x_2(i) * r12^2)
    
    if J(r1,r2) = u12 x v1 x v2
    
    gradu_squared_u_ij_mu = -0.50 x \int r2 \phi_i(2) \phi_j(2) [ v1^2  v2^2 ((grad_1 u12)^2 + (grad_2 u12^2)]) + u12^2  v2^2 (grad_1 v1)^2 + 2 u12 v1 v2^2 (grad_1 u12) . (grad_1 v1) ]
                          = -0.25 x        v1^2      \int r2 \phi_i(2) \phi_j(2) [1 - erf(mu r12)]^2 v2^2
                          + -0.50 x    (grad_1 v1)^2 \int r2 \phi_i(2) \phi_j(2)      u12^2          v2^2
                          + -1.00 x v1 (grad_1 v1)   \int r2 \phi_i(2) \phi_j(2)   (grad_1 u12)      v2^2
                          =                 v1^2        x   int2_grad1u2_grad2u2_j1b2
                          + -0.5 x      (grad_1 v1)^2   x   int2_u2_j1b2
                          + -1.0 X V1 x (grad_1 v1)   \cdot [ int2_u_grad1u_j1b2 x r - int2_u_grad1u_x_j1b ]
    
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`expo_gauss_1_erf_x_2`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1u2_grad2u2_j1b2`
       * :c:data:`int2_u2_j1b2`
       * :c:data:`int2_u_grad1u_j1b2`
       * :c:data:`int2_u_grad1u_x_j1b2`
       * :c:data:`j1b_type`
       * :c:data:`n_max_fit_slat`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`v_1b`
       * :c:data:`v_1b_grad`


 
.. c:var:: int1_grad2_u12_ao


    File : :file:`new_grad_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int1_grad2_u12_ao	(3,ao_num,ao_num,n_points_final_grid)


    
    int1_grad2_u12_ao(:,i,j,ipoint) = \int dr1 [-1 * \grad_r2 J(r1,r2)] \phi_i(r1) \phi_j(r1)
    
    where r1 = r(ipoint)
    
    if J(r1,r2) = u12:
    
    int1_grad2_u12_ao(:,i,j,ipoint) = +0.5 x \int dr1 [-(r1 - r2) (erf(mu * r12)-1)r_12] \phi_i(r1) \phi_j(r1)
                                    = -0.5 * [ v_ij_erf_rk_cst_mu(i,j,ipoint) * r(:) - x_v_ij_erf_rk_cst_mu(i,j,ipoint,:) ]
                                    = -int2_grad1_u12_ao(i,j,ipoint,:)
    
    if J(r1,r2) = u12 x v1 x v2
    
    int1_grad2_u12_ao(:,i,j,ipoint) =      v2    x [ 0.5 x \int dr1 [-(r1 - r2) (erf(mu * r12)-1)r_12] v1 \phi_i(r1) \phi_j(r1) ]
                                    - \grad_2 v2 x [       \int dr1                   u12              v1 \phi_i(r1) \phi_j(r1) ]
                                    =   -0.5 v_1b(ipoint) * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint) * r(:)
                                    +    0.5 v_1b(ipoint) * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,:)
                                    - v_1b_grad[:,ipoint] * v_ij_u_cst_mu_j1b(i,j,ipoint)
    
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`j1b_type`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_1b`
       * :c:data:`v_1b_grad`
       * :c:data:`v_ij_erf_rk_cst_mu_j1b`
       * :c:data:`v_ij_u_cst_mu_j1b`
       * :c:data:`x_v_ij_erf_rk_cst_mu_j1b`


 
.. c:var:: int2_grad1_u12_ao


    File : :file:`new_grad_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao	(ao_num,ao_num,n_points_final_grid,3)


    
    int2_grad1_u12_ao(i,j,ipoint,:) = \int dr2 [-1 * \grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2)
    
    where r1 = r(ipoint)
    
    if J(r1,r2) = u12:
    
    int2_grad1_u12_ao(i,j,ipoint,:) = 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] \phi_i(r2) \phi_j(r2)
                                    = 0.5 * [ v_ij_erf_rk_cst_mu(i,j,ipoint) * r(:) - x_v_ij_erf_rk_cst_mu(i,j,ipoint,:) ]
    
    if J(r1,r2) = u12 x v1 x v2
    
    int2_grad1_u12_ao(i,j,ipoint,:) =      v1    x [ 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] v2 \phi_i(r2) \phi_j(r2) ]
                                    - \grad_1 v1 x [       \int dr2                  u12              v2 \phi_i(r2) \phi_j(r2) ]
                                    =    0.5 v_1b(ipoint) * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint) * r(:)
                                    -    0.5 v_1b(ipoint) * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,:)
                                    - v_1b_grad[:,ipoint] * v_ij_u_cst_mu_j1b(i,j,ipoint)
    
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`j1b_type`
       * :c:data:`n_points_final_grid`
       * :c:data:`read_tc_integ`
       * :c:data:`v_1b`
       * :c:data:`v_1b_grad`
       * :c:data:`v_ij_erf_rk_cst_mu`
       * :c:data:`v_ij_erf_rk_cst_mu_j1b`
       * :c:data:`v_ij_u_cst_mu_j1b`
       * :c:data:`write_tc_integ`
       * :c:data:`x_v_ij_erf_rk_cst_mu_j1b`
       * :c:data:`x_v_ij_erf_rk_cst_mu_transp_bis`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int1_grad2_u12_ao`
       * :c:data:`tc_grad_and_lapl_ao`
       * :c:data:`tc_grad_and_lapl_ao_loop`

 
.. c:var:: int2_grad1_u12_ao_test


    File : :file:`new_grad_tc_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao_test	(ao_num,ao_num,n_points_final_grid,3)


    
    int2_grad1_u12_ao_test(i,j,ipoint,:) = \int dr2 [-1 * \grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2)
    
    where r1 = r(ipoint)
    
    if J(r1,r2) = u12:
    
    int2_grad1_u12_ao_test(i,j,ipoint,:) = 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] \phi_i(r2) \phi_j(r2)
                                    = 0.5 * [ v_ij_erf_rk_cst_mu(i,j,ipoint) * r(:) - x_v_ij_erf_rk_cst_mu(i,j,ipoint,:) ]
    
    if J(r1,r2) = u12 x v1 x v2
    
    int2_grad1_u12_ao_test(i,j,ipoint,:) =      v1    x [ 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] v2 \phi_i(r2) \phi_j(r2) ]
                                    - \grad_1 v1 x [       \int dr2                  u12              v2 \phi_i(r2) \phi_j(r2) ]
                                    =    0.5 v_1b(ipoint) * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint) * r(:)
                                    -    0.5 v_1b(ipoint) * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,:)
                                    - v_1b_grad[:,ipoint] * v_ij_u_cst_mu_j1b(i,j,ipoint)
    
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`j1b_type`
       * :c:data:`n_points_final_grid`
       * :c:data:`read_tc_integ`
       * :c:data:`v_1b`
       * :c:data:`v_1b_grad`
       * :c:data:`v_ij_erf_rk_cst_mu`
       * :c:data:`v_ij_erf_rk_cst_mu_j1b_test`
       * :c:data:`v_ij_u_cst_mu_j1b_test`
       * :c:data:`write_tc_integ`
       * :c:data:`x_v_ij_erf_rk_cst_mu_j1b_test`
       * :c:data:`x_v_ij_erf_rk_cst_mu_tmp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_and_lapl_ao_test`

 
.. c:var:: mo_non_hermit_term


    File : :file:`grad_tc_int.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_non_hermit_term	(mo_num,mo_num,mo_num,mo_num)


                       1 2 1 2      1 2                                1 2
    
    mo_non_hermit_term(k,l,i,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_non_hermit_term_chemist`
       * :c:data:`mo_num`


 
.. c:var:: mo_non_hermit_term_chemist


    File : :file:`grad_tc_int.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_non_hermit_term_chemist	(mo_num,mo_num,mo_num,mo_num)


                               1 1 2 2      1 2                                1 2
    
    mo_non_hermit_term_chemist(k,i,l,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_non_hermit_term_chemist`
       * :c:data:`ao_num`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_non_hermit_term`

 
.. c:var:: tc_grad_and_lapl_ao


    File : :file:`new_grad_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_grad_and_lapl_ao	(ao_num,ao_num,ao_num,ao_num)


    
    tc_grad_and_lapl_ao(k,i,l,j) = < k l | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) . \grad_1 | ij >
    
    = 1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2)
    
    This is obtained by integration by parts.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`n_points_final_grid`
       * :c:data:`read_tc_integ`
       * :c:data:`write_tc_integ`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_tc_int_chemist`

 
.. c:var:: tc_grad_and_lapl_ao_loop


    File : :file:`new_grad_tc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_grad_and_lapl_ao_loop	(ao_num,ao_num,ao_num,ao_num)


    
    tc_grad_and_lapl_ao_loop(k,i,l,j) = < k l | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) . \grad_1 | ij >
    
    = 1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2)
    
    This is obtained by integration by parts.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`n_points_final_grid`


 
.. c:var:: tc_grad_and_lapl_ao_test


    File : :file:`new_grad_tc_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_grad_and_lapl_ao_test	(ao_num,ao_num,ao_num,ao_num)


    
    tc_grad_and_lapl_ao_test(k,i,l,j) = < k l | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) | ij >
    
    = 1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2)
    
    This is obtained by integration by parts.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao_test`
       * :c:data:`n_points_final_grid`
       * :c:data:`read_tc_integ`
       * :c:data:`write_tc_integ`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_tc_int_chemist_test`

 
.. c:var:: tc_grad_square_ao


    File : :file:`grad_squared.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_grad_square_ao	(ao_num,ao_num,ao_num,ao_num)


    
    tc_grad_square_ao(k,i,l,j) = 1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_2 u(r1,r2)|^2 | ij>
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`grad12_j12`
       * :c:data:`n_points_final_grid`
       * :c:data:`read_tc_integ`
       * :c:data:`u12_grad1_u12_j1b_grad1_j1b`
       * :c:data:`u12sq_j1bsq`
       * :c:data:`write_tc_integ`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_tc_int_chemist`
       * :c:data:`ao_vartc_int_chemist`

 
.. c:var:: tc_grad_square_ao_loop


    File : :file:`grad_squared.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_grad_square_ao_loop	(ao_num,ao_num,ao_num,ao_num)


    
    tc_grad_square_ao_loop(k,i,l,j) = 1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_2 u(r1,r2)|^2 | ij>
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`grad12_j12`
       * :c:data:`n_points_final_grid`
       * :c:data:`u12_grad1_u12_j1b_grad1_j1b`
       * :c:data:`u12sq_j1bsq`


 
.. c:var:: tc_grad_square_ao_test


    File : :file:`grad_squared_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_grad_square_ao_test	(ao_num,ao_num,ao_num,ao_num)


    
    tc_grad_square_ao_test(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_1 u(r1,r2)|^2 | ij>
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`grad12_j12_test`
       * :c:data:`n_points_final_grid`
       * :c:data:`read_tc_integ`
       * :c:data:`u12_grad1_u12_j1b_grad1_j1b_test`
       * :c:data:`u12sq_j1bsq_test`
       * :c:data:`write_tc_integ`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_tc_int_chemist_test`
       * :c:data:`ao_vartc_int_chemist`

 
.. c:var:: tc_grad_square_ao_test_ref


    File : :file:`grad_squared_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_grad_square_ao_test_ref	(ao_num,ao_num,ao_num,ao_num)


    
    tc_grad_square_ao_test_ref(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_1 u(r1,r2)|^2 | ij>
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`grad12_j12_test`
       * :c:data:`n_points_final_grid`
       * :c:data:`u12_grad1_u12_j1b_grad1_j1b_test`
       * :c:data:`u12sq_j1bsq_test`


 
.. c:var:: u12_grad1_u12_j1b_grad1_j1b


    File : :file:`grad_squared.irp.f`

    .. code:: fortran

        double precision, allocatable	:: u12_grad1_u12_j1b_grad1_j1b	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u_grad1u_j1b2`
       * :c:data:`int2_u_grad1u_x_j1b2`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_1b`
       * :c:data:`v_1b_grad`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_square_ao`
       * :c:data:`tc_grad_square_ao_loop`

 
.. c:var:: u12_grad1_u12_j1b_grad1_j1b_test


    File : :file:`grad_squared_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: u12_grad1_u12_j1b_grad1_j1b_test	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u_grad1u_j1b2_test`
       * :c:data:`int2_u_grad1u_x_j1b2_test`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_1b`
       * :c:data:`v_1b_grad`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_square_ao_test`
       * :c:data:`tc_grad_square_ao_test_ref`

 
.. c:var:: u12sq_j1bsq


    File : :file:`grad_squared.irp.f`

    .. code:: fortran

        double precision, allocatable	:: u12sq_j1bsq	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_u2_j1b2`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_1b_grad`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_square_ao`
       * :c:data:`tc_grad_square_ao_loop`

 
.. c:var:: u12sq_j1bsq_test


    File : :file:`grad_squared_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: u12sq_j1bsq_test	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_u2_j1b2_test`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_1b_grad`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_square_ao_test`
       * :c:data:`tc_grad_square_ao_test_ref`

 
.. c:var:: v_1b


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: v_1b	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`j1b_pen`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`grad12_j12`
       * :c:data:`grad12_j12_test`
       * :c:data:`gradu_squared_u_ij_mu`
       * :c:data:`int1_grad2_u12_ao`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_ao_test`
       * :c:data:`u12_grad1_u12_j1b_grad1_j1b`
       * :c:data:`u12_grad1_u12_j1b_grad1_j1b_test`

 
.. c:var:: v_1b_grad


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: v_1b_grad	(3,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`j1b_pen`
       * :c:data:`list_all_comb_b2`
       * :c:data:`list_all_comb_b2_size`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradu_squared_u_ij_mu`
       * :c:data:`int1_grad2_u12_ao`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_ao_test`
       * :c:data:`u12_grad1_u12_j1b_grad1_j1b`
       * :c:data:`u12_grad1_u12_j1b_grad1_j1b_test`
       * :c:data:`u12sq_j1bsq`
       * :c:data:`u12sq_j1bsq_test`

 
.. c:var:: v_1b_lapl


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: v_1b_lapl	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`j1b_pen`
       * :c:data:`list_all_comb_b2`
       * :c:data:`list_all_comb_b2_size`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`


 
.. c:var:: v_1b_list_b2


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: v_1b_list_b2	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`list_all_comb_b2_coef`
       * :c:data:`list_all_comb_b2_coef`
       * :c:data:`list_all_comb_b2_coef`
       * :c:data:`list_all_comb_b2_size`
       * :c:data:`n_points_final_grid`


 
.. c:var:: v_1b_list_b3


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: v_1b_list_b3	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`list_all_comb_b3_coef`
       * :c:data:`list_all_comb_b3_coef`
       * :c:data:`list_all_comb_b3_coef`
       * :c:data:`list_all_comb_b3_size`
       * :c:data:`n_points_final_grid`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: debug_fit:


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine debug_fit



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`j1b_pen`
       * :c:data:`mu_erf`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`test_fit_u`
       * :c:func:`test_grad_j1b_nucl`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: debug_integ_jmu_modif:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine debug_integ_jmu_modif



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`j1b_pen`
       * :c:data:`mu_erf`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`test_vect_overlap_gauss_r12_ao_with1s`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: grad1_j12_mu_exc:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        subroutine grad1_j12_mu_exc(r1, r2, grad)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_jmu_modif_num`
       * :c:func:`num_grad12_j12`
       * :c:func:`num_gradu_squared_u_ij_mu`
       * :c:func:`num_u12_grad1_u12_j1b_grad1_j1b`
       * :c:func:`num_u12sq_j1bsq`
       * :c:func:`test_fit_ugradu`

 
.. c:function:: grad1_jmu_modif_num:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        subroutine grad1_jmu_modif_num(r1, r2, grad)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`num_int2_grad1_u12_ao`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_mu_exc`

 
.. c:function:: grad1_x_j12_mu_num:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_x_j12_mu_num(r1, r2)



 
.. c:function:: grad1_x_jmu_modif:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_x_jmu_modif(r1, r2)



 
.. c:function:: grad1_y_j12_mu_num:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_y_j12_mu_num(r1, r2)



 
.. c:function:: grad1_y_jmu_modif:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_y_jmu_modif(r1, r2)



 
.. c:function:: grad1_z_j12_mu_num:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_z_j12_mu_num(r1, r2)



 
.. c:function:: grad1_z_jmu_modif:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_z_jmu_modif(r1, r2)



 
.. c:function:: grad_x_j1b_nucl:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad_x_j1b_nucl(r)



 
.. c:function:: grad_y_j1b_nucl:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad_y_j1b_nucl(r)



 
.. c:function:: grad_z_j1b_nucl:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad_z_j1b_nucl(r)



 
.. c:function:: j12_mu:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function j12_mu(r1, r2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

 
.. c:function:: j12_mu_gauss:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function j12_mu_gauss(r1, r2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`expo_gauss_j_mu_x`
       * :c:data:`n_max_fit_slat`

 
.. c:function:: j12_mu_r12:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function j12_mu_r12(r12)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

 
.. c:function:: j12_nucl:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function j12_nucl(r1, r2)



 
.. c:function:: j1b_nucl:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function j1b_nucl(r)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`j1b_pen`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

 
.. c:function:: jmu_modif:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function jmu_modif(r1, r2)



 
.. c:function:: lapl_j1b_nucl:


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function lapl_j1b_nucl(r)



 
.. c:function:: num_grad12_j12:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_grad12_j12(i, j, ipoint)


    
    -0.50 x \int r2 \phi_i(2) \phi_j(2) x v2^2 [v1^2 ((grad_1 u12)^2 + (grad_2 u12^2)]) ]
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_mu_exc`

 
.. c:function:: num_gradu_squared_u_ij_mu:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_gradu_squared_u_ij_mu(i, j, ipoint)


    
    -0.50 x \int r2 \phi_i(2) \phi_j(2) x v2^2
         [  v1^2    ((grad_1 u12)^2 + (grad_2 u12^2)])
         + u12^2         (grad_1 v1)^2
         + 2 u12 v1 (grad_1 u12) . (grad_1 v1)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_mu_exc`

 
.. c:function:: num_int2_grad1_u12_ao:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        subroutine num_int2_grad1_u12_ao(i, j, ipoint, integ)


    
    \int dr2 [-grad_1 u12] \phi_i(r2) \phi_j(r2) x v12_1b(r1, r2)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_int2_grad1_u12_ao`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_jmu_modif_num`

 
.. c:function:: num_int2_grad1u2_grad2u2_j1b2:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_int2_grad1u2_grad2u2_j1b2(i, j, ipoint)


    
    \int dr2 \frac{-[erf(mu r12) -1]^2}{4} \phi_i(r2) \phi_j(r2) x v_1b(r2)^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

 
.. c:function:: num_int2_u2_j1b2:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_int2_u2_j1b2(i, j, ipoint)


    
    \int dr2 u12^2 \phi_i(r2) \phi_j(r2) x v_1b(r2)^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

 
.. c:function:: num_int2_u_grad1u_total_j1b2:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        subroutine num_int2_u_grad1u_total_j1b2(i, j, ipoint, integ)


    
    \int dr2 u12 (grad_1 u12) \phi_i(r2) \phi_j(r2) x v_1b(r2)^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_int2_u_grad1u_total_j1b2`

 
.. c:function:: num_u12_grad1_u12_j1b_grad1_j1b:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_u12_grad1_u12_j1b_grad1_j1b(i, j, ipoint)


    
    -0.50 x \int r2 \phi_i(2) \phi_j(2) x v2^2 [ 2 u12 v1 (grad_1 u12) . (grad_1 v1) ]
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_mu_exc`

 
.. c:function:: num_u12sq_j1bsq:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_u12sq_j1bsq(i, j, ipoint)


    
    -0.50 x \int r2 \phi_i(2) \phi_j(2) x v2^2 [ u12^2 (grad_1 v1)^2 ]
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_mu_exc`

 
.. c:function:: num_v_ij_erf_rk_cst_mu_j1b:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_v_ij_erf_rk_cst_mu_j1b(i, j, ipoint)


    
    \int dr2 [erf(mu r12) -1]/r12 \phi_i(r2) \phi_j(r2) x v_1b(r2)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

 
.. c:function:: num_v_ij_u_cst_mu_j1b:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_v_ij_u_cst_mu_j1b(i, j, ipoint)


    
    \int dr2 u12 \phi_i(r2) \phi_j(r2) x v_1b(r2)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

 
.. c:function:: num_x_v_ij_erf_rk_cst_mu_j1b:


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        subroutine num_x_v_ij_erf_rk_cst_mu_j1b(i, j, ipoint, integ)


    
    \int dr2 [erf(mu r12) -1]/r12 \phi_i(r2) \phi_j(r2) x v_1b(r2) x r2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_x_v_ij_erf_rk_cst_mu_j1b`

 
.. c:function:: routine_fit:


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine routine_fit



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_non_h`

 
.. c:function:: routine_grad_squared:


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine routine_grad_squared



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_tc_sym_two_e_pot_map`
       * :c:data:`tc_grad_square_ao`

 
.. c:function:: routine_lapl_grad:


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine routine_lapl_grad



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_non_hermit_term_chemist`
       * :c:data:`ao_num`
       * :c:data:`ao_tc_sym_two_e_pot_map`
       * :c:data:`tc_grad_and_lapl_ao`

 
.. c:function:: test_fit_u:


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_fit_u()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`expo_gauss_j_mu_x`
       * :c:data:`final_grid_points`
       * :c:data:`n_max_fit_slat`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_fit`

 
.. c:function:: test_fit_u2:


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_fit_u2()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`expo_gauss_j_mu_x_2`
       * :c:data:`final_grid_points`
       * :c:data:`n_max_fit_slat`
       * :c:data:`n_points_final_grid`

 
.. c:function:: test_fit_ugradu:


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_fit_ugradu()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`expo_gauss_j_mu_1_erf`
       * :c:data:`final_grid_points`
       * :c:data:`n_max_fit_slat`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_mu_exc`

 
.. c:function:: test_grad12_j12:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_grad12_j12()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`grad12_j12`
       * :c:data:`n_points_final_grid`

 
.. c:function:: test_grad_j1b_nucl:


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_grad_j1b_nucl()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_1b_grad`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_fit`

 
.. c:function:: test_gradu_squared_u_ij_mu:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_gradu_squared_u_ij_mu()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`gradu_squared_u_ij_mu`
       * :c:data:`n_points_final_grid`

 
.. c:function:: test_int2_grad1_u12_ao:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_int2_grad1_u12_ao()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`num_int2_grad1_u12_ao`

 
.. c:function:: test_int2_grad1u2_grad2u2_j1b2:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_int2_grad1u2_grad2u2_j1b2()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_grad1u2_grad2u2_j1b2`
       * :c:data:`n_points_final_grid`

 
.. c:function:: test_int2_u2_j1b2:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_int2_u2_j1b2()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_u2_j1b2`
       * :c:data:`n_points_final_grid`

 
.. c:function:: test_int2_u_grad1u_total_j1b2:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_int2_u_grad1u_total_j1b2()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u_grad1u_j1b2`
       * :c:data:`int2_u_grad1u_x_j1b2`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`num_int2_u_grad1u_total_j1b2`

 
.. c:function:: test_j1b_nucl:


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_j1b_nucl()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_1b`

 
.. c:function:: test_lapl_j1b_nucl:


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_lapl_j1b_nucl()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_1b_lapl`

 
.. c:function:: test_list_b2:


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_list_b2()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_1b_list_b2`

 
.. c:function:: test_list_b3:


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_list_b3()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_1b_list_b3`

 
.. c:function:: test_non_h:


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_non_h



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`routine_fit`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: test_u12_grad1_u12_j1b_grad1_j1b:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_u12_grad1_u12_j1b_grad1_j1b()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`u12_grad1_u12_j1b_grad1_j1b`

 
.. c:function:: test_u12sq_j1bsq:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_u12sq_j1bsq()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`u12sq_j1bsq`

 
.. c:function:: test_v_ij_erf_rk_cst_mu_j1b:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_v_ij_erf_rk_cst_mu_j1b()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_erf_rk_cst_mu_j1b`

 
.. c:function:: test_v_ij_u_cst_mu_j1b:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_v_ij_u_cst_mu_j1b()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_u_cst_mu_j1b`

 
.. c:function:: test_vect_overlap_gauss_r12_ao:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_vect_overlap_gauss_r12_ao()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`expo_gauss_j_mu_x_2`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_transp`
       * :c:data:`j1b_pen`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_gauss_r12_ao_v`

 
.. c:function:: test_vect_overlap_gauss_r12_ao_with1s:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_vect_overlap_gauss_r12_ao_with1s()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`expo_gauss_j_mu_x_2`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_transp`
       * :c:data:`j1b_pen`
       * :c:data:`list_all_comb_b3_coef`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_integ_jmu_modif`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_gauss_r12_ao_with1s_v`

 
.. c:function:: test_x_v_ij_erf_rk_cst_mu_j1b:


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_x_v_ij_erf_rk_cst_mu_j1b()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_v_ij_erf_rk_cst_mu_j1b`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`num_x_v_ij_erf_rk_cst_mu_j1b`

