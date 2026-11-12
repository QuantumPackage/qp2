.. _module_non_h_ints_mu: 
 
.. program:: non_h_ints_mu 
 
.. default-role:: option 
 
=============
non_h_ints_mu
=============

The non-Hermitian part of the mu-transcorrelated Hamiltonian.

This module computes the non-Hermitian potential of the mu-TC
Hamiltonian on the |AO| and the BI-ORTHO |MO| basis. The operator is
defined in Eq. 33 of JCP 154, 084119 (2021).

The two providers are:

* `ao_non_hermit_term_chemist`: the non-Hermitian part of the
  two-electron TC Hamiltonian in the |AO| basis
* `mo_non_hermit_term_chemist`: the same operator in the BI-ORTHO
  |MO| basis

The integrals are built from the Jastrow derivative routines
(`jastrow_psi.irp.f`, `jast_deriv*.irp.f`), the mean field
:math:`\mu(r)` and its gradient (`mu_of_r.irp.f`), and the one- and
two-body J-mu integral kernels, in the bump (`j_bump.irp.f`) and
Gaussian (`deb_j_gauss.irp.f`) forms. The `deb_*.irp.f` files are
debug programs for the individual pieces.
 
 
 
Programs 
-------- 
 
 * :ref:`deb_aos` 
 * :ref:`deb_deriv_mu` 
 * :ref:`deb_j_bump` 
 * :ref:`deb_j_gauss` 
 * :ref:`deb_j_mu_of_r` 
 * :ref:`deb_j_psi` 
 * :ref:`deb_mos` 
 * :ref:`debug_fit` 
 * :ref:`debug_integ_jmu_modif` 
 * :ref:`plot_j_gauss` 
 * :ref:`plot_mo` 
 * :ref:`plot_mu_of_r` 
 * :ref:`print_j1ecoef_info` 
 * :ref:`print_jastrow_psi` 
 * :ref:`test_mu_of_r_tc` 
 * :ref:`test_non_h_ints` 
 
Providers 
--------- 
 
.. c:var:: double ao_non_hermit_term_chemist


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

 
.. c:var:: double ao_two_e_tc_tot


    File : :file:`total_tc_int.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_tc_tot	(ao_num,ao_num,ao_num,ao_num)


    
    CHEMIST NOTATION IS USED
    
    ao_two_e_tc_tot(k,i,l,j) = (ki|V^TC(r_12)|lj)
                             = <lk| V^TC(r_12) |ji> where V^TC(r_12) is the total TC operator
                             = tc_grad_and_lapl_ao(k,i,l,j) + tc_grad_square_ao(k,i,l,j) + ao_two_e_coul(k,i,l,j)
    AND IF(var_tc):
    
    ao_two_e_tot(k,i,l,j) = (ki|V^TC(r_12) + [(V^TC)(r_12)]^\dagger|lj) / 2.0
                          = tc_grad_square_ao(k,i,l,j) + ao_two_e_coul(k,i,l,j)
    
    
    where:
    
    tc_grad_and_lapl_ao(k,i,l,j) = < k l | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) . \grad_1 | ij >
                                 = -1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2      \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2)
                                 =  1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 (-1) \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2)
    
    tc_grad_square_ao(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_2 u(r1,r2)|^2 | ij>
    
    ao_two_e_coul(k,i,l,j) = < l k | 1/r12 | j i > = ( k i | 1/r12 | l j )
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_ao_transp`
       * :c:data:`env_square_grad`
       * :c:data:`env_type`
       * :c:data:`ezfio_filename`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_square_ao`
       * :c:data:`int2_u2_env2`
       * :c:data:`is_periodic`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`mpi_master`
       * :c:data:`n_points_final_grid`
       * :c:data:`read_ao_two_e_integrals`
       * :c:data:`read_tc_integ`
       * :c:data:`tc_integ_type`
       * :c:data:`tc_save_mem`
       * :c:data:`tc_save_mem_loops`
       * :c:data:`use_cgtos`
       * :c:data:`use_ipp`
       * :c:data:`use_libcint_effective`
       * :c:data:`use_only_lr`
       * :c:data:`var_tc`


 
.. c:var:: double aos_grad_in_r_array_qmckl


    File : :file:`qmckl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_in_r_array_qmckl	(ao_num,n_points_final_grid)
        double precision, allocatable	:: aos_grad_in_r_array_qmckl	(ao_num,n_points_final_grid,3)
        double precision, allocatable	:: aos_lapl_in_r_array_qmckl	(ao_num,n_points_final_grid)


    AOS computed with qmckl

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`trexio_file`
       * :c:data:`trexio_filename`


 
.. c:var:: double aos_in_r_array_qmckl


    File : :file:`qmckl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_in_r_array_qmckl	(ao_num,n_points_final_grid)
        double precision, allocatable	:: aos_grad_in_r_array_qmckl	(ao_num,n_points_final_grid,3)
        double precision, allocatable	:: aos_lapl_in_r_array_qmckl	(ao_num,n_points_final_grid)


    AOS computed with qmckl

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`trexio_file`
       * :c:data:`trexio_filename`


 
.. c:var:: double aos_lapl_in_r_array_qmckl


    File : :file:`qmckl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_in_r_array_qmckl	(ao_num,n_points_final_grid)
        double precision, allocatable	:: aos_grad_in_r_array_qmckl	(ao_num,n_points_final_grid,3)
        double precision, allocatable	:: aos_lapl_in_r_array_qmckl	(ao_num,n_points_final_grid)


    AOS computed with qmckl

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`trexio_file`
       * :c:data:`trexio_filename`


 
.. c:var:: double c_ij_ab_jastrow


    File : :file:`jastrow_psi.irp.f`

    .. code:: fortran

        double precision, allocatable	:: c_ij_ab_jastrow	(mo_num,mo_num,elec_alpha_num,elec_beta_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`ezfio_work_dir`
       * :c:data:`mo_num`


 
.. c:var:: double env_grad


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: env_grad	(3,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_expo`
       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`list_env1s`
       * :c:data:`list_env1s_size`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_ao_test`
       * :c:data:`int2_grad1_u12_square_ao`
       * :c:data:`int2_grad1_u2e_ao`
       * :c:data:`u12_grad1_u12_env_grad1_env`
       * :c:data:`u12_grad1_u12_env_grad1_env_test`
       * :c:data:`u12sq_envsq`
       * :c:data:`u12sq_envsq_test`

 
.. c:var:: double env_square_grad


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: env_square_grad	(n_points_final_grid,3)
        double precision, allocatable	:: env_square_lapl	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`list_env1s_square_coef`
       * :c:data:`list_env1s_square_coef`
       * :c:data:`list_env1s_square_coef`
       * :c:data:`list_env1s_square_size`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`tc_grad_square_ao_new`

 
.. c:var:: double env_square_lapl


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: env_square_grad	(n_points_final_grid,3)
        double precision, allocatable	:: env_square_lapl	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`list_env1s_square_coef`
       * :c:data:`list_env1s_square_coef`
       * :c:data:`list_env1s_square_coef`
       * :c:data:`list_env1s_square_size`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`tc_grad_square_ao_new`

 
.. c:var:: double env_val


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: env_val	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_expo`
       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`grad12_j12`
       * :c:data:`grad12_j12_test`
       * :c:data:`int2_grad1_u12_ao_test`
       * :c:data:`int2_grad1_u12_square_ao`
       * :c:data:`int2_grad1_u2e_ao`
       * :c:data:`int2_u2e_ao`
       * :c:data:`u12_grad1_u12_env_grad1_env`
       * :c:data:`u12_grad1_u12_env_grad1_env_test`

 
.. c:var:: double grad12_j12


    File : :file:`grad_squared.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grad12_j12	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`env_val`
       * :c:data:`int2_grad1u2_grad2u2_env2`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_square_ao`

 
.. c:var:: double grad12_j12_test


    File : :file:`grad_squared_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grad12_j12_test	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`env_type`
       * :c:data:`env_val`
       * :c:data:`int2_grad1u2_grad2u2_env2_test`
       * :c:data:`j2e_type`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_square_ao_test`
       * :c:data:`tc_grad_square_ao_test_ref`

 
.. c:var:: double grad1_u12_num


    File : :file:`jast_deriv.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grad1_u12_num	(n_points_extra_final_grid,n_points_final_grid,3)
        double precision, allocatable	:: grad1_u12_squared_num	(n_points_extra_final_grid,n_points_final_grid)


    
    
    grad_1 u(r1,r2)
    numerical integration over r1 & r2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`env_expo`
       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`j2e_type`
       * :c:data:`list_env1s`
       * :c:data:`list_env1s_size`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`qmckl_ctx_jastrow`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_ao_num_1shot`

 
.. c:var:: double grad1_u12_squared_num


    File : :file:`jast_deriv.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grad1_u12_num	(n_points_extra_final_grid,n_points_final_grid,3)
        double precision, allocatable	:: grad1_u12_squared_num	(n_points_extra_final_grid,n_points_final_grid)


    
    
    grad_1 u(r1,r2)
    numerical integration over r1 & r2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`env_expo`
       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`j2e_type`
       * :c:data:`list_env1s`
       * :c:data:`list_env1s_size`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`qmckl_ctx_jastrow`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_ao_num_1shot`

 
.. c:var:: double int2_grad1_u12_ao


    File : :file:`tc_integ.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao	(ao_num,ao_num,n_points_final_grid,3)


    
    int2_grad1_u12_ao(i,j,ipoint,:) = \int dr2 [\grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2)
    
    where r1 = r(ipoint)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`elec_num`
       * :c:data:`env_type`
       * :c:data:`ezfio_filename`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`final_grid_points_extra`
       * :c:data:`int2_grad1_u12_ao_num_1shot`
       * :c:data:`int2_grad1_u2e_ao`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_size`
       * :c:data:`mpi_master`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`
       * :c:data:`nu_erf`
       * :c:data:`qp_max_mem`
       * :c:data:`read_tc_integ`
       * :c:data:`tc_integ_type`
       * :c:data:`tc_save_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`tc_grad_and_lapl_ao_new`

 
.. c:var:: double int2_grad1_u12_ao_num


    File : :file:`tc_integ_num.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao_num	(ao_num,ao_num,n_points_final_grid,3)
        double precision, allocatable	:: int2_grad1_u12_square_ao_num	(ao_num,ao_num,n_points_final_grid)


    
    int2_grad1_u12_ao_num(i,j,ipoint,:) = \int dr2 [\grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2)
    
    int2_grad1_u12_square_ao_num = -(1/2) x int dr2 chi_l(r2) chi_j(r2) [grad_1 u(r1,r2)]^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`elec_num`
       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`final_grid_points_extra`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_size`
       * :c:data:`mpi_master`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`
       * :c:data:`nu_erf`
       * :c:data:`qp_max_mem`


 
.. c:var:: double int2_grad1_u12_ao_num_1shot


    File : :file:`tc_integ_num.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao_num_1shot	(ao_num,ao_num,n_points_final_grid,3)
        double precision, allocatable	:: int2_grad1_u12_square_ao_num_1shot	(ao_num,ao_num,n_points_final_grid)


    
    int2_grad1_u12_ao_num_1shot(i,j,ipoint,:) = \int dr2 [\grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2)
    
    int2_grad1_u12_square_ao_num_1shot = -(1/2) x int dr2 chi_l(r2) chi_j(r2) [grad_1 u(r1,r2)]^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`final_grid_points_extra`
       * :c:data:`final_grid_points_extra`
       * :c:data:`grad1_u12_num`
       * :c:data:`grad1_u12_num`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_square_ao`

 
.. c:var:: double int2_grad1_u12_ao_test


    File : :file:`new_grad_tc_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao_test	(ao_num,ao_num,n_points_final_grid,3)


    
                                           !!!!!! WARNING !!!!!!!!!
    
                                             DEFINED WITH - SIGN
    
                                      FOR 3e-iontegrals this doesn't matter
    
                                           !!!!!! WARNING !!!!!!!!!
    
    
    int2_grad1_u12_ao_test(i,j,ipoint,:) = \int dr2 [-1 * \grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2)
    
    where r1 = r(ipoint)
    
    if J(r1,r2) = u12:
    
    int2_grad1_u12_ao_test(i,j,ipoint,:) = 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] \phi_i(r2) \phi_j(r2)
                                    = 0.5 * [ v_ij_erf_rk_cst_mu(i,j,ipoint) * r(:) - x_v_ij_erf_rk_cst_mu(i,j,ipoint,:) ]
    
    if J(r1,r2) = u12 x v1 x v2
    
    int2_grad1_u12_ao_test(i,j,ipoint,:) =      v1    x [ 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] v2 \phi_i(r2) \phi_j(r2) ]
                                    - \grad_1 v1 x [       \int dr2                  u12              v2 \phi_i(r2) \phi_j(r2) ]
                                    =  0.5 env_val(ipoint) * v_ij_erf_rk_cst_mu_env(i,j,ipoint) * r(:)
                                    -  0.5 env_val(ipoint) * x_v_ij_erf_rk_cst_mu_env(i,j,ipoint,:)
                                    -  env_grad[:,ipoint] * v_ij_u_cst_mu_env(i,j,ipoint)
    
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`env_grad`
       * :c:data:`env_type`
       * :c:data:`env_val`
       * :c:data:`ezfio_filename`
       * :c:data:`final_grid_points`
       * :c:data:`j2e_type`
       * :c:data:`mpi_master`
       * :c:data:`n_points_final_grid`
       * :c:data:`read_tc_integ`
       * :c:data:`v_ij_erf_rk_cst_mu_env_test`
       * :c:data:`v_ij_u_cst_mu_env_test`
       * :c:data:`x_v_ij_erf_rk_cst_mu_env_test`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_and_lapl_ao_test`

 
.. c:var:: double int2_grad1_u12_square_ao


    File : :file:`tc_integ.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_square_ao	(ao_num,ao_num,n_points_final_grid)


    
    int2_grad1_u12_square_ao = -(1/2) x int dr2 chi_l(r2) chi_j(r2) [grad_1 u(r1,r2)]^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`elec_num`
       * :c:data:`env_grad`
       * :c:data:`env_type`
       * :c:data:`env_val`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`final_grid_points_extra`
       * :c:data:`grad12_j12`
       * :c:data:`int2_grad1_u12_ao_num_1shot`
       * :c:data:`int2_grad1_u12_ao_num_1shot`
       * :c:data:`int2_grad1_u2e_ao`
       * :c:data:`int2_grad1u2_grad2u2`
       * :c:data:`int2_grad1u2_grad2u2_env2`
       * :c:data:`int2_u2_env2`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_size`
       * :c:data:`mpi_master`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`
       * :c:data:`nu_erf`
       * :c:data:`qp_max_mem`
       * :c:data:`tc_integ_type`
       * :c:data:`tc_save_mem`
       * :c:data:`u12_grad1_u12_env_grad1_env`
       * :c:data:`u12sq_envsq`
       * :c:data:`use_ipp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`tc_grad_square_ao_new`

 
.. c:var:: double int2_grad1_u12_square_ao_num


    File : :file:`tc_integ_num.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao_num	(ao_num,ao_num,n_points_final_grid,3)
        double precision, allocatable	:: int2_grad1_u12_square_ao_num	(ao_num,ao_num,n_points_final_grid)


    
    int2_grad1_u12_ao_num(i,j,ipoint,:) = \int dr2 [\grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2)
    
    int2_grad1_u12_square_ao_num = -(1/2) x int dr2 chi_l(r2) chi_j(r2) [grad_1 u(r1,r2)]^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`elec_num`
       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`final_grid_points_extra`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_size`
       * :c:data:`mpi_master`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`
       * :c:data:`nu_erf`
       * :c:data:`qp_max_mem`


 
.. c:var:: double int2_grad1_u12_square_ao_num_1shot


    File : :file:`tc_integ_num.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao_num_1shot	(ao_num,ao_num,n_points_final_grid,3)
        double precision, allocatable	:: int2_grad1_u12_square_ao_num_1shot	(ao_num,ao_num,n_points_final_grid)


    
    int2_grad1_u12_ao_num_1shot(i,j,ipoint,:) = \int dr2 [\grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2)
    
    int2_grad1_u12_square_ao_num_1shot = -(1/2) x int dr2 chi_l(r2) chi_j(r2) [grad_1 u(r1,r2)]^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`final_grid_points_extra`
       * :c:data:`final_grid_points_extra`
       * :c:data:`grad1_u12_num`
       * :c:data:`grad1_u12_num`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_square_ao`

 
.. c:var:: double int2_grad1_u2e_ao


    File : :file:`jast_2e_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u2e_ao	(ao_num,ao_num,n_points_final_grid,3)


    
    int2_grad1_u2e_ao(i,j,ipoint,:) = \int dr2 [\grad_r1 J_2e(r1,r2)] \phi_i(r2) \phi_j(r2)
    
    where r1 = r(ipoint)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`env_grad`
       * :c:data:`env_type`
       * :c:data:`env_val`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`final_grid_points_extra`
       * :c:data:`ir2_mu_gauss_du`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`mpi_master`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`
       * :c:data:`nu_erf`
       * :c:data:`qp_max_mem`
       * :c:data:`tc_integ_type`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_square_ao`
       * :c:data:`j1e_gradx`

 
.. c:var:: double int2_u2e_ao


    File : :file:`jast_2e_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_u2e_ao	(ao_num,ao_num,n_points_final_grid)


    
    int2_u2e_ao(i,j,ipoint,:) = \int dr2 J_2e(r1,r2) \phi_i(r2) \phi_j(r2)
    
    where r1 = r(ipoint)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`env_type`
       * :c:data:`env_val`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`final_grid_points_extra`
       * :c:data:`ir2_mu_gauss_du`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`mpi_master`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`
       * :c:data:`nu_erf`
       * :c:data:`qp_max_mem`
       * :c:data:`tc_integ_type`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`j1e_gradx`

 
.. c:var:: double j1e_gradx


    File : :file:`jast_1e.irp.f`

    .. code:: fortran

        double precision, allocatable	:: j1e_gradx	(n_points_final_grid)
        double precision, allocatable	:: j1e_grady	(n_points_final_grid)
        double precision, allocatable	:: j1e_gradz	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u2e_ao`
       * :c:data:`int2_u2e_ao`
       * :c:data:`j1e_expo`
       * :c:data:`j1e_expo`
       * :c:data:`j1e_size`
       * :c:data:`j1e_type`
       * :c:data:`mo_coef`
       * :c:data:`mpi_master`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_ao_num`
       * :c:data:`int2_grad1_u12_square_ao`

 
.. c:var:: double j1e_grady


    File : :file:`jast_1e.irp.f`

    .. code:: fortran

        double precision, allocatable	:: j1e_gradx	(n_points_final_grid)
        double precision, allocatable	:: j1e_grady	(n_points_final_grid)
        double precision, allocatable	:: j1e_gradz	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u2e_ao`
       * :c:data:`int2_u2e_ao`
       * :c:data:`j1e_expo`
       * :c:data:`j1e_expo`
       * :c:data:`j1e_size`
       * :c:data:`j1e_type`
       * :c:data:`mo_coef`
       * :c:data:`mpi_master`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_ao_num`
       * :c:data:`int2_grad1_u12_square_ao`

 
.. c:var:: double j1e_gradz


    File : :file:`jast_1e.irp.f`

    .. code:: fortran

        double precision, allocatable	:: j1e_gradx	(n_points_final_grid)
        double precision, allocatable	:: j1e_grady	(n_points_final_grid)
        double precision, allocatable	:: j1e_gradz	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u2e_ao`
       * :c:data:`int2_u2e_ao`
       * :c:data:`j1e_expo`
       * :c:data:`j1e_expo`
       * :c:data:`j1e_size`
       * :c:data:`j1e_type`
       * :c:data:`mo_coef`
       * :c:data:`mpi_master`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_ao_num`
       * :c:data:`int2_grad1_u12_square_ao`

 
.. c:var:: double j1e_lapl


    File : :file:`jast_1e.irp.f`

    .. code:: fortran

        double precision, allocatable	:: j1e_lapl	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`j1e_expo`
       * :c:data:`j1e_expo`
       * :c:data:`j1e_size`
       * :c:data:`j1e_type`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`


 
.. c:var:: double j1e_val


    File : :file:`jast_1e.irp.f`

    .. code:: fortran

        double precision, allocatable	:: j1e_val	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`j1e_expo`
       * :c:data:`j1e_expo`
       * :c:data:`j1e_size`
       * :c:data:`j1e_type`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`


 
.. c:var:: double mo_non_hermit_term


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


 
.. c:var:: double mo_non_hermit_term_chemist


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

 
.. c:var:: double mos_grad_in_r_array_qmckl


    File : :file:`qmckl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_in_r_array_qmckl	(mo_num,n_points_final_grid)
        double precision, allocatable	:: mos_grad_in_r_array_qmckl	(mo_num,n_points_final_grid,3)
        double precision, allocatable	:: mos_lapl_in_r_array_qmckl	(mo_num,n_points_final_grid)


    moS computed with qmckl

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`trexio_file`
       * :c:data:`trexio_filename`


 
.. c:var:: double mos_in_r_array_qmckl


    File : :file:`qmckl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_in_r_array_qmckl	(mo_num,n_points_final_grid)
        double precision, allocatable	:: mos_grad_in_r_array_qmckl	(mo_num,n_points_final_grid,3)
        double precision, allocatable	:: mos_lapl_in_r_array_qmckl	(mo_num,n_points_final_grid)


    moS computed with qmckl

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`trexio_file`
       * :c:data:`trexio_filename`


 
.. c:var:: double mos_lapl_in_r_array_qmckl


    File : :file:`qmckl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_in_r_array_qmckl	(mo_num,n_points_final_grid)
        double precision, allocatable	:: mos_grad_in_r_array_qmckl	(mo_num,n_points_final_grid,3)
        double precision, allocatable	:: mos_lapl_in_r_array_qmckl	(mo_num,n_points_final_grid)


    moS computed with qmckl

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`trexio_file`
       * :c:data:`trexio_filename`


 
.. c:var:: long qmckl_ctx_jastrow


    File : :file:`qmckl.irp.f`

    .. code:: fortran

        integer*8	:: qmckl_ctx_jastrow	


    Context for the QMCKL library

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`jast_qmckl_a_vector`
       * :c:data:`jast_qmckl_aord_num`
       * :c:data:`jast_qmckl_b_vector`
       * :c:data:`jast_qmckl_bord_num`
       * :c:data:`jast_qmckl_c_vector`
       * :c:data:`jast_qmckl_c_vector_size`
       * :c:data:`jast_qmckl_cord_num`
       * :c:data:`jast_qmckl_rescale_ee`
       * :c:data:`jast_qmckl_rescale_en`
       * :c:data:`jast_qmckl_spin_independent`
       * :c:data:`jast_qmckl_type_nucl_num`
       * :c:data:`jast_qmckl_type_nucl_vector`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`grad1_u12_num`

 
.. c:var:: double tc_grad_and_lapl_ao_new


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_grad_and_lapl_ao_new	(ao_num,ao_num,ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`n_points_final_grid`
       * :c:data:`tc_integ_type`


 
.. c:var:: double tc_grad_and_lapl_ao_test


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
       * :c:data:`ezfio_filename`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao_test`
       * :c:data:`mpi_master`
       * :c:data:`n_points_final_grid`
       * :c:data:`read_tc_integ`


 
.. c:var:: double tc_grad_square_ao_new


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_grad_square_ao_new	(ao_num,ao_num,ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`env_square_grad`
       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_square_ao`
       * :c:data:`int2_u2_env2`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`n_points_final_grid`
       * :c:data:`tc_integ_type`
       * :c:data:`use_ipp`


 
.. c:var:: double tc_grad_square_ao_test


    File : :file:`grad_squared_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_grad_square_ao_test	(ao_num,ao_num,ao_num,ao_num)


    
    tc_grad_square_ao_test(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_1 u(r1,r2)|^2 | ij>
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`ezfio_filename`
       * :c:data:`final_grid_points`
       * :c:data:`grad12_j12_test`
       * :c:data:`mpi_master`
       * :c:data:`n_points_final_grid`
       * :c:data:`read_tc_integ`
       * :c:data:`u12_grad1_u12_env_grad1_env_test`
       * :c:data:`u12sq_envsq_test`


 
.. c:var:: double tc_grad_square_ao_test_ref


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
       * :c:data:`u12_grad1_u12_env_grad1_env_test`
       * :c:data:`u12sq_envsq_test`


 
.. c:var:: double two_e_int_mf


    File : :file:`mu_of_r_mean_field.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_int_mf	(elec_beta_num,elec_alpha_num,elec_beta_num,elec_alpha_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`


 
.. c:var:: double u12_grad1_u12_env_grad1_env


    File : :file:`grad_squared.irp.f`

    .. code:: fortran

        double precision, allocatable	:: u12_grad1_u12_env_grad1_env	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`env_grad`
       * :c:data:`env_val`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u_grad1u_env2`
       * :c:data:`int2_u_grad1u_x_env2`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_square_ao`

 
.. c:var:: double u12_grad1_u12_env_grad1_env_test


    File : :file:`grad_squared_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: u12_grad1_u12_env_grad1_env_test	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`env_grad`
       * :c:data:`env_val`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u_grad1u_env2_test`
       * :c:data:`int2_u_grad1u_x_env2_test`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_square_ao_test`
       * :c:data:`tc_grad_square_ao_test_ref`

 
.. c:var:: double u12sq_envsq


    File : :file:`grad_squared.irp.f`

    .. code:: fortran

        double precision, allocatable	:: u12sq_envsq	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`env_grad`
       * :c:data:`int2_u2_env2`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_square_ao`

 
.. c:var:: double u12sq_envsq_test


    File : :file:`grad_squared_manu.irp.f`

    .. code:: fortran

        double precision, allocatable	:: u12sq_envsq_test	(ao_num,ao_num,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`env_grad`
       * :c:data:`int2_u2_env2_test`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_grad_square_ao_test`
       * :c:data:`tc_grad_square_ao_test_ref`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void aos_ik_grad1_esquare(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine aos_ik_grad1_esquare(i, k, r1, val)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_aos_ik_grad1_esquare`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_aos_and_grad_at_r`
       * :c:func:`grad1_env_nucl_square_num`

 
.. c:function:: void denom_jpsi(void)


    File : :file:`jastrow_psi.irp.f`

    .. code:: fortran

        subroutine denom_jpsi(i,j,delta,mos_array_r1,mos_grad_array_r1,mos_array_r2,denom, grad_denom)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad_r1_jastrow_psi`
       * :c:func:`routine_deb_denom`

 
.. c:function:: void env_nucl(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        double precision function env_nucl(r)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_expo`
       * :c:data:`env_type`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

 
.. c:function:: void env_nucl_r1_seq(void)


    File : :file:`jast_deriv_utils_vect.irp.f`

    .. code:: fortran

        subroutine env_nucl_r1_seq(n_grid2, res)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_expo`
       * :c:data:`env_type`
       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad1_u12_2e_r1_seq`
       * :c:func:`get_grad1_u12_withsq_r1_seq`
       * :c:func:`get_u12_2e_r1_seq`

 
.. c:function:: void env_nucl_square(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        double precision function env_nucl_square(r)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_expo`
       * :c:data:`env_type`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

 
.. c:function:: void f_mu_and_deriv_mu(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine f_mu_and_deriv_mu(rho, alpha, mu0, beta, f_mu, d_drho_f_mu)


    function giving mu as a function of rho
    
    f_mu = alpha * rho**beta + mu0 * exp(-rho)
    
    and its derivative with respect to rho d_drho_f_mu

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_f_rho`

 
.. c:function:: void f_mu_and_deriv_mu_erf(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine f_mu_and_deriv_mu_erf(rho,alpha,zeta,mu0,beta,f_mu,d_drho_f_mu)


    function giving mu as a function of rho
    
    f_mu = (alpha * rho)**zeta * erf(beta * rho) + mu0 * (1 - erf(beta*rho))
    
    and its derivative with respect to rho d_drho_f_mu
    
    d_drho_f_mu = 2 beta/sqrt(pi) * exp(-(beta*rho)**2) * ( (alpha*rho)**zeta - mu0)
                  + alpha * zeta * (alpha *rho)**(zeta-1)  * erf(beta*rho)

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_f_rho_erf`

 
.. c:function:: void f_mu_and_deriv_mu_simple(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine f_mu_and_deriv_mu_simple(rho, alpha, mu0, beta, f_mu, d_drho_f_mu)


    function giving mu as a function of rho
    
    f_mu = alpha * rho**beta + mu0
    
    and its derivative with respect to rho d_drho_f_mu

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_f_rho_simple`

 
.. c:function:: void get_all_f_rho(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine get_all_f_rho(rho1, rho2, alpha, mu0, beta, f_rho1, d_drho_f_rho1, f_rho2)


    returns the values f(mu(r1)), f(mu(r2)) and d/drho(1) f(mu(r1))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`mu_r_val_and_grad`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`f_mu_and_deriv_mu`

 
.. c:function:: void get_all_f_rho_erf(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine get_all_f_rho_erf(rho1, rho2, alpha, zeta, mu0, beta, f_rho1, d_drho_f_rho1, f_rho2)


    returns the values f(mu(r1)), f(mu(r2)) and d/drho(1) f(mu(r1))
    with f_mu = (alpha * rho)**zeta * erf(beta * rho) + mu0 * (1 - erf(beta*rho))

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`mu_r_val_and_grad`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`f_mu_and_deriv_mu_erf`

 
.. c:function:: void get_all_f_rho_simple(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine get_all_f_rho_simple(rho1,rho2,alpha,mu0,beta,f_rho1,d_drho_f_rho1,f_rho2)


    returns the values f(mu(r1)), f(mu(r2)) and d/drho(1) f(mu(r1))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`mu_r_val_and_grad`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`f_mu_and_deriv_mu_simple`

 
.. c:function:: void get_all_rho_grad_rho(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine get_all_rho_grad_rho(r1, r2, rho1, rho2, grad_rho1)


    returns the density in r1,r2 and grad_rho at r1

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`mu_r_val_and_grad`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`density_and_grad_alpha_beta`

 
.. c:function:: void get_d_da_j_bump(void)


    File : :file:`j_bump.irp.f`

    .. code:: fortran

        subroutine get_d_da_j_bump(x,a,d_da)


    Derivative with respect by to the parameter "a" of the Jastrow with a bump
    
    j(x,a) = 1/2 * (x-a)* exp[-(x/(a*sqrt(pi)))^2]
    
    d/da j(x,a) = - 1/(pi*a^3) * exp[-(x/(a*sqrt(pi)))^2] * (-2 x^3 + 2 a x^2 + pi a^x3)

 
.. c:function:: void get_deriv_mu_j12(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine get_deriv_mu_j12(x,mu,d_d_mu)


    d/dmu j(mu,r12)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`j2e_type`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`deb_deriv_mu`
       * :c:func:`grad_j_sum_mu_of_r`

 
.. c:function:: void get_deriv_r12_j12(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine get_deriv_r12_j12(x,mu,d_dx_j)


    d/dr12 j(mu,r12)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`j2e_type`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_mu`
       * :c:func:`grad_j_sum_mu_of_r`

 
.. c:function:: void get_f_mf_ab(void)


    File : :file:`mu_of_r_mean_field.irp.f`

    .. code:: fortran

        subroutine get_f_mf_ab(r,f_mf_ab,two_bod_dens, dm_a, dm_b)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`two_e_int_mf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`mu_of_r_mean_field`
       * :c:func:`mu_of_r_mean_field_tc`
       * :c:func:`test_grad_f_mean_field`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:function:: void get_grad1_u12_2e_r1_seq(void)


    File : :file:`jast_deriv_utils_vect.irp.f`

    .. code:: fortran

        subroutine get_grad1_u12_2e_r1_seq(ipoint, n_grid2, resx, resy, resz)


    
    grad_1 u_2e(r1,r2)
    
    we use grid for r1 and extra_grid for r2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`nu_erf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u2e_ao`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`env_nucl_r1_seq`
       * :c:func:`grad1_env_nucl`
       * :c:func:`grad1_j12_r1_seq`
       * :c:func:`grad1_jmu_r1_seq`
       * :c:func:`j12_r1_seq`
       * :c:func:`jmu_r1_seq`

 
.. c:function:: void get_grad1_u12_withsq_r1_seq(void)


    File : :file:`jast_deriv_utils_vect.irp.f`

    .. code:: fortran

        subroutine get_grad1_u12_withsq_r1_seq(ipoint, n_grid2, resx, resy, resz, res)


    
    grad_1 u(r1,r2)
    
    we use grid for r1 and extra_grid for r2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`elec_num`
       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_gradx`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_size`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`nu_erf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_ao_num`
       * :c:data:`int2_grad1_u12_square_ao`
       * :c:func:`test_grad1_u12_withsq_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`env_nucl_r1_seq`
       * :c:func:`grad1_env_nucl`
       * :c:func:`grad1_j12_r1_seq`
       * :c:func:`grad1_jmu_r1_seq`
       * :c:func:`j12_r1_seq`
       * :c:func:`jmu_r1_seq`

 
.. c:function:: void get_grad_f_mf_ab(void)


    File : :file:`mu_of_r_mean_field.irp.f`

    .. code:: fortran

        subroutine get_grad_f_mf_ab(r,grad_f_mf_ab, grad_two_bod_dens,f_mf_ab,two_bod_dens, dm_a, dm_b,grad_dm_a, grad_dm_b)


    gradient of mu(r) mean field, together with the gradient of the one- and two-body HF density.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`two_e_int_mf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad_mu_of_r_mean_field`
       * :c:func:`test_grad_f_mean_field`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_and_grad_at_r`

 
.. c:function:: void get_grad_j_bump(void)


    File : :file:`j_bump.irp.f`

    .. code:: fortran

        subroutine get_grad_j_bump(x,a,grad)


    gradient of the Jastrow with a bump
    
    j(x,a) = 1/2 * (x-a)* exp[-(x/(a*sqrt(pi)))^2]
    
    d/dx j(x,a) = 1/(2 pi a^2) * exp[-(x/(a*sqrt(pi)))^2] * (pi a^2 + 2 a x - 2x^2)

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad_j_bump_mu_of_r`

 
.. c:function:: void get_grad_j_bump_mu_of_r(void)


    File : :file:`j_bump.irp.f`

    .. code:: fortran

        subroutine get_grad_j_bump_mu_of_r(r1,r2,grad_j_bump)


    d/dx1 j(x,a(r1,r2)) where j(x,a) is the Jastrow with a bump
    
    j(x,a) = 1/2 * (x-a)* exp[-(x/(a*sqrt(pi)))^2]
    
    a(r1,r2) = [rho(r1) a(r1) + rho(r2) a(r2)]/[rho(r1) + rho(r2)]
    
    d/dx1 j(x,a) = d/dx1 j(x,a(r1,r2))

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_mu`
       * :c:func:`grad1_j12_r1_seq`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad_j_bump`

 
.. c:function:: void get_grad_r1_jastrow_psi(void)


    File : :file:`jastrow_psi.irp.f`

    .. code:: fortran

        subroutine get_grad_r1_jastrow_psi(r1,r2,grad_j_psi_r1,jast)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`c_ij_ab_jastrow`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`log_jpsi`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_r1_seq`
       * :c:func:`routine_deb_j_psi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`denom_jpsi`
       * :c:func:`give_all_mos_and_grad_at_r`
       * :c:func:`numerator_psi`

 
.. c:function:: void get_j1e_coef_fit_ao(void)


    File : :file:`jast_1e_utils.irp.f`

    .. code:: fortran

        subroutine get_j1e_coef_fit_ao(dim_fit, coef_fit)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u2e_ao`
       * :c:data:`j1e_type`
       * :c:data:`mo_coef`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_j1e_fit_ao`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dgemv`
       * :c:func:`get_inverse`
       * :c:func:`wall_time`

 
.. c:function:: void get_j1e_coef_fit_ao2(void)


    File : :file:`jast_1e_utils.irp.f`

    .. code:: fortran

        subroutine get_j1e_coef_fit_ao2(dim_fit, coef_fit)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u2e_ao`
       * :c:data:`j1e_type`
       * :c:data:`mo_coef`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`j1e_gradx`
       * :c:func:`print_j1ecoef`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dgemv`
       * :c:func:`dgesvd`
       * :c:func:`wall_time`

 
.. c:function:: void get_j1e_coef_fit_ao3(void)


    File : :file:`jast_1e_utils.irp.f`

    .. code:: fortran

        subroutine get_j1e_coef_fit_ao3(dim_fit, coef_fit)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u2e_ao`
       * :c:data:`j1e_type`
       * :c:data:`mo_coef`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dgemv`
       * :c:func:`get_inverse`
       * :c:func:`wall_time`

 
.. c:function:: void get_j_sum_mu_of_r(void)


    File : :file:`jast_deriv_mu_of_r.irp.f`

    .. code:: fortran

        subroutine get_j_sum_mu_of_r(r1,r2,jast)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`murho_type`
       * :c:data:`nu_erf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`routine_test_mu_of_r_tot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad_mu_of_r_mean_field`

 
.. c:function:: void get_u12_2e_r1_seq(void)


    File : :file:`jast_deriv_utils_vect.irp.f`

    .. code:: fortran

        subroutine get_u12_2e_r1_seq(ipoint, n_grid2, res)


    
    u_2e(r1,r2)
    
    we use grid for r1 and extra_grid for r2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_type`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`nu_erf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_u2e_ao`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`env_nucl_r1_seq`
       * :c:func:`j12_r1_seq`
       * :c:func:`jmu_r1_seq`

 
.. c:function:: void grad1_aos_ik_grad1_esquare(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine grad1_aos_ik_grad1_esquare(i, k, r1, grad)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_grade_gradu_naive3`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`aos_ik_grad1_esquare`

 
.. c:function:: void grad1_env_nucl(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine grad1_env_nucl(r, grad)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_expo`
       * :c:data:`env_type`
       * :c:data:`list_env1s`
       * :c:data:`list_env1s_size`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad1_u12_2e_r1_seq`
       * :c:func:`get_grad1_u12_withsq_r1_seq`
       * :c:data:`grad1_u12_num`
       * :c:func:`i_grade_gradu_naive1`

 
.. c:function:: void grad1_env_nucl_square_num(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine grad1_env_nucl_square_num(r1, grad)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`aos_ik_grad1_esquare`
       * :c:func:`i_grade_gradu_naive2`

 
.. c:function:: void grad1_j12_mu(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine grad1_j12_mu(r1, r2, grad)


    
     gradient of j(mu(r1,r2),r12) form of jastrow.
    
    if mu(r1,r2) = cst --->
    
     d/dx1 j(mu,r12) = 0.5 * (1 - erf(mu *r12))/r12 * (x1 - x2)
    
    if mu(r1,r2) /= cst --->
    
    d/dx1 j(mu(r1,r2),r12) = exp(-(mu(r1,r2)*r12)**2) /(2 *sqrt(pi) * mu(r1,r2)**2 ) d/dx1 mu(r1,r2)
                           + 0.5 * (1 - erf(mu(r1,r2) *r12))/r12 * (x1 - x2)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`j2e_type`
       * :c:data:`mu_erf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_jmu_modif_num`
       * :c:data:`grad1_u12_num`
       * :c:func:`i_grade_gradu_naive1`
       * :c:func:`num_grad12_j12`
       * :c:func:`num_u12_grad1_u12_env_grad1_env`
       * :c:func:`num_u12sq_envsq`
       * :c:func:`test_fit_ugradu`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_deriv_r12_j12`
       * :c:func:`get_grad_j_bump_mu_of_r`
       * :c:func:`grad_j_sum_mu_of_r`

 
.. c:function:: void grad1_j12_mu_square_num(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine grad1_j12_mu_square_num(r1, r2, grad)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_grade_gradu_naive2`

 
.. c:function:: void grad1_j12_r1_seq(void)


    File : :file:`jast_deriv_utils_vect.irp.f`

    .. code:: fortran

        subroutine grad1_j12_r1_seq(r1, n_grid2, gradx, grady, gradz)


    
    d/dx1 j_2e(1,2)
    d/dy1 j_2e(1,2)
    d/dz1 j_2e(1,2)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`final_grid_points_extra`
       * :c:data:`j2e_type`
       * :c:data:`jbh_ee`
       * :c:data:`jbh_size`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad1_u12_2e_r1_seq`
       * :c:func:`get_grad1_u12_withsq_r1_seq`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad_j_bump_mu_of_r`
       * :c:func:`get_grad_r1_jastrow_psi`
       * :c:func:`grad_j_sum_mu_of_r`
       * :c:func:`jbh_elem_fct_grad`

 
.. c:function:: void grad1_jmu_modif_num(void)


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

       * :c:func:`grad1_j12_mu`

 
.. c:function:: void grad1_jmu_r1_seq(void)


    File : :file:`jast_deriv_utils_vect.irp.f`

    .. code:: fortran

        subroutine grad1_jmu_r1_seq(mu, r1, n_grid2, gradx, grady, gradz)


    
    d/dx1 jmu(r12) = 0.5 * [(1 - erf(mu * r12)) / r12] * (x1 - x2)
    d/dy1 jmu(r12) = 0.5 * [(1 - erf(mu * r12)) / r12] * (y1 - y2)
    d/dz1 jmu(r12) = 0.5 * [(1 - erf(mu * r12)) / r12] * (z1 - z2)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad1_u12_2e_r1_seq`
       * :c:func:`get_grad1_u12_withsq_r1_seq`

 
.. c:function:: void grad1_x_j12_mu_num(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_x_j12_mu_num(r1, r2)



 
.. c:function:: void grad1_x_jmu_modif(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_x_jmu_modif(r1, r2)



 
.. c:function:: void grad1_y_j12_mu_num(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_y_j12_mu_num(r1, r2)



 
.. c:function:: void grad1_y_jmu_modif(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_y_jmu_modif(r1, r2)



 
.. c:function:: void grad1_z_j12_mu_num(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_z_j12_mu_num(r1, r2)



 
.. c:function:: void grad1_z_jmu_modif(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad1_z_jmu_modif(r1, r2)



 
.. c:function:: void grad_j12_mu_input(void)


    File : :file:`jast_deriv_mu_of_r.irp.f`

    .. code:: fortran

        subroutine grad_j12_mu_input(r1, r2, mu, grad_jmu,r12)


    grad_jmu = d/dx1 j(mu,r12) assuming mu=cst(r1)
    
             = 0.5/r_12 * (x_1 - x_2) * [1 - erf(mu*r12)]

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad_j_sum_mu_of_r`

 
.. c:function:: void grad_j_mu_of_r_1(void)


    File : :file:`jast_deriv_mu_of_r.irp.f`

    .. code:: fortran

        subroutine grad_j_mu_of_r_1(r1,r2,jast, grad_jast, dm_r1, grad_dm_r1)


    grad_r1 of j(mu(r1),r12)
    
    
    d/dx1 j(mu(r1),r12) = exp(-(mu(r1)*r12)**2) /(2 *sqrt(pi) * mu(r1)**2 ) d/dx1 mu(r1)
                        + d/dx1 j(mu(r1),r12)
    
    with
    
              j(mu,r12) = 1/2 r12 (1 - erf(mu r12)) - 1/2 (sqrt(pi) * mu) e^{-(mu*r12)^2}
    
    and d/dx1 j(mu,r12) = 0.5 * (1 - erf(mu *r12))/r12 * (x1 - x2)
    
        d/d mu j(mu,r12) = e^{-(r12*mu)^2}/(2 sqrt(pi) * mu^2)
    
    here mu(r1) is obtained by MU MEAN FIELD

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad_j_sum_mu_of_r`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad_mu_of_r_mean_field`

 
.. c:function:: void grad_j_sum_mu_of_r(void)


    File : :file:`jast_deriv_mu_of_r.irp.f`

    .. code:: fortran

        subroutine grad_j_sum_mu_of_r(r1,r2,jast,grad_jast)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`murho_type`
       * :c:data:`nu_erf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_mu`
       * :c:func:`grad1_j12_r1_seq`
       * :c:func:`routine_test_mu_of_r_tot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_deriv_mu_j12`
       * :c:func:`get_deriv_r12_j12`
       * :c:func:`grad_j12_mu_input`
       * :c:func:`grad_j_mu_of_r_1`
       * :c:func:`grad_mu_of_r_mean_field`

 
.. c:function:: void grad_mu_of_r_mean_field(void)


    File : :file:`mu_of_r.irp.f`

    .. code:: fortran

        subroutine grad_mu_of_r_mean_field(r,mu_mf, dm, grad_mu_mf, grad_dm)


    returns the value and gradients of the mu(r) mean field, together with the HF density and its gradients.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_tc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_j_sum_mu_of_r`
       * :c:func:`grad_j_mu_of_r_1`
       * :c:func:`grad_j_sum_mu_of_r`
       * :c:func:`mu_r_val_and_grad`
       * :c:func:`routine_print`
       * :c:func:`routine_test_mu_of_r`
       * :c:func:`test_grad_mu_mf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad_f_mf_ab`

 
.. c:function:: void grad_x_env_nucl_num(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad_x_env_nucl_num(r)



 
.. c:function:: void grad_y_env_nucl_num(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad_y_env_nucl_num(r)



 
.. c:function:: void grad_z_env_nucl_num(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function grad_z_env_nucl_num(r)



 
.. c:function:: void i_grade_gradu_naive1(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine I_grade_gradu_naive1(i, j, k, l, int)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_ipp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_env_nucl`
       * :c:func:`grad1_j12_mu`

 
.. c:function:: void i_grade_gradu_naive2(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine I_grade_gradu_naive2(i, j, k, l, int)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_ipp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_env_nucl_square_num`
       * :c:func:`grad1_j12_mu_square_num`

 
.. c:function:: void i_grade_gradu_naive3(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine I_grade_gradu_naive3(i, j, k, l, int)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_ipp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_aos_ik_grad1_esquare`

 
.. c:function:: void i_grade_gradu_naive4(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine I_grade_gradu_naive4(i, j, k, l, int)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`env_square_grad`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_ipp`

 
.. c:function:: void i_grade_gradu_seminaive(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine I_grade_gradu_seminaive(i, j, k, l, int)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`env_square_grad`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u2_env2`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_ipp`

 
.. c:function:: void j12_and_grad_j12_mu_input(void)


    File : :file:`jast_deriv_mu_of_r.irp.f`

    .. code:: fortran

        subroutine j12_and_grad_j12_mu_input(r1, r2, mu, jmu, grad_jmu)


    jmu = j(mu,r12)
    grad_jmu = d/dx1 j(mu,r12) assuming mu=cst(r1)
    
             = 0.5/r_12 * (x_1 - x_2) * [1 - erf(mu*r12)]

 
.. c:function:: void j12_mu(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        double precision function j12_mu(r1, r2)


    j(mu,r12) = 1/2 r12 (1 - erf(mu r12)) - 1/2 (sqrt(pi) * mu) e^{-(mu*r12)^2}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`j2e_type`
       * :c:data:`mu_erf`

 
.. c:function:: void j12_mu_gauss(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function j12_mu_gauss(r1, r2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`expo_gauss_j_mu_x`
       * :c:data:`n_max_fit_slat`

 
.. c:function:: void j12_mu_input(void)


    File : :file:`jast_deriv_mu_of_r.irp.f`

    .. code:: fortran

        double precision function j12_mu_input(r1, r2, mu)


    j(mu,r12) = 1/2 r12 (1 - erf(mu r12)) - 1/2 (sqrt(pi) * mu) e^{-(mu*r12)^2}

 
.. c:function:: void j12_mu_r12(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function j12_mu_r12(r12)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

 
.. c:function:: void j12_mu_square(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        double precision function j12_mu_square(r1, r2)



 
.. c:function:: void j12_nucl(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function j12_nucl(r1, r2)



 
.. c:function:: void j12_r1_seq(void)


    File : :file:`jast_deriv_utils_vect.irp.f`

    .. code:: fortran

        subroutine j12_r1_seq(r1, n_grid2, res)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`final_grid_points_extra`
       * :c:data:`j2e_type`
       * :c:data:`mu_erf`
       * :c:data:`n_points_extra_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad1_u12_2e_r1_seq`
       * :c:func:`get_grad1_u12_withsq_r1_seq`
       * :c:func:`get_u12_2e_r1_seq`

 
.. c:function:: void j_bump(void)


    File : :file:`j_bump.irp.f`

    .. code:: fortran

        double precision function j_bump(r1,r2,a)



 
.. c:function:: void j_erf_gauss_deriv(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function j_erf_gauss_deriv(x,mu)


    d/dx (j_gauss(mu,x)+j_mu(mu,x))

 
.. c:function:: void j_erf_gauss_deriv_2(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function j_erf_gauss_deriv_2(x,mu)


    d^2/dx^2 (j_gauss(mu,x)+j_mu(mu,x))

 
.. c:function:: void j_gauss_deriv(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function j_gauss_deriv(x,mu)


    d/dx j_gauss(mu,x) = d/dx c/[4 alpha^2 mu] exp(-(alpha * mu * x)^2)
         with c = 27/(8 sqrt(pi)), alpha=3/2
                       = -0.5 * mu *  c * x * exp(-(alpha * mu * x)^2)

 
.. c:function:: void j_gauss_deriv_2(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function j_gauss_deriv_2(x,mu)


    d/dx j_gauss(mu,x) = d/dx c/[4 alpha^2 mu] exp(-(alpha * mu * x)^2)
         with c = 27/(8 sqrt(pi)), alpha=3/2
                       = 0.5 * mu *  c * exp(-(alpha * mu * x)^2) * (2 (alpha*mu*x)^2 - 1)

 
.. c:function:: void j_gauss_simple(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function j_gauss_simple(x,mu)


    j_mu(mu,x) = c/[4 alpha^2 mu] exp(-(alpha * mu * x)^2)
         with c = 27/(8 sqrt(pi)), alpha=3/2

 
.. c:function:: void j_mu_deriv(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function j_mu_deriv(x,mu)


    d/dx j_mu(mu,x) = d/dx 0.5 x (1 - erf(mu x)) - 1/[2 sqrt(pi)mu] exp(-(x*mu)^2)
                 = 0.5*(1 - erf(mu x))

 
.. c:function:: void j_mu_deriv_2(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function j_mu_deriv_2(x,mu)


    d^2/dx^2 j_mu(mu,x) = d^2/dx^2 0.5 x (1 - erf(mu x)) - 1/[2 sqrt(pi)mu] exp(-(x*mu)^2)
                     = -mu/sqrt(pi) * exp(-(mu x)^2)

 
.. c:function:: void j_mu_simple(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function j_mu_simple(x,mu)


    j_mu(mu,x) = 0.5 x (1 - erf(mu x)) - 1/[2 sqrt(pi)mu] exp(-(x*mu)^2)

 
.. c:function:: void j_simple(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function j_simple(x,mu)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`j2e_type`

 
.. c:function:: void jastrow_psi(void)


    File : :file:`jastrow_psi.irp.f`

    .. code:: fortran

        double precision function jastrow_psi(r1,r2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`c_ij_ab_jastrow`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:function:: void jbh_elem_fct_grad(void)


    File : :file:`jast_deriv_utils_vect.irp.f`

    .. code:: fortran

        subroutine jBH_elem_fct_grad(alpha, r1, r2, fct, grad1_fct)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_j12_r1_seq`

 
.. c:function:: void jmu_modif(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function jmu_modif(r1, r2)



 
.. c:function:: void jmu_r1_seq(void)


    File : :file:`jast_deriv_utils_vect.irp.f`

    .. code:: fortran

        subroutine jmu_r1_seq(mu, r1, n_grid2, res)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad1_u12_2e_r1_seq`
       * :c:func:`get_grad1_u12_withsq_r1_seq`
       * :c:func:`get_u12_2e_r1_seq`

 
.. c:function:: void lapl_env_nucl(void)


    File : :file:`j12_nucl_utils.irp.f`

    .. code:: fortran

        double precision function lapl_env_nucl(r)



 
.. c:function:: void mu_of_r_mean_field(void)


    File : :file:`mu_of_r_mean_field.irp.f`

    .. code:: fortran

        subroutine mu_of_r_mean_field(r,mu_mf, dm)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`plot_mu_of_r_mf`
       * :c:func:`test_grad_mu_mf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_f_mf_ab`

 
.. c:function:: void mu_of_r_mean_field_tc(void)


    File : :file:`mu_of_r_mean_field.irp.f`

    .. code:: fortran

        subroutine mu_of_r_mean_field_tc(r,mu_mf, dm)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`plot_mu_of_r_mf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_f_mf_ab`

 
.. c:function:: void mu_r_val_and_grad(void)


    File : :file:`jast_deriv_utils.irp.f`

    .. code:: fortran

        subroutine mu_r_val_and_grad(r1, r2, mu_val, mu_der)


    various flavours of mu(r1,r2)
    depends on essentially the density and other related quantities
    
    change the variable "murho_type" to change type
    
     murho_type == -1 :: mu(r1,r2) = (rho(r1) mu_mf(r1) + rho(r2) mu_mf(r2))/[rho(r1)+rho(r2)]
    
                ==  0 :: mu(r1,r2) = (sqrt(rho(r1)) mu_mf(r1) + sqrt(rho(r2)) mu_mf(r2))/[sqrt(rho(r1))+sqrt(rho(r2))]
    
                == -2 :: mu(r1,r2) = 0.5(mu_mf(r1) + mu_mf(r2))

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`beta_rho_power`
       * :c:data:`env_type`
       * :c:data:`mu_erf`
       * :c:data:`mu_r_ct`
       * :c:data:`murho_type`
       * :c:data:`zeta_erf_mu_of_r`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`density_and_grad_alpha_beta`
       * :c:func:`get_all_f_rho`
       * :c:func:`get_all_f_rho_erf`
       * :c:func:`get_all_f_rho_simple`
       * :c:func:`get_all_rho_grad_rho`
       * :c:func:`grad_mu_of_r_mean_field`

 
.. c:function:: void num_grad12_j12(void)


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

       * :c:func:`grad1_j12_mu`

 
.. c:function:: void num_int2_grad1_u12_ao(void)


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        subroutine num_int2_grad1_u12_ao(i, j, ipoint, integ)


    
    \int dr2 [-grad_1 u12] \phi_i(r2) \phi_j(r2) x v12_env(r1, r2)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_int2_grad1_u12_ao_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad1_jmu_modif_num`

 
.. c:function:: void num_int2_grad1u2_grad2u2_env2(void)


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_int2_grad1u2_grad2u2_env2(i, j, ipoint)


    
    \int dr2 \frac{-[erf(mu r12) -1]^2}{4} \phi_i(r2) \phi_j(r2) x v_env(r2)^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void num_int2_u2_env2(void)


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_int2_u2_env2(i, j, ipoint)


    
    \int dr2 u12^2 \phi_i(r2) \phi_j(r2) x v_env(r2)^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void num_int2_u_grad1u_total_env2(void)


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        subroutine num_int2_u_grad1u_total_env2(i, j, ipoint, integ)


    
    \int dr2 u12 (grad_1 u12) \phi_i(r2) \phi_j(r2) x v_env(r2)^2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_int2_u_grad1u_total_env2`

 
.. c:function:: void num_u12_grad1_u12_env_grad1_env(void)


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_u12_grad1_u12_env_grad1_env(i, j, ipoint)


    
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

       * :c:func:`grad1_j12_mu`

 
.. c:function:: void num_u12sq_envsq(void)


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_u12sq_envsq(i, j, ipoint)


    
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

       * :c:func:`grad1_j12_mu`

 
.. c:function:: void num_v_ij_erf_rk_cst_mu_env(void)


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_v_ij_erf_rk_cst_mu_env(i, j, ipoint)


    
    \int dr2 [erf(mu r12) -1]/r12 \phi_i(r2) \phi_j(r2) x v_env(r2)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void num_v_ij_u_cst_mu_env(void)


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        double precision function num_v_ij_u_cst_mu_env(i, j, ipoint)


    
    \int dr2 u12 \phi_i(r2) \phi_j(r2) x v_env(r2)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void num_x_v_ij_erf_rk_cst_mu_env(void)


    File : :file:`numerical_integ.irp.f`

    .. code:: fortran

        subroutine num_x_v_ij_erf_rk_cst_mu_env(i, j, ipoint, integ)


    
    \int dr2 [erf(mu r12) -1]/r12 \phi_i(r2) \phi_j(r2) x v_env(r2) x r2
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_x_v_ij_erf_rk_cst_mu_env`

 
.. c:function:: void numerator_psi(void)


    File : :file:`jastrow_psi.irp.f`

    .. code:: fortran

        subroutine numerator_psi(a,b,mos_array_r1,mos_grad_array_r1,mos_array_r2,num, grad_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad_r1_jastrow_psi`

 
.. c:function:: void plot_mu_of_r_mf(void)


    File : :file:`test_mu_of_r_tc.irp.f`

    .. code:: fortran

        subroutine plot_mu_of_r_mf



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_mu_of_r_tc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`mu_of_r_mean_field`
       * :c:func:`mu_of_r_mean_field_tc`

 
.. c:function:: void pot_j_gauss(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function pot_j_gauss(x,mu)


    effective scalar potential associated with the erf_gauss correlation factor
    
    1/x( 1 - 2 * d/dx j_erf_gauss(x,mu)) - d^2/dx^2 j_erf_gauss(x,mu)) - d/dx d/dx (j_erf_gauss(x,mu))^2

 
.. c:function:: void pot_j_mu(void)


    File : :file:`pot_j_gauss.irp.f`

    .. code:: fortran

        double precision function pot_j_mu(x,mu)


    effective scalar potential associated with the correlation factor
    
    1/x( 1 - 2 * d/dx j_erf(x,mu)) - d^2/dx^2 j_erf(x,mu)) - d/dx d/dx (j_erf(x,mu))^2

 
.. c:function:: void print_aos(void)


    File : :file:`deb_aos.irp.f`

    .. code:: fortran

        subroutine print_aos()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_in_r_array_qmckl`
       * :c:data:`aos_lapl_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`deb_aos`

 
.. c:function:: void print_j1ecoef(void)


    File : :file:`print_j1ecoef_info.irp.f`

    .. code:: fortran

        subroutine print_j1ecoef()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_l_char_space`
       * :c:data:`ao_num`
       * :c:data:`mpi_master`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_j1ecoef_info`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_get_jastrow_j1e_coef_ao2`
       * :c:func:`ezfio_has_jastrow_j1e_coef_ao2`
       * :c:func:`ezfio_set_jastrow_j1e_coef_ao2`
       * :c:func:`get_j1e_coef_fit_ao2`

 
.. c:function:: void print_mos(void)


    File : :file:`deb_mos.irp.f`

    .. code:: fortran

        subroutine print_mos()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`mos_in_r_array`
       * :c:data:`mos_in_r_array_qmckl`
       * :c:data:`mos_lapl_in_r_array`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`deb_mos`

 
.. c:function:: void routine_deb_denom(void)


    File : :file:`deb_j_psi.irp.f`

    .. code:: fortran

        subroutine routine_deb_denom



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`a_boys`
       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`denom_jpsi`
       * :c:func:`give_all_mos_and_grad_at_r`

 
.. c:function:: void routine_deb_j_psi(void)


    File : :file:`deb_j_psi.irp.f`

    .. code:: fortran

        subroutine routine_deb_j_psi



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`deb_j_psi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad_r1_jastrow_psi`

 
.. c:function:: void routine_fit(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine routine_fit



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

 
.. c:function:: void routine_print(void)


    File : :file:`plot_mu_of_r.irp.f`

    .. code:: fortran

        subroutine routine_print



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`plot_mu_of_r`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad_mu_of_r_mean_field`

 
.. c:function:: void routine_test_mu_of_r(void)


    File : :file:`deb_j_mu_of_r.irp.f`

    .. code:: fortran

        subroutine routine_test_mu_of_r



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad_mu_of_r_mean_field`

 
.. c:function:: void routine_test_mu_of_r_tot(void)


    File : :file:`deb_j_mu_of_r.irp.f`

    .. code:: fortran

        subroutine routine_test_mu_of_r_tot



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`deb_j_bump`
       * :c:func:`deb_j_gauss`
       * :c:func:`deb_j_mu_of_r`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_j_sum_mu_of_r`
       * :c:func:`grad_j_sum_mu_of_r`

 
.. c:function:: void test_env_nucl(void)


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_env_nucl()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_val`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void test_fit_coef_a1(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_fit_coef_A1()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`wall_time`

 
.. c:function:: void test_fit_coef_inv(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_fit_coef_inv()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dgesvd`
       * :c:func:`get_pseudo_inverse`
       * :c:func:`wall_time`

 
.. c:function:: void test_fit_coef_testinva(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_fit_coef_testinvA()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_non_h_ints`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_pseudo_inverse`
       * :c:func:`wall_time`

 
.. c:function:: void test_fit_u(void)


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

 
.. c:function:: void test_fit_u2(void)


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

 
.. c:function:: void test_fit_ugradu(void)


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

       * :c:func:`grad1_j12_mu`

 
.. c:function:: void test_grad12_j12(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_grad12_j12()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`grad12_j12`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void test_grad1_u12_withsq_num(void)


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_grad1_u12_withsq_num()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`grad1_u12_num`
       * :c:data:`grad1_u12_num`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_fit`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_grad1_u12_withsq_r1_seq`

 
.. c:function:: void test_grad_env_nucl(void)


    File : :file:`debug_fit.irp.f`

    .. code:: fortran

        subroutine test_grad_env_nucl()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`env_grad`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void test_grad_f_mean_field(void)


    File : :file:`test_mu_of_r_tc.irp.f`

    .. code:: fortran

        subroutine test_grad_f_mean_field



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_mu_of_r_tc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_f_mf_ab`
       * :c:func:`get_grad_f_mf_ab`

 
.. c:function:: void test_grad_mu_mf(void)


    File : :file:`test_mu_of_r_tc.irp.f`

    .. code:: fortran

        subroutine test_grad_mu_mf



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_mu_of_r_tc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`grad_mu_of_r_mean_field`
       * :c:func:`mu_of_r_mean_field`

 
.. c:function:: void test_int2_grad1_u12_ao(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_int2_grad1_u12_ao()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_ao_num_1shot`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void test_int2_grad1_u12_ao_num(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_int2_grad1_u12_ao_num()



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

 
.. c:function:: void test_int2_grad1_u12_square_ao(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_int2_grad1_u12_square_ao()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_grad1_u12_ao_num_1shot`
       * :c:data:`int2_grad1_u12_square_ao`
       * :c:data:`int2_grad1_u12_ao_num_1shot`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void test_int2_grad1u2_grad2u2_env2(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_int2_grad1u2_grad2u2_env2()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_grad1u2_grad2u2_env2`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void test_int2_u2_env2(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_int2_u2_env2()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_u2_env2`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void test_int2_u_grad1u_total_env2(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_int2_u_grad1u_total_env2()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u_grad1u_env2`
       * :c:data:`int2_u_grad1u_x_env2`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`num_int2_u_grad1u_total_env2`

 
.. c:function:: void test_ipp(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_ipp()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`env_square_grad`
       * :c:data:`final_grid_points`
       * :c:data:`int2_u2_env2`
       * :c:data:`n_points_final_grid`
       * :c:data:`u12_grad1_u12_env_grad1_env`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`i_grade_gradu_naive1`
       * :c:func:`i_grade_gradu_naive2`
       * :c:func:`i_grade_gradu_naive3`
       * :c:func:`i_grade_gradu_naive4`
       * :c:func:`i_grade_gradu_seminaive`

 
.. c:function:: void test_ir2_mu_long_du_0(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_Ir2_Mu_long_Du_0()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ir2_mu_long_du_0`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_erf_rk_cst_mu_env`

 
.. c:function:: void test_j1e_fit_ao(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_j1e_fit_ao()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`int2_grad1_u2e_ao`
       * :c:data:`mo_coef`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dgemv`
       * :c:func:`get_j1e_coef_fit_ao`

 
.. c:function:: void test_j1e_grad(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_j1e_grad()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`int2_grad1_u2e_ao`
       * :c:data:`j1e_gradx`
       * :c:data:`mo_coef`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void test_tc_grad_and_lapl_ao_new(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_tc_grad_and_lapl_ao_new()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ezfio_filename`
       * :c:data:`tc_grad_and_lapl_ao_new`

 
.. c:function:: void test_tc_grad_square_ao_new(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_tc_grad_square_ao_new()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`env_square_grad`
       * :c:data:`env_type`
       * :c:data:`ezfio_filename`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_square_ao`
       * :c:data:`int2_u2_env2`
       * :c:data:`j1e_type`
       * :c:data:`j2e_type`
       * :c:data:`n_points_final_grid`
       * :c:data:`tc_grad_square_ao_new`
       * :c:data:`tc_integ_type`
       * :c:data:`use_ipp`

 
.. c:function:: void test_u12_grad1_u12_env_grad1_env(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_u12_grad1_u12_env_grad1_env()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`u12_grad1_u12_env_grad1_env`

 
.. c:function:: void test_u12sq_envsq(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_u12sq_envsq()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`u12sq_envsq`

 
.. c:function:: void test_v_ij_erf_rk_cst_mu_env(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_v_ij_erf_rk_cst_mu_env()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_erf_rk_cst_mu_env`

 
.. c:function:: void test_v_ij_u_cst_mu_env(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_v_ij_u_cst_mu_env()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_u_cst_mu_env_fit`

 
.. c:function:: void test_v_ij_u_cst_mu_env_an(void)


    File : :file:`test_non_h_ints.irp.f`

    .. code:: fortran

        subroutine test_v_ij_u_cst_mu_env_an()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_u_cst_mu_env_an`
       * :c:data:`v_ij_u_cst_mu_env_an_old`

 
.. c:function:: void test_vect_overlap_gauss_r12_ao(void)


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
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_gauss_r12_ao_v`

 
.. c:function:: void test_vect_overlap_gauss_r12_ao_with1s(void)


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
       * :c:data:`list_env1s_square_coef`
       * :c:data:`mu_erf`
       * :c:data:`n_points_final_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_gauss_r12_ao_with1s_v`

 
.. c:function:: void test_x_v_ij_erf_rk_cst_mu_env(void)


    File : :file:`debug_integ_jmu_modif.irp.f`

    .. code:: fortran

        subroutine test_x_v_ij_erf_rk_cst_mu_env()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_v_ij_erf_rk_cst_mu_env`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`num_x_v_ij_erf_rk_cst_mu_env`

 
.. c:function:: void wigner_radius(void)


    File : :file:`j_bump.irp.f`

    .. code:: fortran

        double precision function wigner_radius(rho)



