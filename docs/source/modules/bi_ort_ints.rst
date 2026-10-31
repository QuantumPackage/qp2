.. _module_bi_ort_ints: 
 
.. program:: bi_ort_ints 
 
.. default-role:: option 
 
bi_ort_ints
===========

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Providers 
--------- 
 
.. c:var:: ao_one_e_integrals_tc_tot


    File : :file:`one_e_bi_ort.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_one_e_integrals_tc_tot	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`j1b_gauss_hermi`
       * :c:data:`j1b_gauss_hermii`
       * :c:data:`j1b_gauss_nonherm`
       * :c:data:`j1b_type`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_one_e`

 
.. c:var:: ao_two_e_tc_tot


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_tc_tot	(ao_num,ao_num,ao_num,ao_num)


    
    ao_two_e_tc_tot(k,i,l,j) = (ki|V^TC(r_12)|lj) = <lk| V^TC(r_12) |ji> where V^TC(r_12) is the total TC operator
    
    including both hermitian and non hermitian parts. THIS IS IN CHEMIST NOTATION.
    
    WARNING :: non hermitian ! acts on "the right functions" (i,j)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_non_hermit_term_chemist`
       * :c:data:`ao_num`
       * :c:data:`ao_tc_int_chemist`
       * :c:data:`ao_tc_sym_two_e_pot_cache`
       * :c:data:`ao_tc_sym_two_e_pot_cache_min`
       * :c:data:`ao_tc_sym_two_e_pot_in_map`
       * :c:data:`ao_tc_sym_two_e_pot_map`
       * :c:data:`j1b_type`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e_chemist`

 
.. c:var:: ao_two_e_vartc_tot


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_vartc_tot	(ao_num,ao_num,ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_vartc_int_chemist`
       * :c:data:`j1b_type`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`


 
.. c:var:: int2_grad1_u12_ao_t


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao_t	(n_points_final_grid,3,ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`n_points_final_grid`


 
.. c:var:: int2_grad1_u12_ao_transp


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao_transp	(ao_num,ao_num,3,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`int2_grad1_u12_ao_test`
       * :c:data:`n_points_final_grid`
       * :c:data:`test_cycle_tc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_bimo_transp`

 
.. c:var:: int2_grad1_u12_bimo_t


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_bimo_t	(n_points_final_grid,3,mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_bimo_transp`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`three_body_ints_bi_ort`
       * :c:data:`three_e_3_idx_cycle_1_bi_ort`
       * :c:data:`three_e_3_idx_cycle_2_bi_ort`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort_new`
       * :c:data:`three_e_3_idx_exch13_bi_ort`
       * :c:data:`three_e_3_idx_exch23_bi_ort`
       * :c:data:`three_e_4_idx_cycle_1_bi_ort`
       * :c:data:`three_e_4_idx_cycle_2_bi_ort`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_exch12_bi_ort`
       * :c:data:`three_e_4_idx_exch13_bi_ort`
       * :c:data:`three_e_4_idx_exch23_bi_ort`
       * :c:data:`three_e_5_idx_cycle_1_bi_ort`
       * :c:data:`three_e_5_idx_cycle_2_bi_ort`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_5_idx_exch13_bi_ort`
       * :c:data:`three_e_5_idx_exch23_bi_ort`

 
.. c:var:: int2_grad1_u12_bimo_transp


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_bimo_transp	(mo_num,mo_num,3,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_grad1_u12_ao_transp`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_bimo_t`

 
.. c:var:: mo_bi_orth_bipole_x


    File : :file:`one_e_bi_ort.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_orth_bipole_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_orth_bipole_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_orth_bipole_z	(mo_num,mo_num)


    array of the integrals of Left MO_i * x Right MO_j
    array of the integrals of Left MO_i * y Right MO_j
    array of the integrals of Left MO_i * z Right MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_dipole_x`
       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: mo_bi_orth_bipole_y


    File : :file:`one_e_bi_ort.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_orth_bipole_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_orth_bipole_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_orth_bipole_z	(mo_num,mo_num)


    array of the integrals of Left MO_i * x Right MO_j
    array of the integrals of Left MO_i * y Right MO_j
    array of the integrals of Left MO_i * z Right MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_dipole_x`
       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: mo_bi_orth_bipole_z


    File : :file:`one_e_bi_ort.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_orth_bipole_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_orth_bipole_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_orth_bipole_z	(mo_num,mo_num)


    array of the integrals of Left MO_i * x Right MO_j
    array of the integrals of Left MO_i * y Right MO_j
    array of the integrals of Left MO_i * z Right MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_dipole_x`
       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: mo_bi_ortho_coul_e


    File : :file:`biorthog_mo_for_h.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_coul_e	(mo_num,mo_num,mo_num,mo_num)


    
    mo_bi_ortho_coul_e(k,l,i,j) = < k l | 1/r12 | i j > where i,j are right MOs and k,l are left MOs
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_coul_e_chemist`
       * :c:data:`mo_num`


 
.. c:var:: mo_bi_ortho_coul_e_chemist


    File : :file:`biorthog_mo_for_h.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_coul_e_chemist	(mo_num,mo_num,mo_num,mo_num)


    
    mo_bi_ortho_coul_e_chemist(k,i,l,j) = < k l | 1/r12 | i j > where i,j are right MOs and k,l are left MOs
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_coul`
       * :c:data:`mo_l_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef_transp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_coul_e`

 
.. c:var:: mo_bi_ortho_one_e


    File : :file:`biorthog_mo_for_h.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_one_e	(mo_num,mo_num)


    
    mo_bi_ortho_one_e(k,i) = < MO^L_k | h_c | MO^R_i >
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: mo_bi_ortho_tc_one_e


    File : :file:`one_e_bi_ort.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_one_e	(mo_num,mo_num)


    
    mo_bi_ortho_tc_one_e(k,i) = <MO^L_k | h_c | MO^R_i>
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals_tc_tot`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: mo_bi_ortho_tc_two_e


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_two_e	(mo_num,mo_num,mo_num,mo_num)


    
    mo_bi_ortho_tc_two_e(k,l,i,j) = <k l| V(r_12) |i j> where i,j are right MOs and k,l are left MOs
    
    the potential V(r_12) contains ALL TWO-E CONTRIBUTION OF THE TC-HAMILTONIAN
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e_chemist`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e_jj`
       * :c:data:`tc_2e_3idx_coulomb_integrals`

 
.. c:var:: mo_bi_ortho_tc_two_e_chemist


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_two_e_chemist	(mo_num,mo_num,mo_num,mo_num)


    
    mo_bi_ortho_tc_two_e_chemist(k,i,l,j) = <k l|V(r_12)|i j> where i,j are right MOs and k,l are left MOs
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`mo_l_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef_transp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`

 
.. c:var:: mo_bi_ortho_tc_two_e_jj


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_anti	(mo_num,mo_num)


    mo_bi_ortho_tc_two_e_jj(i,j) = J_ij = <ji|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_exchange(i,j) = K_ij = <ij|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_anti(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`


 
.. c:var:: mo_bi_ortho_tc_two_e_jj_anti


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_anti	(mo_num,mo_num)


    mo_bi_ortho_tc_two_e_jj(i,j) = J_ij = <ji|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_exchange(i,j) = K_ij = <ij|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_anti(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`


 
.. c:var:: mo_bi_ortho_tc_two_e_jj_exchange


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_anti	(mo_num,mo_num)


    mo_bi_ortho_tc_two_e_jj(i,j) = J_ij = <ji|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_exchange(i,j) = K_ij = <ij|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_anti(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`


 
.. c:var:: mo_v_ki_bi_ortho_erf_rk_cst_mu


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_v_ki_bi_ortho_erf_rk_cst_mu	(mo_num,mo_num,n_points_final_grid)


    
    mo_v_ki_bi_ortho_erf_rk_cst_mu(k,i,ip) = int dr chi_k(r) phi_i(r) (erf(mu |r - R_ip|) - 1 )/(2|r - R_ip|) on the BI-ORTHO MO basis
    
    where phi_k(r) is a LEFT MOs and phi_i(r) is a RIGHT MO
    
    R_ip = the "ip"-th point of the DFT Grid
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_points_final_grid`
       * :c:data:`v_ij_erf_rk_cst_mu`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_v_ki_bi_ortho_erf_rk_cst_mu_transp`

 
.. c:var:: mo_v_ki_bi_ortho_erf_rk_cst_mu_transp


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_v_ki_bi_ortho_erf_rk_cst_mu_transp	(n_points_final_grid,mo_num,mo_num)


    
    int dr phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1)/(2|r - R|) on the BI-ORTHO MO basis
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mo_v_ki_bi_ortho_erf_rk_cst_mu`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`x_w_ki_bi_ortho_erf_rk`
       * :c:data:`x_w_ki_bi_ortho_erf_rk_diag`

 
.. c:var:: mo_x_v_ki_bi_ortho_erf_rk_cst_mu


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_x_v_ki_bi_ortho_erf_rk_cst_mu	(mo_num,mo_num,3,n_points_final_grid)


    
    mo_x_v_ki_bi_ortho_erf_rk_cst_mu(k,i,m,ip) = int dr x(m) * chi_k(r) phi_i(r) (erf(mu |r - R_ip|) - 1)/2|r - R_ip| on the BI-ORTHO MO basis
    
    where chi_k(r)/phi_i(r) are left/right MOs, m=1 => x(m) = x, m=2 => x(m) = y, m=3 => x(m) = z,
    
    R_ip = the "ip"-th point of the DFT Grid
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_points_final_grid`
       * :c:data:`x_v_ij_erf_rk_cst_mu_transp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp`

 
.. c:var:: mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp	(n_points_final_grid,3,mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mo_x_v_ki_bi_ortho_erf_rk_cst_mu`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`x_w_ki_bi_ortho_erf_rk`
       * :c:data:`x_w_ki_bi_ortho_erf_rk_diag`

 
.. c:var:: tc_2e_3idx_coulomb_integrals


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_2e_3idx_coulomb_integrals	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: tc_2e_3idx_exchange_integrals	(mo_num,mo_num,mo_num)


    tc_2e_3idx_coulomb_integrals(j,k,i)  = <jk|ji>
    
    tc_2e_3idx_exchange_integrals(j,k,i) = <kj|ji>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`


 
.. c:var:: tc_2e_3idx_exchange_integrals


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_2e_3idx_coulomb_integrals	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: tc_2e_3idx_exchange_integrals	(mo_num,mo_num,mo_num)


    tc_2e_3idx_coulomb_integrals(j,k,i)  = <jk|ji>
    
    tc_2e_3idx_exchange_integrals(j,k,i) = <kj|ji>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`


 
.. c:var:: three_body_ints_bi_ort


    File : :file:`three_body_ints_bi_ort.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_body_ints_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num,mo_num)


    matrix element of the -L  three-body operator
    
    notice the -1 sign: in this way three_body_ints_bi_ort can be directly used to compute Slater rules :)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_3_idx_cycle_1_bi_ort


    File : :file:`three_body_ijm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_3_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator ON A BI ORTHONORMAL BASIS for the first cyclic permutation
    
    three_e_3_idx_cycle_1_bi_ort(m,j,i) = <mji|-L|jim>
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_3_idx_cycle_2_bi_ort


    File : :file:`three_body_ijm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_3_idx_cycle_2_bi_ort	(mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator ON A BI ORTHONORMAL BASIS for the second cyclic permutation
    
    three_e_3_idx_direct_bi_ort(m,j,i) = <mji|-L|imj>
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_3_idx_direct_bi_ort


    File : :file:`three_body_ijm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_3_idx_direct_bi_ort	(mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator ON A BI ORTHONORMAL BASIS for the direct terms
    
    three_e_3_idx_direct_bi_ort(m,j,i) = <mji|-L|mji>
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_3_idx_exch12_bi_ort


    File : :file:`three_body_ijm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_3_idx_exch12_bi_ort	(mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator ON A BI ORTHONORMAL BASIS for the permutations of particle 1 and 2
    
    three_e_3_idx_exch12_bi_ort(m,j,i) = <mji|-L|mij>
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_3_idx_exch12_bi_ort_new


    File : :file:`three_body_ijm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_3_idx_exch12_bi_ort_new	(mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator ON A BI ORTHONORMAL BASIS for the permutations of particle 1 and 2
    
    three_e_3_idx_exch12_bi_ort_new(m,j,i) = <mji|-L|mij>
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_3_idx_exch13_bi_ort


    File : :file:`three_body_ijm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_3_idx_exch13_bi_ort	(mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator ON A BI ORTHONORMAL BASIS for the permutations of particle 1 and 3
    
    three_e_3_idx_exch13_bi_ort(m,j,i) = <mji|-L|ijm>
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_3_idx_exch23_bi_ort


    File : :file:`three_body_ijm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_3_idx_exch23_bi_ort	(mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator ON A BI ORTHONORMAL BASIS for the permutations of particle 2 and 3
    
    three_e_3_idx_exch23_bi_ort(m,j,i) = <mji|-L|jmi>
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_4_idx_cycle_1_bi_ort


    File : :file:`three_body_ijmk.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_cycle_1_bi_ort(m,j,k,i) = <mjk|-L|jim> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_4_idx_cycle_2_bi_ort


    File : :file:`three_body_ijmk.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_cycle_2_bi_ort	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_cycle_2_bi_ort(m,j,k,i) = <mjk|-L|imj> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_4_idx_direct_bi_ort


    File : :file:`three_body_ijmk.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_direct_bi_ort(m,j,k,i) = <mjk|-L|mji> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_4_idx_exch12_bi_ort


    File : :file:`three_body_ijmk.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_exch12_bi_ort	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_exch12_bi_ort(m,j,k,i) = <mjk|-L|mij> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_4_idx_exch13_bi_ort


    File : :file:`three_body_ijmk.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_exch13_bi_ort(m,j,k,i) = <mjk|-L|ijm> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_4_idx_exch23_bi_ort


    File : :file:`three_body_ijmk.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_exch23_bi_ort(m,j,k,i) = <mjk|-L|jmi> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_5_idx_cycle_1_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) = <mlk|-L|jim> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_5_idx_cycle_2_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_cycle_2_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i) = <mlk|-L|imj> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_5_idx_direct_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_direct_bi_ort(m,l,j,k,i) = <mlk|-L|mji> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_5_idx_exch12_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_exch12_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_exch12_bi_ort(m,l,j,k,i) = <mlk|-L|mij> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_5_idx_exch13_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_exch13_bi_ort(m,l,j,k,i) = <mlk|-L|ijm> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: three_e_5_idx_exch23_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_exch23_bi_ort(m,l,j,k,i) = <mlk|-L|jmi> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: x_w_ki_bi_ortho_erf_rk


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: x_w_ki_bi_ortho_erf_rk	(n_points_final_grid,3,mo_num,mo_num)


    
    x_W_ki_bi_ortho_erf_rk(ip,m,k,i) = \int dr chi_k(r) \frac{(1 - erf(mu |r-R_ip|))}{2|r-R_ip|} (x(m)-R_ip(m)) phi_i(r) ON THE BI-ORTHO MO BASIS
    
    where chi_k(r)/phi_i(r) are left/right MOs, m=1 => X(m) = x, m=2 => X(m) = y, m=3 => X(m) = z,
    
    R_ip = the "ip"-th point of the DFT Grid

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`mo_v_ki_bi_ortho_erf_rk_cst_mu_transp`
       * :c:data:`mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: x_w_ki_bi_ortho_erf_rk_diag


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: x_w_ki_bi_ortho_erf_rk_diag	(n_points_final_grid,3,mo_num)


    x_W_ki_bi_ortho_erf_rk_diag(ip,m,i) = \int dr chi_i(r) (1 - erf(mu |r-R_ip|)) (x(m)-X(m)_ip) phi_i(r) ON THE BI-ORTHO MO BASIS
    
    where chi_k(r)/phi_i(r) are left/right MOs, m=1 => X(m) = x, m=2 => X(m) = y, m=3 => X(m) = z,
    
    R_ip = the "ip"-th point of the DFT Grid

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`mo_v_ki_bi_ortho_erf_rk_cst_mu_transp`
       * :c:data:`mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp`
       * :c:data:`n_points_final_grid`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: bi_ort_ints:


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine bi_ort_ints


    TODO : Put the documentation of the program here

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`test_3e`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

 
.. c:function:: bi_ortho_mo_coul_ints:


    File : :file:`biorthog_mo_for_h.irp.f`

    .. code:: fortran

        double precision function bi_ortho_mo_coul_ints(l, k, j, i)


    
    < mo^L_k mo^L_l | 1/r12 | mo^R_i mo^R_j >
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_coul`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

 
.. c:function:: bi_ortho_mo_ints:


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision function bi_ortho_mo_ints(l, k, j, i)


    
    <mo^L_k mo^L_l | V^TC(r_12) | mo^R_i mo^R_j>
    
    WARNING :: very naive, super slow, only used to DEBUG.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

 
.. c:function:: give_integrals_3_body_bi_ort:


    File : :file:`three_body_ints_bi_ort.irp.f`

    .. code:: fortran

        subroutine give_integrals_3_body_bi_ort(n, l, k, m, j, i, integral)


    
    < n l k | -L | m j i > with a BI-ORTHONORMAL MOLECULAR ORBITALS
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_3e`
       * :c:data:`three_body_ints_bi_ort`
       * :c:data:`three_e_3_idx_cycle_1_bi_ort`
       * :c:data:`three_e_3_idx_cycle_2_bi_ort`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort_new`
       * :c:data:`three_e_3_idx_exch13_bi_ort`
       * :c:data:`three_e_3_idx_exch23_bi_ort`
       * :c:data:`three_e_4_idx_cycle_1_bi_ort`
       * :c:data:`three_e_4_idx_cycle_2_bi_ort`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_exch12_bi_ort`
       * :c:data:`three_e_4_idx_exch13_bi_ort`
       * :c:data:`three_e_4_idx_exch23_bi_ort`
       * :c:data:`three_e_5_idx_cycle_1_bi_ort`
       * :c:data:`three_e_5_idx_cycle_2_bi_ort`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_exch12_bi_ort`
       * :c:data:`three_e_5_idx_exch13_bi_ort`
       * :c:data:`three_e_5_idx_exch23_bi_ort`

 
.. c:function:: give_integrals_3_body_bi_ort_ao:


    File : :file:`three_body_ints_bi_ort.irp.f`

    .. code:: fortran

        subroutine give_integrals_3_body_bi_ort_ao(n, l, k, m, j, i, integral)


    
    < n l k | -L | m j i > with a BI-ORTHONORMAL ATOMIC ORBITALS
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao_t`
       * :c:data:`n_points_final_grid`

 
.. c:function:: give_integrals_3_body_bi_ort_old:


    File : :file:`three_body_ints_bi_ort.irp.f`

    .. code:: fortran

        subroutine give_integrals_3_body_bi_ort_old(n, l, k, m, j, i, integral)


    
    < n l k | -L | m j i > with a BI-ORTHONORMAL MOLECULAR ORBITALS
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_3e`

 
.. c:function:: test_3e:


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_3e



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`bi_ort_ints`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`
       * :c:func:`give_integrals_3_body_bi_ort_old`

