.. _module_bi_ort_ints: 
 
.. program:: bi_ort_ints 
 
.. default-role:: option 
 
===========
bi_ort_ints
===========

This module contains all necessary integrals for the TC Hamiltonian in a bi-orthonormal (BO) MO Basis.
See in bi_ortho_basis for more information. 
The main providers are : 

One-electron integrals 
----------------------
+) ao_one_e_integrals_tc_tot : total one-electron Hamiltonian which might include non hermitian part coming from one-e correlation factor. 
+) mo_bi_ortho_tc_one_e : one-electron Hamiltonian (h_core+one-J terms) on the BO-MO basis. 
+) mo_bi_orth_bipole_x  : x-component of the dipole operator on the BO-MO basis. (Same for y,z) 

Two-electron integrals 
----------------------
+) ao_two_e_tc_tot : Total two-electron operator (including the non-hermitian term of the TC Hamiltonian) on the AO basis
+) mo_bi_ortho_tc_two_e : Total two-electron operator on the BO-MO basis

Three-electron integrals 
------------------------
+) three_body_ints_bi_ort : 6-indices three-electron tensor (-L) on the BO-MO basis. WARNING :: N^6 storage !
+) three_e_3_idx_direct_bi_ort : DIRECT term with 3 different indices of the -L operator. These terms appear in the DIAGONAL matrix element of the -L operator. The 5 other permutations needed to compute matrix elements can be found in three_body_ijm.irp.f 
+) three_e_4_idx_direct_bi_ort : DIRECT term with 4 different indices of the -L operator. These terms appear in the OFF-DIAGONAL matrix element of the -L operator including SINGLE EXCITATIONS. The 5 other permutations needed to compute matrix elements can be found in three_body_ijmk.irp.f 
+) three_e_5_idx_direct_bi_ort : DIRECT term with 5 different indices of the -L operator. These terms appear in the OFF-DIAGONAL matrix element of the -L operator including DOUBLE EXCITATIONS. The 5 other permutations needed to compute matrix elements can be found in three_body_ijmkl.irp.f 
 
 
 
Programs 
-------- 
 
 * :ref:`bi_ort_ints` 
 
Providers 
--------- 
 
.. c:var:: double ao_one_e_integrals_tc_tot


    File : :file:`one_e_bi_ort.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_one_e_integrals_tc_tot	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_one_e`

 
.. c:var:: double ao_two_e_coul


    File : :file:`biorthog_mo_for_h.irp.f`

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

       * :c:data:`mo_bi_ortho_coul_e_chemist`

 
.. c:var:: double energy_1e_nol_hf


    File : :file:`no_dressing_energy.irp.f`

    .. code:: fortran

        double precision	:: energy_1e_nol_hf	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_bi_ortho_tc_one_e`


 
.. c:var:: double energy_2e_nol_hf


    File : :file:`no_dressing_energy.irp.f`

    .. code:: fortran

        double precision	:: energy_2e_nol_hf	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_bi_ortho_tc_two_e`


 
.. c:var:: double int2_grad1_u12_ao_t


    File : :file:`semi_num_ints_mo_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_ao_t	(n_points_final_grid,3,ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double int2_grad1_u12_ao_transp


    File : :file:`semi_num_ints_mo_old.irp.f`

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
       * :c:data:`three_body_ints_bi_ort`

 
.. c:var:: double int2_grad1_u12_bimo_t


    File : :file:`semi_num_ints_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_bimo_t	(n_points_final_grid,3,mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`int2_grad1_u12_ao`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`nol_0e`
       * :c:data:`nol_0e_naive`
       * :c:data:`nol_0e_v0`
       * :c:data:`nol_1e`
       * :c:data:`nol_1e_naive`
       * :c:data:`nol_1e_v0`
       * :c:data:`nol_2e`
       * :c:data:`nol_2e_naive`
       * :c:data:`nol_2e_v0`
       * :c:data:`three_body_ints_bi_ort`
       * :c:data:`three_e_3_idx_cycle_1_bi_ort`
       * :c:data:`three_e_3_idx_cycle_2_bi_ort`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort_new`
       * :c:data:`three_e_3_idx_exch13_bi_ort`
       * :c:data:`three_e_3_idx_exch23_bi_ort`
       * :c:data:`three_e_4_idx_cycle_1_bi_ort_old`
       * :c:data:`three_e_4_idx_cycle_2_bi_ort_old`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_direct_bi_ort_n4`
       * :c:data:`three_e_4_idx_direct_bi_ort_old`
       * :c:data:`three_e_4_idx_exch12_bi_ort_old`
       * :c:data:`three_e_4_idx_exch13_bi_ort_old`
       * :c:data:`three_e_4_idx_exch23_bi_ort_n4`
       * :c:data:`three_e_4_idx_exch23_bi_ort_old`
       * :c:data:`three_e_5_idx_cycle_1_bi_ort_old`
       * :c:data:`three_e_5_idx_cycle_2_bi_ort_old`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_direct_bi_ort_old`
       * :c:data:`three_e_5_idx_exch12_bi_ort_old`
       * :c:data:`three_e_5_idx_exch13_bi_ort_old`
       * :c:data:`three_e_5_idx_exch23_bi_ort_old`

 
.. c:var:: double int2_grad1_u12_bimo_t_old


    File : :file:`semi_num_ints_mo_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int2_grad1_u12_bimo_t_old	(n_points_final_grid,3,mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_bimo_transp`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double int2_grad1_u12_bimo_transp


    File : :file:`semi_num_ints_mo_old.irp.f`

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

       * :c:data:`int2_grad1_u12_bimo_t_old`

 
.. c:var:: double mo_bi_orth_bipole_x


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


 
.. c:var:: double mo_bi_orth_bipole_y


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


 
.. c:var:: double mo_bi_orth_bipole_z


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


 
.. c:var:: double mo_bi_ortho_coul_e


    File : :file:`biorthog_mo_for_h.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_coul_e	(mo_num,mo_num,mo_num,mo_num)


    
    mo_bi_ortho_coul_e(k,l,i,j) = < k l | 1/r12 | i j > where i,j are right MOs and k,l are left MOs
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_coul_e_chemist`
       * :c:data:`mo_num`


 
.. c:var:: double mo_bi_ortho_coul_e_chemist


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

 
.. c:var:: double mo_bi_ortho_eri_two_e_chemist


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_eri_two_e_chemist	(mo_num,mo_num,mo_num,mo_num)


    
    mo_bi_ortho_eri_two_e_chemist(k,i,l,j) = <k l|V(r_12)|i j> where i,j are right MOs and k,l are left MOs
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_ao_transp`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e_chemist`

 
.. c:var:: double mo_bi_ortho_one_e


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


 
.. c:var:: double mo_bi_ortho_tc_one_e


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
       * :c:data:`nol_1e`
       * :c:data:`nol_standard`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_1e_nol_hf`

 
.. c:var:: double mo_bi_ortho_tc_two_e


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
       * :c:data:`nol_2e`
       * :c:data:`nol_standard`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_2e_nol_hf`
       * :c:data:`mo_bi_ortho_tc_two_e_jj`
       * :c:data:`mo_bi_ortho_tc_two_e_transp`
       * :c:data:`tc_2e_3idx_coulomb_integrals`

 
.. c:var:: double mo_bi_ortho_tc_two_e_chemist


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_two_e_chemist	(mo_num,mo_num,mo_num,mo_num)


    
    mo_bi_ortho_tc_two_e_chemist(k,i,l,j) = <k l|V(r_12)|i j> where i,j are right MOs and k,l are left MOs
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_to_mo_tc_n3`
       * :c:data:`ao_two_e_tc_tot`
       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_bi_ortho_eri_two_e_chemist`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`n_int`
       * :c:data:`read_tc_integ`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`

 
.. c:var:: double mo_bi_ortho_tc_two_e_jj


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_anti	(mo_num,mo_num)


    
    mo_bi_ortho_tc_two_e_jj         (i,j) = J_ij = <ji|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_exchange(i,j) = K_ij = <ij|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_anti    (i,j) = J_ij - K_ij
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`


 
.. c:var:: double mo_bi_ortho_tc_two_e_jj_anti


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_anti	(mo_num,mo_num)


    
    mo_bi_ortho_tc_two_e_jj         (i,j) = J_ij = <ji|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_exchange(i,j) = K_ij = <ij|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_anti    (i,j) = J_ij - K_ij
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`


 
.. c:var:: double mo_bi_ortho_tc_two_e_jj_exchange


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_bi_ortho_tc_two_e_jj_anti	(mo_num,mo_num)


    
    mo_bi_ortho_tc_two_e_jj         (i,j) = J_ij = <ji|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_exchange(i,j) = K_ij = <ij|W-K|ji>
    mo_bi_ortho_tc_two_e_jj_anti    (i,j) = J_ij - K_ij
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`


 
.. c:var:: double mo_bi_ortho_tc_two_e_transp


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_bi_ortho_tc_two_e_transp	(mo_num,mo_num,mo_num,mo_num)


    
    mo_bi_ortho_tc_two_e_transp(i,j,k,l) = <k l| V(r_12) |i j> = transpose of mo_bi_ortho_tc_two_e
    
    the potential V(r_12) contains ALL TWO-E CONTRIBUTION OF THE TC-HAMILTONIAN
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`tc_2e_3idx_coulomb_integrals_transp`

 
.. c:var:: double mo_v_ki_bi_ortho_erf_rk_cst_mu


    File : :file:`semi_num_ints_mo_old.irp.f`

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

 
.. c:var:: double mo_v_ki_bi_ortho_erf_rk_cst_mu_transp


    File : :file:`semi_num_ints_mo_old.irp.f`

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

 
.. c:var:: double mo_x_v_ki_bi_ortho_erf_rk_cst_mu


    File : :file:`semi_num_ints_mo_old.irp.f`

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

 
.. c:var:: double mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp


    File : :file:`semi_num_ints_mo_old.irp.f`

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

 
.. c:var:: double nol_0e


    File : :file:`no_dressing.irp.f`

    .. code:: fortran

        double precision	:: nol_0e	


    
    < Phi_left | L | Phi_right >
    
    ZERO-body part of the normal ordering (on HF reference) of 3-e operator
    
    WARNING: if core_tc_op
    
    then does not take into account any contribution from core orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_tc_op`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double nol_0e_naive


    File : :file:`no_dressing_naive.irp.f`

    .. code:: fortran

        double precision	:: nol_0e_naive	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double nol_0e_v0


    File : :file:`no_dressing.irp.f`

    .. code:: fortran

        double precision	:: nol_0e_v0	


    
    < Phi_left | L | Phi_right > COMPUTED IN A NAIVE WAY
    
    ZERO-body part of the normal ordering (on HF reference) of 3-e operator
    
    WARNING: if core_tc_op
    
    then does not take into account any contribution from core orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_tc_op`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double nol_1e


    File : :file:`no_dressing.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nol_1e	(mo_num,mo_num)


    one-body part of the normal ordering (on HF reference) of 3-e operator
    
    WARNING: if core_tc_op
    
    then does not take into account any contribution from core orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_one_e`

 
.. c:var:: double nol_1e_naive


    File : :file:`no_dressing_naive.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nol_1e_naive	(mo_num,mo_num)


    
    < p | H(1) | s > is dressed with noL_1e_naive(p,s)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double nol_1e_v0


    File : :file:`no_dressing.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nol_1e_v0	(mo_num,mo_num)


    one-body part of the normal ordering (on HF reference) of 3-e operator COMPUTED IN A NAIVE WAY
    
    WARNING: if core_tc_op
    
    then does not take into account any contribution from core orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double nol_2e


    File : :file:`no_dressing.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nol_2e	(mo_num,mo_num,mo_num,mo_num)


    two-body part of the normal ordering (on HF reference) of 3-e operator
    
    WARNING: if core_tc_op
    
    then does not take into account any contribution from core orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`

 
.. c:var:: double nol_2e_naive


    File : :file:`no_dressing_naive.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nol_2e_naive	(mo_num,mo_num,mo_num,mo_num)


    
    < p q | H(2) | s t > is dressed with noL_2e_naive(p,q,s,t)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double nol_2e_v0


    File : :file:`no_dressing.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nol_2e_v0	(mo_num,mo_num,mo_num,mo_num)


    two-body part of the normal ordering (on HF reference) of 3-e operator COMPUTED IN A NAIVE WAY
    
    WARNING: if core_tc_op
    
    then does not take into account any contribution from core orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_tc_op`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double tc_2e_3idx_coulomb_integrals


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_2e_3idx_coulomb_integrals	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: tc_2e_3idx_exchange_integrals	(mo_num,mo_num,mo_num)


    tc_2e_3idx_coulomb_integrals (j,k,i) = <jk|ji>
    tc_2e_3idx_exchange_integrals(j,k,i) = <kj|ji>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`


 
.. c:var:: double tc_2e_3idx_coulomb_integrals_transp


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_2e_3idx_coulomb_integrals_transp	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: tc_2e_3idx_exchange_integrals_transp	(mo_num,mo_num,mo_num)


    tc_2e_3idx_coulomb_integrals_transp (j,k,i) = <jk|ji>
    tc_2e_3idx_exchange_integrals_transp(j,k,i) = <kj|ji>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e_transp`
       * :c:data:`mo_num`


 
.. c:var:: double tc_2e_3idx_exchange_integrals


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_2e_3idx_coulomb_integrals	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: tc_2e_3idx_exchange_integrals	(mo_num,mo_num,mo_num)


    tc_2e_3idx_coulomb_integrals (j,k,i) = <jk|ji>
    tc_2e_3idx_exchange_integrals(j,k,i) = <kj|ji>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e`
       * :c:data:`mo_num`


 
.. c:var:: double tc_2e_3idx_exchange_integrals_transp


    File : :file:`total_twoe_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tc_2e_3idx_coulomb_integrals_transp	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: tc_2e_3idx_exchange_integrals_transp	(mo_num,mo_num,mo_num)


    tc_2e_3idx_coulomb_integrals_transp (j,k,i) = <jk|ji>
    tc_2e_3idx_exchange_integrals_transp(j,k,i) = <kj|ji>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e_transp`
       * :c:data:`mo_num`


 
.. c:var:: double three_body_ints_bi_ort


    File : :file:`three_body_ints_bi_ort.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_body_ints_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num,mo_num)


    matrix element of the -L  three-body operator
    
    notice the -1 sign: in this way three_body_ints_bi_ort can be directly used to compute Slater rules :)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao_transp`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_3_idx_cycle_1_bi_ort


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
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_3_idx_cycle_2_bi_ort


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
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_3_idx_direct_bi_ort


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
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_3_idx_exch12_bi_ort


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
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_3_idx_exch12_bi_ort_new


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
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_3_idx_exch13_bi_ort


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
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_3_idx_exch23_bi_ort


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
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_cycle_1_bi_ort


    File : :file:`three_body_ijmk.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_direct_bi_ort (m,j,k,i) = < m j k | -L | m j i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch13_bi_ort (m,j,k,i) = < m j k | -L | i j m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch23_bi_ort (m,j,k,i) = < m j k | -L | j m i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_cycle_1_bi_ort(m,j,k,i) = < m j k | -L | j i m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_4_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    
    three_e_4_idx_direct_bi_ort (m,j,k,i) : Lk Ri Imm Ijj + Lj Rj Imm Iki + Lm Rm Ijj Iki
    three_e_4_idx_exch13_bi_ort (m,j,k,i) : Lk Rm Imi Ijj + Lj Rj Imi Ikm + Lm Ri Ijj Ikm
    three_e_4_idx_exch23_bi_ort (m,j,k,i) : Lk Ri Imj Ijm + Lj Rm Imj Iki + Lm Rj Ijm Iki
    three_e_4_idx_cycle_1_bi_ort(m,j,k,i) : Lk Rm Imj Iji + Lj Ri Imj Ikm + Lm Rj Iji Ikm
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_cycle_1_bi_ort_n4


    File : :file:`three_body_ijmk_n4.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_direct_bi_ort_n4	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch13_bi_ort_n4	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_cycle_1_bi_ort_n4	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_direct_bi_ort_n4 (m,j,k,i) = < m j k | -L | m j i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch13_bi_ort_n4 (m,j,k,i) = < m j k | -L | i j m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch12_bi_ort_n4 (m,j,k,i) = < m j k | -L | m i j > ::: notice that i is the RIGHT MO and k is the LEFT MO
                                             = three_e_4_idx_exch13_bi_ort_n4 (j,m,k,i)
    three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) = < m j k | -L | j i m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_cycle_2_bi_ort_n4(m,j,k,i) = < m j k | -L | i m j > ::: notice that i is the RIGHT MO and k is the LEFT MO
                                             = three_e_4_idx_cycle_1_bi_ort_n4(j,m,k,i)
    
    notice the -1 sign: in this way three_e_4_idx_direct_bi_ort_n4 can be directly used to compute Slater rules with a + sign
    
    three_e_4_idx_direct_bi_ort_n4 (m,j,k,i) : Lk Ri Imm Ijj + Lj Rj Imm Iki + Lm Rm Ijj Iki
    three_e_4_idx_exch13_bi_ort_n4 (m,j,k,i) : Lk Rm Imi Ijj + Lj Rj Imi Ikm + Lm Ri Ijj Ikm
    three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) : Lk Rm Imj Iji + Lj Ri Imj Ikm + Lm Rj Iji Ikm
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_cycle_1_bi_ort_old


    File : :file:`three_body_ijmk_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_cycle_1_bi_ort_old	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_cycle_1_bi_ort_old(m,j,k,i) = <mjk|-L|jim> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_cycle_2_bi_ort_old


    File : :file:`three_body_ijmk_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_cycle_2_bi_ort_old	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_cycle_2_bi_ort_old(m,j,k,i) = <mjk|-L|imj> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_direct_bi_ort


    File : :file:`three_body_ijmk.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_direct_bi_ort (m,j,k,i) = < m j k | -L | m j i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch13_bi_ort (m,j,k,i) = < m j k | -L | i j m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch23_bi_ort (m,j,k,i) = < m j k | -L | j m i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_cycle_1_bi_ort(m,j,k,i) = < m j k | -L | j i m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_4_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    
    three_e_4_idx_direct_bi_ort (m,j,k,i) : Lk Ri Imm Ijj + Lj Rj Imm Iki + Lm Rm Ijj Iki
    three_e_4_idx_exch13_bi_ort (m,j,k,i) : Lk Rm Imi Ijj + Lj Rj Imi Ikm + Lm Ri Ijj Ikm
    three_e_4_idx_exch23_bi_ort (m,j,k,i) : Lk Ri Imj Ijm + Lj Rm Imj Iki + Lm Rj Ijm Iki
    three_e_4_idx_cycle_1_bi_ort(m,j,k,i) : Lk Rm Imj Iji + Lj Ri Imj Ikm + Lm Rj Iji Ikm
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_direct_bi_ort_n4


    File : :file:`three_body_ijmk_n4.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_direct_bi_ort_n4	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch13_bi_ort_n4	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_cycle_1_bi_ort_n4	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_direct_bi_ort_n4 (m,j,k,i) = < m j k | -L | m j i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch13_bi_ort_n4 (m,j,k,i) = < m j k | -L | i j m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch12_bi_ort_n4 (m,j,k,i) = < m j k | -L | m i j > ::: notice that i is the RIGHT MO and k is the LEFT MO
                                             = three_e_4_idx_exch13_bi_ort_n4 (j,m,k,i)
    three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) = < m j k | -L | j i m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_cycle_2_bi_ort_n4(m,j,k,i) = < m j k | -L | i m j > ::: notice that i is the RIGHT MO and k is the LEFT MO
                                             = three_e_4_idx_cycle_1_bi_ort_n4(j,m,k,i)
    
    notice the -1 sign: in this way three_e_4_idx_direct_bi_ort_n4 can be directly used to compute Slater rules with a + sign
    
    three_e_4_idx_direct_bi_ort_n4 (m,j,k,i) : Lk Ri Imm Ijj + Lj Rj Imm Iki + Lm Rm Ijj Iki
    three_e_4_idx_exch13_bi_ort_n4 (m,j,k,i) : Lk Rm Imi Ijj + Lj Rj Imi Ikm + Lm Ri Ijj Ikm
    three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) : Lk Rm Imj Iji + Lj Ri Imj Ikm + Lm Rj Iji Ikm
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_direct_bi_ort_old


    File : :file:`three_body_ijmk_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_direct_bi_ort_old	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_direct_bi_ort_old(m,j,k,i) = <mjk|-L|mji> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_exch12_bi_ort_old


    File : :file:`three_body_ijmk_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_exch12_bi_ort_old	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_exch12_bi_ort_old(m,j,k,i) = <mjk|-L|mij> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_exch13_bi_ort


    File : :file:`three_body_ijmk.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_direct_bi_ort (m,j,k,i) = < m j k | -L | m j i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch13_bi_ort (m,j,k,i) = < m j k | -L | i j m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch23_bi_ort (m,j,k,i) = < m j k | -L | j m i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_cycle_1_bi_ort(m,j,k,i) = < m j k | -L | j i m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_4_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    
    three_e_4_idx_direct_bi_ort (m,j,k,i) : Lk Ri Imm Ijj + Lj Rj Imm Iki + Lm Rm Ijj Iki
    three_e_4_idx_exch13_bi_ort (m,j,k,i) : Lk Rm Imi Ijj + Lj Rj Imi Ikm + Lm Ri Ijj Ikm
    three_e_4_idx_exch23_bi_ort (m,j,k,i) : Lk Ri Imj Ijm + Lj Rm Imj Iki + Lm Rj Ijm Iki
    three_e_4_idx_cycle_1_bi_ort(m,j,k,i) : Lk Rm Imj Iji + Lj Ri Imj Ikm + Lm Rj Iji Ikm
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_exch13_bi_ort_n4


    File : :file:`three_body_ijmk_n4.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_direct_bi_ort_n4	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch13_bi_ort_n4	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_cycle_1_bi_ort_n4	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_direct_bi_ort_n4 (m,j,k,i) = < m j k | -L | m j i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch13_bi_ort_n4 (m,j,k,i) = < m j k | -L | i j m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch12_bi_ort_n4 (m,j,k,i) = < m j k | -L | m i j > ::: notice that i is the RIGHT MO and k is the LEFT MO
                                             = three_e_4_idx_exch13_bi_ort_n4 (j,m,k,i)
    three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) = < m j k | -L | j i m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_cycle_2_bi_ort_n4(m,j,k,i) = < m j k | -L | i m j > ::: notice that i is the RIGHT MO and k is the LEFT MO
                                             = three_e_4_idx_cycle_1_bi_ort_n4(j,m,k,i)
    
    notice the -1 sign: in this way three_e_4_idx_direct_bi_ort_n4 can be directly used to compute Slater rules with a + sign
    
    three_e_4_idx_direct_bi_ort_n4 (m,j,k,i) : Lk Ri Imm Ijj + Lj Rj Imm Iki + Lm Rm Ijj Iki
    three_e_4_idx_exch13_bi_ort_n4 (m,j,k,i) : Lk Rm Imi Ijj + Lj Rj Imi Ikm + Lm Ri Ijj Ikm
    three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) : Lk Rm Imj Iji + Lj Ri Imj Ikm + Lm Rj Iji Ikm
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_exch13_bi_ort_old


    File : :file:`three_body_ijmk_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_exch13_bi_ort_old	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_exch13_bi_ort_old(m,j,k,i) = <mjk|-L|ijm> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_exch23_bi_ort


    File : :file:`three_body_ijmk.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_4_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_direct_bi_ort (m,j,k,i) = < m j k | -L | m j i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch13_bi_ort (m,j,k,i) = < m j k | -L | i j m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_exch23_bi_ort (m,j,k,i) = < m j k | -L | j m i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    three_e_4_idx_cycle_1_bi_ort(m,j,k,i) = < m j k | -L | j i m > ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_4_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    
    three_e_4_idx_direct_bi_ort (m,j,k,i) : Lk Ri Imm Ijj + Lj Rj Imm Iki + Lm Rm Ijj Iki
    three_e_4_idx_exch13_bi_ort (m,j,k,i) : Lk Rm Imi Ijj + Lj Rj Imi Ikm + Lm Ri Ijj Ikm
    three_e_4_idx_exch23_bi_ort (m,j,k,i) : Lk Ri Imj Ijm + Lj Rm Imj Iki + Lm Rj Ijm Iki
    three_e_4_idx_cycle_1_bi_ort(m,j,k,i) : Lk Rm Imj Iji + Lj Ri Imj Ikm + Lm Rj Iji Ikm
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_exch23_bi_ort_n4


    File : :file:`three_body_ijmk_n4.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_exch23_bi_ort_n4	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_exch23_bi_ort_n4 (m,j,k,i) = < m j k | -L | j m i > ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_4_idx_direct_bi_ort_n4 can be directly used to compute Slater rules with a + sign
    
    three_e_4_idx_exch23_bi_ort_n4 (m,j,k,i) : Lk Ri Imj Ijm + Lj Rm Imj Iki + Lm Rj Ijm Iki
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_4_idx_exch23_bi_ort_old


    File : :file:`three_body_ijmk_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_4_idx_exch23_bi_ort_old	(mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_4_idx_exch23_bi_ort_old(m,j,k,i) = <mjk|-L|jmi> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_cycle_1_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_cycle_2_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_direct_bi_ort(m,l,j,k,i) = <mlk|-L|mji> :: : notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_cycle_1_bi_ort_old


    File : :file:`three_body_ijmkl_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_cycle_1_bi_ort_old	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_cycle_1_bi_ort_old(m,l,j,k,i) = <mlk|-L|jim> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_cycle_2_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_cycle_2_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_direct_bi_ort(m,l,j,k,i) = <mlk|-L|mji> :: : notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_cycle_2_bi_ort_old


    File : :file:`three_body_ijmkl_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_cycle_2_bi_ort_old	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_cycle_2_bi_ort_old(m,l,j,k,i) = <mlk|-L|imj> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_direct_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_cycle_2_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_direct_bi_ort(m,l,j,k,i) = <mlk|-L|mji> :: : notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_direct_bi_ort_old


    File : :file:`three_body_ijmkl_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_direct_bi_ort_old	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_direct_bi_ort_old(m,l,j,k,i) = <mlk|-L|mji> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_exch12_bi_ort_old


    File : :file:`three_body_ijmkl_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_exch12_bi_ort_old	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_exch12_bi_ort_old(m,l,j,k,i) = <mlk|-L|mij> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_exch13_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_cycle_2_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_direct_bi_ort(m,l,j,k,i) = <mlk|-L|mji> :: : notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_exch13_bi_ort_old


    File : :file:`three_body_ijmkl_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_exch13_bi_ort_old	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_exch13_bi_ort_old(m,l,j,k,i) = <mlk|-L|ijm> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_exch23_bi_ort


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_direct_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_exch23_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_exch13_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_cycle_1_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)
        double precision, allocatable	:: three_e_5_idx_cycle_2_bi_ort	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_direct_bi_ort(m,l,j,k,i) = <mlk|-L|mji> :: : notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double three_e_5_idx_exch23_bi_ort_old


    File : :file:`three_body_ijmkl_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: three_e_5_idx_exch23_bi_ort_old	(mo_num,mo_num,mo_num,mo_num,mo_num)


    
    matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
    
    three_e_5_idx_exch23_bi_ort_old(m,l,j,k,i) = <mlk|-L|jmi> ::: notice that i is the RIGHT MO and k is the LEFT MO
    
    notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double x_w_ki_bi_ortho_erf_rk


    File : :file:`semi_num_ints_mo_old.irp.f`

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


 
.. c:var:: double x_w_ki_bi_ortho_erf_rk_diag


    File : :file:`semi_num_ints_mo_old.irp.f`

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
 
.. c:function:: void bi_ortho_mo_coul_ints(void)


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

 
.. c:function:: void bi_ortho_mo_ints(void)


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

 
.. c:function:: void give_integrals_3_body_bi_ort(void)


    File : :file:`three_body_ints_bi_ort.irp.f`

    .. code:: fortran

        subroutine give_integrals_3_body_bi_ort(n, l, k, m, j, i, integral)


    
    < n l k | L | m j i > with a BI-ORTHONORMAL MOLECULAR ORBITALS
    
    /!\ L is defined without the 1/6 factor
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`nol_0e_v0`
       * :c:data:`nol_1e_v0`
       * :c:data:`nol_2e_v0`
       * :c:func:`test_3e`
       * :c:data:`three_body_ints_bi_ort`
       * :c:data:`three_e_3_idx_cycle_1_bi_ort`
       * :c:data:`three_e_3_idx_cycle_2_bi_ort`
       * :c:data:`three_e_3_idx_direct_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort`
       * :c:data:`three_e_3_idx_exch12_bi_ort_new`
       * :c:data:`three_e_3_idx_exch13_bi_ort`
       * :c:data:`three_e_3_idx_exch23_bi_ort`
       * :c:data:`three_e_4_idx_cycle_1_bi_ort_old`
       * :c:data:`three_e_4_idx_cycle_2_bi_ort_old`
       * :c:data:`three_e_4_idx_direct_bi_ort_old`
       * :c:data:`three_e_4_idx_exch12_bi_ort_old`
       * :c:data:`three_e_4_idx_exch13_bi_ort_old`
       * :c:data:`three_e_4_idx_exch23_bi_ort_old`
       * :c:data:`three_e_5_idx_cycle_1_bi_ort_old`
       * :c:data:`three_e_5_idx_cycle_2_bi_ort_old`
       * :c:data:`three_e_5_idx_direct_bi_ort_old`
       * :c:data:`three_e_5_idx_exch12_bi_ort_old`
       * :c:data:`three_e_5_idx_exch13_bi_ort_old`
       * :c:data:`three_e_5_idx_exch23_bi_ort_old`

 
.. c:function:: void give_integrals_3_body_bi_ort_ao(void)


    File : :file:`three_body_ints_bi_ort.irp.f`

    .. code:: fortran

        subroutine give_integrals_3_body_bi_ort_ao(n, l, k, m, j, i, integral)


    
    < n l k | L | m j i > with a BI-ORTHONORMAL ATOMIC ORBITALS
    
    /!\ L is defined without the 1/6 factor
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`aos_in_r_array_transp`
       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_ao_t`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void give_integrals_3_body_bi_ort_old(void)


    File : :file:`three_body_ints_bi_ort.irp.f`

    .. code:: fortran

        subroutine give_integrals_3_body_bi_ort_old(n, l, k, m, j, i, integral)


    
    < n l k | L | m j i > with a BI-ORTHONORMAL MOLECULAR ORBITALS
    
    /!\ L is defined without the 1/6 factor
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_transp`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_3e`

 
.. c:function:: void give_integrals_3_body_bi_ort_spin(void)


    File : :file:`three_body_ints_bi_ort.irp.f`

    .. code:: fortran

        subroutine give_integrals_3_body_bi_ort_spin( n, sigma_n, l, sigma_l, k, sigma_k &
                                                    , m, sigma_m, j, sigma_j, i, sigma_i &
                                                    , integral)


    
    < n l k | L | m j i > with a BI-ORTHONORMAL SPIN-ORBITALS
    
    /!\ L is defined without the 1/6 factor
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`nol_0e_naive`
       * :c:data:`nol_1e_naive`
       * :c:data:`nol_2e_naive`

 
.. c:function:: void test_3e(void)


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_3e



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_integrals_3_body_bi_ort`
       * :c:func:`give_integrals_3_body_bi_ort_old`

 
.. c:function:: void test_4idx(void)


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_4idx()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_cycle_1_bi_ort_old`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_direct_bi_ort_old`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_exch13_bi_ort_old`
       * :c:data:`three_e_4_idx_direct_bi_ort`
       * :c:data:`three_e_4_idx_exch23_bi_ort_old`

 
.. c:function:: void test_4idx2(void)


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_4idx2()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`three_e_4_idx_direct_bi_ort`

 
.. c:function:: void test_4idx_n4(void)


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_4idx_n4()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`three_e_4_idx_direct_bi_ort_n4`
       * :c:data:`three_e_4_idx_cycle_1_bi_ort_old`
       * :c:data:`three_e_4_idx_direct_bi_ort_n4`
       * :c:data:`three_e_4_idx_direct_bi_ort_old`
       * :c:data:`three_e_4_idx_direct_bi_ort_n4`
       * :c:data:`three_e_4_idx_exch13_bi_ort_old`
       * :c:data:`three_e_4_idx_exch23_bi_ort_n4`
       * :c:data:`three_e_4_idx_exch23_bi_ort_old`

 
.. c:function:: void test_5idx(void)


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_5idx



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`three_e_5_idx_direct_bi_ort`
       * :c:data:`three_e_5_idx_direct_bi_ort_old`

 
.. c:function:: void test_5idx2(void)


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_5idx2



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`three_e_5_idx_direct_bi_ort`

 
.. c:function:: void test_int2_grad1_u12_bimo_t(void)


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_int2_grad1_u12_bimo_t()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`int2_grad1_u12_bimo_t`
       * :c:data:`int2_grad1_u12_bimo_t_old`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void test_mos_in_r(void)


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_mos_in_r()



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_l_in_r_array_transp_old`
       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`mos_r_in_r_array_transp_old`
       * :c:data:`n_points_final_grid`

 
.. c:function:: void test_nol_fc(void)


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_nol_fc



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`nol_0e`
       * :c:data:`nol_0e_v0`
       * :c:data:`nol_1e`
       * :c:data:`nol_1e_v0`
       * :c:data:`nol_2e`
       * :c:data:`nol_2e_v0`

 
.. c:function:: void test_total_fc(void)


    File : :file:`bi_ort_ints.irp.f`

    .. code:: fortran

        subroutine test_total_fc



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_tc_two_e_chemist`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`bi_ort_ints`

 
.. c:function:: void three_e_5_idx_exch12_bi_ort(void)


    File : :file:`three_body_ijmkl.irp.f`

    .. code:: fortran

        double precision function three_e_5_idx_exch12_bi_ort(m,l,i,k,j) result(integral)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`three_e_5_idx_direct_bi_ort`

