.. _module_bi_ortho_mos: 
 
.. program:: bi_ortho_mos 
 
.. default-role:: option 
 
bi_ortho_mos
============

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Providers 
--------- 
 
.. c:var:: angle_left_right


    File : :file:`bi_ortho_mos/overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: angle_left_right	(mo_num)
        double precision	:: max_angle_left_right	


    angle_left_right(i) = angle between the left-eigenvector chi_i and the right-eigenvector phi_i

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`overlap_bi_ortho`
       * :c:data:`overlap_mo_r`


 
.. c:var:: max_angle_left_right


    File : :file:`bi_ortho_mos/overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: angle_left_right	(mo_num)
        double precision	:: max_angle_left_right	


    angle_left_right(i) = angle between the left-eigenvector chi_i and the right-eigenvector phi_i

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`overlap_bi_ortho`
       * :c:data:`overlap_mo_r`


 
.. c:var:: mo_l_coef


    File : :file:`bi_ortho_mos/mos_rl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_l_coef	(ao_num,mo_num)


    
    Molecular left-orbital coefficients on |AO| basis set
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_vartc_tot`
       * :c:data:`diag_three_elem_hf`
       * :c:data:`fock_3e_uhf_mo_a`
       * :c:data:`fock_3e_uhf_mo_b`
       * :c:data:`fock_3e_uhf_mo_cs`
       * :c:data:`fock_a_tmp1_bi_ortho`
       * :c:data:`fock_a_tmp2_bi_ortho`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tmp1_bi_ortho`
       * :c:data:`fock_b_tmp2_bi_ortho`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`fock_cs_3e_bi_orth`
       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`fock_matrix_tc_mo_alpha`
       * :c:data:`fock_matrix_tc_mo_beta`
       * :c:data:`fock_matrix_vartc_ao_tot`
       * :c:data:`fock_matrix_vartc_mo_alpha`
       * :c:data:`fock_matrix_vartc_mo_beta`
       * :c:data:`fock_tc_reigvec_ao`
       * :c:data:`fqs_sqf_mo`
       * :c:data:`int2_grad1_u12_bimo_transp`
       * :c:data:`mo_bi_orth_bipole_x`
       * :c:data:`mo_bi_ortho_one_e`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_l_coef_transp`
       * :c:data:`mo_v_ki_bi_ortho_erf_rk_cst_mu`
       * :c:data:`mo_x_v_ki_bi_ortho_erf_rk_cst_mu`
       * :c:data:`overlap_bi_ortho`
       * :c:data:`overlap_mo_r`
       * :c:data:`overlap_mo_r_mo`
       * :c:data:`q_alpha`
       * :c:data:`q_beta`
       * :c:data:`tc_hf_energy`
       * :c:data:`tcscf_bi_ort_dm_ao`
       * :c:data:`tcscf_bi_ort_dm_ao_alpha`
       * :c:data:`tcscf_bi_ort_dm_ao_beta`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`

 
.. c:var:: mo_l_coef_transp


    File : :file:`bi_ortho_mos/mos_rl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_l_coef_transp	(mo_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_coul_e_chemist`
       * :c:data:`mo_bi_ortho_tc_two_e_chemist`
       * :c:data:`mos_l_in_r_array`

 
.. c:var:: mo_r_coef


    File : :file:`bi_ortho_mos/mos_rl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_r_coef	(ao_num,mo_num)


    
    Molecular right-orbital coefficients on |AO| basis set
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ezfio_filename`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_vartc_tot`
       * :c:data:`diag_three_elem_hf`
       * :c:data:`fock_3e_uhf_mo_a`
       * :c:data:`fock_3e_uhf_mo_b`
       * :c:data:`fock_3e_uhf_mo_cs`
       * :c:data:`fock_a_tmp1_bi_ortho`
       * :c:data:`fock_a_tmp2_bi_ortho`
       * :c:data:`fock_a_tot_3e_bi_orth`
       * :c:data:`fock_b_tmp1_bi_ortho`
       * :c:data:`fock_b_tmp2_bi_ortho`
       * :c:data:`fock_b_tot_3e_bi_orth`
       * :c:data:`fock_cs_3e_bi_orth`
       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`fock_matrix_tc_mo_alpha`
       * :c:data:`fock_matrix_tc_mo_beta`
       * :c:data:`fock_matrix_vartc_ao_tot`
       * :c:data:`fock_matrix_vartc_mo_alpha`
       * :c:data:`fock_matrix_vartc_mo_beta`
       * :c:data:`fock_tc_reigvec_ao`
       * :c:data:`fock_vartc_eigvec_ao`
       * :c:data:`fock_vartc_eigvec_mo`
       * :c:data:`fqs_sqf_mo`
       * :c:data:`int2_grad1_u12_bimo_transp`
       * :c:data:`mo_bi_orth_bipole_x`
       * :c:data:`mo_bi_ortho_one_e`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_r_coef_transp`
       * :c:data:`mo_v_ki_bi_ortho_erf_rk_cst_mu`
       * :c:data:`mo_x_v_ki_bi_ortho_erf_rk_cst_mu`
       * :c:data:`overlap_bi_ortho`
       * :c:data:`overlap_mo_r`
       * :c:data:`overlap_mo_r_mo`
       * :c:data:`q_alpha`
       * :c:data:`q_beta`
       * :c:data:`tc_hf_energy`
       * :c:data:`tcscf_bi_ort_dm_ao`
       * :c:data:`tcscf_bi_ort_dm_ao_alpha`
       * :c:data:`tcscf_bi_ort_dm_ao_beta`
       * :c:data:`tcscf_density_matrix_ao_alpha`
       * :c:data:`tcscf_density_matrix_ao_beta`
       * :c:data:`vartc_hf_energy`

 
.. c:var:: mo_r_coef_transp


    File : :file:`bi_ortho_mos/mos_rl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_r_coef_transp	(mo_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_bi_ortho_coul_e_chemist`
       * :c:data:`mo_bi_ortho_tc_two_e_chemist`
       * :c:data:`mos_l_grad_in_r_array`
       * :c:data:`mos_r_grad_in_r_array`
       * :c:data:`mos_r_in_r_array`

 
.. c:var:: mos_l_grad_in_r_array


    File : :file:`bi_ortho_mos/grad_bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_l_grad_in_r_array	(mo_num,n_points_final_grid,3)


    mos_l_grad_in_r_array(i,j,k)          = value of the kth component of the gradient of ith RIGHT mo on the jth grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef_transp`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mos_l_grad_in_r_array_transp`
       * :c:data:`mos_l_grad_in_r_array_transp_bis`

 
.. c:var:: mos_l_grad_in_r_array_transp


    File : :file:`bi_ortho_mos/grad_bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_l_grad_in_r_array_transp	(3,mo_num,n_points_final_grid)


    mos_l_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth RIGHT mo on the ith grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_l_grad_in_r_array`
       * :c:data:`n_points_final_grid`


 
.. c:var:: mos_l_grad_in_r_array_transp_bis


    File : :file:`bi_ortho_mos/grad_bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_l_grad_in_r_array_transp_bis	(3,n_points_final_grid,mo_num)


    mos_l_grad_in_r_array_transp(i,j,k)   = value of the ith component of the gradient on the jth grid point of jth RIGHT MO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_l_grad_in_r_array`
       * :c:data:`n_points_final_grid`


 
.. c:var:: mos_l_in_r_array


    File : :file:`bi_ortho_mos/bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_l_in_r_array	(mo_num,n_points_final_grid)


    mos_in_r_array(i,j) = value of the ith LEFT mo on the jth grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`mo_l_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mos_l_in_r_array_transp`
       * :c:data:`tc_scf_dm_in_r`

 
.. c:var:: mos_l_in_r_array_transp


    File : :file:`bi_ortho_mos/bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_l_in_r_array_transp	(n_points_final_grid,mo_num)


    mos_l_in_r_array_transp(i,j) = value of the jth mo on the ith grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_l_in_r_array`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3e_uhf_mo_a`
       * :c:data:`fock_3e_uhf_mo_b`
       * :c:data:`fock_3e_uhf_mo_cs`
       * :c:data:`fock_a_tmp1_bi_ortho`
       * :c:data:`fock_b_tmp1_bi_ortho`
       * :c:data:`fock_cs_3e_bi_orth`
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
       * :c:data:`w1_diag_in_r`
       * :c:data:`w1_l_in_r`

 
.. c:var:: mos_r_grad_in_r_array


    File : :file:`bi_ortho_mos/grad_bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_r_grad_in_r_array	(mo_num,n_points_final_grid,3)


    mos_r_grad_in_r_array(i,j,k)          = value of the kth component of the gradient of ith RIGHT mo on the jth grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef_transp`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mos_r_grad_in_r_array_transp`
       * :c:data:`mos_r_grad_in_r_array_transp_bis`

 
.. c:var:: mos_r_grad_in_r_array_transp


    File : :file:`bi_ortho_mos/grad_bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_r_grad_in_r_array_transp	(3,mo_num,n_points_final_grid)


    mos_r_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth RIGHT mo on the ith grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_r_grad_in_r_array`
       * :c:data:`n_points_final_grid`


 
.. c:var:: mos_r_grad_in_r_array_transp_bis


    File : :file:`bi_ortho_mos/grad_bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_r_grad_in_r_array_transp_bis	(3,n_points_final_grid,mo_num)


    mos_r_grad_in_r_array_transp(i,j,k)   = value of the ith component of the gradient on the jth grid point of jth RIGHT MO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_r_grad_in_r_array`
       * :c:data:`n_points_final_grid`


 
.. c:var:: mos_r_in_r_array


    File : :file:`bi_ortho_mos/bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_r_in_r_array	(mo_num,n_points_final_grid)


    mos_in_r_array(i,j) = value of the ith RIGHT mo on the jth grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef_transp`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mos_r_in_r_array_transp`
       * :c:data:`tc_scf_dm_in_r`

 
.. c:var:: mos_r_in_r_array_transp


    File : :file:`bi_ortho_mos/bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_r_in_r_array_transp	(n_points_final_grid,mo_num)


    mos_r_in_r_array_transp(i,j) = value of the jth mo on the ith grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_r_in_r_array`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3e_uhf_mo_a`
       * :c:data:`fock_3e_uhf_mo_b`
       * :c:data:`fock_3e_uhf_mo_cs`
       * :c:data:`fock_a_tmp1_bi_ortho`
       * :c:data:`fock_b_tmp1_bi_ortho`
       * :c:data:`fock_cs_3e_bi_orth`
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
       * :c:data:`w1_diag_in_r`
       * :c:data:`w1_in_r`
       * :c:data:`w1_r_in_r`

 
.. c:var:: overlap_bi_ortho


    File : :file:`bi_ortho_mos/overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: overlap_bi_ortho	(mo_num,mo_num)
        double precision, allocatable	:: overlap_diag_bi_ortho	(mo_num)


    Overlap matrix between the RIGHT and LEFT MOs. Should be the identity matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`angle_left_right`

 
.. c:var:: overlap_diag_bi_ortho


    File : :file:`bi_ortho_mos/overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: overlap_bi_ortho	(mo_num,mo_num)
        double precision, allocatable	:: overlap_diag_bi_ortho	(mo_num)


    Overlap matrix between the RIGHT and LEFT MOs. Should be the identity matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`angle_left_right`

 
.. c:var:: overlap_mo_l


    File : :file:`bi_ortho_mos/overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: overlap_mo_r	(mo_num,mo_num)
        double precision, allocatable	:: overlap_mo_l	(mo_num,mo_num)


    overlap_mo_r_mo(j,i) = <MO_i|MO_R_j>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`angle_left_right`

 
.. c:var:: overlap_mo_l_mo


    File : :file:`bi_ortho_mos/overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: overlap_mo_r_mo	(mo_num,mo_num)
        double precision, allocatable	:: overlap_mo_l_mo	(mo_num,mo_num)


    overlap_mo_r_mo(j,i) = <MO_j|MO_R_i>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: overlap_mo_r


    File : :file:`bi_ortho_mos/overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: overlap_mo_r	(mo_num,mo_num)
        double precision, allocatable	:: overlap_mo_l	(mo_num,mo_num)


    overlap_mo_r_mo(j,i) = <MO_i|MO_R_j>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`angle_left_right`

 
.. c:var:: overlap_mo_r_mo


    File : :file:`bi_ortho_mos/overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: overlap_mo_r_mo	(mo_num,mo_num)
        double precision, allocatable	:: overlap_mo_l_mo	(mo_num,mo_num)


    overlap_mo_r_mo(j,i) = <MO_j|MO_R_i>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`


 
.. c:var:: tcscf_bi_ort_dm_ao


    File : :file:`bi_ortho_mos/bi_density.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tcscf_bi_ort_dm_ao	(ao_num,ao_num)


    TCSCF_bi_ort_dm_ao(i,j) = <Chi_0| a^dagger_i,beta+alpha a_j,beta+alpha |Phi_0> where i,j are AO basis.
    
    This is the equivalent of the total electronic density of the HF Slater determinant, but with a couple of bi-orthonormal Slater determinant |Chi_0> and |Phi_0>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`
       * :c:data:`tcscf_bi_ort_dm_ao_alpha`
       * :c:data:`tcscf_bi_ort_dm_ao_beta`


 
.. c:var:: tcscf_bi_ort_dm_ao_alpha


    File : :file:`bi_ortho_mos/bi_density.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tcscf_bi_ort_dm_ao_alpha	(ao_num,ao_num)


    TCSCF_bi_ort_dm_ao_alpha(i,j) = <Chi_0| a^dagger_i,alpha a_j,alpha |Phi_0> where i,j are AO basis.
    
    This is the equivalent of the alpha density of the HF Slater determinant, but with a couple of bi-orthonormal Slater determinant |Chi_0> and |Phi_0>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3e_uhf_ao_a`
       * :c:data:`fock_3e_uhf_ao_b`
       * :c:data:`tcscf_bi_ort_dm_ao`
       * :c:data:`tcscf_density_matrix_ao_alpha`

 
.. c:var:: tcscf_bi_ort_dm_ao_beta


    File : :file:`bi_ortho_mos/bi_density.irp.f`

    .. code:: fortran

        double precision, allocatable	:: tcscf_bi_ort_dm_ao_beta	(ao_num,ao_num)


    TCSCF_bi_ort_dm_ao_beta(i,j) = <Chi_0| a^dagger_i,beta a_j,beta |Phi_0> where i,j are AO basis.
    
    This is the equivalent of the beta density of the HF Slater determinant, but with a couple of bi-orthonormal Slater determinant |Chi_0> and |Phi_0>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_r_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_3e_uhf_ao_a`
       * :c:data:`fock_3e_uhf_ao_b`
       * :c:data:`tcscf_bi_ort_dm_ao`
       * :c:data:`tcscf_density_matrix_ao_beta`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: ao_to_mo_bi_ortho:


    File : :file:`bi_ortho_mos/mos_rl.irp.f`

    .. code:: fortran

        subroutine ao_to_mo_bi_ortho(A_ao, LDA_ao, A_mo, LDA_mo)


    
    Transform A from the |AO| basis to the BI ORTHONORMAL MOS
    
    $C_L^\dagger.A_{ao}.C_R$ where C_L and C_R are the LEFT and RIGHT MO coefs
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_mo_alpha`
       * :c:data:`fock_matrix_tc_mo_beta`
       * :c:data:`fock_matrix_vartc_mo_alpha`
       * :c:data:`fock_matrix_vartc_mo_beta`
       * :c:data:`fqs_sqf_mo`
       * :c:data:`int2_grad1_u12_bimo_transp`
       * :c:data:`mo_bi_orth_bipole_x`
       * :c:data:`mo_bi_ortho_one_e`
       * :c:data:`mo_bi_ortho_tc_one_e`
       * :c:data:`mo_v_ki_bi_ortho_erf_rk_cst_mu`
       * :c:data:`mo_x_v_ki_bi_ortho_erf_rk_cst_mu`
       * :c:func:`rh_tcscf_diis`
       * :c:func:`test_fock_3e_uhf_ao`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: give_all_mos_l_at_r:


    File : :file:`bi_ortho_mos/bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        subroutine give_all_mos_l_at_r(r, mos_l_array)


    mos_l_array(i) = ith LEFT MO function evaluated at "r"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_l_coef_transp`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`mos_l_in_r_array`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`give_all_aos_at_r`

 
.. c:function:: give_all_mos_r_at_r:


    File : :file:`bi_ortho_mos/bi_ort_mos_in_r.irp.f`

    .. code:: fortran

        subroutine give_all_mos_r_at_r(r, mos_r_array)


    mos_r_array(i) = ith RIGHT MO function evaluated at "r"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef_transp`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`mos_r_in_r_array`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`give_all_aos_at_r`

 
.. c:function:: mo_to_ao_bi_ortho:


    File : :file:`bi_ortho_mos/mos_rl.irp.f`

    .. code:: fortran

        subroutine mo_to_ao_bi_ortho(A_mo, LDA_mo, A_ao, LDA_ao)


    
    mo_l_coef.T x     A_ao   x mo_r_coef = A_mo
    mo_l_coef.T x ao_overlap x mo_r_coef =  I
    
    ==> A_ao = (ao_overlap x mo_r_coef) x A_mo x (ao_overlap x mo_l_coef).T
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_l_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_r_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_tc_ao_tot`
       * :c:data:`fock_matrix_vartc_ao_tot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

