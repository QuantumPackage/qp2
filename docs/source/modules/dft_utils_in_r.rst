.. _module_dft_utils_in_r: 
 
.. program:: dft_utils_in_r 
 
.. default-role:: option 
 
==============
dft_utils_in_r
==============

This module contains most of the fundamental quantities (AOs, MOs or density derivatives) evaluated in real-space representation that are needed for the various DFT modules.

As these quantities might be used and re-used, the values at each point of the grid are stored (see ``becke_numerical_grid`` for more information on the grid).

The main providers for this module are:

* `aos_in_r_array`: values of the |AO| basis on the grid point.
* `mos_in_r_array`: values of the |MO| basis on the grid point.
* `one_e_dm_and_grad_alpha_in_r`: values of the density and its gradienst on the grid points.

 
 
 
Providers 
--------- 
 
.. c:var:: alpha_dens_kin_in_r


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: alpha_dens_kin_in_r	(n_points_final_grid)
        double precision, allocatable	:: beta_dens_kin_in_r	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mos_grad_in_r_array_tranp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: ao_abs_int_grid


    File : :file:`dft_utils_in_r/ao_prod_mlti_pl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_abs_int_grid	(ao_num)


    ao_abs_int_grid(i) = \int dr |phi_i(r) |

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`


 
.. c:var:: ao_overlap_abs_grid


    File : :file:`dft_utils_in_r/ao_prod_mlti_pl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_abs_grid	(ao_num,ao_num)


    ao_overlap_abs_grid(j,i) = \int dr |phi_i(r) phi_j(r)|

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_prod_center`
       * :c:data:`ao_prod_sigma`

 
.. c:var:: ao_prod_abs_r


    File : :file:`dft_utils_in_r/ao_prod_mlti_pl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_prod_abs_r	(ao_num,ao_num)


    ao_prod_abs_r(i,j) = \int |phi_i(r) phi_j(r)| dsqrt((x - <|i|x|j|>)^2 + (y - <|i|y|j|>)^2 +(z - <|i|z|j|>)^2) / \int |phi_i(r) phi_j(r)|
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_prod_center`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_prod_sigma`

 
.. c:var:: ao_prod_center


    File : :file:`dft_utils_in_r/ao_prod_mlti_pl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_prod_center	(3,ao_num,ao_num)


    ao_prod_center(1:3,j,i) = \int dr |phi_i(r) phi_j(r)| x/y/z / \int |phi_i(r) phi_j(r)|
    
    if \int |phi_i(r) phi_j(r)| < 1.d-10 then ao_prod_center = 10000.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs_grid`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_prod_abs_r`
       * :c:data:`ao_prod_dist_grid`

 
.. c:var:: ao_prod_dist_grid


    File : :file:`dft_utils_in_r/ao_prod_mlti_pl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_prod_dist_grid	(ao_num,ao_num,n_points_final_grid)


    ao_prod_dist_grid(j,i,ipoint) = distance between the center of |phi_i(r) phi_j(r)| and the grid point r(ipoint)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_prod_center`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`


 
.. c:var:: ao_prod_sigma


    File : :file:`dft_utils_in_r/ao_prod_mlti_pl.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_prod_sigma	(ao_num,ao_num)


    Gaussian exponent reproducing the product |chi_i(r) chi_j(r)|
    
    Therefore |chi_i(r) chi_j(r)|  \approx e^{-ao_prod_sigma(j,i) (r - ao_prod_center(1:3,j,i))**2}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs_grid`
       * :c:data:`ao_prod_abs_r`


 
.. c:var:: aos_grad_in_r_array


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_grad_in_r_array	(ao_num,n_points_final_grid,3)


    
    aos_grad_in_r_array(i,j,k) = value of the kth component of the gradient of ith ao on the jth grid point
    
    k = 1 : x, k= 2, y, k  3, z
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power_ordered_transp_per_nucl`
       * :c:data:`ao_prim_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_grad_in_r_array_transp_3`
       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`mos_grad_in_r_array`

 
.. c:var:: aos_grad_in_r_array_extra


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_grad_in_r_array_extra	(ao_num,n_points_extra_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power_ordered_transp_per_nucl`
       * :c:data:`ao_prim_num`
       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`


 
.. c:var:: aos_grad_in_r_array_transp


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_grad_in_r_array_transp	(3,ao_num,n_points_final_grid)


    aos_grad_in_r_array_transp(k,i,j)   = value of the kth component of the gradient of jth ao on the ith grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`

 
.. c:var:: aos_grad_in_r_array_transp_3


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_grad_in_r_array_transp_3	(3,n_points_final_grid,ao_num)


    Transposed gradients
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`n_points_final_grid`


 
.. c:var:: aos_grad_in_r_array_transp_bis


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_grad_in_r_array_transp_bis	(n_points_final_grid,ao_num,3)


    Transposed gradients
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`n_points_final_grid`


 
.. c:var:: aos_in_r_array


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_in_r_array	(ao_num,n_points_final_grid)


    aos_in_r_array(i,j) = value of the ith ao on the jth grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power_ordered_transp_per_nucl`
       * :c:data:`ao_prim_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_abs_int_grid`
       * :c:data:`ao_overlap_abs_grid`
       * :c:data:`ao_prod_abs_r`
       * :c:data:`ao_prod_center`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`f_hf_cholesky_sparse_bis`
       * :c:data:`pot_scal_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_sr_pbe`
       * :c:data:`potential_c_alpha_ao_lda`
       * :c:data:`potential_c_alpha_ao_sr_lda`
       * :c:data:`potential_x_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_sr_lda`
       * :c:data:`potential_xc_alpha_ao_lda`
       * :c:data:`potential_xc_alpha_ao_sr_lda`

 
.. c:var:: aos_in_r_array_extra


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_in_r_array_extra	(ao_num,n_points_extra_final_grid)


    aos_in_r_array_extra(i,j)        = value of the ith ao on the jth grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power_ordered_transp_per_nucl`
       * :c:data:`ao_prim_num`
       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_in_r_array_extra_transp`

 
.. c:var:: aos_in_r_array_extra_transp


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_in_r_array_extra_transp	(n_points_extra_final_grid,ao_num)


    aos_in_r_array_extra_transp(i,j) = value of the jth ao on the ith grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`n_points_extra_final_grid`


 
.. c:var:: aos_in_r_array_transp


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_in_r_array_transp	(n_points_final_grid,ao_num)


    aos_in_r_array_transp(i,j) = value of the jth ao on the ith grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_grad_xc_alpha_ao_pbe`
       * :c:data:`pot_grad_xc_alpha_ao_sr_pbe`

 
.. c:var:: aos_lapl_in_r_array


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_lapl_in_r_array	(3,ao_num,n_points_final_grid)


    aos_lapl_in_r_array(i,j,k)   = value of the kth component of the laplacian of jth ao on the ith grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power_ordered_transp_per_nucl`
       * :c:data:`ao_prim_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_lapl_in_r_array_transp`

 
.. c:var:: aos_lapl_in_r_array_transp


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_lapl_in_r_array_transp	(ao_num,n_points_final_grid,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_lapl_in_r_array`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mos_lapl_in_r_array`

 
.. c:var:: beta_dens_kin_in_r


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: alpha_dens_kin_in_r	(n_points_final_grid)
        double precision, allocatable	:: beta_dens_kin_in_r	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mos_grad_in_r_array_tranp`
       * :c:data:`n_points_final_grid`


 
.. c:var:: elec_alpha_num_grid_becke


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: elec_beta_num_grid_becke	(N_states)
        double precision, allocatable	:: elec_alpha_num_grid_becke	(N_states)
        double precision, allocatable	:: elec_num_grid_becke	(N_states)


    number of electrons when the one-e alpha/beta densities are numerically integrated on the DFT grid
    
    !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_average_prov`

 
.. c:var:: elec_beta_num_grid_becke


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: elec_beta_num_grid_becke	(N_states)
        double precision, allocatable	:: elec_alpha_num_grid_becke	(N_states)
        double precision, allocatable	:: elec_num_grid_becke	(N_states)


    number of electrons when the one-e alpha/beta densities are numerically integrated on the DFT grid
    
    !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_average_prov`

 
.. c:var:: elec_num_grid_becke


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: elec_beta_num_grid_becke	(N_states)
        double precision, allocatable	:: elec_alpha_num_grid_becke	(N_states)
        double precision, allocatable	:: elec_num_grid_becke	(N_states)


    number of electrons when the one-e alpha/beta densities are numerically integrated on the DFT grid
    
    !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_average_prov`

 
.. c:var:: kinetic_density_generalized


    File : :file:`dft_utils_in_r/kin_dens.irp.f`

    .. code:: fortran

        double precision, allocatable	:: kinetic_density_generalized	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_grad_in_r_array_tranp`
       * :c:data:`n_points_final_grid`
       * :c:data:`one_e_dm_mo_for_dft`


 
.. c:var:: mo_grad_ints


    File : :file:`dft_utils_in_r/ints_grad.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_grad_ints	(mo_num,mo_num,3)


    mo_grad_ints(i,j,m) = <phi_i^MO | d/dx | phi_j^MO>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mo_num`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`mos_in_r_array`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_grad_ints_transp`

 
.. c:var:: mo_grad_ints_transp


    File : :file:`dft_utils_in_r/ints_grad.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_grad_ints_transp	(3,mo_num,mo_num)


    mo_grad_ints(i,j,m) = <phi_i^MO | d/dx | phi_j^MO>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_grad_ints`
       * :c:data:`mo_num`


 
.. c:var:: mos_grad_in_r_array


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_grad_in_r_array	(mo_num,n_points_final_grid,3)


    mos_grad_in_r_array(i,j,k)          = value of the kth component of the gradient of ith mo on the jth grid point
    
    mos_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth mo on the ith grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_inact_act_mos_grad_in_r_array`
       * :c:data:`mo_grad_ints`
       * :c:data:`mos_grad_in_r_array_tranp`
       * :c:data:`mos_grad_in_r_array_transp_3`
       * :c:data:`mos_grad_in_r_array_transp_bis`

 
.. c:var:: mos_grad_in_r_array_tranp


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_grad_in_r_array_tranp	(3,mo_num,n_points_final_grid)


    mos_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth mo on the ith grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`alpha_dens_kin_in_r`
       * :c:data:`kinetic_density_generalized`

 
.. c:var:: mos_grad_in_r_array_transp_3


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_grad_in_r_array_transp_3	(3,n_points_final_grid,mo_num)


    Transposed gradients
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`n_points_final_grid`


 
.. c:var:: mos_grad_in_r_array_transp_bis


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_grad_in_r_array_transp_bis	(n_points_final_grid,mo_num,3)


    Transposed gradients
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`n_points_final_grid`


 
.. c:var:: mos_in_r_array


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_in_r_array	(mo_num,n_points_final_grid)


    mos_in_r_array(i,j)        = value of the ith mo on the jth grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`basis_mos_in_r_array`
       * :c:data:`mo_grad_ints`

 
.. c:var:: mos_in_r_array_omp


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_in_r_array_omp	(mo_num,n_points_final_grid)


    mos_in_r_array(i,j)        = value of the ith mo on the jth grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_hf_cholesky_sparse`
       * :c:data:`f_hf_cholesky_sparse_bis`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`mos_times_cholesky_r1`
       * :c:data:`mos_times_cholesky_r2`
       * :c:data:`on_top_hf_grid`

 
.. c:var:: mos_in_r_array_transp


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_in_r_array_transp	(n_points_final_grid,mo_num)


    mos_in_r_array_transp(i,j) = value of the jth mo on the ith grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_mos_in_r_array`
       * :c:data:`core_inact_act_mos_in_r_array`
       * :c:data:`core_mos_in_r_array`
       * :c:data:`inact_mos_in_r_array`
       * :c:data:`virt_mos_in_r_array`

 
.. c:var:: mos_lapl_in_r_array


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_lapl_in_r_array	(mo_num,n_points_final_grid,3)


    mos_lapl_in_r_array(i,j,k)          = value of the kth component of the laplacian of ith mo on the jth grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_lapl_in_r_array_transp`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mos_lapl_in_r_array_tranp`

 
.. c:var:: mos_lapl_in_r_array_tranp


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_lapl_in_r_array_tranp	(3,mo_num,n_points_final_grid)


    mos_lapl_in_r_array_transp(i,j,k)   = value of the kth component of the laplient of jth mo on the ith grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mos_lapl_in_r_array`
       * :c:data:`n_points_final_grid`


 
.. c:var:: one_e_dm_and_grad_alpha_in_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_and_grad_alpha_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_and_grad_beta_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_beta_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: scal_prod_grad_one_e_dm_ab	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_stuff_for_pbe	(3,n_points_final_grid,N_states)


    one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
    
    one_e_grad_2_dm_alpha_at_r(i,istate)      = (d\dx n_alpha(r_i,istate))^2 + (d\dy n_alpha(r_i,istate))^2 + (d\dz n_alpha(r_i,istate))^2
    
    scal_prod_grad_one_e_dm_ab(i,istate)      = grad n_alpha(r_i) . grad n_beta(r_i)
    
    where r_i is the ith point of the grid and istate is the state number
    
    !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`effective_alpha_dm`
       * :c:data:`effective_spin_dm`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`energy_c_lda`
       * :c:data:`energy_c_sr_lda`
       * :c:data:`energy_x_lda`
       * :c:data:`energy_x_pbe`
       * :c:data:`energy_x_sr_lda`
       * :c:data:`energy_x_sr_pbe`
       * :c:data:`mu_average_prov`

 
.. c:var:: one_e_dm_and_grad_beta_in_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_and_grad_alpha_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_and_grad_beta_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_beta_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: scal_prod_grad_one_e_dm_ab	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_stuff_for_pbe	(3,n_points_final_grid,N_states)


    one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
    
    one_e_grad_2_dm_alpha_at_r(i,istate)      = (d\dx n_alpha(r_i,istate))^2 + (d\dy n_alpha(r_i,istate))^2 + (d\dz n_alpha(r_i,istate))^2
    
    scal_prod_grad_one_e_dm_ab(i,istate)      = grad n_alpha(r_i) . grad n_beta(r_i)
    
    where r_i is the ith point of the grid and istate is the state number
    
    !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`effective_alpha_dm`
       * :c:data:`effective_spin_dm`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`energy_c_lda`
       * :c:data:`energy_c_sr_lda`
       * :c:data:`energy_x_lda`
       * :c:data:`energy_x_pbe`
       * :c:data:`energy_x_sr_lda`
       * :c:data:`energy_x_sr_pbe`
       * :c:data:`mu_average_prov`

 
.. c:var:: one_e_grad_2_dm_alpha_at_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_and_grad_alpha_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_and_grad_beta_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_beta_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: scal_prod_grad_one_e_dm_ab	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_stuff_for_pbe	(3,n_points_final_grid,N_states)


    one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
    
    one_e_grad_2_dm_alpha_at_r(i,istate)      = (d\dx n_alpha(r_i,istate))^2 + (d\dy n_alpha(r_i,istate))^2 + (d\dz n_alpha(r_i,istate))^2
    
    scal_prod_grad_one_e_dm_ab(i,istate)      = grad n_alpha(r_i) . grad n_beta(r_i)
    
    where r_i is the ith point of the grid and istate is the state number
    
    !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`effective_alpha_dm`
       * :c:data:`effective_spin_dm`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`energy_c_lda`
       * :c:data:`energy_c_sr_lda`
       * :c:data:`energy_x_lda`
       * :c:data:`energy_x_pbe`
       * :c:data:`energy_x_sr_lda`
       * :c:data:`energy_x_sr_pbe`
       * :c:data:`mu_average_prov`

 
.. c:var:: one_e_grad_2_dm_beta_at_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_and_grad_alpha_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_and_grad_beta_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_beta_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: scal_prod_grad_one_e_dm_ab	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_stuff_for_pbe	(3,n_points_final_grid,N_states)


    one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
    
    one_e_grad_2_dm_alpha_at_r(i,istate)      = (d\dx n_alpha(r_i,istate))^2 + (d\dy n_alpha(r_i,istate))^2 + (d\dz n_alpha(r_i,istate))^2
    
    scal_prod_grad_one_e_dm_ab(i,istate)      = grad n_alpha(r_i) . grad n_beta(r_i)
    
    where r_i is the ith point of the grid and istate is the state number
    
    !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`effective_alpha_dm`
       * :c:data:`effective_spin_dm`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`energy_c_lda`
       * :c:data:`energy_c_sr_lda`
       * :c:data:`energy_x_lda`
       * :c:data:`energy_x_pbe`
       * :c:data:`energy_x_sr_lda`
       * :c:data:`energy_x_sr_pbe`
       * :c:data:`mu_average_prov`

 
.. c:var:: one_e_stuff_for_pbe


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_and_grad_alpha_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_and_grad_beta_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_beta_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: scal_prod_grad_one_e_dm_ab	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_stuff_for_pbe	(3,n_points_final_grid,N_states)


    one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
    
    one_e_grad_2_dm_alpha_at_r(i,istate)      = (d\dx n_alpha(r_i,istate))^2 + (d\dy n_alpha(r_i,istate))^2 + (d\dz n_alpha(r_i,istate))^2
    
    scal_prod_grad_one_e_dm_ab(i,istate)      = grad n_alpha(r_i) . grad n_beta(r_i)
    
    where r_i is the ith point of the grid and istate is the state number
    
    !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`effective_alpha_dm`
       * :c:data:`effective_spin_dm`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`energy_c_lda`
       * :c:data:`energy_c_sr_lda`
       * :c:data:`energy_x_lda`
       * :c:data:`energy_x_pbe`
       * :c:data:`energy_x_sr_lda`
       * :c:data:`energy_x_sr_pbe`
       * :c:data:`mu_average_prov`

 
.. c:var:: scal_prod_grad_one_e_dm_ab


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_and_grad_alpha_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_and_grad_beta_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_beta_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: scal_prod_grad_one_e_dm_ab	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_stuff_for_pbe	(3,n_points_final_grid,N_states)


    one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
    
    one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
    
    one_e_grad_2_dm_alpha_at_r(i,istate)      = (d\dx n_alpha(r_i,istate))^2 + (d\dy n_alpha(r_i,istate))^2 + (d\dz n_alpha(r_i,istate))^2
    
    scal_prod_grad_one_e_dm_ab(i,istate)      = grad n_alpha(r_i) . grad n_beta(r_i)
    
    where r_i is the ith point of the grid and istate is the state number
    
    !!!!! WARNING !!!! if no_core_density = .True. then all core electrons are removed

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`effective_alpha_dm`
       * :c:data:`effective_spin_dm`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`energy_c_lda`
       * :c:data:`energy_c_sr_lda`
       * :c:data:`energy_x_lda`
       * :c:data:`energy_x_pbe`
       * :c:data:`energy_x_sr_lda`
       * :c:data:`energy_x_sr_pbe`
       * :c:data:`mu_average_prov`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: dens_grad_a_b_no_core_and_aos_grad_aos_at_r:


    File : :file:`dft_utils_in_r/dm_in_r_routines.irp.f`

    .. code:: fortran

        subroutine dens_grad_a_b_no_core_and_aos_grad_aos_at_r(r,dm_a,dm_b, grad_dm_a, grad_dm_b, aos_array, grad_aos_array)


    input:
    
    * r(1) ==> r(1) = x, r(2) = y, r(3) = z
    
    output:
    
    * dm_a = alpha density evaluated at r without the core orbitals
    * dm_b = beta  density evaluated at r without the core orbitals
    * aos_array(i) = ao(i) evaluated at r without the core orbitals
    * grad_dm_a(1) = X gradient of the alpha density evaluated in r without the core orbitals
    * grad_dm_a(1) = X gradient of the beta  density evaluated in r without the core orbitals
    * grad_aos_array(1) = X gradient of the aos(i) evaluated at r
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft_no_core`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsymv`
       * :c:func:`give_all_aos_and_grad_at_r`

 
.. c:function:: density_and_grad_alpha_beta:


    File : :file:`dft_utils_in_r/dm_in_r_routines.irp.f`

    .. code:: fortran

        subroutine density_and_grad_alpha_beta(r,dm_a,dm_b, grad_dm_a, grad_dm_b)


    input:
    
    * r(1) ==> r(1) = x, r(2) = y, r(3) = z
    
    output:
    
    * dm_a = alpha density evaluated at r
    * dm_b = beta  density evaluated at r
    * grad_dm_a(1) = X gradient of the alpha density evaluated in r
    * grad_dm_a(1) = X gradient of the beta  density evaluated in r
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`mu_grad_rho_func`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsymv`
       * :c:func:`give_all_aos_and_grad_at_r`

 
.. c:function:: density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r:


    File : :file:`dft_utils_in_r/dm_in_r_routines.irp.f`

    .. code:: fortran

        subroutine density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a,dm_b, grad_dm_a, grad_dm_b, aos_array, grad_aos_array)


    input:
    
    * r(1) ==> r(1) = x, r(2) = y, r(3) = z
    
    output:
    
    * dm_a = alpha density evaluated at r
    * dm_b = beta  density evaluated at r
    * aos_array(i) = ao(i) evaluated at r
    * grad_dm_a(1) = X gradient of the alpha density evaluated in r
    * grad_dm_a(1) = X gradient of the beta  density evaluated in r
    * grad_aos_array(1) = X gradient of the aos(i) evaluated at r
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_md_on_top_pbe_mu_corrected`
       * :c:func:`ecmd_pbe_ueg_at_r`
       * :c:func:`give_all_stuffs_in_r_for_lyp_88`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsymv`
       * :c:func:`give_all_aos_and_grad_at_r`

 
.. c:function:: density_and_grad_lapl_alpha_beta_and_all_aos_and_grad_aos_at_r:


    File : :file:`dft_utils_in_r/dm_in_r_routines.irp.f`

    .. code:: fortran

        subroutine density_and_grad_lapl_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a,dm_b, grad_dm_a, grad_dm_b, lapl_dm_a, lapl_dm_b, aos_array, grad_aos_array, lapl_aos_array)


    input:
    
    * r(1) ==> r(1) = x, r(2) = y, r(3) = z
    
    output:
    
    * dm_a = alpha density evaluated at r
    * dm_b = beta  density evaluated at r
    * aos_array(i) = ao(i) evaluated at r
    * grad_dm_a(1) = X gradient of the alpha density evaluated in r
    * grad_dm_a(1) = X gradient of the beta  density evaluated in r
    * grad_aos_array(1) = X gradient of the aos(i) evaluated at r
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsymv`
       * :c:func:`give_all_aos_and_grad_and_lapl_at_r`

 
.. c:function:: dm_dft_alpha_beta_and_all_aos_at_r:


    File : :file:`dft_utils_in_r/dm_in_r_routines.irp.f`

    .. code:: fortran

        subroutine dm_dft_alpha_beta_and_all_aos_at_r(r,dm_a,dm_b,aos_array)


    input: r(1) ==> r(1) = x, r(2) = y, r(3) = z
    output : dm_a = alpha density evaluated at r
    output : dm_b = beta  density evaluated at r
    output : aos_array(i) = ao(i) evaluated at r

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsymv`
       * :c:func:`give_all_aos_at_r`

 
.. c:function:: dm_dft_alpha_beta_at_r:


    File : :file:`dft_utils_in_r/dm_in_r_routines.irp.f`

    .. code:: fortran

        subroutine dm_dft_alpha_beta_at_r(r,dm_a,dm_b)


    input: r(1) ==> r(1) = x, r(2) = y, r(3) = z
    output : dm_a = alpha density evaluated at r(3)
    output : dm_b = beta  density evaluated at r(3)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`correction_to_on_top_from_ueg`
       * :c:data:`mu_of_r_dft_average`
       * :c:data:`mu_rsc_of_r`
       * :c:func:`print_mos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`give_all_aos_at_r`

 
.. c:function:: dm_dft_alpha_beta_no_core_at_r:


    File : :file:`dft_utils_in_r/dm_in_r_routines.irp.f`

    .. code:: fortran

        subroutine dm_dft_alpha_beta_no_core_at_r(r,dm_a,dm_b)


    input: r(1) ==> r(1) = x, r(2) = y, r(3) = z
    output : dm_a = alpha density evaluated at r(3) without the core orbitals
    output : dm_b = beta  density evaluated at r(3) without the core orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_ao_for_dft_no_core`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`give_all_aos_at_r`

