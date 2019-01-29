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
 
.. c:var:: aos_grad_in_r_array


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_grad_in_r_array	(ao_num,n_points_final_grid,3)
        double precision, allocatable	:: aos_grad_in_r_array_transp	(n_points_final_grid,ao_num,3)
        double precision, allocatable	:: aos_grad_in_r_array_transp_xyz	(3,n_points_final_grid,ao_num)


    aos_grad_in_r_array(i,j,k)          = value of the kth component of the gradient of ith ao on the jth grid point
    
    aos_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth ao on the ith grid point
    
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

       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`potential_sr_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_grad_in_r_array_transp


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_grad_in_r_array	(ao_num,n_points_final_grid,3)
        double precision, allocatable	:: aos_grad_in_r_array_transp	(n_points_final_grid,ao_num,3)
        double precision, allocatable	:: aos_grad_in_r_array_transp_xyz	(3,n_points_final_grid,ao_num)


    aos_grad_in_r_array(i,j,k)          = value of the kth component of the gradient of ith ao on the jth grid point
    
    aos_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth ao on the ith grid point
    
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

       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`potential_sr_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_grad_in_r_array_transp_xyz


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_grad_in_r_array	(ao_num,n_points_final_grid,3)
        double precision, allocatable	:: aos_grad_in_r_array_transp	(n_points_final_grid,ao_num,3)
        double precision, allocatable	:: aos_grad_in_r_array_transp_xyz	(3,n_points_final_grid,ao_num)


    aos_grad_in_r_array(i,j,k)          = value of the kth component of the gradient of ith ao on the jth grid point
    
    aos_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth ao on the ith grid point
    
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

       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`potential_sr_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_in_r_array


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_in_r_array	(ao_num,n_points_final_grid)
        double precision, allocatable	:: aos_in_r_array_transp	(n_points_final_grid,ao_num)


    aos_in_r_array(i,j)        = value of the ith ao on the jth grid point
    
    aos_in_r_array_transp(i,j) = value of the jth ao on the ith grid point

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

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`potential_sr_c_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_in_r_array_transp


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_in_r_array	(ao_num,n_points_final_grid)
        double precision, allocatable	:: aos_in_r_array_transp	(n_points_final_grid,ao_num)


    aos_in_r_array(i,j)        = value of the ith ao on the jth grid point
    
    aos_in_r_array_transp(i,j) = value of the jth ao on the ith grid point

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

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`potential_sr_c_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_lapl_in_r_array


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_lapl_in_r_array	(ao_num,n_points_final_grid,3)
        double precision, allocatable	:: aos_lapl_in_r_array_transp	(n_points_final_grid,ao_num,3)


    aos_lapl_in_r_array(i,j,k)          = value of the kth component of the laplacian of ith ao on the jth grid point
    
    aos_lapl_in_r_array_transp(i,j,k)   = value of the kth component of the laplacian of jth ao on the ith grid point
    
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

       * :c:data:`mos_lapl_in_r_array`

 
.. c:var:: aos_lapl_in_r_array_transp


    File : :file:`dft_utils_in_r/ao_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_lapl_in_r_array	(ao_num,n_points_final_grid,3)
        double precision, allocatable	:: aos_lapl_in_r_array_transp	(n_points_final_grid,ao_num,3)


    aos_lapl_in_r_array(i,j,k)          = value of the kth component of the laplacian of ith ao on the jth grid point
    
    aos_lapl_in_r_array_transp(i,j,k)   = value of the kth component of the laplacian of jth ao on the ith grid point
    
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

       * :c:data:`mos_lapl_in_r_array`

 
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


 
.. c:var:: mos_in_r_array


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_in_r_array	(mo_num,n_points_final_grid)
        double precision, allocatable	:: mos_in_r_array_transp	(n_points_final_grid,mo_num)


    mos_in_r_array(i,j)        = value of the ith mo on the jth grid point
    
    mos_in_r_array_transp(i,j) = value of the jth mo on the ith grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`


 
.. c:var:: mos_in_r_array_transp


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_in_r_array	(mo_num,n_points_final_grid)
        double precision, allocatable	:: mos_in_r_array_transp	(n_points_final_grid,mo_num)


    mos_in_r_array(i,j)        = value of the ith mo on the jth grid point
    
    mos_in_r_array_transp(i,j) = value of the jth mo on the ith grid point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`final_grid_points`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`


 
.. c:var:: mos_lapl_in_r_array


    File : :file:`dft_utils_in_r/mo_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_lapl_in_r_array	(mo_num,n_points_final_grid,3)


    mos_lapl_in_r_array(i,j,k)          = value of the kth component of the laplacian of ith mo on the jth grid point
    
    mos_lapl_in_r_array_transp(i,j,k)   = value of the kth component of the laplacian of jth mo on the ith grid point
    
    k = 1 : x, k= 2, y, k  3, z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_lapl_in_r_array`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_points_final_grid`


 
.. c:var:: one_e_dm_alpha_at_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_beta_at_r	(n_points_final_grid,N_states)


    one_e_dm_alpha_at_r(i,istate) = n_alpha(r_i,istate)
    one_e_dm_beta_at_r(i,istate) =  n_beta(r_i,istate)
    where r_i is the ith point of the grid and istate is the state number

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
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`energy_sr_x_lda`
       * :c:data:`energy_x_lda`

 
.. c:var:: one_e_dm_alpha_in_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_alpha_in_r	(n_points_integration_angular,n_points_radial_grid,nucl_num,N_states)
        double precision, allocatable	:: one_e_dm_beta_in_r	(n_points_integration_angular,n_points_radial_grid,nucl_num,N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`grid_points_per_atom`
       * :c:data:`mo_num`
       * :c:data:`n_points_radial_grid`
       * :c:data:`n_states`
       * :c:data:`nucl_num`
       * :c:data:`one_e_dm_alpha_ao_for_dft`


 
.. c:var:: one_e_dm_and_grad_alpha_in_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_and_grad_alpha_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_and_grad_beta_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_beta_at_r	(n_points_final_grid,N_states)


    one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
    one_e_grad_2_dm_alpha_at_r(i,istate)      = d\dx n_alpha(r_i,istate)^2 + d\dy n_alpha(r_i,istate)^2 + d\dz n_alpha(r_i,istate)^2
    where r_i is the ith point of the grid and istate is the state number

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

       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`energy_sr_x_pbe`
       * :c:data:`energy_x_pbe`

 
.. c:var:: one_e_dm_and_grad_beta_in_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_and_grad_alpha_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_and_grad_beta_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_beta_at_r	(n_points_final_grid,N_states)


    one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
    one_e_grad_2_dm_alpha_at_r(i,istate)      = d\dx n_alpha(r_i,istate)^2 + d\dy n_alpha(r_i,istate)^2 + d\dz n_alpha(r_i,istate)^2
    where r_i is the ith point of the grid and istate is the state number

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

       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`energy_sr_x_pbe`
       * :c:data:`energy_x_pbe`

 
.. c:var:: one_e_dm_beta_at_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_beta_at_r	(n_points_final_grid,N_states)


    one_e_dm_alpha_at_r(i,istate) = n_alpha(r_i,istate)
    one_e_dm_beta_at_r(i,istate) =  n_beta(r_i,istate)
    where r_i is the ith point of the grid and istate is the state number

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
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`energy_sr_x_lda`
       * :c:data:`energy_x_lda`

 
.. c:var:: one_e_dm_beta_in_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_alpha_in_r	(n_points_integration_angular,n_points_radial_grid,nucl_num,N_states)
        double precision, allocatable	:: one_e_dm_beta_in_r	(n_points_integration_angular,n_points_radial_grid,nucl_num,N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`grid_points_per_atom`
       * :c:data:`mo_num`
       * :c:data:`n_points_radial_grid`
       * :c:data:`n_states`
       * :c:data:`nucl_num`
       * :c:data:`one_e_dm_alpha_ao_for_dft`


 
.. c:var:: one_e_grad_2_dm_alpha_at_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_and_grad_alpha_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_and_grad_beta_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_beta_at_r	(n_points_final_grid,N_states)


    one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
    one_e_grad_2_dm_alpha_at_r(i,istate)      = d\dx n_alpha(r_i,istate)^2 + d\dy n_alpha(r_i,istate)^2 + d\dz n_alpha(r_i,istate)^2
    where r_i is the ith point of the grid and istate is the state number

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

       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`energy_sr_x_pbe`
       * :c:data:`energy_x_pbe`

 
.. c:var:: one_e_grad_2_dm_beta_at_r


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_and_grad_alpha_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_dm_and_grad_beta_in_r	(4,n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_alpha_at_r	(n_points_final_grid,N_states)
        double precision, allocatable	:: one_e_grad_2_dm_beta_at_r	(n_points_final_grid,N_states)


    one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
    one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
    one_e_grad_2_dm_alpha_at_r(i,istate)      = d\dx n_alpha(r_i,istate)^2 + d\dy n_alpha(r_i,istate)^2 + d\dz n_alpha(r_i,istate)^2
    where r_i is the ith point of the grid and istate is the state number

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

       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`energy_sr_x_pbe`
       * :c:data:`energy_x_pbe`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r:


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

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
       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsymv`
       * :c:func:`give_all_aos_and_grad_at_r`

 
.. c:function:: dm_dft_alpha_beta_and_all_aos_at_r:


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

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
       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsymv`
       * :c:func:`give_all_aos_at_r`

 
.. c:function:: dm_dft_alpha_beta_at_r:


    File : :file:`dft_utils_in_r/dm_in_r.irp.f`

    .. code:: fortran

        subroutine dm_dft_alpha_beta_at_r(r,dm_a,dm_b)


    input: r(1) ==> r(1) = x, r(2) = y, r(3) = z
    output : dm_a = alpha density evaluated at r(3)
    output : dm_b = beta  density evaluated at r(3)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_alpha_at_r`
       * :c:data:`one_e_dm_alpha_in_r`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`give_all_aos_at_r`

