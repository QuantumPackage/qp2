.. _module_functionals: 
 
.. program:: functionals 
 
.. default-role:: option 
 
functionals
===========

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Providers 
--------- 
 
.. c:var:: double aos_d_vc_alpha_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

 
.. c:var:: double aos_d_vc_alpha_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

 
.. c:var:: double aos_d_vc_beta_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

 
.. c:var:: double aos_d_vc_beta_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

 
.. c:var:: double aos_d_vx_alpha_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

 
.. c:var:: double aos_d_vx_alpha_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

 
.. c:var:: double aos_d_vx_beta_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

 
.. c:var:: double aos_d_vx_beta_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

 
.. c:var:: double aos_d_vxc_alpha_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vxc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vxc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_xc_alpha_ao_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_pbe`

 
.. c:var:: double aos_d_vxc_alpha_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vxc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vxc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_xc_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_sr_pbe`

 
.. c:var:: double aos_d_vxc_beta_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vxc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vxc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_xc_alpha_ao_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_pbe`

 
.. c:var:: double aos_d_vxc_beta_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vxc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vxc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_xc_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_sr_pbe`

 
.. c:var:: double aos_sr_vc_alpha_lda_w


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_lda_w	(ao_num,n_points_final_grid,N_states)


    aos_sr_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (sr_v^x_alpha(r_j) + sr_v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao_sr_lda`
       * :c:data:`potential_x_alpha_ao_sr_lda`

 
.. c:var:: double aos_sr_vc_beta_lda_w


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_lda_w	(ao_num,n_points_final_grid,N_states)


    aos_sr_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (sr_v^x_alpha(r_j) + sr_v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao_sr_lda`
       * :c:data:`potential_x_alpha_ao_sr_lda`

 
.. c:var:: double aos_sr_vx_alpha_lda_w


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_lda_w	(ao_num,n_points_final_grid,N_states)


    aos_sr_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (sr_v^x_alpha(r_j) + sr_v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao_sr_lda`
       * :c:data:`potential_x_alpha_ao_sr_lda`

 
.. c:var:: double aos_sr_vx_beta_lda_w


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_lda_w	(ao_num,n_points_final_grid,N_states)


    aos_sr_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (sr_v^x_alpha(r_j) + sr_v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao_sr_lda`
       * :c:data:`potential_x_alpha_ao_sr_lda`

 
.. c:var:: double aos_sr_vxc_alpha_lda_w


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vxc_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vxc_beta_lda_w	(ao_num,n_points_final_grid,N_states)


    aos_sr_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_sr_lda`

 
.. c:var:: double aos_sr_vxc_beta_lda_w


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vxc_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vxc_beta_lda_w	(ao_num,n_points_final_grid,N_states)


    aos_sr_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_sr_lda`

 
.. c:var:: double aos_vc_alpha_lda_w


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_lda_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_lda`

 
.. c:var:: double aos_vc_alpha_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

 
.. c:var:: double aos_vc_alpha_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

 
.. c:var:: double aos_vc_beta_lda_w


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_lda_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_lda`

 
.. c:var:: double aos_vc_beta_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

 
.. c:var:: double aos_vc_beta_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

 
.. c:var:: double aos_vx_alpha_lda_w


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_lda_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_lda`

 
.. c:var:: double aos_vx_alpha_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

 
.. c:var:: double aos_vx_alpha_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

 
.. c:var:: double aos_vx_beta_lda_w


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_lda_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_lda_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_lda`

 
.. c:var:: double aos_vx_beta_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

 
.. c:var:: double aos_vx_beta_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vx_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    intermediates to compute the sr_pbe potentials
    
    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

 
.. c:var:: double aos_vxc_alpha_lda_w


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vxc_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vxc_beta_lda_w	(n_points_final_grid,ao_num,N_states)


    aos_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_lda`

 
.. c:var:: double aos_vxc_alpha_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vxc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vxc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_xc_alpha_ao_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_pbe`

 
.. c:var:: double aos_vxc_alpha_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vxc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vxc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_xc_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_sr_pbe`

 
.. c:var:: double aos_vxc_beta_lda_w


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vxc_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vxc_beta_lda_w	(n_points_final_grid,ao_num,N_states)


    aos_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_lda`

 
.. c:var:: double aos_vxc_beta_pbe_w


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vxc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vxc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_xc_alpha_ao_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_pbe`

 
.. c:var:: double aos_vxc_beta_sr_pbe_w


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vxc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vxc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_alpha_sr_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_d_vxc_beta_sr_pbe_w	(ao_num,n_points_final_grid,N_states)


    aos_vxc_alpha_sr_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pot_grad_xc_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_sr_pbe`

 
.. c:var:: double energy_c_lda


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_c_lda	(N_states)


    correlation energy with the lda functional

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

       * :c:data:`energy_c`

 
.. c:var:: double energy_c_none


    File : :file:`functionals/none.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_c_none	(N_states)


    null correlation energy

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_c`

 
.. c:var:: double energy_c_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_pbe	(N_states)
        double precision, allocatable	:: energy_c_pbe	(N_states)


    exchange / correlation energies  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)
    exchange/correlation energy with the short range pbe functional

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

       * :c:data:`energy_c`
       * :c:data:`energy_x`

 
.. c:var:: double energy_c_sr_lda


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_c_sr_lda	(N_states)


    exchange             energy with the short range lda functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_c`

 
.. c:var:: double energy_c_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_sr_pbe	(N_states)
        double precision, allocatable	:: energy_c_sr_pbe	(N_states)


    exchange / correlation energies  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)
    exchange/correlation energy with the short range pbe functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_c`
       * :c:data:`energy_x`

 
.. c:var:: double energy_x_lda


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_lda	(N_states)


    exchange energy with the lda functional

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

       * :c:data:`energy_x`

 
.. c:var:: double energy_x_none


    File : :file:`functionals/none.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_none	(N_states)


    null exchange energy

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x`

 
.. c:var:: double energy_x_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_pbe	(N_states)
        double precision, allocatable	:: energy_c_pbe	(N_states)


    exchange / correlation energies  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)
    exchange/correlation energy with the short range pbe functional

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

       * :c:data:`energy_c`
       * :c:data:`energy_x`

 
.. c:var:: double energy_x_sr_lda


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_sr_lda	(N_states)


    exchange             energy with the short range lda functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x`

 
.. c:var:: double energy_x_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_sr_pbe	(N_states)
        double precision, allocatable	:: energy_c_sr_pbe	(N_states)


    exchange / correlation energies  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)
    exchange/correlation energy with the short range pbe functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_c`
       * :c:data:`energy_x`

 
.. c:var:: double pot_grad_c_alpha_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: double pot_grad_c_alpha_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_sr_pbe`

 
.. c:var:: double pot_grad_c_beta_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: double pot_grad_c_beta_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_sr_pbe`

 
.. c:var:: double pot_grad_x_alpha_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: double pot_grad_x_alpha_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_sr_pbe`

 
.. c:var:: double pot_grad_x_beta_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: double pot_grad_x_beta_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_sr_pbe`

 
.. c:var:: double pot_grad_xc_alpha_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_xc_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_xc_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_pbe`

 
.. c:var:: double pot_grad_xc_alpha_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_xc_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_xc_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_sr_pbe`

 
.. c:var:: double pot_grad_xc_beta_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_xc_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_xc_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_pbe`

 
.. c:var:: double pot_grad_xc_beta_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_grad_xc_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_grad_xc_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_sr_pbe`

 
.. c:var:: double pot_scal_c_alpha_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: double pot_scal_c_alpha_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_sr_pbe`

 
.. c:var:: double pot_scal_c_beta_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: double pot_scal_c_beta_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_sr_pbe`

 
.. c:var:: double pot_scal_x_alpha_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: double pot_scal_x_alpha_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_sr_pbe`

 
.. c:var:: double pot_scal_x_beta_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: double pot_scal_x_beta_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_sr_pbe`

 
.. c:var:: double pot_scal_xc_alpha_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_xc_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_xc_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_pbe`

 
.. c:var:: double pot_scal_xc_alpha_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_xc_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_xc_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_sr_pbe`

 
.. c:var:: double pot_scal_xc_beta_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_xc_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_xc_beta_ao_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vxc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_pbe`

 
.. c:var:: double pot_scal_xc_beta_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pot_scal_xc_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: pot_scal_xc_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vxc_alpha_sr_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao_sr_pbe`

 
.. c:var:: double potential_c_alpha_ao_lda


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_c_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_lda	(ao_num,ao_num,N_states)


    short range correlation alpha/beta potentials with lda functional on the |AO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`

 
.. c:var:: double potential_c_alpha_ao_none


    File : :file:`functionals/none.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_c_alpha_ao_none	(ao_num,ao_num,N_states)


    Potential for a null exchange-correlation functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`

 
.. c:var:: double potential_c_alpha_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`
       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_c_alpha_ao_sr_lda


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_c_alpha_ao_sr_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_sr_lda	(ao_num,ao_num,N_states)


    short range correlation alpha/beta potentials with lda functional on the |AO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`

 
.. c:var:: double potential_c_alpha_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`
       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_c_beta_ao_lda


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_c_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_lda	(ao_num,ao_num,N_states)


    short range correlation alpha/beta potentials with lda functional on the |AO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`

 
.. c:var:: double potential_c_beta_ao_none


    File : :file:`functionals/none.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_c_beta_ao_none	(ao_num,ao_num,N_states)


    Potential for a null exchange-correlation functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`

 
.. c:var:: double potential_c_beta_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`
       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_c_beta_ao_sr_lda


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_c_alpha_ao_sr_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_sr_lda	(ao_num,ao_num,N_states)


    short range correlation alpha/beta potentials with lda functional on the |AO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`

 
.. c:var:: double potential_c_beta_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`
       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_x_alpha_ao_lda


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_lda	(ao_num,ao_num,N_states)


    short range exchange alpha/beta potentials with lda functional on the |AO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_x_alpha_ao_none


    File : :file:`functionals/none.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_none	(ao_num,ao_num,N_states)


    Potential for a null exchange-correlation functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_x_alpha_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`
       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_x_alpha_ao_sr_lda


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_sr_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_sr_lda	(ao_num,ao_num,N_states)


    short range exchange alpha/beta potentials with lda functional on the |AO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_x_alpha_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`
       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_x_beta_ao_lda


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_lda	(ao_num,ao_num,N_states)


    short range exchange alpha/beta potentials with lda functional on the |AO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_x_beta_ao_none


    File : :file:`functionals/none.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_beta_ao_none	(ao_num,ao_num,N_states)


    Potential for a null exchange-correlation functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_x_beta_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_x_alpha_ao_pbe`
       * :c:data:`pot_scal_x_alpha_ao_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`
       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_x_beta_ao_sr_lda


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_sr_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_sr_lda	(ao_num,ao_num,N_states)


    short range exchange alpha/beta potentials with lda functional on the |AO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_x_beta_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional
    
    defined in Chem. Phys.329, 276 (2006)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_x_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_x_alpha_ao_sr_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_c_alpha_ao`
       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: double potential_xc_alpha_ao_lda


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_xc_beta_ao_lda	(ao_num,ao_num,N_states)


    short range exchange/correlation alpha/beta potentials with lda functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao`

 
.. c:var:: double potential_xc_alpha_ao_none


    File : :file:`functionals/none.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao_none	(ao_num,ao_num,N_states)


    Potential for a null exchange-correlation functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao`

 
.. c:var:: double potential_xc_alpha_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_xc_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the Perdew-Burke-Ernzerhof GGA functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_xc_alpha_ao_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao`

 
.. c:var:: double potential_xc_alpha_ao_sr_lda


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao_sr_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_xc_beta_ao_sr_lda	(ao_num,ao_num,N_states)


    short range exchange/correlation alpha/beta potentials with lda functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao`

 
.. c:var:: double potential_xc_alpha_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_xc_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the Perdew-Burke-Ernzerhof GGA functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_xc_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_sr_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao`

 
.. c:var:: double potential_xc_beta_ao_lda


    File : :file:`functionals/lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_xc_beta_ao_lda	(ao_num,ao_num,N_states)


    short range exchange/correlation alpha/beta potentials with lda functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vxc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao`

 
.. c:var:: double potential_xc_beta_ao_none


    File : :file:`functionals/none.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_beta_ao_none	(ao_num,ao_num,N_states)


    Potential for a null exchange-correlation functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao`

 
.. c:var:: double potential_xc_beta_ao_pbe


    File : :file:`functionals/pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_xc_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the Perdew-Burke-Ernzerhof GGA functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_xc_alpha_ao_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao`

 
.. c:var:: double potential_xc_beta_ao_sr_lda


    File : :file:`functionals/sr_lda.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao_sr_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_xc_beta_ao_sr_lda	(ao_num,ao_num,N_states)


    short range exchange/correlation alpha/beta potentials with lda functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_sr_vxc_alpha_lda_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao`

 
.. c:var:: double potential_xc_beta_ao_sr_pbe


    File : :file:`functionals/sr_pbe.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_xc_alpha_ao_sr_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_xc_beta_ao_sr_pbe	(ao_num,ao_num,N_states)


    exchange / correlation potential for alpha / beta electrons  with the Perdew-Burke-Ernzerhof GGA functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`n_states`
       * :c:data:`pot_grad_xc_alpha_ao_sr_pbe`
       * :c:data:`pot_scal_xc_alpha_ao_sr_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_xc_alpha_ao`

