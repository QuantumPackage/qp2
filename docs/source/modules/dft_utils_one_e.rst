.. _module_dft_utils_one_e: 
 
.. program:: dft_utils_one_e 
 
.. default-role:: option 
 
===============
dft_utils_one_e
===============

This module contains all the one-body related quantities needed to perform DFT or RS-DFT calculations.
Therefore, it contains most of the properties which depends on the one-body density and density matrix.

The most important files and variables are:

* The general *providers* for the x/c energies in :file:`e_xc_general.irp.f`
* The general *providers* for the x/c potentials in :file:`pot_general.irp.f`
* The short-range hartree operator and all related quantities in :file:`sr_coulomb.irp.f`

These *providers* will be used in many DFT-related programs, such as :file:`ks_scf.irp.f` or :file:`rs_ks_scf.irp.f`.
It is also needed to compute the effective one-body operator needed in multi-determinant RS-DFT (see plugins by eginer).

Some other interesting quantities:

* The LDA and PBE *providers* for the x/c energies in :file:`e_xc.irp.f` and :file:`sr_exc.irp.f`
* The LDA and PBE *providers* for the x/c potentials on the AO basis in :file:`pot_ao.irp.f` and  :file:`sr_pot_ao.irp.f`
* The :math:`h_{core}` energy computed directly with the one-body density matrix in :file:`one_e_energy_dft.irp.f`
* LDA and PBE short-range functionals *subroutines* in :file:`exc_sr_lda.irp.f` and :file:`exc_sr_pbe.irp.f`


 
 
 
Providers 
--------- 
 
.. c:var:: ao_effective_one_e_potential


    File : :file:`dft_utils_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_effective_one_e_potential	(ao_num,ao_num,N_states)
        double precision, allocatable	:: ao_effective_one_e_potential_without_kin	(ao_num,ao_num,N_states)


    ao_effective_one_e_potential(i,j) = :math:`\rangle i_{AO}| v_{H}^{sr} |j_{AO}\rangle  + \rangle i_{AO}| h_{core} |j_{AO}\rangle  + \rangle i_{AO}|v_{xc} |j_{AO}\rangle` 
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`effective_one_e_potential`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`


 
.. c:var:: ao_effective_one_e_potential_without_kin


    File : :file:`dft_utils_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_effective_one_e_potential	(ao_num,ao_num,N_states)
        double precision, allocatable	:: ao_effective_one_e_potential_without_kin	(ao_num,ao_num,N_states)


    ao_effective_one_e_potential(i,j) = :math:`\rangle i_{AO}| v_{H}^{sr} |j_{AO}\rangle  + \rangle i_{AO}| h_{core} |j_{AO}\rangle  + \rangle i_{AO}|v_{xc} |j_{AO}\rangle` 
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`effective_one_e_potential`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`


 
.. c:var:: aos_dsr_vc_alpha_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: aos_dsr_vc_beta_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: aos_dsr_vx_alpha_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: aos_dsr_vx_beta_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: aos_dvc_alpha_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_dvc_beta_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_dvx_alpha_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_dvx_beta_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_sr_vc_alpha_lda_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_lda_w	(n_points_final_grid,ao_num,N_states)


    aos_sr_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (sr_v^x_alpha(r_j) + sr_v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_c_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_lda`

 
.. c:var:: aos_sr_vc_alpha_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: aos_sr_vc_beta_lda_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_lda_w	(n_points_final_grid,ao_num,N_states)


    aos_sr_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (sr_v^x_alpha(r_j) + sr_v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_c_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_lda`

 
.. c:var:: aos_sr_vc_beta_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: aos_sr_vx_alpha_lda_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_lda_w	(n_points_final_grid,ao_num,N_states)


    aos_sr_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (sr_v^x_alpha(r_j) + sr_v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_c_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_lda`

 
.. c:var:: aos_sr_vx_alpha_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: aos_sr_vx_beta_lda_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_lda_w	(n_points_final_grid,ao_num,N_states)


    aos_sr_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (sr_v^x_alpha(r_j) + sr_v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_c_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_lda`

 
.. c:var:: aos_sr_vx_beta_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: aos_vc_alpha_lda_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vc_beta_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vx_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vx_beta_lda_w	(n_points_final_grid,ao_num,N_states)


    aos_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_lda`

 
.. c:var:: aos_vc_alpha_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_vc_beta_lda_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vc_beta_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vx_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vx_beta_lda_w	(n_points_final_grid,ao_num,N_states)


    aos_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_lda`

 
.. c:var:: aos_vc_beta_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_vx_alpha_lda_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vc_beta_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vx_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vx_beta_lda_w	(n_points_final_grid,ao_num,N_states)


    aos_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_lda`

 
.. c:var:: aos_vx_alpha_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: aos_vx_beta_lda_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vc_beta_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vx_alpha_lda_w	(n_points_final_grid,ao_num,N_states)
        double precision, allocatable	:: aos_vx_beta_lda_w	(n_points_final_grid,ao_num,N_states)


    aos_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_lda`

 
.. c:var:: aos_vx_beta_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: effective_one_e_potential


    File : :file:`dft_utils_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_one_e_potential	(mo_num,mo_num,N_states)
        double precision, allocatable	:: effective_one_e_potential_without_kin	(mo_num,mo_num,N_states)


    Effective_one_e_potential(i,j) = :math:`\rangle i_{MO}| v_{H}^{sr} |j_{MO}\rangle  + \rangle i_{MO}| h_{core} |j_{MO}\rangle  + \rangle i_{MO}|v_{xc} |j_{MO}\rangle` 
    
    on the |MO| basis
    Taking the expectation value does not provide any energy, but
    effective_one_e_potential(i,j) is the potential coupling DFT and WFT part to
    be used in any WFT calculation.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_x_alpha_mo`
       * :c:data:`short_range_hartree_operator`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_effective_one_e_potential`

 
.. c:var:: effective_one_e_potential_without_kin


    File : :file:`dft_utils_one_e/effective_pot.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_one_e_potential	(mo_num,mo_num,N_states)
        double precision, allocatable	:: effective_one_e_potential_without_kin	(mo_num,mo_num,N_states)


    Effective_one_e_potential(i,j) = :math:`\rangle i_{MO}| v_{H}^{sr} |j_{MO}\rangle  + \rangle i_{MO}| h_{core} |j_{MO}\rangle  + \rangle i_{MO}|v_{xc} |j_{MO}\rangle` 
    
    on the |MO| basis
    Taking the expectation value does not provide any energy, but
    effective_one_e_potential(i,j) is the potential coupling DFT and WFT part to
    be used in any WFT calculation.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_x_alpha_mo`
       * :c:data:`short_range_hartree_operator`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_effective_one_e_potential`

 
.. c:var:: energy_c


    File : :file:`dft_utils_one_e/e_xc_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x	(N_states)
        double precision, allocatable	:: energy_c	(N_states)


    correlation and exchange energies general providers.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`correlation_functional`
       * :c:data:`energy_sr_x_lda`
       * :c:data:`energy_sr_x_pbe`
       * :c:data:`exchange_functional`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_correlation_dft`
       * :c:data:`e_exchange_dft`
       * :c:data:`shifting_constant`

 
.. c:var:: energy_c_lda


    File : :file:`dft_utils_one_e/e_xc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_lda	(N_states)
        double precision, allocatable	:: energy_c_lda	(N_states)


    exchange/correlation energy with the short range LDA functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`


 
.. c:var:: energy_c_pbe


    File : :file:`dft_utils_one_e/e_xc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_pbe	(N_states)
        double precision, allocatable	:: energy_c_pbe	(N_states)


    exchange/correlation energy with the short range PBE functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: energy_sr_c_lda


    File : :file:`dft_utils_one_e/sr_exc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_sr_x_lda	(N_states)
        double precision, allocatable	:: energy_sr_c_lda	(N_states)


    exchange/correlation energy with the short range LDA functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x`

 
.. c:var:: energy_sr_c_pbe


    File : :file:`dft_utils_one_e/sr_exc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_sr_x_pbe	(N_states)
        double precision, allocatable	:: energy_sr_c_pbe	(N_states)


    exchange/correlation energy with the short range PBE functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x`

 
.. c:var:: energy_sr_x_lda


    File : :file:`dft_utils_one_e/sr_exc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_sr_x_lda	(N_states)
        double precision, allocatable	:: energy_sr_c_lda	(N_states)


    exchange/correlation energy with the short range LDA functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x`

 
.. c:var:: energy_sr_x_pbe


    File : :file:`dft_utils_one_e/sr_exc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_sr_x_pbe	(N_states)
        double precision, allocatable	:: energy_sr_c_pbe	(N_states)


    exchange/correlation energy with the short range PBE functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x`

 
.. c:var:: energy_x


    File : :file:`dft_utils_one_e/e_xc_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x	(N_states)
        double precision, allocatable	:: energy_c	(N_states)


    correlation and exchange energies general providers.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`correlation_functional`
       * :c:data:`energy_sr_x_lda`
       * :c:data:`energy_sr_x_pbe`
       * :c:data:`exchange_functional`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_correlation_dft`
       * :c:data:`e_exchange_dft`
       * :c:data:`shifting_constant`

 
.. c:var:: energy_x_lda


    File : :file:`dft_utils_one_e/e_xc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_lda	(N_states)
        double precision, allocatable	:: energy_c_lda	(N_states)


    exchange/correlation energy with the short range LDA functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_alpha_at_r`


 
.. c:var:: energy_x_pbe


    File : :file:`dft_utils_one_e/e_xc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_x_pbe	(N_states)
        double precision, allocatable	:: energy_c_pbe	(N_states)


    exchange/correlation energy with the short range PBE functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:function:: gga_sr_type_functionals:


    File : :file:`dft_utils_one_e/utils.irp.f`

    .. code:: fortran

        subroutine GGA_sr_type_functionals(r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
                                ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &
                                ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b )


    routine that helps in building the x/c potentials on the AO basis for a GGA functional with a short-range interaction

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf_dft`
       * :c:data:`exchange_functional`
       * :c:data:`correlation_functional`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`energy_sr_x_pbe`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_pbe_sr`
       * :c:func:`ex_pbe_sr`
       * :c:func:`grad_rho_ab_to_grad_rho_oc`
       * :c:func:`rho_ab_to_rho_oc`
       * :c:func:`v_grad_rho_oc_to_v_grad_rho_ab`
       * :c:func:`v_rho_oc_to_v_rho_ab`

 
.. c:function:: gga_type_functionals:


    File : :file:`dft_utils_one_e/utils.irp.f`

    .. code:: fortran

        subroutine GGA_type_functionals(r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
                                ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &
                                ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b )


    routine that helps in building the x/c potentials on the AO basis for a GGA functional

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`exchange_functional`
       * :c:data:`correlation_functional`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`energy_x_pbe`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_pbe_sr`
       * :c:func:`ex_pbe_sr`
       * :c:func:`grad_rho_ab_to_grad_rho_oc`
       * :c:func:`rho_ab_to_rho_oc`
       * :c:func:`v_grad_rho_oc_to_v_grad_rho_ab`
       * :c:func:`v_rho_oc_to_v_rho_ab`

 
.. c:var:: grad_aos_dsr_vc_alpha_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: grad_aos_dsr_vc_beta_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: grad_aos_dsr_vx_alpha_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: grad_aos_dsr_vx_beta_pbe_w


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_sr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_sr_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dsr_vx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`mu_erf_dft`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_sr_x_alpha_ao_pbe`

 
.. c:var:: grad_aos_dvc_alpha_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: grad_aos_dvc_beta_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: grad_aos_dvx_alpha_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: grad_aos_dvx_beta_pbe_w


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: aos_vc_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vc_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_alpha_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_vx_beta_pbe_w	(ao_num,n_points_final_grid,N_states)
        double precision, allocatable	:: aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvc_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_alpha_pbe_w	(ao_num,n_points_final_grid,3,N_states)
        double precision, allocatable	:: grad_aos_dvx_beta_pbe_w	(ao_num,n_points_final_grid,3,N_states)


    aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao_pbe`

 
.. c:var:: mu_erf_dft


    File : :file:`dft_utils_one_e/mu_erf_dft.irp.f`

    .. code:: fortran

        double precision	:: mu_erf_dft	


    range separation parameter used in RS-DFT. It is set to mu_erf in order to be consistent with the two electrons integrals erf

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`energy_sr_x_lda`
       * :c:data:`energy_sr_x_pbe`

 
.. c:var:: potential_c_alpha_ao


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao	(ao_num,ao_num,N_states)


    general providers for the alpha/beta exchange/correlation potentials on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`n_states`
       * :c:data:`potential_sr_c_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`potential_x_alpha_mo`

 
.. c:var:: potential_c_alpha_ao_lda


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_lda	(ao_num,ao_num,N_states)


    short range exchange/correlation alpha/beta potentials with LDA functional on the AO basis

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

 
.. c:var:: potential_c_alpha_ao_pbe


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange/correlation alpha/beta potentials with the short range PBE functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: potential_c_alpha_mo


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_x_beta_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_c_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_c_beta_mo	(mo_num,mo_num,N_states)


    general providers for the alpha/beta exchange/correlation potentials on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_x_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`trace_v_xc`

 
.. c:var:: potential_c_beta_ao


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao	(ao_num,ao_num,N_states)


    general providers for the alpha/beta exchange/correlation potentials on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`n_states`
       * :c:data:`potential_sr_c_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`potential_x_alpha_mo`

 
.. c:var:: potential_c_beta_ao_lda


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_lda	(ao_num,ao_num,N_states)


    short range exchange/correlation alpha/beta potentials with LDA functional on the AO basis

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

 
.. c:var:: potential_c_beta_ao_pbe


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange/correlation alpha/beta potentials with the short range PBE functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: potential_c_beta_mo


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_x_beta_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_c_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_c_beta_mo	(mo_num,mo_num,N_states)


    general providers for the alpha/beta exchange/correlation potentials on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_x_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`trace_v_xc`

 
.. c:var:: potential_sr_c_alpha_ao_lda


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_sr_c_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_c_beta_ao_lda	(ao_num,ao_num,N_states)


    short range correlation alpha/beta potentials with LDA functional on the |AO| basis

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

 
.. c:var:: potential_sr_c_alpha_ao_pbe


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_sr_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange/correlation alpha/beta potentials with the short range PBE functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: potential_sr_c_beta_ao_lda


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_sr_c_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_c_beta_ao_lda	(ao_num,ao_num,N_states)


    short range correlation alpha/beta potentials with LDA functional on the |AO| basis

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

 
.. c:var:: potential_sr_c_beta_ao_pbe


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_sr_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange/correlation alpha/beta potentials with the short range PBE functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: potential_sr_x_alpha_ao_lda


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_sr_x_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_x_beta_ao_lda	(ao_num,ao_num,N_states)


    short range exchange alpha/beta potentials with LDA functional on the |AO| basis

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

 
.. c:var:: potential_sr_x_alpha_ao_pbe


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_sr_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange/correlation alpha/beta potentials with the short range PBE functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: potential_sr_x_beta_ao_lda


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_sr_x_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_x_beta_ao_lda	(ao_num,ao_num,N_states)


    short range exchange alpha/beta potentials with LDA functional on the |AO| basis

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

 
.. c:var:: potential_sr_x_beta_ao_pbe


    File : :file:`dft_utils_one_e/sr_pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_sr_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_sr_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange/correlation alpha/beta potentials with the short range PBE functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_sr_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: potential_x_alpha_ao


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao	(ao_num,ao_num,N_states)


    general providers for the alpha/beta exchange/correlation potentials on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`n_states`
       * :c:data:`potential_sr_c_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`potential_x_alpha_mo`

 
.. c:var:: potential_x_alpha_ao_lda


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_lda	(ao_num,ao_num,N_states)


    short range exchange/correlation alpha/beta potentials with LDA functional on the AO basis

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

 
.. c:var:: potential_x_alpha_ao_pbe


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange/correlation alpha/beta potentials with the short range PBE functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: potential_x_alpha_mo


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_x_beta_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_c_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_c_beta_mo	(mo_num,mo_num,N_states)


    general providers for the alpha/beta exchange/correlation potentials on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_x_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`trace_v_xc`

 
.. c:var:: potential_x_beta_ao


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao	(ao_num,ao_num,N_states)


    general providers for the alpha/beta exchange/correlation potentials on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`correlation_functional`
       * :c:data:`exchange_functional`
       * :c:data:`n_states`
       * :c:data:`potential_sr_c_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_lda`
       * :c:data:`potential_sr_x_alpha_ao_pbe`
       * :c:data:`potential_x_alpha_ao_lda`
       * :c:data:`potential_x_alpha_ao_pbe`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_potential_alpha_xc`
       * :c:data:`potential_x_alpha_mo`

 
.. c:var:: potential_x_beta_ao_lda


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_lda	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_lda	(ao_num,ao_num,N_states)


    short range exchange/correlation alpha/beta potentials with LDA functional on the AO basis

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

 
.. c:var:: potential_x_beta_ao_pbe


    File : :file:`dft_utils_one_e/pot_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_x_beta_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_alpha_ao_pbe	(ao_num,ao_num,N_states)
        double precision, allocatable	:: potential_c_beta_ao_pbe	(ao_num,ao_num,N_states)


    exchange/correlation alpha/beta potentials with the short range PBE functional on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_vc_alpha_pbe_w`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`potential_x_alpha_ao`

 
.. c:var:: potential_x_beta_mo


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: potential_x_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_x_beta_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_c_alpha_mo	(mo_num,mo_num,N_states)
        double precision, allocatable	:: potential_c_beta_mo	(mo_num,mo_num,N_states)


    general providers for the alpha/beta exchange/correlation potentials on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`potential_x_alpha_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`trace_v_xc`

 
.. c:var:: psi_dft_energy_h_core


    File : :file:`dft_utils_one_e/one_e_energy_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_dft_energy_kinetic	(N_states)
        double precision, allocatable	:: psi_dft_energy_nuclear_elec	(N_states)
        double precision, allocatable	:: psi_dft_energy_h_core	(N_states)


    kinetic, electron-nuclear and total h_core energy computed with the density matrix one_e_dm_mo_beta_for_dft+one_e_dm_mo_alpha_for_dft

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`


 
.. c:var:: psi_dft_energy_kinetic


    File : :file:`dft_utils_one_e/one_e_energy_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_dft_energy_kinetic	(N_states)
        double precision, allocatable	:: psi_dft_energy_nuclear_elec	(N_states)
        double precision, allocatable	:: psi_dft_energy_h_core	(N_states)


    kinetic, electron-nuclear and total h_core energy computed with the density matrix one_e_dm_mo_beta_for_dft+one_e_dm_mo_alpha_for_dft

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`


 
.. c:var:: psi_dft_energy_nuclear_elec


    File : :file:`dft_utils_one_e/one_e_energy_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_dft_energy_kinetic	(N_states)
        double precision, allocatable	:: psi_dft_energy_nuclear_elec	(N_states)
        double precision, allocatable	:: psi_dft_energy_h_core	(N_states)


    kinetic, electron-nuclear and total h_core energy computed with the density matrix one_e_dm_mo_beta_for_dft+one_e_dm_mo_alpha_for_dft

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`


 
.. c:var:: shifting_constant


    File : :file:`dft_utils_one_e/shifted_potential.irp.f`

    .. code:: fortran

        double precision, allocatable	:: shifting_constant	(N_states)


    shifting_constant = (E_{Hxc} - <\Psi | V_{Hxc} | \Psi>) / N_elec
    constant to add to the potential in order to obtain the variational energy as
    the eigenvalue of the effective long-range Hamiltonian
    (see original paper of Levy PRL 113, 113002 (2014), equation (17) )

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`energy_x`
       * :c:data:`n_states`
       * :c:data:`short_range_hartree_operator`
       * :c:data:`trace_v_xc`


 
.. c:var:: short_range_hartree


    File : :file:`dft_utils_one_e/sr_coulomb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: short_range_hartree_operator	(mo_num,mo_num,N_states)
        double precision, allocatable	:: short_range_hartree	(N_states)


    short_range_Hartree_operator(i,j) = :math:`\int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}` 
    
    short_range_Hartree = :math:`1/2  \sum_{i,j} \rho_{ij} \mathtt{short_range_Hartree_operator}(i,j)` 
    
                        = :math:`1/2  \int dr \int r' \rho(r) \rho(r') W_{ee}^{sr}` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_average_mo_for_dft`
       * :c:data:`one_e_dm_mo_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`shifting_constant`
       * :c:data:`trace_v_xc`

 
.. c:var:: short_range_hartree_operator


    File : :file:`dft_utils_one_e/sr_coulomb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: short_range_hartree_operator	(mo_num,mo_num,N_states)
        double precision, allocatable	:: short_range_hartree	(N_states)


    short_range_Hartree_operator(i,j) = :math:`\int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}` 
    
    short_range_Hartree = :math:`1/2  \sum_{i,j} \rho_{ij} \mathtt{short_range_Hartree_operator}(i,j)` 
    
                        = :math:`1/2  \int dr \int r' \rho(r) \rho(r') W_{ee}^{sr}` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_average_mo_for_dft`
       * :c:data:`one_e_dm_mo_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_one_e_potential`
       * :c:data:`shifting_constant`
       * :c:data:`trace_v_xc`

 
.. c:var:: trace_v_h


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: trace_v_xc	(N_states)
        double precision, allocatable	:: trace_v_h	(N_states)
        double precision, allocatable	:: trace_v_hxc	(N_states)


    Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
    Trace_v_Hxc = \sum_{i,j} v^{H}_{ij} (rho_{ij}_\alpha + rho_{ij}_\beta)
    Trace_v_Hxc = \sum_{i,j} rho_{ij} v^{Hxc}_{ij}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`potential_x_alpha_mo`
       * :c:data:`short_range_hartree_operator`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`shifting_constant`

 
.. c:var:: trace_v_hxc


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: trace_v_xc	(N_states)
        double precision, allocatable	:: trace_v_h	(N_states)
        double precision, allocatable	:: trace_v_hxc	(N_states)


    Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
    Trace_v_Hxc = \sum_{i,j} v^{H}_{ij} (rho_{ij}_\alpha + rho_{ij}_\beta)
    Trace_v_Hxc = \sum_{i,j} rho_{ij} v^{Hxc}_{ij}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`potential_x_alpha_mo`
       * :c:data:`short_range_hartree_operator`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`shifting_constant`

 
.. c:var:: trace_v_xc


    File : :file:`dft_utils_one_e/pot_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: trace_v_xc	(N_states)
        double precision, allocatable	:: trace_v_h	(N_states)
        double precision, allocatable	:: trace_v_hxc	(N_states)


    Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
    Trace_v_Hxc = \sum_{i,j} v^{H}_{ij} (rho_{ij}_\alpha + rho_{ij}_\beta)
    Trace_v_Hxc = \sum_{i,j} rho_{ij} v^{Hxc}_{ij}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`potential_x_alpha_mo`
       * :c:data:`short_range_hartree_operator`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`shifting_constant`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: berf:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        function berf(a)



 
.. c:function:: dberfda:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        function dberfda(a)



 
.. c:function:: dpol:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function dpol(rs)



 
.. c:function:: dpold:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function dpold(rs)



 
.. c:function:: dpoldd:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function dpoldd(rs)



 
.. c:function:: ec_lda:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ec_lda(rho_a,rho_b,ec,vc_a,vc_b)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ec_pbe_only`
       * :c:func:`ec_pbe_sr`
       * :c:data:`energy_x_lda`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ecpw`

 
.. c:function:: ec_lda_sr:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ec_lda_sr(mu,rho_a,rho_b,ec,vc_a,vc_b)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:func:`ec_pbe_only`
       * :c:func:`ec_pbe_sr`
       * :c:data:`energy_sr_x_lda`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ecorrlr`
       * :c:func:`ecpw`
       * :c:func:`vcorrlr`

 
.. c:function:: ec_only_lda_sr:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ec_only_lda_sr(mu,rho_a,rho_b,ec)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ecorrlr`
       * :c:func:`ecpw`

 
.. c:function:: ec_pbe_only:


    File : :file:`dft_utils_one_e/exc_sr_pbe.irp.f`

    .. code:: fortran

        subroutine ec_pbe_only(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec)


    Short-range PBE correlation energy functional for erf interaction
    
    input : ==========
    
    mu = range separated parameter
    
    rhoc, rhoo = total density and spin density
    
    sigmacc    = square of the gradient of the total density
    
    sigmaco    = square of the gradient of the spin density
    
    sigmaoo    = scalar product between the gradient of the total density and the one of the spin density
    
    output: ==========
    
    ec         = correlation energy
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_lda`
       * :c:func:`ec_lda_sr`

 
.. c:function:: ec_pbe_sr:


    File : :file:`dft_utils_one_e/exc_sr_pbe.irp.f`

    .. code:: fortran

        subroutine ec_pbe_sr(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)


    Short-range PBE correlation energy functional for erf interaction
    
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
       * :c:func:`gga_type_functionals`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ec_lda`
       * :c:func:`ec_lda_sr`

 
.. c:function:: ecorrlr:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

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

 
.. c:function:: ecpw:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

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

 
.. c:function:: ex_lda:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ex_lda(rho_a,rho_b,ex,vx_a,vx_b)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_x_lda`

 
.. c:function:: ex_lda_sr:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine ex_lda_sr(mu,rho_a,rho_b,ex,vx_a,vx_b)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_sr_vc_alpha_lda_w`
       * :c:data:`aos_vc_alpha_lda_w`
       * :c:data:`energy_sr_x_lda`
       * :c:func:`ex_pbe_sr`
       * :c:func:`ex_pbe_sr_only`

 
.. c:function:: ex_pbe_sr:


    File : :file:`dft_utils_one_e/exc_sr_pbe.irp.f`

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
       * :c:func:`gga_type_functionals`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ex_lda_sr`

 
.. c:function:: ex_pbe_sr_only:


    File : :file:`dft_utils_one_e/exc_sr_pbe.irp.f`

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

 
.. c:function:: g0d:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function g0d(rs)



 
.. c:function:: g0dd:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function g0dd(rs)



 
.. c:function:: g0f:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function g0f(x)



 
.. c:function:: gpw:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        subroutine GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G,Gd,Gdd)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ecpw`

 
.. c:function:: grad_rho_ab_to_grad_rho_oc:


    File : :file:`dft_utils_one_e/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_o_2,grad_rho_c_2,grad_rho_o_c)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gga_sr_type_functionals`
       * :c:func:`gga_type_functionals`

 
.. c:function:: qrpa:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function Qrpa(x)



 
.. c:function:: qrpad:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function Qrpad(x)



 
.. c:function:: qrpadd:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

    .. code:: fortran

        double precision function Qrpadd(x)



 
.. c:function:: rho_ab_to_rho_oc:


    File : :file:`dft_utils_one_e/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gga_sr_type_functionals`
       * :c:func:`gga_type_functionals`

 
.. c:function:: rho_oc_to_rho_ab:


    File : :file:`dft_utils_one_e/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine rho_oc_to_rho_ab(rho_o,rho_c,rho_a,rho_b)



 
.. c:function:: v_grad_rho_oc_to_v_grad_rho_ab:


    File : :file:`dft_utils_one_e/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine v_grad_rho_oc_to_v_grad_rho_ab(v_grad_rho_o_2,v_grad_rho_c_2,v_grad_rho_o_c,v_grad_rho_a_2,v_grad_rho_b_2,v_grad_rho_a_b)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gga_sr_type_functionals`
       * :c:func:`gga_type_functionals`

 
.. c:function:: v_rho_ab_to_v_rho_oc:


    File : :file:`dft_utils_one_e/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine v_rho_ab_to_v_rho_oc(v_rho_a,v_rho_b,v_rho_o,v_rho_c)



 
.. c:function:: v_rho_oc_to_v_rho_ab:


    File : :file:`dft_utils_one_e/rho_ab_to_rho_tot.irp.f`

    .. code:: fortran

        subroutine v_rho_oc_to_v_rho_ab(v_rho_o,v_rho_c,v_rho_a,v_rho_b)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gga_sr_type_functionals`
       * :c:func:`gga_type_functionals`

 
.. c:function:: vcorrlr:


    File : :file:`dft_utils_one_e/exc_sr_lda.irp.f`

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

