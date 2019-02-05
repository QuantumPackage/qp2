.. _module_density_for_dft: 
 
.. program:: density_for_dft 
 
.. default-role:: option 
 
===============
density_for_dft
===============


This module defines the *provider* of the density used for the |DFT| related
calculations.  This definition is done through the keyword
:option:`density_for_dft density_for_dft`.  The density can be:

* `WFT`: the density is computed with a potentially multi determinant wave
  function (see variables `psi_det` and `psi_det`)# input_density: the density
  is set to a density previously stored in the |EZFIO| directory (see
  ``aux_quantities``)
* `damping_rs_dft`: the density is damped between the input_density and the WFT
  density, with a damping factor of :option:`density_for_dft damping_for_rs_dft`

 
 
 
EZFIO parameters 
---------------- 
 
.. option:: density_for_dft
 
    Type of density used for DFT calculation. If set to WFT , it uses the density of the wave function stored in (psi_det,psi_coef). If set to input_density it uses the one-body dm stored in aux_quantities/ . If set to damping_rs_dft it uses the damped density between WFT and input_density. In the ks_scf and rs_ks_scf programs, it is set to WFT.
 
    Default: WFT
 
.. option:: damping_for_rs_dft
 
    damping factor for the density used in RSFT.
 
    Default: 0.5
 
 
Providers 
--------- 
 
.. c:var:: one_body_dm_mo_alpha_one_det


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_body_dm_mo_alpha_one_det	(mo_num,mo_num,N_states)
        double precision, allocatable	:: one_body_dm_mo_beta_one_det	(mo_num,mo_num,N_states)


    One body density matrix on the |MO| basis for a single determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`

 
.. c:var:: one_body_dm_mo_beta_one_det


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_body_dm_mo_alpha_one_det	(mo_num,mo_num,N_states)
        double precision, allocatable	:: one_body_dm_mo_beta_one_det	(mo_num,mo_num,N_states)


    One body density matrix on the |MO| basis for a single determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`

 
.. c:var:: one_e_dm_alpha_ao_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_alpha_ao_for_dft	(ao_num,ao_num,N_states)
        double precision, allocatable	:: one_e_dm_beta_ao_for_dft	(ao_num,ao_num,N_states)


    one body density matrix on the AO basis based on one_e_dm_mo_alpha_for_dft

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_alpha_at_r`
       * :c:data:`one_e_dm_alpha_in_r`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

 
.. c:var:: one_e_dm_average_mo_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_average_mo_for_dft	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_for_dft`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`short_range_hartree_operator`

 
.. c:var:: one_e_dm_beta_ao_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_alpha_ao_for_dft	(ao_num,ao_num,N_states)
        double precision, allocatable	:: one_e_dm_beta_ao_for_dft	(ao_num,ao_num,N_states)


    one body density matrix on the AO basis based on one_e_dm_mo_alpha_for_dft

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_alpha_at_r`
       * :c:data:`one_e_dm_alpha_in_r`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

 
.. c:var:: one_e_dm_mo_alpha_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_alpha_for_dft	(mo_num,mo_num,N_states)


    density matrix for alpha electrons in the MO basis used for all DFT calculations based on the density

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`damping_for_rs_dft`
       * :c:data:`data_one_e_dm_alpha_mo`
       * :c:data:`density_for_dft`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_body_dm_mo_alpha_one_det`
       * :c:data:`one_e_dm_mo_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`one_e_dm_mo_for_dft`
       * :c:data:`psi_dft_energy_kinetic`
       * :c:data:`trace_v_xc`

 
.. c:var:: one_e_dm_mo_beta_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_beta_for_dft	(mo_num,mo_num,N_states)


    density matrix for beta  electrons in the MO basis used for all DFT calculations based on the density

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`damping_for_rs_dft`
       * :c:data:`data_one_e_dm_beta_mo`
       * :c:data:`density_for_dft`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_body_dm_mo_alpha_one_det`
       * :c:data:`one_e_dm_mo_alpha`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`one_e_dm_mo_for_dft`
       * :c:data:`psi_dft_energy_kinetic`
       * :c:data:`trace_v_xc`

 
.. c:var:: one_e_dm_mo_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_for_dft	(mo_num,mo_num,N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_average_mo_for_dft`
       * :c:data:`short_range_hartree_operator`

