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
 
.. option:: no_core_density
 
    if [no_core_dm] then all elements of the density matrix involving at least one orbital set as core are set to zero
 
    Default: full_density
 
.. option:: normalize_dm
 
    if .True., then you normalize the no_core_dm to elec_alpha_num - n_core_orb  and elec_beta_num - n_core_orb
 
    Default: True
 
 
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
       * :c:data:`data_one_e_dm_alpha_ao`
       * :c:data:`data_one_e_dm_beta_ao`
       * :c:data:`density_for_dft`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`


 
.. c:var:: one_e_dm_alpha_ao_for_dft_no_core


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_alpha_ao_for_dft_no_core	(ao_num,ao_num,N_states)
        double precision, allocatable	:: one_e_dm_beta_ao_for_dft_no_core	(ao_num,ao_num,N_states)


    one body density matrix on the AO basis based on one_e_dm_mo_alpha_for_dft_no_core

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft_no_core`
       * :c:data:`one_e_dm_mo_beta_for_dft_no_core`


 
.. c:var:: one_e_dm_average_alpha_mo_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_average_alpha_mo_for_dft	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_average_mo_for_dft`

 
.. c:var:: one_e_dm_average_beta_mo_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_average_beta_mo_for_dft	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_average_mo_for_dft`

 
.. c:var:: one_e_dm_average_mo_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_average_mo_for_dft	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`one_e_dm_average_alpha_mo_for_dft`
       * :c:data:`one_e_dm_average_beta_mo_for_dft`


 
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
       * :c:data:`data_one_e_dm_alpha_ao`
       * :c:data:`data_one_e_dm_beta_ao`
       * :c:data:`density_for_dft`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft`


 
.. c:var:: one_e_dm_beta_ao_for_dft_no_core


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_alpha_ao_for_dft_no_core	(ao_num,ao_num,N_states)
        double precision, allocatable	:: one_e_dm_beta_ao_for_dft_no_core	(ao_num,ao_num,N_states)


    one body density matrix on the AO basis based on one_e_dm_mo_alpha_for_dft_no_core

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft_no_core`
       * :c:data:`one_e_dm_mo_beta_for_dft_no_core`


 
.. c:var:: one_e_dm_mo_alpha_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_alpha_for_dft	(mo_num,mo_num,N_states)


    density matrix for alpha electrons in the MO basis used for all DFT calculations based on the density

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`damping_for_rs_dft`
       * :c:data:`data_one_e_dm_alpha_mo`
       * :c:data:`density_for_dft`
       * :c:data:`elec_alpha_num`
       * :c:data:`list_core`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_states`
       * :c:data:`no_core_density`
       * :c:data:`normalize_dm`
       * :c:data:`one_body_dm_mo_alpha_one_det`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`one_e_dm_mo_alpha_average`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`one_e_dm_average_alpha_mo_for_dft`
       * :c:data:`one_e_dm_mo_alpha_for_dft_no_core`
       * :c:data:`one_e_dm_mo_for_dft`

 
.. c:var:: one_e_dm_mo_alpha_for_dft_no_core


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_alpha_for_dft_no_core	(mo_num,mo_num,N_states)


    density matrix for alpha electrons in the MO basis without the core orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_alpha_ao_for_dft_no_core`

 
.. c:var:: one_e_dm_mo_beta_for_dft


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_beta_for_dft	(mo_num,mo_num,N_states)


    density matrix for beta  electrons in the MO basis used for all DFT calculations based on the density

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`damping_for_rs_dft`
       * :c:data:`data_one_e_dm_beta_mo`
       * :c:data:`density_for_dft`
       * :c:data:`elec_beta_num`
       * :c:data:`list_core`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_states`
       * :c:data:`no_core_density`
       * :c:data:`normalize_dm`
       * :c:data:`one_body_dm_mo_alpha_one_det`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`one_e_dm_mo_alpha_average`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`one_e_dm_average_beta_mo_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft_no_core`
       * :c:data:`one_e_dm_mo_for_dft`

 
.. c:var:: one_e_dm_mo_beta_for_dft_no_core


    File : :file:`density_for_dft/density_for_dft.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_dm_mo_beta_for_dft_no_core	(mo_num,mo_num,N_states)


    density matrix for beta  electrons in the MO basis without the core orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_beta_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_alpha_ao_for_dft_no_core`

 
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


