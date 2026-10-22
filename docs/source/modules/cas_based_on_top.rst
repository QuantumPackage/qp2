.. _module_cas_based_on_top: 
 
.. program:: cas_based_on_top 
 
.. default-role:: option 
 
================
cas_based_on_top
================

On-top pair density of a |CAS|-based wave function.

For a wave function in which all the determinants belong to an active
space, this module computes the on-top pair density :math:`n_2(r, r)`
and the related quantities:

* :c:data:`total_cas_on_top_density`: :math:`n_2(r, r)` at each point
  of the Becke grid, with all the core, inactive and active
  contributions
* :c:func:`give_n2_cas`: :math:`\mu(r)`, :math:`f_\Psi` and
  :math:`n_2^\Psi` for a general |CAS| wave function
* the effective density matrices :c:data:`effective_alpha_dm`,
  :c:data:`effective_beta_dm` and :c:data:`effective_spin_dm`, with
  their gradients, used by the range separated |KS-DFT| calculation
  of the :ref:`module_kohn_sham_rs` module
* the one-electron densities of the core, inactive and active parts:
  :c:data:`core_density`, :c:data:`inact_density`,
  :c:data:`one_e_act_density_alpha` and the related providers

The programs :ref:`cas_based_on_top` and :ref:`cas_based_density`
are small examples of the use of these quantities.
 
 
 
Programs 
-------- 
 
 * :ref:`cas_based_density` 
 * :ref:`cas_based_on_top` 
 
Providers 
--------- 
 
.. c:var:: act_mos_in_r_array


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: act_mos_in_r_array	(n_act_orb,n_points_final_grid)
        double precision, allocatable	:: act_mos_in_r_array_transp	(n_points_final_grid,n_act_orb)


    all ACTIVE MOs on the grid points, arranged in two different ways

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_act_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_act_density_alpha`
       * :c:data:`one_e_act_density_beta`
       * :c:data:`total_cas_on_top_density`

 
.. c:var:: act_mos_in_r_array_transp


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: act_mos_in_r_array	(n_act_orb,n_points_final_grid)
        double precision, allocatable	:: act_mos_in_r_array_transp	(n_points_final_grid,n_act_orb)


    all ACTIVE MOs on the grid points, arranged in two different ways

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_act_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_act_density_alpha`
       * :c:data:`one_e_act_density_beta`
       * :c:data:`total_cas_on_top_density`

 
.. c:var:: average_on_top


    File : :file:`on_top_cas_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: average_on_top	(n_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`total_cas_on_top_density`


 
.. c:var:: core_density


    File : :file:`cas_dens_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: core_density	(n_points_final_grid)


    CORE part of the density for alpha/beta.
    
    WARNING :: IF YOU NEED THE TOTAL DENSITY COMING FROM THE CORE,
    
               YOU MUST MULTIPLY BY TWO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_mos_in_r_array`
       * :c:data:`n_core_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_cas_total_density`
       * :c:data:`total_cas_on_top_density`

 
.. c:var:: core_inact_act_mos_grad_in_r_array


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: core_inact_act_mos_grad_in_r_array	(3,n_core_inact_act_orb,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core_inact_act`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`grad_total_cas_on_top_density`

 
.. c:var:: core_inact_act_mos_in_r_array


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: core_inact_act_mos_in_r_array	(n_core_inact_act_orb,n_points_final_grid)
        double precision, allocatable	:: core_inact_act_mos_in_r_array_transp	(n_points_final_grid,n_core_inact_act_orb)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core_inact_act`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`grad_total_cas_on_top_density`

 
.. c:var:: core_inact_act_mos_in_r_array_transp


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: core_inact_act_mos_in_r_array	(n_core_inact_act_orb,n_points_final_grid)
        double precision, allocatable	:: core_inact_act_mos_in_r_array_transp	(n_points_final_grid,n_core_inact_act_orb)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core_inact_act`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`grad_total_cas_on_top_density`

 
.. c:var:: core_mos_in_r_array


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: core_mos_in_r_array	(n_core_orb,n_points_final_grid)
        double precision, allocatable	:: core_mos_in_r_array_transp	(n_points_final_grid,n_core_orb)


    all COREE  MOs on the grid points, arranged in two different ways

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_density`

 
.. c:var:: core_mos_in_r_array_transp


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: core_mos_in_r_array	(n_core_orb,n_points_final_grid)
        double precision, allocatable	:: core_mos_in_r_array_transp	(n_points_final_grid,n_core_orb)


    all COREE  MOs on the grid points, arranged in two different ways

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_core_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_density`

 
.. c:var:: effective_alpha_dm


    File : :file:`eff_spin_dens.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_alpha_dm	(n_points_final_grid,N_states)
        double precision, allocatable	:: effective_beta_dm	(n_points_final_grid,N_states)
        double precision, allocatable	:: grad_effective_alpha_dm	(3,n_points_final_grid,N_states)
        double precision, allocatable	:: grad_effective_beta_dm	(3,n_points_final_grid,N_states)


    effective_alpha_dm(r_i) = 1/2 * (effective_spin_dm(r_i) + n(r_i))
    effective_beta_dm(r_i) = 1/2 * (-effective_spin_dm(r_i) + n(r_i))

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`effective_spin_dm`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`
       * :c:data:`total_cas_on_top_density`


 
.. c:var:: effective_beta_dm


    File : :file:`eff_spin_dens.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_alpha_dm	(n_points_final_grid,N_states)
        double precision, allocatable	:: effective_beta_dm	(n_points_final_grid,N_states)
        double precision, allocatable	:: grad_effective_alpha_dm	(3,n_points_final_grid,N_states)
        double precision, allocatable	:: grad_effective_beta_dm	(3,n_points_final_grid,N_states)


    effective_alpha_dm(r_i) = 1/2 * (effective_spin_dm(r_i) + n(r_i))
    effective_beta_dm(r_i) = 1/2 * (-effective_spin_dm(r_i) + n(r_i))

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`effective_spin_dm`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`
       * :c:data:`total_cas_on_top_density`


 
.. c:var:: effective_spin_dm


    File : :file:`eff_spin_dens.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_spin_dm	(n_points_final_grid,N_states)
        double precision, allocatable	:: grad_effective_spin_dm	(3,n_points_final_grid,N_states)


    effective_spin_dm(r_i) = \sqrt( n(r)^2 - 4 * ontop(r) )
    effective spin density obtained from the total density and on-top pair density
    see equation (6) of Phys. Chem. Chem. Phys., 2015, 17, 22412--22422 | 22413

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`grad_total_cas_on_top_density`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`
       * :c:data:`total_cas_on_top_density`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_alpha_dm`

 
.. c:var:: grad_effective_alpha_dm


    File : :file:`eff_spin_dens.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_alpha_dm	(n_points_final_grid,N_states)
        double precision, allocatable	:: effective_beta_dm	(n_points_final_grid,N_states)
        double precision, allocatable	:: grad_effective_alpha_dm	(3,n_points_final_grid,N_states)
        double precision, allocatable	:: grad_effective_beta_dm	(3,n_points_final_grid,N_states)


    effective_alpha_dm(r_i) = 1/2 * (effective_spin_dm(r_i) + n(r_i))
    effective_beta_dm(r_i) = 1/2 * (-effective_spin_dm(r_i) + n(r_i))

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`effective_spin_dm`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`
       * :c:data:`total_cas_on_top_density`


 
.. c:var:: grad_effective_beta_dm


    File : :file:`eff_spin_dens.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_alpha_dm	(n_points_final_grid,N_states)
        double precision, allocatable	:: effective_beta_dm	(n_points_final_grid,N_states)
        double precision, allocatable	:: grad_effective_alpha_dm	(3,n_points_final_grid,N_states)
        double precision, allocatable	:: grad_effective_beta_dm	(3,n_points_final_grid,N_states)


    effective_alpha_dm(r_i) = 1/2 * (effective_spin_dm(r_i) + n(r_i))
    effective_beta_dm(r_i) = 1/2 * (-effective_spin_dm(r_i) + n(r_i))

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`effective_spin_dm`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`
       * :c:data:`total_cas_on_top_density`


 
.. c:var:: grad_effective_spin_dm


    File : :file:`eff_spin_dens.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_spin_dm	(n_points_final_grid,N_states)
        double precision, allocatable	:: grad_effective_spin_dm	(3,n_points_final_grid,N_states)


    effective_spin_dm(r_i) = \sqrt( n(r)^2 - 4 * ontop(r) )
    effective spin density obtained from the total density and on-top pair density
    see equation (6) of Phys. Chem. Chem. Phys., 2015, 17, 22412--22422 | 22413

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`grad_total_cas_on_top_density`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`
       * :c:data:`total_cas_on_top_density`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_alpha_dm`

 
.. c:var:: grad_total_cas_on_top_density


    File : :file:`on_top_grad.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grad_total_cas_on_top_density	(4,n_points_final_grid,N_states)
        double precision	:: wall_time_core_inact_act_on_top_of_r	


    grad_total_cas_on_top_density(1:3,ipoint,istate) : provider for the on top pair density gradient (x,y,z) for the point 'ipoint' and state 'istate'
    
    grad_total_cas_on_top_density(4,ipoint,istate)   : on top pair density for the point 'ipoint' and state 'istate'

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_inact_act_mos_grad_in_r_array`
       * :c:data:`core_inact_act_mos_in_r_array`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_spin_dm`

 
.. c:var:: inact_density


    File : :file:`cas_dens_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: inact_density	(n_points_final_grid)


    INACTIVE part of the density for alpha/beta.
    
    WARNING :: IF YOU NEED THE TOTAL DENSITY COMING FROM THE INACTIVE,
    
               YOU MUST MULTIPLY BY TWO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`inact_mos_in_r_array`
       * :c:data:`n_inact_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_cas_total_density`
       * :c:data:`total_cas_on_top_density`

 
.. c:var:: inact_mos_in_r_array


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: inact_mos_in_r_array	(n_inact_orb,n_points_final_grid)
        double precision, allocatable	:: inact_mos_in_r_array_transp	(n_points_final_grid,n_inact_orb)


    all INACTIVE MOs on the grid points, arranged in two different ways

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_inact_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`inact_density`

 
.. c:var:: inact_mos_in_r_array_transp


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: inact_mos_in_r_array	(n_inact_orb,n_points_final_grid)
        double precision, allocatable	:: inact_mos_in_r_array_transp	(n_points_final_grid,n_inact_orb)


    all INACTIVE MOs on the grid points, arranged in two different ways

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_inact_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`inact_density`

 
.. c:var:: one_e_act_density_alpha


    File : :file:`cas_dens_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_act_density_alpha	(n_points_final_grid,N_states)


    one_e_act_density_alpha = pure ACTIVE part of the DENSITY for ALPHA ELECTRONS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_mos_in_r_array`
       * :c:data:`n_act_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_act_dm_alpha_mo_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_cas_total_density`
       * :c:data:`total_cas_on_top_density`

 
.. c:var:: one_e_act_density_beta


    File : :file:`cas_dens_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_act_density_beta	(n_points_final_grid,N_states)


    one_e_act_density_beta = pure ACTIVE part of the DENSITY for BETA ELECTRONS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_mos_in_r_array`
       * :c:data:`n_act_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_act_dm_beta_mo_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_cas_total_density`
       * :c:data:`total_cas_on_top_density`

 
.. c:var:: one_e_act_dm_alpha_mo_for_dft


    File : :file:`cas_one_e_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_act_dm_alpha_mo_for_dft	(n_act_orb,n_act_orb,N_states)


    one_e_act_dm_alpha_mo_for_dft = pure ACTIVE part of the ONE ELECTRON REDUCED DENSITY MATRIX for the ALPHA ELECTRONS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_alpha_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_act_density_alpha`

 
.. c:var:: one_e_act_dm_beta_mo_for_dft


    File : :file:`cas_one_e_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_act_dm_beta_mo_for_dft	(n_act_orb,n_act_orb,N_states)


    one_e_act_dm_beta_mo_for_dft = pure ACTIVE part of the ONE ELECTRON REDUCED DENSITY MATRIX for the BETA ELECTRONS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_mo_beta_for_dft`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_act_density_beta`

 
.. c:var:: one_e_cas_total_density


    File : :file:`cas_dens_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: one_e_cas_total_density	(n_points_final_grid,N_states)


    one_e_cas_total_density = TOTAL DENSITY FOR a CAS wave function
    
    WARNING : if "no_core_density" == .True. then the core part of density is ignored

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_density`
       * :c:data:`inact_density`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`no_core_density`
       * :c:data:`one_e_act_density_alpha`
       * :c:data:`one_e_act_density_beta`


 
.. c:var:: total_cas_on_top_density


    File : :file:`on_top_cas_prov.irp.f`

    .. code:: fortran

        double precision, allocatable	:: total_cas_on_top_density	(n_points_final_grid,N_states)


    on top pair density :: n2(r,r) at each of the Becke's grid point of a CAS-BASED wf
    
    Contains all core/inact/act contribution.
    
    !!!!! WARNING !!!!! If no_core_density then you REMOVE ALL CONTRIBUTIONS COMING FROM THE CORE ORBITALS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_mos_in_r_array`
       * :c:data:`core_density`
       * :c:data:`inact_density`
       * :c:data:`n_act_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`no_core_density`
       * :c:data:`one_e_act_density_alpha`
       * :c:data:`one_e_act_density_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`average_on_top`
       * :c:data:`effective_alpha_dm`
       * :c:data:`effective_spin_dm`

 
.. c:var:: virt_mos_in_r_array


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: virt_mos_in_r_array	(n_virt_orb,n_points_final_grid)
        double precision, allocatable	:: virt_mos_in_r_array_transp	(n_points_final_grid,n_virt_orb)


    all VIRTUAL MOs on the grid points, arranged in two different ways

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_virt`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_virt_orb`


 
.. c:var:: virt_mos_in_r_array_transp


    File : :file:`c_i_a_v_mos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: virt_mos_in_r_array	(n_virt_orb,n_points_final_grid)
        double precision, allocatable	:: virt_mos_in_r_array_transp	(n_points_final_grid,n_virt_orb)


    all VIRTUAL MOs on the grid points, arranged in two different ways

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_virt`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_virt_orb`


 
.. c:var:: wall_time_core_inact_act_on_top_of_r


    File : :file:`on_top_grad.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grad_total_cas_on_top_density	(4,n_points_final_grid,N_states)
        double precision	:: wall_time_core_inact_act_on_top_of_r	


    grad_total_cas_on_top_density(1:3,ipoint,istate) : provider for the on top pair density gradient (x,y,z) for the point 'ipoint' and state 'istate'
    
    grad_total_cas_on_top_density(4,ipoint,istate)   : on top pair density for the point 'ipoint' and state 'istate'

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_inact_act_mos_grad_in_r_array`
       * :c:data:`core_inact_act_mos_in_r_array`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_spin_dm`

 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: act_on_top_on_grid_pt:


    File : :file:`on_top_cas_prov.irp.f`

    .. code:: fortran

        subroutine act_on_top_on_grid_pt(ipoint,istate,pure_act_on_top_of_r)


    act_on_top_on_grid_pt returns the purely ACTIVE part of the on top pair density
    
    at the grid point ipoint, for the state istate

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_mos_in_r_array`
       * :c:data:`n_act_orb`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`total_cas_on_top_density`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`dger`

 
.. c:macro:: give_active_on_top_in_r_one_state:


    File : :file:`on_top_cas_rout.irp.f`

    .. code:: fortran

        subroutine give_active_on_top_in_r_one_state(r,istate,mos_array,act_on_top)


    gives the purely active on-top pair density for a given state

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_on_top_in_r_one_state`

 
.. c:macro:: give_cas_density_in_r:


    File : :file:`cas_dens_rout.irp.f`

    .. code:: fortran

        subroutine give_cas_density_in_r(core_dens,inact_dens,act_dens,total_cas_dens,r)


    returns the different component of the density at grid point r(3) for a CAS wave function
    
    core_dens      : density coming from the CORE orbitals
    
    inact_dens     : density coming from the INACT orbitals
    
    act_dens(1/2,1:N_states)  : active part of the alpha/beta electrons for all states
    
    total_cas_dens : total density of the cas wave function
    
    WARNING : if "no_core_density" == .True. then the core part of density is ignored in total_cas_dens

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_states`
       * :c:data:`no_core_density`
       * :c:data:`one_e_act_dm_alpha_mo_for_dft`
       * :c:data:`one_e_act_dm_beta_mo_for_dft`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`routine_test_cas_based_density`
       * :c:func:`write_on_top_in_real_space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:macro:: give_core_inact_act_density_in_r:


    File : :file:`on_top_cas_rout.irp.f`

    .. code:: fortran

        subroutine give_core_inact_act_density_in_r(r,mos_array,core_density_in_r,inact_density_in_r,act_density_in_r, total_density)


    core, inactive and active part of the density for alpha/beta electrons
    
    the density coming from the core and inactive are the same for alpha/beta electrons
    
    act_density(1/2, i) = alpha/beta density for the ith state
    
    total_density(i) = 2 * (core_density_in_r+inact_density_in_r) + act_density_in_r(1,i) + act_density_in_r(2,i)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_states`
       * :c:data:`one_e_act_dm_alpha_mo_for_dft`
       * :c:data:`one_e_act_dm_beta_mo_for_dft`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_on_top_in_r_one_state`

 
.. c:macro:: give_n2_aa_val_ab:


    File : :file:`two_body_dens_rout.irp.f`

    .. code:: fortran

        subroutine give_n2_aa_val_ab(r1,r2,two_bod_dens,istate)


    contribution from purely active orbitals to n2_{\Psi^B}(r_1,r_2) for the "istate" state of a CAS wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_n2_cas`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:macro:: give_n2_cas:


    File : :file:`two_body_dens_rout.irp.f`

    .. code:: fortran

        subroutine give_n2_cas(r1,r2,istate,n2_psi)


    returns mu(r), f_psi, n2_psi for a general cas wave function

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_n2_aa_val_ab`
       * :c:func:`give_n2_ia_val_ab`
       * :c:func:`give_n2_ii_val_ab`

 
.. c:macro:: give_n2_ia_val_ab:


    File : :file:`two_body_dens_rout.irp.f`

    .. code:: fortran

        subroutine give_n2_ia_val_ab(r1,r2,two_bod_dens,istate)


    contribution from inactive and active orbitals to n2_{\Psi^B}(r_1,r_2) for the "istate" state of a CAS wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_inact`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`one_e_act_dm_alpha_mo_for_dft`
       * :c:data:`one_e_act_dm_beta_mo_for_dft`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_n2_cas`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:macro:: give_n2_ii_val_ab:


    File : :file:`two_body_dens_rout.irp.f`

    .. code:: fortran

        subroutine give_n2_ii_val_ab(r1,r2,two_bod_dens)


    contribution from purely inactive orbitals to n2_{\Psi^B}(r_1,r_2) for a CAS wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`mo_num`
       * :c:data:`n_inact_orb`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_n2_cas`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:macro:: give_on_top_gradient:


    File : :file:`on_top_grad.irp.f`

    .. code:: fortran

        subroutine give_on_top_gradient(ipoint,istate,ontop_grad)


    on top pair density and its gradient evaluated at a given point of the grid
    ontop_grad(1:3) :: gradients of the on-top pair density
    ontop_grad(4)   :: on-top pair density

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_inact_act_mos_grad_in_r_array`
       * :c:data:`core_inact_act_mos_in_r_array`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`n_core_inact_act_orb`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`grad_total_cas_on_top_density`

 
.. c:macro:: give_on_top_in_r_one_state:


    File : :file:`on_top_cas_rout.irp.f`

    .. code:: fortran

        subroutine give_on_top_in_r_one_state(r,istate,on_top_in_r)


    on top pair density in r for the state istate a CAS-BASED wf
    
    note that if no_core_density .EQ. .True., all core contributions are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`no_core_density`
       * :c:data:`one_e_act_dm_alpha_mo_for_dft`
       * :c:data:`one_e_act_dm_beta_mo_for_dft`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`write_on_top_in_real_space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_active_on_top_in_r_one_state`
       * :c:func:`give_all_mos_at_r`
       * :c:func:`give_core_inact_act_density_in_r`

 
.. c:macro:: routine:


    File : :file:`cas_based_on_top.irp.f`

    .. code:: fortran

        subroutine routine



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`average_on_top`
       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`total_cas_on_top_density`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`cas_based_on_top`

 
.. c:macro:: routine_test_cas_based_density:


    File : :file:`cas_based_density.irp.f`

    .. code:: fortran

        subroutine routine_test_cas_based_density



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_cas_total_density`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`cas_based_density`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_cas_density_in_r`

 
.. c:macro:: write_on_top_in_real_space:


    File : :file:`example.irp.f`

    .. code:: fortran

        subroutine write_on_top_in_real_space


    This routines is a simple example of how to plot the on-top pair density on a simple 1D grid

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`nucl_coord`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_cas_density_in_r`
       * :c:func:`give_on_top_in_r_one_state`

