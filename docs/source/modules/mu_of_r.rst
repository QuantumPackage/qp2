.. _module_mu_of_r: 
 
.. program:: mu_of_r 
 
.. default-role:: option 
 
=======
mu_of_r
=======

The on-top pair density :math:`\mu(r)` of the wave function.

For a |HF| wave function, :c:data:`mu_of_r_hf` is :math:`\mu(r)`, as
in Eq. (37) of :cite:`Giner_2018` for a |HF| reference. For a wave
function expanded in an active space, :c:data:`mu_of_r_psi_cas` is
the same quantity for a |CAS| wave function, and
:c:func:`give_mu_of_r_cas` returns :math:`\mu(r)`, :math:`f_\Psi`
and :math:`n_2^\Psi`.

The Fock space components of the two-body density matrix that enter
the expression of :math:`\mu(r)`, with the Cholesky decomposition of
the two-electron integrals, are provided by
:c:data:`f_hf_cholesky`, :c:data:`mos_times_cholesky_r1` and
:c:data:`mos_times_cholesky_r2`, and the :c:data:`f_psi_hf_ab` /
:c:data:`f_psi_cas_ab` variants.

The `read_mu_of_r` and `write_mu_of_r` routines store and restore
:math:`\mu(r)` in the |EZFIO| database, and the program
:ref:`test_proj_op` is a small example.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: mu_of_r_disk
 
    array of the values of mu(r)
 
 
.. option:: mu_of_r_potential
 
    type of potential for the mu(r) interaction: can be [ hf| hf_sparse | cas_full | cas_truncated | pure_act]
 
    Default: hf_sparse
 
.. option:: io_mu_of_r
 
    Read/Write the mu(r) from/to disk [ Write | Read | None ]
 
    Default: None
 
 
Programs 
-------- 
 
 * :ref:`test_proj_op` 
 
Providers 
--------- 
 
.. c:var:: double basis_mos_in_r_array


    File : :file:`basis_def.irp.f`

    .. code:: fortran

        double precision, allocatable	:: basis_mos_in_r_array	(n_basis_orb,n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_basis`
       * :c:data:`mos_in_r_array`
       * :c:data:`n_basis_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`full_occ_v_kl_cntrctd`

 
.. c:var:: double core_inact_act_v_kl_contracted


    File : :file:`f_psi_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: core_inact_act_v_kl_contracted	(n_core_inact_act_orb,n_core_inact_act_orb,n_points_final_grid)


    core_inact_act_V_kl_contracted(k,l,ipoint) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
    
    This is needed to build the function f_{\Psi^B}(X_1,X_2) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`full_occ_v_kl_cntrctd`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab_old`

 
.. c:var:: double f_hf_cholesky


    File : :file:`f_hf_cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: f_hf_cholesky	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mos_times_cholesky_r1`
       * :c:data:`mos_times_cholesky_r2`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_hf`

 
.. c:var:: double f_hf_cholesky_sparse


    File : :file:`f_hf_cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: f_hf_cholesky_sparse	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cholesky_threshold`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`n_basis_orb`
       * :c:data:`n_occ_val_orb_for_hf`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_hf_sparse`

 
.. c:var:: double f_hf_cholesky_sparse_bis


    File : :file:`f_hf_cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: f_hf_cholesky_sparse_bis	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cholesky_threshold`
       * :c:data:`ao_num`
       * :c:data:`aos_in_r_array`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`n_basis_orb`
       * :c:data:`n_occ_val_orb_for_hf`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double f_psi_cas_ab


    File : :file:`f_val_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: f_psi_cas_ab	(n_points_final_grid,N_states)
        double precision, allocatable	:: on_top_cas_mu_r	(n_points_final_grid,N_states)


    
    Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid and for all states and for a CAS wave function
    
    Assumes that the wave function in psi_det is developped within an active space defined
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`final_grid_points`
       * :c:data:`list_act`
       * :c:data:`list_basis`
       * :c:data:`list_inact`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_basis_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_act_dm_alpha_mo_for_dft`
       * :c:data:`one_e_act_dm_beta_mo_for_dft`
       * :c:data:`two_e_int_aa_f`
       * :c:data:`two_e_int_ia_f`
       * :c:data:`two_e_int_ii_f`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_average_prov_ot`
       * :c:data:`mu_of_r_psi_cas`

 
.. c:var:: double f_psi_cas_ab_old


    File : :file:`f_psi_old.irp.f`

    .. code:: fortran

        double precision, allocatable	:: f_psi_cas_ab_old	(n_points_final_grid,N_states)


    
    Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid and for all states
    
    Assumes that the wave function in psi_det is developped within an active space defined
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_inact_act_v_kl_contracted`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_cntrctd`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`


 
.. c:var:: double f_psi_hf_ab


    File : :file:`f_val_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: f_psi_hf_ab	(n_points_final_grid)
        double precision, allocatable	:: on_top_hf_mu_r	(n_points_final_grid)


    
    Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid for a HF wave function
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mo_num`
       * :c:data:`n_basis_orb`
       * :c:data:`n_occ_val_orb_for_hf`
       * :c:data:`n_points_final_grid`
       * :c:data:`two_e_int_hf_f`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_hf_old`

 
.. c:var:: double full_occ_2_rdm_cntrctd


    File : :file:`f_psi_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: full_occ_2_rdm_cntrctd	(n_core_inact_act_orb,n_core_inact_act_orb,n_points_final_grid,N_states)


    full_occ_2_rdm_cntrctd(k,l,ipoint,istate) = \sum_{ij} \Gamma_{ij}^{kl}  phi_i(r_ipoint) phi_j(r_ipoint)
    
    where \Gamma_{ij}^{kl}(istate)  = <Psi_{istate}\| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} \|Psi_{istate}>
    
    This is needed to build the function f_{\Psi^B}(X_1,X_2) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`full_occ_2_rdm_cntrctd_trans`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab_old`

 
.. c:var:: double full_occ_2_rdm_cntrctd_trans


    File : :file:`f_psi_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: full_occ_2_rdm_cntrctd_trans	(n_points_final_grid,n_core_inact_act_orb,n_core_inact_act_orb,N_states)


    full_occ_2_rdm_cntrctd_trans(ipoint,k,l,istate) = \sum_{ij} \Gamma_{ij}^{kl}  phi_i(r_ipoint) phi_j(r_ipoint)
    
    where \Gamma_{ij}^{kl}(istate)  = <Psi_{istate}\| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} \|Psi_{istate}>
    
    This is needed to build the function f_{\Psi^B}(X_1,X_2) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_inact_act_mos_in_r_array`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`full_occ_2_rdm_cntrctd`

 
.. c:var:: double full_occ_v_kl_cntrctd


    File : :file:`f_psi_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: full_occ_v_kl_cntrctd	(n_points_final_grid,n_core_inact_act_orb,n_core_inact_act_orb)


    full_occ_v_kl_cntrctd(ipoint,k,l) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
    
    This is needed to build the function f_{\Psi^B}(X_1,X_2) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018)
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`basis_mos_in_r_array`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`gpu_num`
       * :c:data:`list_basis`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`n_basis_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`nproc_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_inact_act_v_kl_contracted`

 
.. c:var:: int list_basis


    File : :file:`basis_def.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_basis	(n_basis_orb)


    Defines the set of orbitals you will use to explore the basis set
    
    This determines the space \mathcal{B} of Eqs. (15-16) of Phys.Chem.Lett.2019, 10, 2931   2937
    
    It corresponds to all MOs except those defined as "deleted"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_all_but_del_orb`
       * :c:data:`mu_of_r_potential`
       * :c:data:`n_act_orb`
       * :c:data:`n_all_but_del_orb`
       * :c:data:`n_basis_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`basis_mos_in_r_array`
       * :c:data:`f_hf_cholesky_sparse`
       * :c:data:`f_hf_cholesky_sparse_bis`
       * :c:data:`f_psi_cas_ab`
       * :c:data:`f_psi_hf_ab`
       * :c:data:`full_occ_v_kl_cntrctd`
       * :c:data:`list_couple_hf_orb_r1`
       * :c:data:`list_couple_hf_orb_r2`
       * :c:data:`mos_times_cholesky_r1`
       * :c:data:`mos_times_cholesky_r2`
       * :c:data:`two_e_int_aa_f`
       * :c:data:`two_e_int_hf_f`
       * :c:data:`two_e_int_ia_f`
       * :c:data:`two_e_int_ii_f`

 
.. c:var:: int list_couple_hf_orb_r1


    File : :file:`f_hf_cholesky.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_couple_hf_orb_r1	(2,n_couple_orb_r1)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`n_basis_orb`
       * :c:data:`n_couple_orb_r1`
       * :c:data:`n_occ_val_orb_for_hf`


 
.. c:var:: int list_couple_hf_orb_r2


    File : :file:`f_hf_cholesky.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_couple_hf_orb_r2	(2,n_couple_orb_r2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`n_basis_orb`
       * :c:data:`n_couple_orb_r2`
       * :c:data:`n_occ_val_orb_for_hf`


 
.. c:var:: int list_valence_orb_for_hf


    File : :file:`basis_def.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_valence_orb_for_hf	(n_max_occ_val_orb_for_hf,2)


    List of OCCUPIED valence orbitals for each spin to build the f_{HF}(r_1,r_2) function
    
    This corresponds to ALL OCCUPIED orbitals in the HF wave function, except those defined as "core"
    
    This determines the space \mathcal{A} of Eqs. (15-16) of Phys.Chem.Lett.2019, 10, 2931   2937

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_class`
       * :c:data:`n_occ_val_orb_for_hf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_hf_cholesky_sparse`
       * :c:data:`f_hf_cholesky_sparse_bis`
       * :c:data:`f_psi_hf_ab`
       * :c:data:`list_couple_hf_orb_r1`
       * :c:data:`list_couple_hf_orb_r2`
       * :c:data:`mos_times_cholesky_r1`
       * :c:data:`mos_times_cholesky_r2`
       * :c:data:`on_top_hf_grid`
       * :c:data:`two_e_int_hf_f`

 
.. c:var:: double mos_times_cholesky_r1


    File : :file:`f_hf_cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_times_cholesky_r1	(cholesky_mo_num,n_points_final_grid)


    V1_AR = \sum_{I}V_AI Phi_IR where "R" specifies the index of the grid point and A the number of cholesky point
    
    here Phi_IR is phi_i(R)xphi_b(R) for r1 and V_AI = (ib|A) chollesky vector

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`n_basis_orb`
       * :c:data:`n_couple_orb_r1`
       * :c:data:`n_occ_val_orb_for_hf`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_hf_cholesky`

 
.. c:var:: double mos_times_cholesky_r2


    File : :file:`f_hf_cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mos_times_cholesky_r2	(cholesky_mo_num,n_points_final_grid)


    V1_AR = \sum_{I}V_AI Phi_IR where "R" specifies the index of the grid point and A the number of cholesky point
    
    here Phi_IR is phi_i(R)xphi_b(R) for r2 and V_AI = (ib|A) chollesky vector

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`n_basis_orb`
       * :c:data:`n_couple_orb_r2`
       * :c:data:`n_occ_val_orb_for_hf`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_hf_cholesky`

 
.. c:var:: double mu_average_proj


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_average_proj	(N_states)


    average value of mu(r) weighted with the total one-e density and divided by the number of electrons
    
    !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
    
    in the one- and two-body density matrix are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_projector_mo`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double mu_average_prov


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_average_prov	(N_states)


    average value of mu(r) weighted with the total one-e density and divided by the number of electrons
    
    !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
    
    in the one- and two-body density matrix are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double mu_average_prov2


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_average_prov2	(N_states)


    average value of mu(r) weighted with square of the total one-e density
    
    !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
    
    in the one- and two-body density matrix are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double mu_average_prov3


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_average_prov3	(N_states)


    average value of mu(r) weighted with square of the total one-e density
    
    !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
    
    in the one- and two-body density matrix are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double mu_average_prov_ot


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_average_prov_ot	(N_states)


    average value of mu(r) weighted with on-top pair density
    
    !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
    
    in the one- and two-body density matrix are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`final_grid_points`
       * :c:data:`mu_of_r_prov`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`


 
.. c:var:: double mu_of_r_hf


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_of_r_hf	(n_points_final_grid)


    mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
    
    corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
    
    !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
    
    in the two-body density matrix are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_hf_cholesky`
       * :c:data:`n_points_final_grid`
       * :c:data:`on_top_hf_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_prov`

 
.. c:var:: double mu_of_r_hf_old


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_of_r_hf_old	(n_points_final_grid)


    mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
    
    corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
    
    !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
    
    in the two-body density matrix are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`f_psi_hf_ab`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_prov`

 
.. c:var:: double mu_of_r_hf_sparse


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_of_r_hf_sparse	(n_points_final_grid)


    mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
    
    corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
    
    !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
    
    in the two-body density matrix are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_hf_cholesky_sparse`
       * :c:data:`n_points_final_grid`
       * :c:data:`on_top_hf_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_prov`

 
.. c:var:: double mu_of_r_projector_ao_prod


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_of_r_projector_ao_prod	(n_points_final_grid)


    mu(r) computed with the projector onto the atomic basis
     P_B(\mathbf{r},\mathbf{r}') = \sum_{ij} |
     \chi_{i} \rangle \left[S^{-1}\right]_{ij} \langle \chi_{j} |
     \] where :math:`i`  and :math:`j`  denote all atomic orbitals.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`aos_in_r_array`
       * :c:data:`basis_prod_idx`
       * :c:data:`basis_prod_num`
       * :c:data:`basis_prod_overlap`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_prov`

 
.. c:var:: double mu_of_r_projector_mo


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_of_r_projector_mo	(n_points_final_grid)


    mu(r) computed with the projector onto the atomic basis
     P_B(\mathbf{r},\mathbf{r}') = \sum_{ij} |
     \chi_{i} \rangle \left[S^{-1}\right]_{ij} \langle \chi_{j} |
     \] where :math:`i`  and :math:`j`  denote all atomic orbitals.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact_act`
       * :c:data:`list_virt`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`n_inact_act_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_average_proj`

 
.. c:var:: double mu_of_r_prov


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_of_r_prov	(n_points_final_grid,N_states)


    general variable for mu(r)
    
    corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018)
    
    !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
    
    in the two-body density matrix are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_hf`
       * :c:data:`mu_of_r_hf_old`
       * :c:data:`mu_of_r_hf_sparse`
       * :c:data:`mu_of_r_potential`
       * :c:data:`mu_of_r_projector_ao_prod`
       * :c:data:`mu_of_r_psi_cas`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`read_mu_of_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_average_prov`
       * :c:data:`mu_average_prov2`
       * :c:data:`mu_average_prov3`
       * :c:data:`mu_average_prov_ot`

 
.. c:var:: double mu_of_r_psi_cas


    File : :file:`mu_of_r_conditions.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mu_of_r_psi_cas	(n_points_final_grid,N_states)


    mu(r) computed with a wave function developped in an active space
    
    corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018)
    
    !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
    
    in the one- and two-body density matrix are excluded

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_prov`

 
.. c:var:: int n_basis_orb


    File : :file:`basis_def.irp.f`

    .. code:: fortran

        integer	:: n_basis_orb	


    Defines the number of orbitals you will use to explore the basis set
    
    This determines the size of the space \mathcal{B} of Eqs. (15-16) of Phys.Chem.Lett.2019, 10, 2931   2937
    
    It corresponds to all MOs except those defined as "deleted"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_potential`
       * :c:data:`n_act_orb`
       * :c:data:`n_all_but_del_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`basis_mos_in_r_array`
       * :c:data:`f_hf_cholesky_sparse`
       * :c:data:`f_hf_cholesky_sparse_bis`
       * :c:data:`f_psi_cas_ab`
       * :c:data:`f_psi_hf_ab`
       * :c:data:`full_occ_v_kl_cntrctd`
       * :c:data:`list_basis`
       * :c:data:`list_couple_hf_orb_r1`
       * :c:data:`list_couple_hf_orb_r2`
       * :c:data:`mos_times_cholesky_r1`
       * :c:data:`mos_times_cholesky_r2`
       * :c:data:`n_couple_orb_r1`
       * :c:data:`n_couple_orb_r2`
       * :c:data:`two_e_int_aa_f`
       * :c:data:`two_e_int_hf_f`
       * :c:data:`two_e_int_ia_f`
       * :c:data:`two_e_int_ii_f`

 
.. c:var:: int n_couple_orb_r1


    File : :file:`f_hf_cholesky.irp.f`

    .. code:: fortran

        integer	:: n_couple_orb_r1	


    number of couples of alpha occupied times any basis orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_basis_orb`
       * :c:data:`n_occ_val_orb_for_hf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_couple_hf_orb_r1`
       * :c:data:`mos_times_cholesky_r1`

 
.. c:var:: int n_couple_orb_r2


    File : :file:`f_hf_cholesky.irp.f`

    .. code:: fortran

        integer	:: n_couple_orb_r2	


    number of couples of beta occupied times any basis orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_basis_orb`
       * :c:data:`n_occ_val_orb_for_hf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_couple_hf_orb_r2`
       * :c:data:`mos_times_cholesky_r2`

 
.. c:var:: int n_max_occ_val_orb_for_hf


    File : :file:`basis_def.irp.f`

    .. code:: fortran

        integer, allocatable	:: n_occ_val_orb_for_hf	(2)
        integer	:: n_max_occ_val_orb_for_hf	


    Number of OCCUPIED VALENCE ORBITALS for each spin to build the f_{HF}(r_1,r_2) function
    
    This is typically elec_alpha_num - n_core_orb for alpha electrons and elec_beta_num - n_core_orb for beta electrons
    
    This determines the size of the space \mathcal{A} of Eqs. (15-16) of Phys.Chem.Lett.2019, 10, 2931   2937

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_class`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_hf_cholesky_sparse`
       * :c:data:`f_hf_cholesky_sparse_bis`
       * :c:data:`f_psi_hf_ab`
       * :c:data:`list_couple_hf_orb_r1`
       * :c:data:`list_couple_hf_orb_r2`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mos_times_cholesky_r1`
       * :c:data:`mos_times_cholesky_r2`
       * :c:data:`n_couple_orb_r1`
       * :c:data:`n_couple_orb_r2`
       * :c:data:`on_top_hf_grid`
       * :c:data:`two_e_int_hf_f`

 
.. c:var:: int n_occ_val_orb_for_hf


    File : :file:`basis_def.irp.f`

    .. code:: fortran

        integer, allocatable	:: n_occ_val_orb_for_hf	(2)
        integer	:: n_max_occ_val_orb_for_hf	


    Number of OCCUPIED VALENCE ORBITALS for each spin to build the f_{HF}(r_1,r_2) function
    
    This is typically elec_alpha_num - n_core_orb for alpha electrons and elec_beta_num - n_core_orb for beta electrons
    
    This determines the size of the space \mathcal{A} of Eqs. (15-16) of Phys.Chem.Lett.2019, 10, 2931   2937

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_class`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_hf_cholesky_sparse`
       * :c:data:`f_hf_cholesky_sparse_bis`
       * :c:data:`f_psi_hf_ab`
       * :c:data:`list_couple_hf_orb_r1`
       * :c:data:`list_couple_hf_orb_r2`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mos_times_cholesky_r1`
       * :c:data:`mos_times_cholesky_r2`
       * :c:data:`n_couple_orb_r1`
       * :c:data:`n_couple_orb_r2`
       * :c:data:`on_top_hf_grid`
       * :c:data:`two_e_int_hf_f`

 
.. c:var:: double on_top_cas_mu_r


    File : :file:`f_val_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: f_psi_cas_ab	(n_points_final_grid,N_states)
        double precision, allocatable	:: on_top_cas_mu_r	(n_points_final_grid,N_states)


    
    Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid and for all states and for a CAS wave function
    
    Assumes that the wave function in psi_det is developped within an active space defined
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`final_grid_points`
       * :c:data:`list_act`
       * :c:data:`list_basis`
       * :c:data:`list_inact`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_basis_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_states`
       * :c:data:`one_e_act_dm_alpha_mo_for_dft`
       * :c:data:`one_e_act_dm_beta_mo_for_dft`
       * :c:data:`two_e_int_aa_f`
       * :c:data:`two_e_int_ia_f`
       * :c:data:`two_e_int_ii_f`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_average_prov_ot`
       * :c:data:`mu_of_r_psi_cas`

 
.. c:var:: double on_top_hf_grid


    File : :file:`f_hf_cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: on_top_hf_grid	(n_points_final_grid)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`n_occ_val_orb_for_hf`
       * :c:data:`n_points_final_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_hf`
       * :c:data:`mu_of_r_hf_sparse`

 
.. c:var:: double on_top_hf_mu_r


    File : :file:`f_val_general.irp.f`

    .. code:: fortran

        double precision, allocatable	:: f_psi_hf_ab	(n_points_final_grid)
        double precision, allocatable	:: on_top_hf_mu_r	(n_points_final_grid)


    
    Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid for a HF wave function
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mo_num`
       * :c:data:`n_basis_orb`
       * :c:data:`n_occ_val_orb_for_hf`
       * :c:data:`n_points_final_grid`
       * :c:data:`two_e_int_hf_f`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_hf_old`

 
.. c:var:: double two_e_int_aa_f


    File : :file:`f_psi_i_a_v_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_int_aa_f	(n_basis_orb,n_basis_orb,n_act_orb,n_act_orb)


    list of two-electron integrals (built with the MOs belonging to the \mathcal{B} space)
    
    needed to compute the function f_{ii}(r_1,r_2)
    
    two_e_int_aa_f(j,i,n,m) = < j i | n m > where all orbitals belong to "list_basis"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`banned_excitation`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`gpu_num`
       * :c:data:`list_act`
       * :c:data:`list_basis`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_act_orb`
       * :c:data:`n_basis_orb`
       * :c:data:`nproc_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`

 
.. c:var:: double two_e_int_hf_f


    File : :file:`f_hf_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_int_hf_f	(n_basis_orb,n_basis_orb,n_max_occ_val_orb_for_hf,n_max_occ_val_orb_for_hf)


    list of two-electron integrals (built with the MOs belonging to the \mathcal{B} space)
    
    needed to compute the function f_{HF}(r_1,r_2)
    
    two_e_int_hf_f(j,i,n,m) = < j i | n m > where all orbitals belong to "list_basis"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`banned_excitation`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_basis_orb`
       * :c:data:`n_occ_val_orb_for_hf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_hf_ab`

 
.. c:var:: double two_e_int_ia_f


    File : :file:`f_psi_i_a_v_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_int_ia_f	(n_basis_orb,n_basis_orb,n_inact_orb,n_act_orb)


    list of two-electron integrals (built with the MOs belonging to the \mathcal{B} space)
    
    needed to compute the function f_{ia}(r_1,r_2)
    
    two_e_int_aa_f(j,i,n,m) = < j i | n m > where all orbitals belong to "list_basis"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`gpu_num`
       * :c:data:`list_act`
       * :c:data:`list_basis`
       * :c:data:`list_inact`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_basis_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`nproc_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`

 
.. c:var:: double two_e_int_ii_f


    File : :file:`f_psi_i_a_v_utils.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_int_ii_f	(n_basis_orb,n_basis_orb,n_inact_orb,n_inact_orb)


    list of two-electron integrals (built with the MOs belonging to the \mathcal{B} space)
    
    needed to compute the function f_{ii}(r_1,r_2)
    
    two_e_int_ii_f(j,i,n,m) = < j i | n m > where all orbitals belong to "list_basis"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`gpu_num`
       * :c:data:`list_basis`
       * :c:data:`list_inact`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`n_basis_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`nproc_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void f_hf_valence_ab(void)


    File : :file:`f_hf_utils.irp.f`

    .. code:: fortran

        subroutine f_HF_valence_ab(r1,r2,f_HF_val_ab,two_bod_dens)


    f_HF_val_ab(r1,r2) = function f_{\Psi^B}(X_1,X_2) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018)
    
    for alpha beta spins and an HF wave function and excluding the "core" orbitals (see Eq. 16a of Phys.Chem.Lett.2019, 10, 2931   2937)
    
    two_bod_dens = on-top pair density of the HF wave function
    
    < HF | wee_{\alpha\beta} | HF > =  \int (r1,r2) f_HF_ab(r1,r2) excluding all contributions from "core" "electrons"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mo_num`
       * :c:data:`n_basis_orb`
       * :c:data:`n_occ_val_orb_for_hf`
       * :c:data:`two_e_int_hf_f`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_hf_ab`
       * :c:func:`test_f_hf_valence_ab`
       * :c:func:`test_f_ii_valence_ab`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:function:: void give_f_aa_val_ab(void)


    File : :file:`f_psi_i_a_v_utils.irp.f`

    .. code:: fortran

        subroutine give_f_aa_val_ab(r1,r2,f_aa_val_ab,two_bod_dens,istate)


    contribution from purely active orbitals to f_{\Psi^B}(r_1,r_2) for the "istate" state of a CAS wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`list_act`
       * :c:data:`list_basis`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_basis_orb`
       * :c:data:`two_e_int_aa_f`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:func:`give_mu_of_r_cas`
       * :c:func:`test_f_ia_valence_ab`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`give_all_mos_at_r`

 
.. c:function:: void give_f_ia_val_ab(void)


    File : :file:`f_psi_i_a_v_utils.irp.f`

    .. code:: fortran

        subroutine give_f_ia_val_ab(r1,r2,f_ia_val_ab,two_bod_dens,istate)


    contribution from inactive and active orbitals to f_{\Psi^B}(r_1,r_2) for the "istate" state of a CAS wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_basis`
       * :c:data:`list_inact`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_basis_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`one_e_act_dm_alpha_mo_for_dft`
       * :c:data:`one_e_act_dm_beta_mo_for_dft`
       * :c:data:`two_e_int_ia_f`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:func:`give_mu_of_r_cas`
       * :c:func:`test_f_ia_valence_ab`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:function:: void give_f_ii_val_ab(void)


    File : :file:`f_psi_i_a_v_utils.irp.f`

    .. code:: fortran

        subroutine give_f_ii_val_ab(r1,r2,f_ii_val_ab,two_bod_dens)


    contribution from purely inactive orbitals to f_{\Psi^B}(r_1,r_2) for a CAS wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_basis`
       * :c:data:`list_inact`
       * :c:data:`mo_num`
       * :c:data:`n_basis_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`two_e_int_ii_f`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:func:`give_mu_of_r_cas`
       * :c:func:`test_f_ia_valence_ab`
       * :c:func:`test_f_ii_valence_ab`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:function:: void give_mu_of_r_cas(void)


    File : :file:`f_psi_i_a_v_utils.irp.f`

    .. code:: fortran

        subroutine give_mu_of_r_cas(r,istate,mu_of_r,f_psi,n2_psi)


    returns mu(r), f_psi, n2_psi for a general cas wave function

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_f_aa_val_ab`
       * :c:func:`give_f_ia_val_ab`
       * :c:func:`give_f_ii_val_ab`

 
.. c:function:: void integral_f_hf_valence_ab(void)


    File : :file:`f_hf_utils.irp.f`

    .. code:: fortran

        subroutine integral_f_HF_valence_ab(r1,int_f_HF_val_ab)


    in_f_HF_val_ab(r_1) = \int dr_2 f_{\Psi^B}(r_1,r_2)
    
    where f_{\Psi^B}(r_1,r_2) is defined by Eq. (22) of J. Chem. Phys. 149, 194301 (2018)
    
    for alpha beta spins and an HF wave function and excluding the "core" orbitals (see Eq. 16a of Phys.Chem.Lett.2019, 10, 2931   2937)
    
    Such function can be used to test if the f_HF_val_ab(r_1,r_2) is correctly built.
    
    < HF | wee_{\alpha\beta} | HF > =  \int (r1) int_f_HF_val_ab(r_1)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_basis`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mo_num`
       * :c:data:`n_basis_orb`
       * :c:data:`n_occ_val_orb_for_hf`
       * :c:data:`two_e_int_hf_f`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_f_hf_valence_ab`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:function:: void test(void)


    File : :file:`test_proj_op.irp.f`

    .. code:: fortran

        subroutine test



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_hf_cholesky_sparse`
       * :c:data:`f_hf_cholesky_sparse_bis`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_proj_op`

 
.. c:function:: void test_f_hf_valence_ab(void)


    File : :file:`example.irp.f`

    .. code:: fortran

        subroutine test_f_HF_valence_ab


    routine to test the function f_HF(r1,r2)
    
    the integral over r1,r2 should be equal to the alpha/beta interaction of HF determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`list_valence_orb_for_hf`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_occ_val_orb_for_hf`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_proj_op`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`f_hf_valence_ab`
       * :c:func:`integral_f_hf_valence_ab`

 
.. c:function:: void test_f_ia_valence_ab(void)


    File : :file:`example.irp.f`

    .. code:: fortran

        subroutine test_f_ia_valence_ab


    routine to test the function f_ii(r1,r2), f_ia(r1,r2) and f_aa(r1,r2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab_old`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`total_cas_on_top_density`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_proj_op`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_f_aa_val_ab`
       * :c:func:`give_f_ia_val_ab`
       * :c:func:`give_f_ii_val_ab`

 
.. c:function:: void test_f_ii_ia_aa_valence_ab(void)


    File : :file:`example.irp.f`

    .. code:: fortran

        subroutine test_f_ii_ia_aa_valence_ab


    routine to test the function f_Psi(r1,r2) based on core/inactive/active orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`f_psi_cas_ab`
       * :c:data:`f_psi_cas_ab_old`
       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`
       * :c:data:`total_cas_on_top_density`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_proj_op`

 
.. c:function:: void test_f_ii_valence_ab(void)


    File : :file:`example.irp.f`

    .. code:: fortran

        subroutine test_f_ii_valence_ab


    routine to test the function f_ii(r1,r2)
    
    it should be the same that f_HF(r1,r2) only for inactive orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_proj_op`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`f_hf_valence_ab`
       * :c:func:`give_f_ii_val_ab`

