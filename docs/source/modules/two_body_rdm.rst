.. _module_two_body_rdm: 
 
.. program:: two_body_rdm 
 
.. default-role:: option 
 
============
two_body_rdm
============

Two-electron density matrices (2-RDMs) in the |MO| basis.

The 2-RDMs of the :math:`\alpha\alpha`, :math:`\beta\beta` and
:math:`\alpha\beta` electron pairs are stored as arrays, with the
physicists' notation, consistent with the two-electron integrals in
the |MO| basis. The state specific matrices are
:c:data:`full_occ_2_rdm_aa_mo`, :c:data:`full_occ_2_rdm_ab_mo` and
:c:data:`full_occ_2_rdm_bb_mo`, with the spin traced variant
:c:data:`full_occ_2_rdm_spin_trace_mo`. The versions restricted to
the active orbitals are :c:data:`act_2_rdm_aa_mo`,
:c:data:`act_2_rdm_ab_mo`, :c:data:`act_2_rdm_bb_mo` and
:c:data:`act_2_rdm_spin_trace_mo`.

The `io_two_body_rdm_*` subroutines read and write the 2-RDMs in the
|EZFIO| database.

The program :ref:`test_2_rdm` is a small example of the use of these
providers.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: io_two_body_rdm_ab
 
    Read/Write the active part of the two-body rdm for alpha/beta electrons from/to disk [ Write \| Read \| None ]
 
    Default: None
 
.. option:: io_two_body_rdm_aa
 
    Read/Write the active part of the two-body rdm for alpha/alpha electrons from/to disk [ Write \| Read \| None ]
 
    Default: None
 
.. option:: io_two_body_rdm_bb
 
    Read/Write the active part of the two-body rdm for beta/beta electrons from/to disk [ Write \| Read \| None ]
 
    Default: None
 
.. option:: io_two_body_rdm_spin_trace
 
    Read/Write the active part of the two-body rdm for spin trace electrons from/to disk [ Write \| Read \| None ]
 
    Default: None
 
 
Programs 
-------- 
 
 * :ref:`test_2_rdm` 
 
Providers 
--------- 
 
.. c:var:: double act_2_rdm_aa_mo


    File : :file:`act_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: act_2_rdm_aa_mo	(n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)


    act_2_rdm_aa_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of ALPHA/ALPHA electrons
    
    <Psi_{istate}\| a^{\dagger}_{i \alpha} a^{\dagger}_{j \alpha} a_{l \alpha} a_{k \alpha} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act"
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha}^{act} * (N_{\alpha}^{act} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`ezfio_filename`
       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`read_two_body_rdm_aa`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`full_occ_2_rdm_aa_mo`

 
.. c:var:: double act_2_rdm_ab_mo


    File : :file:`act_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: act_2_rdm_ab_mo	(n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)


    act_2_rdm_ab_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/beta+beta/alpha electrons
    
    ::
    
     +---------------------------+
     |              12 12        |
     | 1 2 1 2 == <ij|kl>        |
     +---------------------------+
    
    (12 12 labels the positions of i j k l in <ij|kl>)
    
      <Psi_{istate}\| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} \|Psi_{istate}>
    
    + <Psi_{istate}\| a^{\dagger}_{i \beta} a^{\dagger}_{j \alpha} a_{l \alpha} a_{k \beta} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act"
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha}^{act} * N_{\beta}^{act} * 2
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`ezfio_filename`
       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`read_two_body_rdm_ab`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`full_occ_2_rdm_ab_mo`

 
.. c:var:: double act_2_rdm_bb_mo


    File : :file:`act_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: act_2_rdm_bb_mo	(n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)


    act_2_rdm_bb_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of BETA/BETA electrons
    
    <Psi_{istate}\| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act"
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\beta}^{act} * (N_{\beta}^{act} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`ezfio_filename`
       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`read_two_body_rdm_bb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`full_occ_2_rdm_bb_mo`

 
.. c:var:: double act_2_rdm_spin_trace_mo


    File : :file:`act_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: act_2_rdm_spin_trace_mo	(n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)


    act_2_rdm_spin_trace_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM
    
    \sum_{\sigma,\sigma'}<Psi_{istate}\| a^{\dagger}_{i \sigma} a^{\dagger}_{j \sigma'} a_{l \sigma'} a_{k \sigma} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act"
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{elec}^{act} * (N_{elec}^{act} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`ezfio_filename`
       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`read_two_body_rdm_spin_trace`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`full_occ_2_rdm_spin_trace_mo`

 
.. c:var:: double act_2_rdm_trans_spin_trace_mo


    File : :file:`act_2_transition_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: act_2_rdm_trans_spin_trace_mo	(n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states,N_states)


    act_2_rdm_trans_spin_trace_mo(i,j,k,l,istate,jstate) =  STATE SPECIFIC physicist notation for 2rdm_trans
    
    \sum_{\sigma,\sigma'}<Psi_{istate}\| a^{\dagger}_{i \sigma} a^{\dagger}_{j \sigma'} a_{l \sigma'} a_{k \sigma} \|Psi_{jstate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act"
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{elec}^{act} * (N_{elec}^{act} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`


 
.. c:var:: double full_occ_2_rdm_aa_mo


    File : :file:`full_orb_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: full_occ_2_rdm_aa_mo	(n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)


    full_occ_2_rdm_aa_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/alpha electrons
    
    <Psi\| a^{\dagger}_{i \alpha} a^{\dagger}_{j \alpha} a_{l \alpha} a_{k \alpha} \|Psi>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO ALL OCCUPIED ORBITALS : core, inactive and active
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha} * (N_{\alpha} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
    
     !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_states`
       * :c:data:`no_core_density`
       * :c:data:`one_e_dm_mo_alpha`


 
.. c:var:: double full_occ_2_rdm_ab_mo


    File : :file:`full_orb_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: full_occ_2_rdm_ab_mo	(n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)


    full_occ_2_rdm_ab_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/beta + beta/alpha electrons
    
      <Psi\| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} \|Psi>
    
    + <Psi\| a^{\dagger}_{i \beta} a^{\dagger}_{j \alpha} a_{l \alpha} a_{k \beta} \|Psi>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO ALL OCCUPIED ORBITALS : core, inactive and active
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha} * N_{\beta} * 2
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
    
     !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO ARE SET TO ZERO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_states`
       * :c:data:`no_core_density`
       * :c:data:`one_e_dm_mo_alpha`


 
.. c:var:: double full_occ_2_rdm_bb_mo


    File : :file:`full_orb_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: full_occ_2_rdm_bb_mo	(n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)


    full_occ_2_rdm_bb_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of beta/beta electrons
    
    <Psi\| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} \|Psi>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO ALL OCCUPIED ORBITALS : core, inactive and active
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\beta} * (N_{\beta} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
    
     !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_states`
       * :c:data:`no_core_density`
       * :c:data:`one_e_dm_mo_alpha`


 
.. c:var:: double full_occ_2_rdm_spin_trace_mo


    File : :file:`full_orb_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: full_occ_2_rdm_spin_trace_mo	(n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)


    full_occ_2_rdm_bb_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of beta/beta electrons
    
    <Psi\| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} \|Psi>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO ALL OCCUPIED ORBITALS : core, inactive and active
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{elec} * (N_{elec} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
    
    !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero
    The two-electron energy of each state can be computed as:
    
      \sum_{i,j,k,l = 1, n_core_inact_act_orb} full_occ_2_rdm_spin_trace_mo(i,j,k,l,istate) * 1/2 * < ii jj \| kk ll >
    
      with ii = list_core_inact_act(i), jj = list_core_inact_act(j), kk = list_core_inact_act(k), ll = list_core_inact_act(l)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_states`
       * :c:data:`no_core_density`
       * :c:data:`one_e_dm_mo_alpha`


 
.. c:var:: double state_av_act_2_rdm_aa_mo


    File : :file:`state_av_act_2rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: state_av_act_2_rdm_aa_mo	(n_act_orb,n_act_orb,n_act_orb,n_act_orb)


    state_av_act_2_rdm_aa_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for alpha-alpha electron pairs
    
                                        = \sum_{istate} w(istate) * <Psi_{istate}\| a^{\dagger}_{i,alpha} a^{\dagger}_{j,alpha} a_{l,alpha} a_{k,alpha} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act"
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha}^{act} * (N_{\alpha}^{act} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`state_av_act_2_rdm_spin_trace_mo`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`

 
.. c:var:: double state_av_act_2_rdm_ab_mo


    File : :file:`state_av_act_2rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: state_av_act_2_rdm_ab_mo	(n_act_orb,n_act_orb,n_act_orb,n_act_orb)


    state_av_act_2_rdm_ab_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for alpha/beta+beta/alpha electron pairs
    
                                        = \sum_{istate} w(istate) * <Psi_{istate}\| a^{\dagger}_{i,alpha} a^{\dagger}_{j,beta} a_{l,beta} a_{k,alpha} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act"
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha}^{act} * N_{\beta}^{act} * 2
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`state_av_act_2_rdm_spin_trace_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`

 
.. c:var:: double state_av_act_2_rdm_bb_mo


    File : :file:`state_av_act_2rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: state_av_act_2_rdm_bb_mo	(n_act_orb,n_act_orb,n_act_orb,n_act_orb)


    state_av_act_2_rdm_bb_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for beta-beta electron pairs
    
                                        = \sum_{istate} w(istate) * <Psi_{istate}\| a^{\dagger}_{i,beta} a^{\dagger}_{j,beta} a_{l,beta} a_{k,beta} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act"
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\beta}^{act} * (N_{\beta}^{act} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_coef`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`state_av_act_2_rdm_spin_trace_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`

 
.. c:var:: double state_av_act_2_rdm_spin_trace_mo


    File : :file:`state_av_act_2rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: state_av_act_2_rdm_spin_trace_mo	(n_act_orb,n_act_orb,n_act_orb,n_act_orb)


    state_av_act_2_rdm_spin_trace_mo(i,j,k,l) = state average physicist spin trace two-body rdm restricted to the ACTIVE indices
    
                                        = \sum_{istate} w(istate) * \sum_{sigma,sigma'} <Psi_{istate}\| a^{\dagger}_{i,sigma} a^{\dagger'}_{j,sigma} a_{l,sigma'} a_{k,sigma} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act"
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{elec} * (N_{elec} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_act_orb`
       * :c:data:`n_states`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`

 
.. c:var:: double state_av_full_occ_2_rdm_aa_mo


    File : :file:`state_av_full_orb_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: state_av_full_occ_2_rdm_aa_mo	(n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb)


    state_av_full_occ_2_rdm_aa_mo(i,j,k,l) =  STATE AVERAGE physicist notation for 2RDM of alpha/alpha electrons
    
                                        = \sum_{istate} w(istate) * <Psi_{istate}\| a^{\dagger}_{i,alpha} a^{\dagger}_{j,alpha} a_{l,alpha} a_{k,alpha} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONG TO ALL OCCUPIED ORBITALS : core, inactive and active
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha} * (N_{\alpha} - 1)
    
     !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
    
     !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`no_core_density`
       * :c:data:`one_e_dm_mo_alpha_average`
       * :c:data:`state_av_act_2_rdm_aa_mo`


 
.. c:var:: double state_av_full_occ_2_rdm_ab_mo


    File : :file:`state_av_full_orb_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: state_av_full_occ_2_rdm_ab_mo	(n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb)


    state_av_full_occ_2_rdm_ab_mo(i,j,k,l) =  STATE AVERAGE physicist notation for 2RDM of alpha/beta + beta/alpha electrons
    
                                        = \sum_{istate} w(istate) * <Psi_{istate}\| a^{\dagger}_{i,alpha} a^{\dagger}_{j,beta} a_{l,beta} a_{k,alpha} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONG TO ALL OCCUPIED ORBITALS : core, inactive and active
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha} * N_{\beta} * 2
    
     !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
    
     !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`no_core_density`
       * :c:data:`one_e_dm_mo_alpha_average`
       * :c:data:`state_av_act_2_rdm_ab_mo`


 
.. c:var:: double state_av_full_occ_2_rdm_bb_mo


    File : :file:`state_av_full_orb_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: state_av_full_occ_2_rdm_bb_mo	(n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb)


    state_av_full_occ_2_rdm_bb_mo(i,j,k,l) =  STATE AVERAGE physicist notation for 2RDM of beta/beta electrons
    
                                        = \sum_{istate} w(istate) * <Psi_{istate}\| a^{\dagger}_{i,beta} a^{\dagger}_{j,beta} a_{l,beta} a_{k,beta} \|Psi_{istate}>
    
    WHERE ALL ORBITALS (i,j,k,l) BELONG TO ALL OCCUPIED ORBITALS : core, inactive and active
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\beta} * (N_{\beta} - 1)
    
    !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
    
     !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`no_core_density`
       * :c:data:`one_e_dm_mo_alpha_average`
       * :c:data:`state_av_act_2_rdm_bb_mo`


 
.. c:var:: double state_av_full_occ_2_rdm_spin_trace_mo


    File : :file:`state_av_full_orb_2_rdm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: state_av_full_occ_2_rdm_spin_trace_mo	(n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb)


    state_av_full_occ_2_rdm_bb_mo(i,j,k,l) =  STATE AVERAGE physicist notation for 2RDM of beta/beta electrons
    
                                        = \sum_{istate} w(istate) * \sum_{sigma,sigma'} <Psi_{istate}\| a^{\dagger}_{i,sigma} a^{\dagger'}_{j,sigma} a_{l,sigma'} a_{k,sigma} \|Psi_{istate}>
    
    
    WHERE ALL ORBITALS (i,j,k,l) BELONG TO ALL OCCUPIED ORBITALS : core, inactive and active
    
    THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{elec} * (N_{elec} - 1)
    
     !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
    
     !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`no_core_density`
       * :c:data:`one_e_dm_mo_alpha_average`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`two_e_dm_mo`

 
.. c:var:: double two_e_dm_mo


    File : :file:`two_e_dm_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: two_e_dm_mo	(mo_num,mo_num,mo_num,mo_num)


    \sum_{\sigma \sigma'}
    <Psi\| a^{\dagger}_{i \sigma} a^{\dagger}_{j \sigma'} a_{l \sigma'} a_{k \sigma} \|Psi>
    
    where the indices (i,j,k,l) belong to all MOs.
    
    The normalization (i.e. sum of diagonal elements) is set to :math:`N_{elec} * (N_{elec} - 1)/2` 
    
    !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO are set to zero
    The state-averaged two-electron energy :
    
      \sum_{i,j,k,l = 1, mo_num} two_e_dm_mo(i,j,k,l) * < kk ll \| ii jj >

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core_inact_act`
       * :c:data:`mo_num`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void read_array_two_rdm(void)


    File : :file:`io_two_rdm.irp.f`

    .. code:: fortran

        subroutine read_array_two_rdm(n_orb,nstates,array_tmp,name_file)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lock_io`
       * :c:func:`unlock_io`

 
.. c:function:: void read_array_two_trans_rdm(void)


    File : :file:`io_two_rdm.irp.f`

    .. code:: fortran

        subroutine read_array_two_trans_rdm(n_orb,nstates,array_tmp,name_file)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lock_io`
       * :c:func:`unlock_io`

 
.. c:function:: void routine_active_only(void)


    File : :file:`example.irp.f`

    .. code:: fortran

        subroutine routine_active_only


    This routine computes the two electron repulsion within the active space using various providers
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`list_act`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_act_orb`
       * :c:data:`n_states`
       * :c:data:`psi_energy_two_e`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`
       * :c:data:`state_average_weight`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_2_rdm`

 
.. c:function:: void routine_active_only_trans(void)


    File : :file:`example.irp.f`

    .. code:: fortran

        subroutine routine_active_only_trans


    This routine computes the two electron repulsion within the active space using various providers
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_trans_spin_trace_mo`
       * :c:data:`list_act`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_act_orb`
       * :c:data:`n_states`
       * :c:data:`psi_energy_two_e_trans`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_2_rdm`

 
.. c:function:: void routine_full_mos(void)


    File : :file:`example.irp.f`

    .. code:: fortran

        subroutine routine_full_mos


    This routine computes the two electron repulsion using various providers
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_states`
       * :c:data:`psi_energy_two_e`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`
       * :c:data:`state_average_weight`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_2_rdm`

 
.. c:function:: void write_array_two_rdm(void)


    File : :file:`io_two_rdm.irp.f`

    .. code:: fortran

        subroutine write_array_two_rdm(n_orb,nstates,array_tmp,name_file)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lock_io`
       * :c:func:`unlock_io`

 
.. c:function:: void write_array_two_trans_rdm(void)


    File : :file:`io_two_rdm.irp.f`

    .. code:: fortran

        subroutine write_array_two_trans_rdm(n_orb,nstates,array_tmp,name_file)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lock_io`
       * :c:func:`unlock_io`

