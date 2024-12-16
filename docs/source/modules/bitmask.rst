.. _module_bitmask: 
 
.. program:: bitmask 
 
.. default-role:: option 
 
==============
bitmask module
==============

The central part of this module is the :file:`bitmasks_module.f90` file. It contains
the constants that will be used to define on which kind of integer the bitmasks
will be defined.

In the program, to represent a determinant as a pair of bitstrings,
the determinant should be defined as

.. code-block:: fortran

  use bitmasks
  integer(bit_kind)  :: determinant(N_int,2)


:file:`bitmasks_routines.irp.f` contains helper routines to manipulate bitmask, like
transforming a bit string to a list of integers for example.


`bit_kind_shift`, `bit_kind_size` and `bit_kind` are supposed to be consistent::

   2**bit_kind_shift = bit_kind_size
   bit_kind = bit_kind_size / 8


For an example of how to use the bitmaks, see the file :file:`example.irp.f`.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: n_act_orb
 
    Number of active |MOs|
 
 
.. option:: do_ormas
 
    if |true| restrict selection based on ORMAS rules
 
    Default: false
 
.. option:: ormas_n_space
 
    Number of active spaces
 
    Default: 1
 
.. option:: ormas_mstart
 
    starting orb for each ORMAS space
 
 
.. option:: ormas_min_e
 
    min number of electrons in each ORMAS space
 
 
.. option:: ormas_max_e
 
    max number of electrons in each ORMAS space
 
 
 
Providers 
--------- 
 
.. c:var:: act_bitmask


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)


    Bitmask identifying the active MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`closed_shell_ref_bitmask`
       * :c:data:`n_det_generators`
       * :c:data:`psi_cas`
       * :c:data:`psi_det_generators`
       * :c:data:`reunion_of_act_virt_bitmask`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_inact_act_bitmask`

 
.. c:var:: closed_shell_ref_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: closed_shell_ref_bitmask	(N_int,2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`


 
.. c:var:: core_bitmask


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)


    Bitmask identifying the core MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`inact_virt_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`

 
.. c:var:: core_inact_act_bitmask_4


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: core_inact_act_bitmask_4	(N_int,4)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_act_bitmask`


 
.. c:var:: core_inact_virt_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: inact_virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: core_inact_virt_bitmask	(N_int,2)


    Reunion of the inactive and virtual bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`inact_bitmask`
       * :c:data:`n_int`
       * :c:data:`virt_bitmask`


 
.. c:var:: del_bitmask


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    Bitmask identifying the deleted MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_del`
       * :c:data:`n_del_orb`
       * :c:data:`n_int`


 
.. c:var:: dim_list_act_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_act_orb	


    dimensions for the allocation of list_act.
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_act_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`

 
.. c:var:: dim_list_core_inact_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_core_inact_orb	


    dimensions for the allocation of list_core.
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_inact_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_core_inact`

 
.. c:var:: dim_list_core_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_core_orb	


    dimensions for the allocation of list_core.
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_core`

 
.. c:var:: dim_list_del_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_del_orb	


    dimensions for the allocation of list_del.
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_del_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_del`

 
.. c:var:: dim_list_inact_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_inact_orb	


    dimensions for the allocation of list_inact.
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_inact_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`

 
.. c:var:: dim_list_virt_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_virt_orb	


    dimensions for the allocation of list_virt.
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_virt`

 
.. c:var:: full_ijkl_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: full_ijkl_bitmask	(N_int)


    Bitmask to include all possible MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`full_ijkl_bitmask_4`
       * :c:data:`generators_bitmask`

 
.. c:var:: full_ijkl_bitmask_4


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: full_ijkl_bitmask_4	(N_int,4)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`full_ijkl_bitmask`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`

 
.. c:var:: generators_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: generators_bitmask	(N_int,2,6)


    Bitmasks for generator determinants.
    (N_int, alpha/beta, hole/particle, generator).
    
    3rd index is :
    
    * 1 : hole     for single exc
    
    * 2 : particle for single exc
    
    * 3 : hole     for 1st exc of double
    
    * 4 : particle for 1st exc of double
    
    * 5 : hole     for 2nd exc of double
    
    * 6 : particle for 2nd exc of double
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`full_ijkl_bitmask`
       * :c:data:`n_int`
       * :c:data:`reunion_of_act_virt_bitmask`
       * :c:data:`reunion_of_inact_act_bitmask`


 
.. c:var:: hf_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: hf_bitmask	(N_int,2)


    Hartree Fock bit mask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`double_exc_bitmask`
       * :c:data:`max_degree_exc`
       * :c:data:`psi_cas`
       * :c:data:`psi_det`
       * :c:data:`ref_bitmask`
       * :c:data:`single_exc_bitmask`
       * :c:data:`unpaired_alpha_electrons`

 
.. c:var:: inact_bitmask


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)


    Bitmask identifying the  inactive MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_inact_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`inact_virt_bitmask`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`reunion_of_inact_act_bitmask`

 
.. c:var:: inact_virt_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: inact_virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: core_inact_virt_bitmask	(N_int,2)


    Reunion of the inactive and virtual bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`inact_bitmask`
       * :c:data:`n_int`
       * :c:data:`virt_bitmask`


 
.. c:var:: index_holes_bitmask


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        integer, allocatable	:: index_holes_bitmask	(3)


    Index of the holes in the generators_bitmasks


 
.. c:var:: index_particl_bitmask


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        integer, allocatable	:: index_particl_bitmask	(3)


    Index of the holes in the generators_bitmasks


 
.. c:var:: list_act


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)


    List of MO indices which are in the active.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_act_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`act_2_rdm_trans_spin_trace_mo`
       * :c:data:`act_bitmask`
       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pqxx_no_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`bielec_pxxq_no_array`
       * :c:data:`bielecci`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`d0tu`
       * :c:data:`d0tu_alpha_ao`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`etwo`
       * :c:data:`excit`
       * :c:data:`fapq`
       * :c:data:`fock_matrix_mo`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`natorbsci_mos`
       * :c:data:`occnum`
       * :c:data:`one_ints_no`
       * :c:data:`p0tuvx_peter`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`
       * :c:data:`super_ci_dm`
       * :c:data:`umat`

 
.. c:var:: list_act_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)


    List of MO indices which are in the active.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_act_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`act_2_rdm_trans_spin_trace_mo`
       * :c:data:`act_bitmask`
       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pqxx_no_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`bielec_pxxq_no_array`
       * :c:data:`bielecci`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`d0tu`
       * :c:data:`d0tu_alpha_ao`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`etwo`
       * :c:data:`excit`
       * :c:data:`fapq`
       * :c:data:`fock_matrix_mo`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`natorbsci_mos`
       * :c:data:`occnum`
       * :c:data:`one_ints_no`
       * :c:data:`p0tuvx_peter`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`
       * :c:data:`super_ci_dm`
       * :c:data:`umat`

 
.. c:var:: list_all_but_del_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_all_but_del_orb	(n_all_but_del_orb)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_all_but_del_orb`


 
.. c:var:: list_core


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)


    List of MO indices which are in the core.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_alpha_for_dft_no_core`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft_no_core`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`

 
.. c:var:: list_core_inact


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_core_inact	(dim_list_core_inact_orb)
        integer, allocatable	:: list_core_inact_reverse	(mo_num)


    List of indices of the core and inactive MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_inact_orb`
       * :c:data:`mo_num`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`d0tu_alpha_ao`
       * :c:data:`etwo`
       * :c:data:`excit`
       * :c:data:`fipq`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`occnum`
       * :c:data:`super_ci_dm`
       * :c:data:`umat`

 
.. c:var:: list_core_inact_act


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_core_inact_act	(n_core_inact_act_orb)
        integer, allocatable	:: list_core_inact_act_reverse	(mo_num)


    List of indices of the core inactive and active MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_act_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`etwo`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`two_e_dm_mo`

 
.. c:var:: list_core_inact_act_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_core_inact_act	(n_core_inact_act_orb)
        integer, allocatable	:: list_core_inact_act_reverse	(mo_num)


    List of indices of the core inactive and active MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_act_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`etwo`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`two_e_dm_mo`

 
.. c:var:: list_core_inact_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_core_inact	(dim_list_core_inact_orb)
        integer, allocatable	:: list_core_inact_reverse	(mo_num)


    List of indices of the core and inactive MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_inact_orb`
       * :c:data:`mo_num`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`d0tu_alpha_ao`
       * :c:data:`etwo`
       * :c:data:`excit`
       * :c:data:`fipq`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`occnum`
       * :c:data:`super_ci_dm`
       * :c:data:`umat`

 
.. c:var:: list_core_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)


    List of MO indices which are in the core.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_alpha_for_dft_no_core`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft_no_core`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`

 
.. c:var:: list_del


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_del	(dim_list_del_orb)
        integer, allocatable	:: list_del_reverse	(mo_num)


    List of MO indices which are deleted.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_del_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_del_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`del_bitmask`

 
.. c:var:: list_del_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_del	(dim_list_del_orb)
        integer, allocatable	:: list_del_reverse	(mo_num)


    List of MO indices which are deleted.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_del_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_del_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`del_bitmask`

 
.. c:var:: list_inact


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)


    List of MO indices which are inactive.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_inact_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_inact_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`inact_bitmask`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`

 
.. c:var:: list_inact_act


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact_act	(n_inact_act_orb)
        integer, allocatable	:: list_inact_act_reverse	(mo_num)


    List of indices of the inactive and active MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_inact_act_orb`
       * :c:data:`n_int`
       * :c:data:`reunion_of_inact_act_bitmask`


 
.. c:var:: list_inact_act_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact_act	(n_inact_act_orb)
        integer, allocatable	:: list_inact_act_reverse	(mo_num)


    List of indices of the inactive and active MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_inact_act_orb`
       * :c:data:`n_int`
       * :c:data:`reunion_of_inact_act_bitmask`


 
.. c:var:: list_inact_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)


    List of MO indices which are inactive.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_inact_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_inact_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`inact_bitmask`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`

 
.. c:var:: list_virt


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_virt_reverse	(mo_num)


    List of MO indices which are virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_virt_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_no_total_transp`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`excit`
       * :c:data:`fock_matrix_mo`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`super_ci_dm`
       * :c:data:`umat`
       * :c:data:`virt_bitmask`

 
.. c:var:: list_virt_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_virt_reverse	(mo_num)


    List of MO indices which are virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_virt_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_no_total_transp`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`excit`
       * :c:data:`fock_matrix_mo`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`super_ci_dm`
       * :c:data:`umat`
       * :c:data:`virt_bitmask`

 
.. c:var:: mo_coef_begin_iteration


    File : :file:`bitmask/track_orb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_coef_begin_iteration	(ao_num,mo_num)


    Void provider to store the coefficients of the |MO| basis at the beginning of the SCF iteration
    
    Useful to track some orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_num`


 
.. c:var:: mpi_bit_kind


    File : :file:`bitmask/mpi.irp.f`

    .. code:: fortran

        integer	:: mpi_bit_kind	


    MPI bit kind type


 
.. c:var:: n_act_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_act_orb	


    Number of active MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`act_2_rdm_trans_spin_trace_mo`
       * :c:data:`act_bitmask`
       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pqxx_no_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`bielec_pxxq_no_array`
       * :c:data:`bielecci`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`d0tu`
       * :c:data:`d0tu_alpha_ao`
       * :c:data:`dim_list_act_orb`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`etwo`
       * :c:data:`excit`
       * :c:data:`fapq`
       * :c:data:`fock_matrix_mo`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`list_act`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`n_c_a_prov`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_inact_act_orb`
       * :c:data:`natorbsci`
       * :c:data:`natorbsci_mos`
       * :c:data:`nmonoex`
       * :c:data:`nsomomax`
       * :c:data:`occnum`
       * :c:data:`one_ints_no`
       * :c:data:`p0tuvx`
       * :c:data:`p0tuvx_no`
       * :c:data:`p0tuvx_peter`
       * :c:data:`state_av_act_2_rdm_aa_mo`
       * :c:data:`state_av_act_2_rdm_ab_mo`
       * :c:data:`state_av_act_2_rdm_bb_mo`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`
       * :c:data:`super_ci_dm`
       * :c:data:`umat`

 
.. c:var:: n_all_but_del_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_all_but_del_orb	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_all_but_del_orb`

 
.. c:var:: n_core_inact_act_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_core_inact_act_orb	


    Number of core inactive and active MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pqxx_no_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`bielec_pxxq_no_array`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`list_core_inact_act`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`
       * :c:data:`two_e_dm_mo`

 
.. c:var:: n_core_inact_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_core_inact_orb	


    n_core + n_inact

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pqxx_no_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`bielec_pxxq_no_array`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`d0tu_alpha_ao`
       * :c:data:`dim_list_core_inact_orb`
       * :c:data:`etwo`
       * :c:data:`excit`
       * :c:data:`fipq`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`list_core_inact`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`n_c_a_prov`
       * :c:data:`nmonoex`
       * :c:data:`occnum`
       * :c:data:`super_ci_dm`
       * :c:data:`umat`

 
.. c:var:: n_core_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_core_orb	


    Number of core MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`dim_list_core_orb`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`list_core`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`one_e_dm_mo_alpha_for_dft`
       * :c:data:`one_e_dm_mo_alpha_for_dft_no_core`
       * :c:data:`one_e_dm_mo_beta_for_dft`
       * :c:data:`one_e_dm_mo_beta_for_dft_no_core`
       * :c:data:`pt2_f`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`

 
.. c:var:: n_del_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_del_orb	


    Number of deleted MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`del_bitmask`
       * :c:data:`dim_list_del_orb`
       * :c:data:`list_del`

 
.. c:var:: n_inact_act_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_inact_act_orb	


    n_inact + n_act

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_act_orb`
       * :c:data:`n_inact_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact_act`

 
.. c:var:: n_inact_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_inact_orb	


    Number of inactive MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_inact_orb`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`full_occ_2_rdm_aa_mo`
       * :c:data:`full_occ_2_rdm_ab_mo`
       * :c:data:`full_occ_2_rdm_bb_mo`
       * :c:data:`full_occ_2_rdm_spin_trace_mo`
       * :c:data:`inact_bitmask`
       * :c:data:`list_inact`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_inact_act_orb`
       * :c:data:`state_av_full_occ_2_rdm_aa_mo`
       * :c:data:`state_av_full_occ_2_rdm_ab_mo`
       * :c:data:`state_av_full_occ_2_rdm_bb_mo`
       * :c:data:`state_av_full_occ_2_rdm_spin_trace_mo`

 
.. c:var:: n_int


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_int	


    Number of 64-bit integers needed to represent determinants as binary strings

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`cfg_seniority_index`
       * :c:data:`ci_electronic_energy`
       * :c:data:`closed_shell_ref_bitmask`
       * :c:data:`coef_hf_selector`
       * :c:data:`core_bitmask`
       * :c:data:`core_inact_act_bitmask_4`
       * :c:data:`del_bitmask`
       * :c:data:`det_to_configuration`
       * :c:data:`dettocsftransformationmatrix`
       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`dominant_dets_of_cfgs`
       * :c:data:`double_exc_bitmask`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`full_ijkl_bitmask`
       * :c:data:`full_ijkl_bitmask_4`
       * :c:data:`generators_bitmask`
       * :c:data:`global_selection_buffer`
       * :c:data:`gradvec_old`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`h_matrix_diag_all_dets`
       * :c:data:`hessmat_old`
       * :c:data:`hf_bitmask`
       * :c:data:`inact_bitmask`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`max_degree_exc`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_det_generators`
       * :c:data:`n_dominant_dets_of_cfgs`
       * :c:data:`n_elec_alpha_for_psi_configuration`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`one_e_tr_dm_mo`
       * :c:data:`one_e_tr_dm_mo_alpha`
       * :c:data:`orb_swap`
       * :c:data:`ormas_bitmask`
       * :c:data:`p0tuvx`
       * :c:data:`p0tuvx_peter`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_cas`
       * :c:data:`psi_cas_sorted_bit`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_sorted`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_csf_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_det_hii`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy_two_e`
       * :c:data:`psi_energy_two_e_trans`
       * :c:data:`psi_non_cas`
       * :c:data:`psi_non_cas_sorted_bit`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_diag_h_mat`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`
       * :c:data:`ref_closed_shell_bitmask`
       * :c:data:`reunion_of_act_virt_bitmask`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`reunion_of_inact_act_bitmask`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`s2_values`
       * :c:data:`single_exc_bitmask`
       * :c:data:`singles_alpha_csc`
       * :c:data:`singles_alpha_csc_idx`
       * :c:data:`singles_alpha_csc_map`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`
       * :c:data:`singles_beta_csc_map`
       * :c:data:`unpaired_alpha_electrons`
       * :c:data:`virt_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: n_virt_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_virt_orb	


    Number of virtual MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_no_total_transp`
       * :c:data:`dim_list_virt_orb`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`excit`
       * :c:data:`fock_matrix_mo`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`list_virt`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`n_c_a_prov`
       * :c:data:`nmonoex`
       * :c:data:`super_ci_dm`
       * :c:data:`umat`
       * :c:data:`virt_bitmask`

 
.. c:var:: ormas_bitmask


    File : :file:`bitmask/bitmasks_ormas.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: ormas_bitmask	(N_int,ormas_n_space)


    bitmask for each ormas space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`ormas_list_orb`
       * :c:data:`ormas_n_orb`
       * :c:data:`ormas_n_space`


 
.. c:var:: ormas_list_orb


    File : :file:`bitmask/bitmasks_ormas.irp.f`

    .. code:: fortran

        integer, allocatable	:: ormas_list_orb	(ormas_max_n_orb,ormas_n_space)


    list of orbitals in each ormas space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ormas_n_orb`
       * :c:data:`ormas_n_space`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ormas_bitmask`

 
.. c:var:: ormas_max_e


    File : :file:`bitmask/bitmasks_ormas.irp.f`

    .. code:: fortran

        integer, allocatable	:: ormas_max_e	(ormas_n_space)


    max nelec in each active space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`ormas_n_space`


 
.. c:var:: ormas_max_n_orb


    File : :file:`bitmask/bitmasks_ormas.irp.f`

    .. code:: fortran

        integer, allocatable	:: ormas_n_orb	(ormas_n_space)
        integer	:: ormas_max_n_orb	


    number of orbitals in each ormas space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`ormas_mstart`
       * :c:data:`ormas_n_space`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ormas_bitmask`
       * :c:data:`ormas_list_orb`

 
.. c:var:: ormas_min_e


    File : :file:`bitmask/bitmasks_ormas.irp.f`

    .. code:: fortran

        integer, allocatable	:: ormas_min_e	(ormas_n_space)


    min nelec in each active space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`ormas_n_space`


 
.. c:var:: ormas_mstart


    File : :file:`bitmask/bitmasks_ormas.irp.f`

    .. code:: fortran

        integer, allocatable	:: ormas_mstart	(ormas_n_space)


    first orbital idx in each active space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`ormas_n_space`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ormas_n_orb`

 
.. c:var:: ormas_n_orb


    File : :file:`bitmask/bitmasks_ormas.irp.f`

    .. code:: fortran

        integer, allocatable	:: ormas_n_orb	(ormas_n_space)
        integer	:: ormas_max_n_orb	


    number of orbitals in each ormas space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`ormas_mstart`
       * :c:data:`ormas_n_space`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ormas_bitmask`
       * :c:data:`ormas_list_orb`

 
.. c:var:: ref_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: ref_bitmask	(N_int,2)


    Reference bit mask, used in Slater rules, chosen as Hartree-Fock bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`closed_shell_ref_bitmask`
       * :c:data:`coef_hf_selector`
       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`psi_det_hii`
       * :c:data:`psi_selectors_diag_h_mat`
       * :c:data:`ref_bitmask_energy`
       * :c:data:`ref_closed_shell_bitmask`

 
.. c:var:: reunion_of_act_virt_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: reunion_of_act_virt_bitmask	(N_int,2)


    Reunion of the  inactive and active bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_int`
       * :c:data:`virt_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`

 
.. c:var:: reunion_of_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: reunion_of_bitmask	(N_int,2)


    Reunion of the inactive, active and virtual bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`inact_bitmask`
       * :c:data:`n_int`
       * :c:data:`virt_bitmask`


 
.. c:var:: reunion_of_core_inact_act_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: reunion_of_core_inact_act_bitmask	(N_int,2)


    Reunion of the core, inactive and active bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_inact_act_bitmask_4`
       * :c:data:`list_core_inact_act`

 
.. c:var:: reunion_of_core_inact_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: reunion_of_core_inact_bitmask	(N_int,2)


    Reunion of the core and inactive and virtual bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_bitmask`
       * :c:data:`inact_bitmask`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_core_inact`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_det_generators`
       * :c:data:`psi_det_generators`
       * :c:data:`reunion_of_core_inact_act_bitmask`

 
.. c:var:: reunion_of_inact_act_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: reunion_of_inact_act_bitmask	(N_int,2)


    Reunion of the  inactive and active bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`inact_bitmask`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`list_inact_act`

 
.. c:var:: unpaired_alpha_electrons


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: unpaired_alpha_electrons	(N_int)


    Bitmask reprenting the unpaired alpha electrons in the HF_bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_int`


 
.. c:var:: virt_bitmask


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)


    Bitmask identifying the virtual MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_virt`
       * :c:data:`n_int`
       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`inact_virt_bitmask`
       * :c:data:`n_det_generators`
       * :c:data:`psi_det_generators`
       * :c:data:`reunion_of_act_virt_bitmask`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: virt_bitmask_4


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: virt_bitmask_4	(N_int,4)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`virt_bitmask`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: bitstring_to_hexa:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine bitstring_to_hexa( output, string, Nint )


    Transform a bit string to a string in hexadecimal format for printing

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_cfg`
       * :c:func:`debug_det`
       * :c:func:`debug_spindet`

 
.. c:function:: bitstring_to_list:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine bitstring_to_list( string, list, n_elements, Nint)


    Gives the indices(+1) of the bits set to 1 in the bit string

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map`
       * :c:func:`add_integrals_to_map_erf`
       * :c:func:`create_microlist`
       * :c:func:`example_bitmask`
       * :c:func:`generate_cas_space`
       * :c:func:`getmobiles`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`ref_bitmask_energy`
       * :c:func:`splash_p`
       * :c:func:`spot_hasbeen`

 
.. c:function:: bitstring_to_str:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine bitstring_to_str( output, string, Nint )


    Transform a bit string to a string for printing

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map_erf`
       * :c:func:`example_bitmask`
       * :c:func:`print_det`
       * :c:func:`print_det_one_dimension`
       * :c:func:`print_spindet`

 
.. c:function:: broadcast_chunks_bit_kind:


    File : :file:`bitmask/mpi.irp.f`

    .. code:: fortran

        subroutine broadcast_chunks_bit_kind(A, LDA)


    Broadcast with chunks of ~2GB

 
.. c:function:: clear_bit_to_integer:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine clear_bit_to_integer(i_physical,key,Nint)


    set to 0 the bit number i_physical in the bitstring key

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`example_bitmask`
       * :c:data:`ref_closed_shell_bitmask`

 
.. c:function:: configuration_to_str:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine configuration_to_str( output, string, Nint )


    Transform the bit string of a configuration to a string for printing

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_cfg`

 
.. c:function:: debug_cfg:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine debug_cfg(string,Nint)


    Subroutine to print the content of a determinant in '+-' notation and
    hexadecimal representation.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_hexa`
       * :c:func:`configuration_to_str`

 
.. c:function:: debug_det:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine debug_det(string,Nint)


    Subroutine to print the content of a determinant in '+-' notation and
    hexadecimal representation.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`build_fock_tmp`
       * :c:func:`example_determinants`
       * :c:func:`get_excitation_degree_vector_single_or_exchange_verbose`
       * :c:func:`get_particles_general`
       * :c:func:`number_of_holes_verbose`
       * :c:func:`number_of_particles_verbose`
       * :c:func:`routine_example_psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_hexa`
       * :c:func:`print_det`

 
.. c:function:: debug_spindet:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine debug_spindet(string,Nint)


    Subroutine to print the content of a determinant in '+-' notation and
    hexadecimal representation.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_hexa`
       * :c:func:`print_spindet`

 
.. c:function:: det_allowed_ormas:


    File : :file:`bitmask/bitmasks_ormas.irp.f`

    .. code:: fortran

        logical function det_allowed_ormas(key_in)


    return true if det has allowable ormas occupations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`ormas_bitmask`
       * :c:data:`ormas_max_e`
       * :c:data:`ormas_min_e`
       * :c:data:`ormas_n_space`

 
.. c:function:: example_bitmask:


    File : :file:`bitmask/example.irp.f`

    .. code:: fortran

        subroutine example_bitmask


    subroutine that illustrates the main features available in bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_int`
       * :c:data:`n_virt_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list`
       * :c:func:`bitstring_to_str`
       * :c:func:`clear_bit_to_integer`
       * :c:func:`set_bit_to_integer`

 
.. c:function:: initialize_mo_coef_begin_iteration:


    File : :file:`bitmask/track_orb.irp.f`

    .. code:: fortran

        subroutine initialize_mo_coef_begin_iteration


    
    Initialize :c:data:`mo_coef_begin_iteration` to the current :c:data:`mo_coef`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_coef`
       * :c:data:`mo_coef_begin_iteration`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`damping_scf`
       * :c:func:`roothaan_hall_scf`

 
.. c:function:: is_a_1h:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_1h(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_1h1p:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_1h1p(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_1h2p:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_1h2p(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_1p:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_1p(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_2h:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_2h(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_2h1p:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_2h1p(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_2p:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_2p(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_two_holes_two_particles:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_two_holes_two_particles(key_in)


    logical function that returns True if the determinant 'key_in'
    belongs to the 2h-2p excitation class of the DDCI space
    this is calculated using the act_bitmask that defines the active
    orbital space, the inact_bitmasl that defines the inactive oribital space
    and the virt_bitmask that defines the virtual orbital space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask`

 
.. c:function:: is_i_in_virtual:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_i_in_virtual(i)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`virt_bitmask`

 
.. c:function:: is_integer_in_string:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        logical function is_integer_in_string(bite,string,Nint)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`set_bit_to_integer`

 
.. c:function:: is_the_hole_in_det:


    File : :file:`bitmask/find_hole.irp.f`

    .. code:: fortran

        logical function is_the_hole_in_det(key_in,ispin,i_hole)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_the_particl_in_det:


    File : :file:`bitmask/find_hole.irp.f`

    .. code:: fortran

        logical function is_the_particl_in_det(key_in,ispin,i_particl)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: list_to_bitstring:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine list_to_bitstring( string, list, n_elements, Nint)


    Returns the physical string "string(N_int,2)" from the array of
    occupations "list(N_int*bit_kind_size,2)

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`core_bitmask`
       * :c:data:`del_bitmask`
       * :c:func:`generate_cas_space`
       * :c:data:`hf_bitmask`
       * :c:data:`inact_bitmask`
       * :c:func:`orb_range_2_rdm_openmp_work_1`
       * :c:func:`orb_range_2_rdm_openmp_work_2`
       * :c:func:`orb_range_2_rdm_openmp_work_3`
       * :c:func:`orb_range_2_rdm_openmp_work_4`
       * :c:func:`orb_range_2_rdm_openmp_work_n_int`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_1`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_2`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_3`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_4`
       * :c:func:`orb_range_2_rdm_state_av_openmp_work_n_int`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_1`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_2`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_3`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_4`
       * :c:func:`orb_range_2_trans_rdm_openmp_work_n_int`
       * :c:data:`ormas_bitmask`
       * :c:data:`virt_bitmask`

 
.. c:function:: modify_bitmasks_for_hole:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine modify_bitmasks_for_hole(i_hole)


    modify the generators_bitmask in order that one can only excite
    the electrons occupying i_hole

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`index_holes_bitmask`
       * :c:data:`n_int`

 
.. c:function:: modify_bitmasks_for_hole_in_out:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine modify_bitmasks_for_hole_in_out(i_hole)


    modify the generators_bitmask in order that one can only excite
    the electrons occupying i_hole

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`index_holes_bitmask`

 
.. c:function:: modify_bitmasks_for_particl:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine modify_bitmasks_for_particl(i_part)


    modify the generators_bitmask in order that one can only excite
    the electrons to the orbital i_part

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`index_particl_bitmask`
       * :c:data:`n_int`

 
.. c:function:: number_of_holes:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        integer function number_of_holes(key_in)


    Function that returns the number of holes in the inact space
    
      popcnt(
         xor(
           iand(
             reunion_of_core_inact_bitmask(1,1),
             xor(
               key_in(1,1),
               iand(
                 key_in(1,1),
                 act_bitmask(1,1))
             )
           ),
           reunion_of_core_inact_bitmask(1,1)) )
    
    (key_in && act_bitmask)
    +---------------------+
       electrons in cas     xor key_in
    +---------------------------------+
           electrons outside of cas     && reunion_of_core_inact_bitmask
    +------------------------------------------------------------------+
               electrons in the core/inact space     xor reunion_of_core_inact_bitmask
    +---------------------------------------------------------------------------------+
                 holes

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_bitmask`

 
.. c:function:: number_of_holes_verbose:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        integer function number_of_holes_verbose(key_in)


    function that returns the number of holes in the inact space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_bitmask`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`

 
.. c:function:: number_of_particles:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        integer function number_of_particles(key_in)


    function that returns the number of particles in the virtual space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_int`
       * :c:data:`virt_bitmask`

 
.. c:function:: number_of_particles_verbose:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        integer function number_of_particles_verbose(key_in)


    function that returns the number of particles in the inact space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_int`
       * :c:data:`virt_bitmask`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`

 
.. c:function:: ormas_occ:


    File : :file:`bitmask/bitmasks_ormas.irp.f`

    .. code:: fortran

        subroutine ormas_occ(key_in, occupancies)


    number of electrons in each ormas space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`ormas_bitmask`
       * :c:data:`ormas_n_space`

 
.. c:function:: print_det:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine print_det(string,Nint)


    Subroutine to print the content of a determinant using the '+-' notation

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`
       * :c:func:`example_determinants`
       * :c:func:`print_generators_bitmasks_particles`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_str`

 
.. c:function:: print_det_one_dimension:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine print_det_one_dimension(string,Nint)


    Subroutine to print the content of a determinant using the '+-' notation

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_str`

 
.. c:function:: print_generators_bitmasks_holes:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine print_generators_bitmasks_holes



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`index_holes_bitmask`
       * :c:data:`n_int`

 
.. c:function:: print_generators_bitmasks_particles:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine print_generators_bitmasks_particles



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`index_particl_bitmask`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_det`

 
.. c:function:: print_spindet:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine print_spindet(string,Nint)


    Subroutine to print the content of a determinant using the '+-' notation

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_spindet`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_str`

 
.. c:function:: reorder_core_orb:


    File : :file:`bitmask/track_orb.irp.f`

    .. code:: fortran

        subroutine reorder_core_orb


    routines that takes the current :c:data:`mo_coef` and reorder the core orbitals (see :c:data:`list_core` and :c:data:`n_core_orb`) according to the overlap with :c:data:`mo_coef_begin_iteration`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`list_core`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_begin_iteration`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`damping_scf`
       * :c:func:`roothaan_hall_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsort`

 
.. c:function:: set_bit_to_integer:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine set_bit_to_integer(i_physical,key,Nint)


    set to 1 the bit number i_physical in the bitstring key

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`example_bitmask`
       * :c:func:`is_integer_in_string`
       * :c:data:`orb_swap`

 
.. c:function:: set_bitmask_hole_as_input:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine set_bitmask_hole_as_input(input_bitmask)


    set the generators_bitmask for the holes
    as the input_bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`index_holes_bitmask`
       * :c:data:`n_int`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`

 
.. c:function:: set_bitmask_particl_as_input:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine set_bitmask_particl_as_input(input_bitmask)


    set the generators_bitmask for the particles
    as the input_bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`index_particl_bitmask`
       * :c:data:`n_int`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`

