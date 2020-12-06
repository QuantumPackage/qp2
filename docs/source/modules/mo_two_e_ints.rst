.. _module_mo_two_e_ints: 
 
.. program:: mo_two_e_ints 
 
.. default-role:: option 
 
==================
mo_two_e_ints
==================

Here, all two-electron integrals (:math:`1/r_{12}`) are computed.
As they have 4 indices and many are zero, they are stored in a map, as defined
in :file:`Utils/map_module.f90`.

To fetch an |AO| integral, use the
`get_ao_two_e_integral(i,j,k,l,ao_integrals_map)` function, and
to fetch an |MO| integral, use
`get_two_e_integral(i,j,k,l,mo_integrals_map)` or
`mo_two_e_integral(i,j,k,l)`.

The conventions are:

* For |AO| integrals : (ik|jl) = (11|22)
* For |MO| integrals : <ij|kl> = <12|12>



 
 
 
EZFIO parameters 
---------------- 
 
.. option:: io_mo_two_e_integrals
 
    Read/Write |MO| integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: mo_integrals_threshold
 
    If | <ij|kl> | < `mo_integrals_threshold` then <ij|kl> is zero
 
    Default: 1.e-15
 
 
Providers 
--------- 
 
.. c:var:: banned_excitation


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        logical, allocatable	:: banned_excitation	(mo_num,mo_num)


    If true, the excitation is banned in the selection. Useful with local MOs.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: big_array_coulomb_integrals


    File : :file:`mo_two_e_ints/integrals_3_index.irp.f`

    .. code:: fortran

        double precision, allocatable	:: big_array_coulomb_integrals	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: big_array_exchange_integrals	(mo_num,mo_num,mo_num)


    big_array_coulomb_integrals(i,j)  = <ij|ij> = (ii|jj)
    
    big_array_exchange_integrals(i,j) = <ij|ji> = (ij|ij)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`

 
.. c:var:: big_array_exchange_integrals


    File : :file:`mo_two_e_ints/integrals_3_index.irp.f`

    .. code:: fortran

        double precision, allocatable	:: big_array_coulomb_integrals	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: big_array_exchange_integrals	(mo_num,mo_num,mo_num)


    big_array_coulomb_integrals(i,j)  = <ij|ij> = (ii|jj)
    
    big_array_exchange_integrals(i,j) = <ij|ji> = (ij|ij)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`

 
.. c:var:: core_energy


    File : :file:`mo_two_e_ints/core_quantities.irp.f`

    .. code:: fortran

        double precision	:: core_energy	


    energy from the core : contains all core-core contributions

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`n_core_orb`
       * :c:data:`nuclear_repulsion`


 
.. c:var:: core_fock_operator


    File : :file:`mo_two_e_ints/core_quantities.irp.f`

    .. code:: fortran

        double precision, allocatable	:: core_fock_operator	(mo_num,mo_num)


    this is the contribution to the Fock operator from the core electrons

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`


 
.. c:function:: insert_into_mo_integrals_map:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine insert_into_mo_integrals_map(n_integrals,                 &
      buffer_i, buffer_values, thr)


    Create new entry into MO map, or accumulate in an existing entry

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map`
       * :c:func:`add_integrals_to_map_no_exit_34`
       * :c:func:`add_integrals_to_map_three_indices`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_update`

 
.. c:var:: mo_coef_novirt


    File : :file:`mo_two_e_ints/four_idx_novvvv.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_coef_novirt	(ao_num,n_core_inact_act_orb)


    MO coefficients without virtual MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_coef`
       * :c:data:`n_core_inact_act_orb`


 
.. c:var:: mo_integrals_cache


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_integrals_cache	(0_8:128_8*128_8*128_8*128_8)


    Cache of MO integrals for fast access

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_max


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer*4	:: mo_integrals_cache_min	
        integer*4	:: mo_integrals_cache_max	
        integer*8	:: mo_integrals_cache_min_8	
        integer*8	:: mo_integrals_cache_max_8	


    Min and max values of the MOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_max_8


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer*4	:: mo_integrals_cache_min	
        integer*4	:: mo_integrals_cache_max	
        integer*8	:: mo_integrals_cache_min_8	
        integer*8	:: mo_integrals_cache_max_8	


    Min and max values of the MOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_min


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer*4	:: mo_integrals_cache_min	
        integer*4	:: mo_integrals_cache_max	
        integer*8	:: mo_integrals_cache_min_8	
        integer*8	:: mo_integrals_cache_max_8	


    Min and max values of the MOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_min_8


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer*4	:: mo_integrals_cache_min	
        integer*4	:: mo_integrals_cache_max	
        integer*8	:: mo_integrals_cache_min_8	
        integer*8	:: mo_integrals_cache_max_8	


    Min and max values of the MOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_map


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        type(map_type)	:: mo_integrals_map	


    MO integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_two_e_integral_jj_from_ao


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_integral_jj_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_exchange_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_anti_from_ao	(mo_num,mo_num)


    mo_two_e_integral_jj_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_exchange_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_anti_from_ao(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`


 
.. c:var:: mo_two_e_integrals_in_map


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        logical	:: mo_two_e_integrals_in_map	


    If True, the map of MO two-electron integrals is provided

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`ezfio_filename`
       * :c:data:`full_ijkl_bitmask_4`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_int`
       * :c:data:`no_vvvv_integrals`
       * :c:data:`read_mo_two_e_integrals`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`ci_electronic_energy`
       * :c:data:`core_fock_operator`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_two_e_integrals_jj


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_integrals_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_anti	(mo_num,mo_num)


    mo_two_e_integrals_jj(i,j) = J_ij
    mo_two_e_integrals_jj_exchange(i,j) = K_ij
    mo_two_e_integrals_jj_anti(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`ref_bitmask_energy`

 
.. c:var:: mo_two_e_integrals_jj_anti


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_integrals_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_anti	(mo_num,mo_num)


    mo_two_e_integrals_jj(i,j) = J_ij
    mo_two_e_integrals_jj_exchange(i,j) = K_ij
    mo_two_e_integrals_jj_anti(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`ref_bitmask_energy`

 
.. c:var:: mo_two_e_integrals_jj_anti_from_ao


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_integral_jj_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_exchange_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_anti_from_ao	(mo_num,mo_num)


    mo_two_e_integral_jj_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_exchange_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_anti_from_ao(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`


 
.. c:var:: mo_two_e_integrals_jj_exchange


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_integrals_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_anti	(mo_num,mo_num)


    mo_two_e_integrals_jj(i,j) = J_ij
    mo_two_e_integrals_jj_exchange(i,j) = K_ij
    mo_two_e_integrals_jj_anti(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`ref_bitmask_energy`

 
.. c:var:: mo_two_e_integrals_jj_exchange_from_ao


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_integral_jj_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_exchange_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_jj_anti_from_ao	(mo_num,mo_num)


    mo_two_e_integral_jj_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_exchange_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_anti_from_ao(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`


 
.. c:var:: mo_two_e_integrals_vv_anti_from_ao


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_integrals_vv_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_vv_exchange_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_vv_anti_from_ao	(mo_num,mo_num)


    mo_two_e_integrals_vv_from_ao(i,j) = J_ij
    mo_two_e_integrals_vv_exchange_from_ao(i,j) = J_ij
    mo_two_e_integrals_vv_anti_from_ao(i,j) = J_ij - K_ij
    but only for the virtual orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`list_virt`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_virt_orb`


 
.. c:var:: mo_two_e_integrals_vv_exchange_from_ao


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_integrals_vv_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_vv_exchange_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_vv_anti_from_ao	(mo_num,mo_num)


    mo_two_e_integrals_vv_from_ao(i,j) = J_ij
    mo_two_e_integrals_vv_exchange_from_ao(i,j) = J_ij
    mo_two_e_integrals_vv_anti_from_ao(i,j) = J_ij - K_ij
    but only for the virtual orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`list_virt`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_virt_orb`


 
.. c:var:: mo_two_e_integrals_vv_from_ao


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_integrals_vv_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_vv_exchange_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_integrals_vv_anti_from_ao	(mo_num,mo_num)


    mo_two_e_integrals_vv_from_ao(i,j) = J_ij
    mo_two_e_integrals_vv_exchange_from_ao(i,j) = J_ij
    mo_two_e_integrals_vv_anti_from_ao(i,j) = J_ij - K_ij
    but only for the virtual orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`list_virt`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`
       * :c:data:`n_virt_orb`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: add_integrals_to_map:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine add_integrals_to_map(mask_ijkl)


    Adds integrals to tha MO map according to some bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`four_idx_novvvv2`
       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list`
       * :c:func:`cpu_time`
       * :c:func:`get_ao_two_e_integrals`
       * :c:func:`insert_into_mo_integrals_map`
       * :c:func:`map_merge`
       * :c:func:`mo_two_e_integrals_index`
       * :c:func:`wall_time`

 
.. c:function:: add_integrals_to_map_no_exit_34:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine add_integrals_to_map_no_exit_34(mask_ijkl)


    Adds integrals to tha MO map according to some bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list`
       * :c:func:`cpu_time`
       * :c:func:`get_ao_two_e_integrals`
       * :c:func:`insert_into_mo_integrals_map`
       * :c:func:`map_merge`
       * :c:func:`mo_two_e_integrals_index`
       * :c:func:`wall_time`

 
.. c:function:: add_integrals_to_map_three_indices:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine add_integrals_to_map_three_indices(mask_ijk)


    Adds integrals to tha MO map according to some bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list`
       * :c:func:`cpu_time`
       * :c:func:`get_ao_two_e_integrals`
       * :c:func:`insert_into_mo_integrals_map`
       * :c:func:`map_merge`
       * :c:func:`mo_two_e_integrals_index`
       * :c:func:`wall_time`

 
.. c:function:: ao_to_mo_novirt:


    File : :file:`mo_two_e_ints/four_idx_novvvv.irp.f`

    .. code:: fortran

        subroutine ao_to_mo_novirt(A_ao,LDA_ao,A_mo,LDA_mo)


    Transform A from the |AO| basis to the |MO| basis excluding virtuals
    
    $C^\dagger.A_{ao}.C$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef_novirt`
       * :c:data:`n_core_inact_act_orb`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`four_idx_novvvv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: clear_mo_map:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine clear_mo_map


    Frees the memory of the MO map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_deinit`

 
.. c:function:: dump_mo_integrals:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine dump_mo_integrals(filename)


    Save to disk the |MO| integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`
       * :c:data:`mpi_master`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_work_empty`

 
.. c:function:: four_idx_novvvv:


    File : :file:`mo_two_e_ints/four_idx_novvvv.irp.f`

    .. code:: fortran

        subroutine four_idx_novvvv


    Retransform MO integrals for next CAS-SCF step

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`n_core_inact_act_orb`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_to_mo`
       * :c:func:`ao_to_mo_novirt`
       * :c:func:`map_append`
       * :c:func:`map_shrink`
       * :c:func:`map_sort`
       * :c:func:`map_unique`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: four_idx_novvvv2:


    File : :file:`mo_two_e_ints/four_idx_novvvv.irp.f`

    .. code:: fortran

        subroutine four_idx_novvvv2



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`core_inact_act_bitmask_4`
       * :c:data:`full_ijkl_bitmask_4`
       * :c:data:`n_int`
       * :c:data:`virt_bitmask`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map`

 
.. c:function:: get_mo_map_size:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer*8 function get_mo_map_size()


    Return the number of elements in the MO map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`

 
.. c:function:: get_mo_two_e_integrals:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals(j,k,l,sze,out_val,map)


    Returns multiple integrals <ij|kl> in the MO basis, all
    i for j,k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_mo_two_e_integrals_i1j1`
       * :c:func:`get_mo_two_e_integrals_ij`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`

 
.. c:function:: get_mo_two_e_integrals_coulomb_ii:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_coulomb_ii(k,l,sze,out_val,map)


    Returns multiple integrals <ki|li>
    k(1)i(2) 1/r12 l(1)i(2) :: out_val(i1)
    for k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`

 
.. c:function:: get_mo_two_e_integrals_exch_ii:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_exch_ii(k,l,sze,out_val,map)


    Returns multiple integrals <ki|il>
    k(1)i(2) 1/r12 i(1)l(2) :: out_val(i1)
    for k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`

 
.. c:function:: get_mo_two_e_integrals_i1j1:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_i1j1(k,l,sze,out_array,map)


    Returns multiple integrals <ik|jl> in the MO basis, all
    i(1)j(1) 1/r12 k(2)l(2)
    i, j for k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_mo_two_e_integrals`

 
.. c:function:: get_mo_two_e_integrals_ij:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_ij(k,l,sze,out_array,map)


    Returns multiple integrals <ij|kl> in the MO basis, all
    i(1)j(2) 1/r12 k(1)l(2)
    i, j for k,l fixed.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_mo_two_e_integrals`

 
.. c:function:: get_two_e_integral:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision function get_two_e_integral(i,j,k,l,map)


    Returns one integral <ij|kl> in the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: load_mo_integrals:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer function load_mo_integrals(filename)


    Read from disk the |MO| integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`cache_map_reallocate`
       * :c:func:`map_deinit`
       * :c:func:`map_sort`

 
.. c:function:: mo_two_e_integral:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision function mo_two_e_integral(i,j,k,l)


    Returns one integral <ij|kl> in the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`

 
.. c:function:: mo_two_e_integrals_index:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine mo_two_e_integrals_index(i,j,k,l,i1)


    Computes an unique index for i,j,k,l integrals

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map`
       * :c:func:`add_integrals_to_map_no_exit_34`
       * :c:func:`add_integrals_to_map_three_indices`

