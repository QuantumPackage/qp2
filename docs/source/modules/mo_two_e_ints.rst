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
 
.. option:: io_mo_cholesky
 
    Read/Write |MO| Cholesky integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: io_mo_two_e_integrals
 
    Read/Write |MO| integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: mo_integrals_cache_shift
 
    Adjusts the size of the MO integrals cache. 2: 2KB, 3: 32KB, 4: 512KB, 5: 8MB, 6: 128MB, 7: 2GB, 8: 32GB, 9: 512GB
 
    Default: 7
 
.. option:: mo_integrals_threshold
 
    If | <ij|kl> | < `mo_integrals_threshold` then <ij|kl> is zero
 
    Default: 1.e-15
 
.. option:: io_mo_two_e_integrals_erf
 
    Read/Write MO integrals with the long range interaction from/to disk [    Write | Read | None ]
 
    Default: None
 
 
Providers 
--------- 
 
.. c:var:: banned_excitation


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        logical, allocatable	:: banned_excitation	(mo_num,mo_num)
        logical	:: use_banned_excitation	


    If true, the excitation is banned in the selection. Useful with local MOs.

    Needs:

    .. hlist::
       :columns: 3

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


    big_array_coulomb_integrals(j,i,k)  = <ij|kj> = (ik|jj)
    
    big_array_exchange_integrals(j,i,k) = <ij|jk> = (ij|kj)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_hf_selector`
       * :c:data:`h_core_ri`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`h_matrix_diag_all_dets`
       * :c:data:`psi_energy_two_e_trans`

 
.. c:var:: big_array_exchange_integrals


    File : :file:`mo_two_e_ints/integrals_3_index.irp.f`

    .. code:: fortran

        double precision, allocatable	:: big_array_coulomb_integrals	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: big_array_exchange_integrals	(mo_num,mo_num,mo_num)


    big_array_coulomb_integrals(j,i,k)  = <ij|kj> = (ik|jj)
    
    big_array_exchange_integrals(j,i,k) = <ij|jk> = (ij|kj)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_hf_selector`
       * :c:data:`h_core_ri`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`h_matrix_diag_all_dets`
       * :c:data:`psi_energy_two_e_trans`

 
.. c:var:: cholesky_mo


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cholesky_mo	(mo_num,mo_num,cholesky_mo_num)


    Cholesky vectors in MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`mo_num`


 
.. c:var:: cholesky_mo_num


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        integer	:: cholesky_mo_num	


    Number of Cholesky vectors in MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_ao_num`
       * :c:data:`ezfio_work_dir`
       * :c:data:`read_mo_cholesky`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:data:`core_fock_operator`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: cholesky_mo_transp


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cholesky_mo_transp	(cholesky_mo_num,mo_num,mo_num)


    Cholesky vectors in MO basis. Warning: it is transposed wrt cholesky_ao:
    
    -  cholesky_ao        is (ao_num^2 x cholesky_ao_num)
    
    - cholesky_mo_transp is (cholesky_mo_num x mo_num^2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_mo_num`
       * :c:data:`ezfio_work_dir`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`read_mo_cholesky`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:data:`core_fock_operator`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: cholesky_semi_mo_transp_simple


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cholesky_semi_mo_transp_simple	(cholesky_mo_num,ao_num,mo_num)


    Cholesky vectors in MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`


 
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


 
.. c:var:: core_energy_erf


    File : :file:`mo_two_e_ints/core_quantities_erf.irp.f`

    .. code:: fortran

        double precision	:: core_energy_erf	


    energy from the core : contains all core-core contributionswith the erf interaction

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core`
       * :c:data:`mo_two_e_int_erf_jj`
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
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`


 
.. c:var:: core_fock_operator_erf


    File : :file:`mo_two_e_ints/core_quantities_erf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: core_fock_operator_erf	(mo_num,mo_num)


    this is the contribution to the Fock operator from the core electrons with the erf interaction

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_integrals_erf_cache_min`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`


 
.. c:var:: do_mo_cholesky


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        logical	:: do_mo_cholesky	


    If True, use Cholesky vectors for MO integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`do_ao_cholesky`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: h_core_ri


    File : :file:`mo_two_e_ints/core_quantities.irp.f`

    .. code:: fortran

        double precision, allocatable	:: h_core_ri	(mo_num,mo_num)


    Core Hamiltonian with 3-index exchange integrals:
    
    :math:`\tilde{h}{pq} = h_{pq} - \frac{1}{2}\sum_{k} g(pk,kq)` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`


 
.. c:function:: insert_into_mo_integrals_erf_map:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine insert_into_mo_integrals_erf_map(n_integrals,                 &
      buffer_i, buffer_values, thr)


    Create new entry into |MO| map, or accumulate in an existing entry

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_update`

 
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

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_update`

 
.. c:var:: int_erf_3_index


    File : :file:`mo_two_e_ints/ints_erf_3_index.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int_erf_3_index	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: int_erf_3_index_exc	(mo_num,mo_num,mo_num)


    int_erf_3_index(i,j)     = <ij|ij> = (ii|jj) with the erf interaction
    
    int_erf_3_index_exc(i,j) = <ij|ji> = (ij|ij) with the erf interaction

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_integrals_erf_cache_min`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`


 
.. c:var:: int_erf_3_index_exc


    File : :file:`mo_two_e_ints/ints_erf_3_index.irp.f`

    .. code:: fortran

        double precision, allocatable	:: int_erf_3_index	(mo_num,mo_num,mo_num)
        double precision, allocatable	:: int_erf_3_index_exc	(mo_num,mo_num,mo_num)


    int_erf_3_index(i,j)     = <ij|ij> = (ii|jj) with the erf interaction
    
    int_erf_3_index_exc(i,j) = <ij|ji> = (ij|ij) with the erf interaction

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_integrals_erf_cache_min`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`


 
.. c:var:: mo_integrals_cache


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_integrals_cache	(0_8:(1_8*mo_integrals_cache_size)**4)


    Cache of MO integrals for fast access

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_max


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer	:: mo_integrals_cache_min	
        integer	:: mo_integrals_cache_max	
        integer	:: mo_integrals_cache_size	


    Min and max values of the MOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_min


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer	:: mo_integrals_cache_min	
        integer	:: mo_integrals_cache_max	
        integer	:: mo_integrals_cache_size	


    Min and max values of the MOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_size


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer	:: mo_integrals_cache_min	
        integer	:: mo_integrals_cache_max	
        integer	:: mo_integrals_cache_size	


    Min and max values of the MOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_erf_cache


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_integrals_erf_cache	(0:64*64*64*64)


    Cache of |MO| integrals for fast access

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_cache_min`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_fock_operator_erf`
       * :c:data:`int_erf_3_index`
       * :c:data:`mo_two_e_int_erf_jj`

 
.. c:var:: mo_integrals_erf_cache_max


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        integer	:: mo_integrals_erf_cache_min	
        integer	:: mo_integrals_erf_cache_max	


    Min and max values of the MOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_fock_operator_erf`
       * :c:data:`int_erf_3_index`
       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_two_e_int_erf_jj`

 
.. c:var:: mo_integrals_erf_cache_min


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        integer	:: mo_integrals_erf_cache_min	
        integer	:: mo_integrals_erf_cache_max	


    Min and max values of the MOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_fock_operator_erf`
       * :c:data:`int_erf_3_index`
       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_two_e_int_erf_jj`

 
.. c:var:: mo_integrals_erf_map


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        type(map_type)	:: mo_integrals_erf_map	


    |MO| integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_fock_operator_erf`
       * :c:data:`int_erf_3_index`
       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_two_e_int_erf_jj`
       * :c:data:`mo_two_e_integrals_erf_in_map`

 
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

       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`bielecci`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`coef_hf_selector`
       * :c:data:`core_fock_operator`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`h_matrix_diag_all_dets`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`psi_energy_two_e_trans`

 
.. c:var:: mo_two_e_int_erf_jj


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_int_erf_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_anti	(mo_num,mo_num)


    mo_two_e_integrals_jj(i,j) = J_ij
    mo_two_e_integrals_jj_exchange(i,j) = K_ij
    mo_two_e_integrals_jj_anti(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_integrals_erf_cache_min`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy_erf`

 
.. c:var:: mo_two_e_int_erf_jj_anti


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_int_erf_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_anti	(mo_num,mo_num)


    mo_two_e_integrals_jj(i,j) = J_ij
    mo_two_e_integrals_jj_exchange(i,j) = K_ij
    mo_two_e_integrals_jj_anti(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_integrals_erf_cache_min`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy_erf`

 
.. c:var:: mo_two_e_int_erf_jj_anti_from_ao


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_int_erf_jj_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_exchange_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_anti_from_ao	(mo_num,mo_num)


    mo_two_e_integral_jj_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_exchange_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_anti_from_ao(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integral_erf_schwartz`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`


 
.. c:var:: mo_two_e_int_erf_jj_exchange


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_int_erf_jj	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_exchange	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_anti	(mo_num,mo_num)


    mo_two_e_integrals_jj(i,j) = J_ij
    mo_two_e_integrals_jj_exchange(i,j) = K_ij
    mo_two_e_integrals_jj_anti(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_integrals_erf_cache_min`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy_erf`

 
.. c:var:: mo_two_e_int_erf_jj_exchange_from_ao


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_int_erf_jj_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_exchange_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_anti_from_ao	(mo_num,mo_num)


    mo_two_e_integral_jj_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_exchange_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_anti_from_ao(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integral_erf_schwartz`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`


 
.. c:var:: mo_two_e_int_erf_jj_from_ao


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_two_e_int_erf_jj_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_exchange_from_ao	(mo_num,mo_num)
        double precision, allocatable	:: mo_two_e_int_erf_jj_anti_from_ao	(mo_num,mo_num)


    mo_two_e_integral_jj_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_exchange_from_ao(i,j) = J_ij
    mo_two_e_integrals_jj_anti_from_ao(i,j) = J_ij - K_ij

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integral_erf_schwartz`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_num`


 
.. c:var:: mo_two_e_integrals_erf_in_map


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        logical	:: mo_two_e_integrals_erf_in_map	


    If True, the map of MO two-electron integrals is provided

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`ezfio_filename`
       * :c:data:`full_ijkl_bitmask_4`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_int`
       * :c:data:`qp_max_mem`
       * :c:data:`read_mo_two_e_integrals_erf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_fock_operator_erf`
       * :c:data:`int_erf_3_index`
       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_two_e_int_erf_jj`

 
.. c:var:: mo_two_e_integrals_in_map


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        logical	:: mo_two_e_integrals_in_map	


    If True, the map of MO two-electron integrals is provided

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_ao_cholesky`
       * :c:data:`do_mo_cholesky`
       * :c:data:`ezfio_filename`
       * :c:data:`full_ijkl_bitmask_4`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_int`
       * :c:data:`qp_max_mem`
       * :c:data:`read_mo_two_e_integrals`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`act_2_rdm_aa_mo`
       * :c:data:`act_2_rdm_ab_mo`
       * :c:data:`act_2_rdm_bb_mo`
       * :c:data:`act_2_rdm_spin_trace_mo`
       * :c:data:`act_2_rdm_trans_spin_trace_mo`
       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`
       * :c:data:`bielecci`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`ci_electronic_energy`
       * :c:data:`coef_hf_selector`
       * :c:data:`core_fock_operator`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`h_matrix_diag_all_dets`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`psi_energy_two_e_trans`

 
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
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
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
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
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
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
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

 
.. c:var:: use_banned_excitation


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        logical, allocatable	:: banned_excitation	(mo_num,mo_num)
        logical	:: use_banned_excitation	


    If true, the excitation is banned in the selection. Useful with local MOs.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_two_e_integrals_jj`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: add_integrals_to_map:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine add_integrals_to_map(mask_ijkl)


    Adds integrals to the MO map according to some bitmask

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

       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list`
       * :c:func:`get_ao_two_e_integrals`
       * :c:func:`insert_into_mo_integrals_map`
       * :c:func:`map_merge`
       * :c:func:`mo_two_e_integrals_index`
       * :c:func:`wall_time`

 
.. c:function:: add_integrals_to_map_cholesky:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine add_integrals_to_map_cholesky


    Adds integrals to the MO map using Cholesky vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`map_append`
       * :c:func:`map_sort`
       * :c:func:`map_unique`
       * :c:func:`mo_two_e_integrals_index`
       * :c:func:`set_multiple_levels_omp`

 
.. c:function:: add_integrals_to_map_erf:


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        subroutine add_integrals_to_map_erf(mask_ijkl)


    Adds integrals to tha MO map according to some bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list`
       * :c:func:`bitstring_to_str`
       * :c:func:`cpu_time`
       * :c:func:`get_ao_two_e_integrals_erf`
       * :c:func:`insert_into_mo_integrals_erf_map`
       * :c:func:`map_merge`
       * :c:func:`mo_two_e_integrals_index`
       * :c:func:`wall_time`

 
.. c:function:: clear_mo_erf_map:


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        subroutine clear_mo_erf_map


    Frees the memory of the MO map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_deinit`

 
.. c:function:: clear_mo_map:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine clear_mo_map


    Frees the memory of the MO map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`
       * :c:func:`update_parameters`

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

 
.. c:function:: four_idx_dgemm:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine four_idx_dgemm



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`mo_coef`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_ao_two_e_integrals`
       * :c:func:`map_append`
       * :c:func:`map_sort`
       * :c:func:`map_unique`
       * :c:func:`mo_two_e_integrals_index`

 
.. c:function:: four_idx_dgemm_erf:


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        subroutine four_idx_dgemm_erf



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_ao_two_e_integrals_erf`
       * :c:func:`map_append`
       * :c:func:`map_sort`
       * :c:func:`map_unique`
       * :c:func:`mo_two_e_integrals_index`

 
.. c:function:: get_mo_erf_map_size:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        integer*8 function get_mo_erf_map_size()


    Returns the number of elements in the |MO| map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`

 
.. c:function:: get_mo_map_size:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer*8 function get_mo_map_size()


    Return the number of elements in the MO map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`

 
.. c:function:: get_mo_two_e_integral_erf:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        double precision function get_mo_two_e_integral_erf(i,j,k,l,map)


    Returns one integral $\langle ij|kl \rangle$ in the |MO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_integrals_erf_cache_min`
       * :c:data:`mo_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index`

 
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
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_d0`
       * :c:func:`get_d1`
       * :c:func:`get_mo_two_e_integrals_i1j1`
       * :c:func:`get_mo_two_e_integrals_ij`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`get_mo_two_e_integrals_cache`
       * :c:func:`map_get`

 
.. c:function:: get_mo_two_e_integrals_cache:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_cache(j,k,l,sze,out_val)


    Returns multiple integrals <ij|kl> in the MO basis, all
    i for j,k,l fixed, all integrals from the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_mo_two_e_integrals`
       * :c:func:`get_mo_two_e_integrals_i1j1`
       * :c:func:`get_mo_two_e_integrals_ij`

 
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

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`

 
.. c:function:: get_mo_two_e_integrals_erf:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_erf(j,k,l,sze,out_val,map)


    Returns multiple integrals $\langle ij|kl \rangle$ in the |MO| basis, all
    i for j,k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get_many`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_mo_two_e_integrals_erf_coulomb_ii:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_erf_coulomb_ii(k,l,sze,out_val,map)


    Returns multiple integrals $\langle ki|li \rangle$
    
    k(1)i(2) 1/r12 l(1)i(2) :: out_val(i1)
    for k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get_many`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_mo_two_e_integrals_erf_exch_ii:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_erf_exch_ii(k,l,sze,out_val,map)


    Returns multiple integrals $\langle ki|il \rangle$
    
    $\int k(1)i(2) \frac{1}{r_{12}} i(1)l(2)$ :: out_val(i1)
    for k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get_many`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_mo_two_e_integrals_erf_i1j1:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_erf_i1j1(k,l,sze,out_array,map)


    Returns multiple integrals $\langle ik|jl \rangle$ in the |MO| basis, all
    $\int i(1)j(1) \frac{\erf(\mu * r_{12})}{r_{12}} k(2)l(2)$
    i, j for k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i2sort`
       * :c:func:`i8sort`
       * :c:func:`isort`
       * :c:func:`map_get_many`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_mo_two_e_integrals_erf_ij:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_erf_ij(k,l,sze,out_array,map)


    Returns multiple integrals $\langle ij|kl \rangle$ in the |MO| basis, all
    $\int i(1)j(2) \frac{1}{r_{12}} k(1)l(2)$
    i, j for k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i2sort`
       * :c:func:`i8sort`
       * :c:func:`isort`
       * :c:func:`map_get_many`
       * :c:func:`two_e_integrals_index`

 
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

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_num`
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

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`get_mo_two_e_integrals`
       * :c:func:`get_mo_two_e_integrals_cache`

 
.. c:function:: get_mo_two_e_integrals_ij:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_mo_two_e_integrals_ij(k,l,sze,out_array,map)


    Returns multiple integrals <ij|kl> in the MO basis, all
    i(1)j(2) 1/r12 k(1)l(2)
    i, j for k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pxxq_array`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_mo_two_e_integrals`
       * :c:func:`get_mo_two_e_integrals_cache`

 
.. c:function:: get_two_e_integral:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision function get_two_e_integral(i,j,k,l,map)


    Returns one integral <ij|kl> in the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_two_e_integral_cache:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision function get_two_e_integral_cache(i,j,k,l)


    Returns one integral <ij|kl> in the MO basis taken from the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`

 
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
       * :c:func:`lock_io`
       * :c:func:`map_deinit`
       * :c:func:`map_sort`
       * :c:func:`unlock_io`

 
.. c:function:: load_mo_integrals_erf:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        integer function load_mo_integrals_erf(filename)


    Read from disk the |MO| erf integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`

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

 
.. c:function:: mo_two_e_integral_erf:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        double precision function mo_two_e_integral_erf(i,j,k,l)


    Returns one integral $\langle ij|kl \rangle$ in the |MO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`

 
.. c:function:: mo_two_e_integrals_erf_index:


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        subroutine mo_two_e_integrals_erf_index(i,j,k,l,i1)


    Computes an unique index for i,j,k,l integrals

 
.. c:function:: mo_two_e_integrals_index:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine mo_two_e_integrals_index(i,j,k,l,i1)


    Computes an unique index for i,j,k,l integrals

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map`
       * :c:func:`add_integrals_to_map_cholesky`
       * :c:func:`add_integrals_to_map_erf`
       * :c:func:`four_idx_dgemm`
       * :c:func:`four_idx_dgemm_erf`

 
.. c:function:: provide_all_mo_integrals_erf:


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        subroutine provide_all_mo_integrals_erf



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_two_e_int_erf_jj`
       * :c:data:`mo_two_e_int_erf_jj`
       * :c:data:`mo_two_e_int_erf_jj`
       * :c:data:`mo_two_e_integrals_erf_in_map`

 
.. c:function:: save_erf_two_e_integrals_mo:


    File : :file:`mo_two_e_ints/routines_save_integrals_erf.irp.f`

    .. code:: fortran

        subroutine save_erf_two_e_integrals_mo



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_mo_two_e_ints_io_mo_two_e_integrals_erf`
       * :c:func:`ezfio_set_work_empty`
       * :c:func:`map_save_to_disk`

 
.. c:function:: save_erf_two_e_ints_mo_into_ints_mo:


    File : :file:`mo_two_e_ints/routines_save_integrals_erf.irp.f`

    .. code:: fortran

        subroutine save_erf_two_e_ints_mo_into_ints_mo



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_mo_two_e_ints_io_mo_two_e_integrals`
       * :c:func:`ezfio_set_work_empty`
       * :c:func:`map_save_to_disk`

