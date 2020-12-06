.. _module_mo_two_e_erf_ints: 
 
.. program:: mo_two_e_erf_ints 
 
.. default-role:: option 
 
======================
mo_two_e_erf_ints
======================

Here, all two-electron integrals (:math:`erf({\mu}_{erf} * r_{12})/r_{12}`) are computed.
As they have 4 indices and many are zero, they are stored in a map, as defined
in :file:`Utils/map_module.f90`.

The range separation parameter :math:`{\mu}_{erf}` is the variable :option:`ao_two_e_erf_ints mu_erf`.

To fetch an |MO| integral, use
`get_mo_two_e_integral_erf(i,j,k,l,mo_integrals_map_erf)`

The conventions are:

* For |MO| integrals : <ij|kl> = <12|12>

Be aware that it might not be the same conventions for |MO| and |AO| integrals.


 
 
 
EZFIO parameters 
---------------- 
 
.. option:: io_mo_two_e_integrals_erf
 
    Read/Write MO integrals with the long range interaction from/to disk [    Write | Read | None ]
 
    Default: None
 
 
Providers 
--------- 
 
.. c:var:: core_energy_erf


    File : :file:`mo_two_e_erf_ints/core_quantities_erf.irp.f`

    .. code:: fortran

        double precision	:: core_energy_erf	


    energy from the core : contains all core-core contributionswith the erf interaction

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_int_erf_jj`
       * :c:data:`n_core_orb`
       * :c:data:`nuclear_repulsion`


 
.. c:var:: core_fock_operator_erf


    File : :file:`mo_two_e_erf_ints/core_quantities_erf.irp.f`

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


 
.. c:function:: insert_into_mo_integrals_erf_map:


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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

 
.. c:var:: int_erf_3_index


    File : :file:`mo_two_e_erf_ints/ints_erf_3_index.irp.f`

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


    File : :file:`mo_two_e_erf_ints/ints_erf_3_index.irp.f`

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


 
.. c:var:: mo_integrals_erf_cache


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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

 
.. c:var:: mo_two_e_int_erf_jj


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

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
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_transp`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`read_mo_two_e_integrals_erf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_fock_operator_erf`
       * :c:data:`int_erf_3_index`
       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_two_e_int_erf_jj`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: add_integrals_to_map_erf:


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

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

 
.. c:function:: get_mo_erf_map_size:


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        integer*8 function get_mo_erf_map_size()


    Returns the number of elements in the |MO| map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`

 
.. c:function:: get_mo_two_e_integral_erf:


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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

 
.. c:function:: get_mo_two_e_integrals_erf:


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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

       * :c:func:`i2radix_sort`
       * :c:func:`i8radix_sort`
       * :c:func:`iradix_sort`
       * :c:func:`map_get_many`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_mo_two_e_integrals_erf_ij:


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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

       * :c:func:`i2radix_sort`
       * :c:func:`i8radix_sort`
       * :c:func:`iradix_sort`
       * :c:func:`map_get_many`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: load_mo_integrals_erf:


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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

 
.. c:function:: mo_two_e_integral_erf:


    File : :file:`mo_two_e_erf_ints/map_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        subroutine mo_two_e_integrals_erf_index(i,j,k,l,i1)


    Computes an unique index for i,j,k,l integrals

 
.. c:function:: provide_all_mo_integrals_erf:


    File : :file:`mo_two_e_erf_ints/mo_bi_integrals_erf.irp.f`

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


    File : :file:`mo_two_e_erf_ints/routines_save_integrals_erf.irp.f`

    .. code:: fortran

        subroutine save_erf_two_e_integrals_mo



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`routine`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_mo_two_e_erf_ints_io_mo_two_e_integrals_erf`
       * :c:func:`ezfio_set_work_empty`
       * :c:func:`map_save_to_disk`

 
.. c:function:: save_erf_two_e_ints_mo_into_ints_mo:


    File : :file:`mo_two_e_erf_ints/routines_save_integrals_erf.irp.f`

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

