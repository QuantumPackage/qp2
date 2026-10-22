.. _module_mo_two_e_ints: 
 
.. program:: mo_two_e_ints 
 
.. default-role:: option 
 
==================
mo_two_e_ints
==================

This module computes and manages molecular orbital (MO) two-electron integrals.

The module implements several algorithms and optimization strategies for handling
the large number of two-electron integrals that arise in quantum chemical calculations.

Key Features
------------

1. **Convention Awareness**: 
   - AO integrals use the chemists' convention: (ik|jl) = (11|22)
   - MO integrals use the physicist's convention: <ij|kl> = <12|12>

2. **Efficient Storage**: Two-electron integrals are stored in a sparse map structure
   to minimize memory usage, as most integrals are zero due to symmetry and
   selection rules.

3. **Specialized Dense Arrays**:
   - For integrals with 2 identical indices such as <ij|ik> and <ij|jk>,
     dedicated dense arrays are maintained to accelerate computations
     involving single excitations
   - These arrays significantly improve performance for CI calculations
     that frequently compute single excitation amplitudes

4. **Multi-Level Caching System**:
   - Implements a sophisticated cache mechanism for frequently accessed integrals
   - Cache size is configurable via ``mo_integrals_cache_shift`` parameter
   - Automatically manages memory usage while maximizing hit rates for common integrals
   - Reduces redundant computations and improves overall performance

5. **Multiple Computation Strategies**:
   - Direct AO-to-MO transformation using DGEMM operations
   - Sparse transformation methods for memory-efficient calculation
   - Cholesky decomposition for large-scale calculations

6. **Hardware Acceleration**:
   - GPU acceleration support for faster computation
   - Multi-threading with OpenMP
   - Optimized memory access patterns

7. **Flexible Configuration**:
   - Configurable integral thresholds for pruning negligible values
   - Adjustable cache sizes for optimal performance
   - Support for both single and double precision

8. **I/O Options**:
   - Read/write integrals from/to disk for memory-constrained environments
   - Integration with Cholesky decompositions for enhanced efficiency

Implementation Details
----------------------

The module provides several key routines:

- ``get_two_e_integral(i,j,k,l,mo_integrals_map)``: Fetches a specific MO integral
- ``mo_two_e_integral(i,j,k,l)``: Direct access to MO integrals
- ``mo_two_e_integrals_jj``, ``mo_two_e_integrals_jj_exchange``, ``mo_two_e_integrals_jj_anti``: Specialized routines for 2-index integrals
- ``big_array_exchange_integrals`` and ``big_array_coulomb_integrals``: 3-index arrays for efficient exchange and Coulomb integrals computation

Storage Format
--------------

Integrals are stored in a sparse map data structure optimized for:
- Memory efficiency by storing only non-zero integrals
- Fast lookup of integrals by their 4-index notation
- Support for both symmetric and antisymmetric integral combinations
- Integration with multi-level caching for frequently accessed integrals

The caching system automatically manages memory usage while maximizing hit rates for common integrals,
reducing redundant computations and improving overall performance. The cache size can be adjusted via
the ``mo_integrals_cache_shift`` parameter in the EZFIO configuration.

The module automatically determines the most efficient computation strategy based on:
- Available memory (``qp_max_mem``)
- System architecture (CPU/GPU)
- Configuration settings (``do_mo_cholesky``, ``mo_integrals_threshold``)

Cholesky Decomposition
----------------------

When ``do_mo_cholesky`` is enabled, the module uses Cholesky decomposition to:
- Reduce memory footprint for large systems
- Enable faster integral evaluation in CI calculations
- Support efficient handling of the large number of integrals in correlated wavefunction methods

Usage Example
-------------

To use the two-electron integrals in your code:

.. code-block:: fortran

   ! Get a specific 4-index MO integral
   double precision :: integral_value
   integral_value = mo_two_e_integral(i, j, k, l)

   ! Or use the map-based approach for more control
   integral_value = get_two_e_integral(i, j, k, l, mo_integrals_map)

The module is automatically initialized during the configuration process and
provides the integrals as part of the standard quantum package environment.



 
 
 
EZFIO parameters 
---------------- 
 
.. option:: do_mo_cholesky
 
    Use Cholesky decomposition of MO integrals in CI calculations
 
    Default: False
 
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
 
.. option:: mo_cholesky_double
 
    Use double precision to build integrals from Cholesky vectors
 
    Default: True
 
 
Providers 
--------- 
 
.. c:var:: all_mo_integrals


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        logical	:: all_mo_integrals	


    Used to provide everything needed before using MO integrals
    PROVIDE all_mo_integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`mo_two_e_integrals_jj`
       * :c:data:`mo_two_e_integrals_jj`

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
       * :c:data:`ci_electronic_energy`
       * :c:data:`coef_hf_selector`
       * :c:data:`energy_by_projection`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`h_matrix_diag_all_dets`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`psi_energy_two_e_trans`

 
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

       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_num`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`coef_hf_selector`
       * :c:data:`core_fock_operator`
       * :c:data:`energy_by_projection`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`h_core_ri`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`h_matrix_diag_all_dets`
       * :c:data:`mo_two_e_integrals_jj`
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

       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_num`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`coef_hf_selector`
       * :c:data:`core_fock_operator`
       * :c:data:`energy_by_projection`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`h_core_ri`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`h_matrix_diag_all_dets`
       * :c:data:`mo_two_e_integrals_jj`
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
        integer, allocatable	:: cholesky_mo_num_split	(1:5)


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
       * :c:data:`cholesky_mo_sp`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:data:`core_fock_operator`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: cholesky_mo_num_split


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        integer	:: cholesky_mo_num	
        integer, allocatable	:: cholesky_mo_num_split	(1:5)


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
       * :c:data:`cholesky_mo_sp`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:data:`core_fock_operator`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: cholesky_mo_sp


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        real, allocatable	:: cholesky_mo_sp	(mo_num,mo_num,cholesky_mo_num)


    Cholesky vectors in MO basis in stored in single precision

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`mo_num`


 
.. c:var:: cholesky_mo_transp


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cholesky_mo_transp	(cholesky_mo_num,mo_num,mo_num)
        type(gpu_double3), allocatable	:: cholesky_mo_transp_d	(0:gpu_num)


    Cholesky vectors in MO basis. Warning: it is transposed wrt cholesky_ao:
    
    - cholesky_ao        is (ao_num^2 x cholesky_ao_num)
    
    - cholesky_mo_transp is (cholesky_mo_num x mo_num^2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`blas_handle`
       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_mo_num`
       * :c:data:`ezfio_work_dir`
       * :c:data:`gpu_num`
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
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: cholesky_mo_transp_d


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cholesky_mo_transp	(cholesky_mo_num,mo_num,mo_num)
        type(gpu_double3), allocatable	:: cholesky_mo_transp_d	(0:gpu_num)


    Cholesky vectors in MO basis. Warning: it is transposed wrt cholesky_ao:
    
    - cholesky_ao        is (ao_num^2 x cholesky_ao_num)
    
    - cholesky_mo_transp is (cholesky_mo_num x mo_num^2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`blas_handle`
       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_mo_num`
       * :c:data:`ezfio_work_dir`
       * :c:data:`gpu_num`
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
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`cholesky_semi_mo_transp_simple`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: cholesky_mo_transp_sp


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        real, allocatable	:: cholesky_mo_transp_sp	(cholesky_mo_num,mo_num,mo_num)


    Cholesky vectors in MO basis in s. Warning: it is transposed wrt cholesky_ao:
    
    - cholesky_ao        is (ao_num^2 x cholesky_ao_num)
    
    - cholesky_mo_transp is (cholesky_mo_num x mo_num^2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pxxq_array`
       * :c:data:`cholesky_mo_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: cholesky_mo_transp_sp_d


    File : :file:`mo_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        type(gpu_real3), allocatable	:: cholesky_mo_transp_sp_d	(0:gpu_num)


    Cholesky vectors in MO basis, on GPU, single precision.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`gpu_num`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pxxq_array`

 
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
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
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


 
.. c:var:: gpu_integral_buffer1


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        type(gpu_double1), allocatable	:: gpu_integral_buffer1	(0:nproc_max)
        type(gpu_real1), allocatable	:: gpu_integral_buffer_sp1	(0:nproc_max)
        type(gpu_double2), allocatable	:: gpu_integral_buffer2	(0:nproc_max)
        type(gpu_real2), allocatable	:: gpu_integral_buffer_sp2	(0:nproc_max)


    Temporary buffers for gpu_dgemv

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_num`
       * :c:data:`mo_num`
       * :c:data:`nproc_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`

 
.. c:var:: gpu_integral_buffer2


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        type(gpu_double1), allocatable	:: gpu_integral_buffer1	(0:nproc_max)
        type(gpu_real1), allocatable	:: gpu_integral_buffer_sp1	(0:nproc_max)
        type(gpu_double2), allocatable	:: gpu_integral_buffer2	(0:nproc_max)
        type(gpu_real2), allocatable	:: gpu_integral_buffer_sp2	(0:nproc_max)


    Temporary buffers for gpu_dgemv

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_num`
       * :c:data:`mo_num`
       * :c:data:`nproc_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`

 
.. c:var:: gpu_integral_buffer_sp1


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        type(gpu_double1), allocatable	:: gpu_integral_buffer1	(0:nproc_max)
        type(gpu_real1), allocatable	:: gpu_integral_buffer_sp1	(0:nproc_max)
        type(gpu_double2), allocatable	:: gpu_integral_buffer2	(0:nproc_max)
        type(gpu_real2), allocatable	:: gpu_integral_buffer_sp2	(0:nproc_max)


    Temporary buffers for gpu_dgemv

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_num`
       * :c:data:`mo_num`
       * :c:data:`nproc_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`

 
.. c:var:: gpu_integral_buffer_sp2


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        type(gpu_double1), allocatable	:: gpu_integral_buffer1	(0:nproc_max)
        type(gpu_real1), allocatable	:: gpu_integral_buffer_sp1	(0:nproc_max)
        type(gpu_double2), allocatable	:: gpu_integral_buffer2	(0:nproc_max)
        type(gpu_real2), allocatable	:: gpu_integral_buffer_sp2	(0:nproc_max)


    Temporary buffers for gpu_dgemv

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_num`
       * :c:data:`mo_num`
       * :c:data:`nproc_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`bielec_pqxx_array`
       * :c:data:`bielec_pxxq_array`

 
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

        double precision, allocatable	:: mo_integrals_cache	(0_8:mo_integrals_cache_size_8)


    Cache of MO integrals for fast access

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_num`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`core_fock_operator`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_max


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer	:: mo_integrals_cache_min	
        integer	:: mo_integrals_cache_max	
        integer	:: mo_integrals_cache_size	
        integer*8	:: mo_integrals_cache_size_8	


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
       * :c:data:`core_fock_operator`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_min


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer	:: mo_integrals_cache_min	
        integer	:: mo_integrals_cache_max	
        integer	:: mo_integrals_cache_size	
        integer*8	:: mo_integrals_cache_size_8	


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
       * :c:data:`core_fock_operator`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_size


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer	:: mo_integrals_cache_min	
        integer	:: mo_integrals_cache_max	
        integer	:: mo_integrals_cache_size	
        integer*8	:: mo_integrals_cache_size_8	


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
       * :c:data:`core_fock_operator`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: mo_integrals_cache_size_8


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer	:: mo_integrals_cache_min	
        integer	:: mo_integrals_cache_max	
        integer	:: mo_integrals_cache_size	
        integer*8	:: mo_integrals_cache_size_8	


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
       * :c:data:`core_fock_operator`
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
       * :c:data:`energy_by_projection`
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
       * :c:data:`gpu_num`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_ao_cholesky`
       * :c:data:`do_mo_cholesky`
       * :c:data:`ezfio_filename`
       * :c:data:`full_ijkl_bitmask_4`
       * :c:data:`gpu_num`
       * :c:data:`mo_cholesky_double`
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

       * :c:data:`all_mo_integrals`
       * :c:data:`bielec_pqxx_array`
       * :c:data:`core_fock_operator`
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
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_num`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`core_energy`
       * :c:data:`ref_bitmask_energy_ab`
       * :c:data:`ref_bitmask_two_e_energy`

 
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
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_num`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`core_energy`
       * :c:data:`ref_bitmask_energy_ab`
       * :c:data:`ref_bitmask_two_e_energy`

 
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
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_num`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`core_energy`
       * :c:data:`ref_bitmask_energy_ab`
       * :c:data:`ref_bitmask_two_e_energy`

 
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

       * :c:data:`core_fock_operator`
       * :c:data:`mo_two_e_integrals_jj`

 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: add_integrals_to_map:


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

 
.. c:macro:: add_integrals_to_map_cholesky:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine add_integrals_to_map_cholesky


    Adds integrals to the MO map using Cholesky vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_num`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`gpu_num`
       * :c:data:`mo_cholesky_double`
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

       * :c:func:`add_integrals_to_map_cholesky_sp`
       * :c:func:`dgemm`
       * :c:func:`gpu_allocate`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_download`
       * :c:func:`gpu_set_device`
       * :c:func:`map_append`
       * :c:func:`map_sort`
       * :c:func:`map_unique`
       * :c:func:`mo_two_e_integrals_index`

 
.. c:macro:: add_integrals_to_map_cholesky_sp:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine add_integrals_to_map_cholesky_sp


    Adds integrals to the MO map using Cholesky vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_num`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`gpu_num`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_integrals_threshold`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map_cholesky`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_allocate`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_download`
       * :c:func:`gpu_set_device`
       * :c:func:`gpu_sgemm`
       * :c:func:`map_append`
       * :c:func:`map_sort`
       * :c:func:`map_unique`
       * :c:func:`mo_two_e_integrals_index`
       * :c:func:`sgemm`

 
.. c:macro:: add_integrals_to_map_erf:


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

 
.. c:macro:: clear_mo_erf_map:


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

 
.. c:macro:: clear_mo_map:


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

 
.. c:macro:: dump_mo_integrals:


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

 
.. c:macro:: four_idx_dgemm:


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

 
.. c:macro:: four_idx_dgemm_erf:


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

 
.. c:macro:: get_mo_erf_map_size:


    File : :file:`mo_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        integer*8 function get_mo_erf_map_size()


    Returns the number of elements in the |MO| map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_erf_map`

 
.. c:macro:: get_mo_map_size:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer*8 function get_mo_map_size()


    Return the number of elements in the MO map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_map`

 
.. c:macro:: get_mo_two_e_integral_erf:


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

 
.. c:macro:: get_mo_two_e_integrals:


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
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`gpu_num`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`nproc_max`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_d1`
       * :c:func:`get_mo_two_e_integrals_i1j1`
       * :c:func:`get_mo_two_e_integrals_ij`
       * :c:func:`i_h_j_double_alpha_beta_s2_sorted`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`get_mo_two_e_integrals_cache`
       * :c:func:`gpu_dgemv`
       * :c:func:`gpu_download`
       * :c:func:`gpu_set_busy`
       * :c:func:`gpu_set_device`
       * :c:func:`gpu_sgemv`
       * :c:func:`gpu_unset_busy`
       * :c:func:`map_get`
       * :c:func:`qp_bug`
       * :c:func:`sgemv`

 
.. c:macro:: get_mo_two_e_integrals_cache:


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

 
.. c:macro:: get_mo_two_e_integrals_erf:


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

 
.. c:macro:: get_mo_two_e_integrals_erf_coulomb_ii:


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

 
.. c:macro:: get_mo_two_e_integrals_erf_exch_ii:


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

 
.. c:macro:: get_mo_two_e_integrals_erf_i1j1:


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

 
.. c:macro:: get_mo_two_e_integrals_erf_ij:


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

 
.. c:macro:: get_mo_two_e_integrals_i1j1:


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
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`gpu_num`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`nproc_max`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:func:`i_h_j_double_alpha_beta_s2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`get_mo_two_e_integrals`
       * :c:func:`get_mo_two_e_integrals_cache`
       * :c:func:`gpu_dgemv`
       * :c:func:`gpu_download`
       * :c:func:`gpu_set_busy`
       * :c:func:`gpu_set_device`
       * :c:func:`gpu_unset_busy`
       * :c:func:`qp_bug`

 
.. c:macro:: get_mo_two_e_integrals_ij:


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
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`gpu_num`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_num`
       * :c:data:`nproc_max`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pxxq_array`
       * :c:func:`get_d0`
       * :c:func:`h_s2_u_0_nstates_openmp_work_1_cholesky`
       * :c:func:`h_s2_u_0_nstates_openmp_work_2_cholesky`
       * :c:func:`h_s2_u_0_nstates_openmp_work_3_cholesky`
       * :c:func:`h_s2_u_0_nstates_openmp_work_4_cholesky`
       * :c:func:`h_s2_u_0_nstates_openmp_work_n_int_cholesky`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_mo_two_e_integrals`
       * :c:func:`get_mo_two_e_integrals_cache`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_download`
       * :c:func:`gpu_set_busy`
       * :c:func:`gpu_set_device`
       * :c:func:`gpu_sgemm`
       * :c:func:`gpu_unset_busy`
       * :c:func:`qp_bug`
       * :c:func:`sgemm`

 
.. c:macro:: get_two_e_integral:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision function get_two_e_integral(i,j,k,l,map)


    Returns one integral <ij|kl> in the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`banned_excitation`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index`

 
.. c:macro:: get_two_e_integral_cache:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision function get_two_e_integral_cache(i,j,k,l,cache,m,s)


    Returns one integral <ij|kl> in the MO basis taken from the cache

 
.. c:macro:: insert_into_mo_integrals_erf_map:


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

 
.. c:macro:: insert_into_mo_integrals_map:


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

 
.. c:macro:: load_mo_integrals:


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

 
.. c:macro:: load_mo_integrals_erf:


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

 
.. c:macro:: mo_two_e_integral:


    File : :file:`mo_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision function mo_two_e_integral(i,j,k,l)


    Returns one integral <ij|kl> in the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`mo_integrals_map`

 
.. c:macro:: mo_two_e_integral_erf:


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

 
.. c:macro:: mo_two_e_integrals_erf_index:


    File : :file:`mo_two_e_ints/mo_bi_integrals_erf.irp.f`

    .. code:: fortran

        subroutine mo_two_e_integrals_erf_index(i,j,k,l,i1)


    Computes an unique index for i,j,k,l integrals

 
.. c:macro:: mo_two_e_integrals_index:


    File : :file:`mo_two_e_ints/mo_bi_integrals.irp.f`

    .. code:: fortran

        subroutine mo_two_e_integrals_index(i,j,k,l,i1)


    Computes an unique index for i,j,k,l integrals

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map`
       * :c:func:`add_integrals_to_map_cholesky`
       * :c:func:`add_integrals_to_map_cholesky_sp`
       * :c:func:`add_integrals_to_map_erf`
       * :c:func:`four_idx_dgemm`
       * :c:func:`four_idx_dgemm_erf`

 
.. c:macro:: provide_all_mo_integrals_erf:


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

 
.. c:macro:: save_erf_two_e_integrals_mo:


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

 
.. c:macro:: save_erf_two_e_ints_mo_into_ints_mo:


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

