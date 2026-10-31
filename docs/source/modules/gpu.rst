.. _module_gpu: 
 
.. program:: gpu 
 
.. default-role:: option 
 
===
gpu
===

Architecture independent bindings to the `simple_gpu` library. The
architecture specific backend is selected by the LIB file of the
`gpu_arch` module, which is a link to `plugins/local/gpu_x86` in this
tree and provides the CPU implementation of `simple_gpu`, used when
no GPU is available.

The provider :c:data:`blas_handle` always creates a handle of the GPU
BLAS library: cuBLAS or RocBLAS when a GPU is present (:c:data:`gpu_num`
greater than 0), or the 1-byte CPU-emulation handle when no GPU is
available, since the simple_gpu CPU backend fully emulates the GPU BLAS
API on host memory. :c:data:`gpu_default_stream` is the default
stream and :c:data:`gpu_mem` tracks the amount of GPU memory in use,
released with :c:func:`gpu_free_memory`.

Several routines may use the same GPU at the same time. The busy state
is counted in :c:data:`gpu_busy`, and the peak values of the most
memory demanding routines are stored in :c:data:`gpu_busy_max_ddot`,
:c:data:`gpu_busy_max_dgemv` and :c:data:`gpu_busy_max_dgemm`. The C
file `gpu_affinity.c` maps each CPU core to its closest GPU.
 
 
 
Providers 
--------- 
 
.. c:var:: int blas_handle


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        type(gpu_blas)	:: blas_handle	


    Handle for cuBLAS/RocBLAS, or the 1-byte CPU-emulation handle when
    gpu_num = 0 (the simple_gpu CPU backend fully emulates the GPU BLAS
    API on host memory)

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: int blas_handle_mt


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        integer	:: gpu_num	
        type(gpu_blas), allocatable	:: blas_handle_mt	(0:nthreads_pt2+1)
        integer, allocatable	:: igpu_mt	(0:nthreads_pt2+1)


    Number of usable GPUs
    Handle for cuBLAS or RocBLAS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_pt2`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: int gpu_busy


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        integer, allocatable	:: gpu_busy	(0:gpu_num)
        integer	:: gpu_busy_max_ddot	
        integer	:: gpu_busy_max_dgemv	
        integer	:: gpu_busy_max_dgemm	


    Number of ongoing GPU operations, one entry per GPU,
    indexed from 0 to gpu_num. gpu_set_busy and gpu_unset_busy
    increment and decrement it (OpenMP atomic). Callers compare
    it with the corresponding cap before a BLAS call. The caps
    gpu_busy_max_ddot, gpu_busy_max_dgemv and gpu_busy_max_dgemm
    are set when gpu_num > 0 to 0, 4 and nthreads_pt2/(gpu_num*2).
    The three caps are set from empirical measurements
    (author, 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_num`
       * :c:data:`nthreads_pt2`


 
.. c:var:: int gpu_busy_max_ddot


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        integer, allocatable	:: gpu_busy	(0:gpu_num)
        integer	:: gpu_busy_max_ddot	
        integer	:: gpu_busy_max_dgemv	
        integer	:: gpu_busy_max_dgemm	


    Number of ongoing GPU operations, one entry per GPU,
    indexed from 0 to gpu_num. gpu_set_busy and gpu_unset_busy
    increment and decrement it (OpenMP atomic). Callers compare
    it with the corresponding cap before a BLAS call. The caps
    gpu_busy_max_ddot, gpu_busy_max_dgemv and gpu_busy_max_dgemm
    are set when gpu_num > 0 to 0, 4 and nthreads_pt2/(gpu_num*2).
    The three caps are set from empirical measurements
    (author, 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_num`
       * :c:data:`nthreads_pt2`


 
.. c:var:: int gpu_busy_max_dgemm


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        integer, allocatable	:: gpu_busy	(0:gpu_num)
        integer	:: gpu_busy_max_ddot	
        integer	:: gpu_busy_max_dgemv	
        integer	:: gpu_busy_max_dgemm	


    Number of ongoing GPU operations, one entry per GPU,
    indexed from 0 to gpu_num. gpu_set_busy and gpu_unset_busy
    increment and decrement it (OpenMP atomic). Callers compare
    it with the corresponding cap before a BLAS call. The caps
    gpu_busy_max_ddot, gpu_busy_max_dgemv and gpu_busy_max_dgemm
    are set when gpu_num > 0 to 0, 4 and nthreads_pt2/(gpu_num*2).
    The three caps are set from empirical measurements
    (author, 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_num`
       * :c:data:`nthreads_pt2`


 
.. c:var:: int gpu_busy_max_dgemv


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        integer, allocatable	:: gpu_busy	(0:gpu_num)
        integer	:: gpu_busy_max_ddot	
        integer	:: gpu_busy_max_dgemv	
        integer	:: gpu_busy_max_dgemm	


    Number of ongoing GPU operations, one entry per GPU,
    indexed from 0 to gpu_num. gpu_set_busy and gpu_unset_busy
    increment and decrement it (OpenMP atomic). Callers compare
    it with the corresponding cap before a BLAS call. The caps
    gpu_busy_max_ddot, gpu_busy_max_dgemv and gpu_busy_max_dgemm
    are set when gpu_num > 0 to 0, 4 and nthreads_pt2/(gpu_num*2).
    The three caps are set from empirical measurements
    (author, 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_num`
       * :c:data:`nthreads_pt2`


 
.. c:var:: int gpu_default_stream


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        type(gpu_stream)	:: gpu_default_stream	


    Default stream


 
.. c:var:: int gpu_mem


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        integer	:: gpu_mem	


    Total GPU memory (GB)


 
.. c:var:: int gpu_num


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        integer	:: gpu_num	
        type(gpu_blas), allocatable	:: blas_handle_mt	(0:nthreads_pt2+1)
        integer, allocatable	:: igpu_mt	(0:nthreads_pt2+1)


    Number of usable GPUs
    Handle for cuBLAS or RocBLAS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_pt2`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
.. c:var:: int igpu_mt


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        integer	:: gpu_num	
        type(gpu_blas), allocatable	:: blas_handle_mt	(0:nthreads_pt2+1)
        integer, allocatable	:: igpu_mt	(0:nthreads_pt2+1)


    Number of usable GPUs
    Handle for cuBLAS or RocBLAS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_pt2`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_jj`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void gpu_free_memory(void)


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        subroutine gpu_free_memory(value)


    Returns the current used memory in gigabytes used by the current process.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_ao_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_get_memory`

 
.. c:function:: void gpu_set_busy(void)


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        subroutine gpu_set_busy(igpu)


    Set the GPU as busy

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_busy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_mo_two_e_integrals`
       * :c:func:`get_mo_two_e_integrals_i1j1`
       * :c:func:`get_mo_two_e_integrals_ij`

 
.. c:function:: void gpu_unset_busy(void)


    File : :file:`gpu/gpu.irp.f`

    .. code:: fortran

        subroutine gpu_unset_busy(igpu)


    Set the GPU as busy

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gpu_busy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_mo_two_e_integrals`
       * :c:func:`get_mo_two_e_integrals_i1j1`
       * :c:func:`get_mo_two_e_integrals_ij`

