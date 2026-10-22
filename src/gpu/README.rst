===
gpu
===

Architecture independent bindings to the `simple_gpu` library. The
architecture specific backend is selected by the LIB file of the
`gpu_arch` module, which is a link to `plugins/local/gpu_x86` in this
tree and provides the CPU implementation of `simple_gpu`, used when
no GPU is available.

When a GPU is present (:c:data:`gpu_num` greater than 0), the provider
:c:data:`blas_handle` creates the handle of the GPU BLAS library
(cuBLAS or RocBLAS). :c:data:`gpu_default_stream` is the default
stream and :c:data:`gpu_mem` tracks the amount of GPU memory in use,
released with :c:func:`gpu_free_memory`.

Several routines may use the same GPU at the same time. The busy state
is counted in :c:data:`gpu_busy`, and the peak values of the most
memory demanding routines are stored in :c:data:`gpu_busy_max_ddot`,
:c:data:`gpu_busy_max_dgemv` and :c:data:`gpu_busy_max_dgemm`. The C
file `gpu_affinity.c` maps each CPU core to its closest GPU.
