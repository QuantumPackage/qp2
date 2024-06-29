use gpu

BEGIN_PROVIDER [ type(gpu_blas), blas_handle ]
 implicit none
 BEGIN_DOC
 ! Handle for cuBLAS or RocBLAS
 END_DOC
 call gpu_blas_create(blas_handle)
END_PROVIDER

BEGIN_PROVIDER [ type(gpu_stream), gpu_default_stream ]
 implicit none
 BEGIN_DOC
 ! Default stream
 END_DOC
 gpu_default_stream%c = C_NULL_PTR
END_PROVIDER

