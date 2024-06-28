use gpu

BEGIN_PROVIDER [ type(gpu_blas), blas_handle ]
 implicit none
 BEGIN_DOC
 ! Handle for cuBLAS or RocBLAS
 END_DOC
 call gpu_blas_create(blas_handle)
END_PROVIDER


