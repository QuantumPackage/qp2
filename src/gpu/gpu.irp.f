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

BEGIN_PROVIDER [ integer, gpu_num ]
 implicit none
 BEGIN_DOC
 ! Number of usable GPUs
 END_DOC
 gpu_num = gpu_ndevices()
END_PROVIDER

BEGIN_PROVIDER [ integer, gpu_mem ]
 implicit none
 BEGIN_DOC
 ! Total GPU memory (GB)
 END_DOC
 integer(c_size_t) :: free, total
 call gpu_get_memory(free, total)
 gpu_mem = int(total/(1024_8)**3,4)

END_PROVIDER


subroutine gpu_free_memory(value)
  use gpu
  implicit none
  BEGIN_DOC
! Returns the current used memory in gigabytes used by the current process.
  END_DOC
  double precision, intent(out) :: value
  integer(c_size_t) :: free, total
  call gpu_get_memory(free, total)

  value = dble(free)
  value = value / (1024.d0**3)
end function
