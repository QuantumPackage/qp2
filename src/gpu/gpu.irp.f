use gpu

BEGIN_PROVIDER [ type(gpu_blas), blas_handle ]
 implicit none
 BEGIN_DOC
 ! Handle for cuBLAS or RocBLAS
 END_DOC
 if (gpu_num > 0) then
   call gpu_blas_create(blas_handle)
 endif
! call gpu_set_stream(blas_handle, gpu_default_stream)
END_PROVIDER

 BEGIN_PROVIDER [ integer, gpu_busy, (0:gpu_num-1) ]
&BEGIN_PROVIDER [ integer, gpu_busy_max_ddot ]
&BEGIN_PROVIDER [ integer, gpu_busy_max_dgemv ]
&BEGIN_PROVIDER [ integer, gpu_busy_max_dgemm ]
 implicit none
 gpu_busy = 0
 if (gpu_num > 0) then
   gpu_busy_max_ddot = 0
   gpu_busy_max_dgemv = 4
   gpu_busy_max_dgemm = nthreads_pt2 / (gpu_num * 2)
 endif
END_PROVIDER

subroutine gpu_set_busy(igpu)
 implicit none
 BEGIN_DOC
! Set the GPU as busy
 END_DOC
 integer :: igpu
 !$OMP ATOMIC
 gpu_busy(igpu) = gpu_busy(igpu) + 1
end

subroutine gpu_unset_busy(igpu)
 implicit none
 BEGIN_DOC
! Set the GPU as busy
 END_DOC
 integer :: igpu
 !$OMP ATOMIC
 gpu_busy(igpu) = gpu_busy(igpu) - 1
end


 BEGIN_PROVIDER [ type(gpu_blas), blas_handle_mt, (0:nthreads_pt2+1) ]
&BEGIN_PROVIDER [ integer, igpu_mt, (0:nthreads_pt2+1) ]
 implicit none
 BEGIN_DOC
 ! Handle for cuBLAS or RocBLAS
 END_DOC
 integer :: i
 if (gpu_num > 0) then
   do i=0,nthreads_pt2+1
     igpu_mt(i) = mod(i, gpu_num)
     call gpu_set_device(igpu_mt(i))
     call gpu_blas_create(blas_handle_mt(i))
   enddo
   call gpu_set_device(0)
 endif
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
