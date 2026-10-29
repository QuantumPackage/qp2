use gpu
use gpu_affinity_mod

BEGIN_PROVIDER [ type(gpu_blas), blas_handle ]
 implicit none
 BEGIN_DOC
 ! Handle for cuBLAS/RocBLAS, or the 1-byte CPU-emulation handle when
 ! gpu_num = 0 (the simple_gpu CPU backend fully emulates the GPU BLAS
 ! API on host memory)
 END_DOC
 call gpu_blas_create(blas_handle)
! call gpu_set_stream(blas_handle, gpu_default_stream)
END_PROVIDER

 BEGIN_PROVIDER [ integer, gpu_busy, (0:gpu_num) ]
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


 BEGIN_PROVIDER [ integer, gpu_num ]
&BEGIN_PROVIDER [ type(gpu_blas), blas_handle_mt, (0:nthreads_pt2+1) ]
&BEGIN_PROVIDER [ integer, igpu_mt, (0:nthreads_pt2+1) ]
 use omp_lib
 implicit none
 BEGIN_DOC
 ! Number of usable GPUs
 ! Handle for cuBLAS or RocBLAS
 END_DOC
 integer :: i
 integer :: tid, igpu

 gpu_num = gpu_ndevices()
 igpu_mt = 0
 if (gpu_num > 0) then


  tid  = omp_get_thread_num()
  if (tid > 0) then
    call qp_bug(irp_here, tid, "blas_handle_mt provided in OpenMP section")
  endif

  ! ── Build the cpu->gpu map once, in serial ────────────────────────────
  call build_gpu_affinity_map()

!$omp parallel private(tid, igpu) num_threads(nthreads_pt2+2)
  tid  = omp_get_thread_num()         ! 0-based thread index
  igpu = get_closest_gpu()            ! closest GPU for THIS thread

  igpu_mt(tid) = igpu
  call gpu_set_device(igpu)
  call gpu_blas_create(blas_handle_mt(tid))
!$omp end parallel
   print *, 'CPU Thread/GPU mapping:'
   print *, int(igpu_mt(:),2)
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
