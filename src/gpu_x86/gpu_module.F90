module gpu
  use, intrinsic :: iso_c_binding, only : c_int32_t, c_int64_t, c_double, c_size_t, c_char
  implicit none

  interface
    integer function gpu_ndevices() bind(C)
    end function

    subroutine gpu_set_device(id) bind(C)
      import
      integer(c_int32_t), value :: id
    end subroutine

    subroutine gpu_allocate_c(ptr, n) bind(C, name='gpu_allocate')
      import
      type(c_ptr) :: ptr
      integer(c_int64_t), value :: n
    end subroutine

    subroutine gpu_free_c(ptr) bind(C, name='gpu_free')
      import
      type(c_ptr) :: ptr
    end subroutine

    subroutine gpu_upload_c(cpu_ptr, gpu_ptr, n) bind(C, name='gpu_upload')
      import
      type(c_ptr), value :: cpu_ptr
      type(c_ptr), value :: gpu_ptr
      integer(c_int64_t), value :: n
    end subroutine

    subroutine gpu_download_c(gpu_ptr, cpu_ptr, n) bind(C, name='gpu_download')
      import
      type(c_ptr), value :: gpu_ptr
      type(c_ptr), value :: cpu_ptr
      integer(c_int64_t), value :: n
    end subroutine

    subroutine gpu_copy_c(gpu_ptr_src, gpu_ptr_dest, n) bind(C, name='gpu_copy')
      import
      type(c_ptr), value :: gpu_ptr_src
      type(c_ptr), value :: gpu_ptr_dest
      integer(c_int64_t), value :: n
    end subroutine

    subroutine gpu_stream_create(stream) bind(C)
      import
      type(c_ptr) :: stream
    end subroutine

    subroutine gpu_stream_destroy(stream) bind(C)
      import
      type(c_ptr) :: stream
    end subroutine

    subroutine gpu_set_stream(handle, stream) bind(C)
      import
      type(c_ptr) :: handle, stream
    end subroutine

    subroutine gpu_synchronize()
    end subroutine

    subroutine gpu_blas_create(handle) bind(C)
      import
      type(c_ptr) :: handle
    end subroutine

    subroutine gpu_blas_destroy(handle) bind(C)
      import
      type(c_ptr) :: handle
    end subroutine

    subroutine gpu_ddot(handle, n, dx, incx, dy, incy, res) bind(C)
      import
      type(c_ptr), intent(in)     :: handle
      integer(c_int64_t), value   :: n, incx, incy
      real(c_double), intent(in)  :: dx(*), dy(*)
      real(c_double), intent(out) :: res
    end subroutine

    subroutine gpu_sdot(handle, n, dx, incx, dy, incy, res) bind(C)
      import
      type(c_ptr), intent(in)     :: handle
      integer(c_int64_t), value   :: n, incx, incy
      real(c_float), intent(in)   :: dx(*), dy(*)
      real(c_float), intent(out)  :: res
    end subroutine

  end interface

end module

subroutine gpu_allocate_double(ptr, s)
  use gpu
  implicit none
  double precision, pointer, intent(inout) :: ptr
  integer*8, intent(in) :: s(*)
  type(c_ptr) :: cptr

  call gpu_allocate_c(cptr, sum(s)*8_8)
  call c_f_pointer(cptr, ptr, s)
end subroutine

subroutine gpu_free_double(ptr)
  use gpu
  implicit none
  double precision, pointer, intent(inout) :: ptr
  type(c_ptr) :: cptr
  cptr = cloc(ptr)
  call gpu_free(cptr)
  NULLIFY(ptr)
end subroutine

subroutine gpu_upload_double(cpu_ptr, gpu_ptr, n)
  use gpu
  implicit none
  double precision, intent(in)   :: cpu_ptr(*)
  double precision, intent(out)  :: gpu_ptr(*)
  integer(c_int64_t), intent(in) :: n
  call gpu_upload_c(cpu_ptr, gpu_ptr, 8_8*n)
end subroutine

subroutine gpu_download_double(gpu_ptr, cpu_ptr, n)
  use gpu
  implicit none
  double precision, intent(in)   :: gpu_ptr(*)
  double precision, intent(out)  :: cpu_ptr(*)
  integer(c_int64_t), intent(in) :: n
  call gpu_download_c(gpu_ptr, cpu_ptr, 8_8*n)
end subroutine

subroutine gpu_copy_double(gpu_ptr_src, gpu_ptr_dest, n)
  use gpu
  implicit none
  double precision, intent(in)   :: gpu_ptr_src(*)
  double precision, intent(out)  :: gpu_ptr_dest(*)
  integer(c_int64_t), intent(in) :: n
  call gpu_copy_c(gpu_ptr_src, gpu_ptr_dest, 8_8*n)
end subroutine

