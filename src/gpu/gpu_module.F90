module gpu
  use, intrinsic :: iso_c_binding
  implicit none

! Data types
! ----------

  type gpu_double1
    type(c_ptr) :: c
    double precision, pointer :: f(:)
  end type

  type gpu_double2
    type(c_ptr) :: c
    double precision, pointer :: f(:,:)
  end type

  type gpu_double3
    type(c_ptr) :: c
    double precision, pointer :: f(:,:,:)
  end type

  type gpu_double4
    type(c_ptr) :: c
    double precision, pointer :: f(:,:,:,:)
  end type

  type gpu_double5
    type(c_ptr) :: c
    double precision, pointer :: f(:,:,:,:,:)
  end type

  type gpu_double6
    type(c_ptr) :: c
    double precision, pointer :: f(:,:,:,:,:,:)
  end type


  type gpu_blas
    type(c_ptr) :: c
  end type

  type gpu_stream
    type(c_ptr) :: c
  end type


! C interfaces
! ------------

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

    subroutine gpu_deallocate_c(ptr) bind(C, name='gpu_deallocate')
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

    subroutine gpu_stream_create_c(stream) bind(C, name='gpu_stream_create')
      import
      type(c_ptr) :: stream
    end subroutine

    subroutine gpu_stream_destroy_c(stream) bind(C, name='gpu_stream_destroy')
      import
      type(c_ptr) :: stream
    end subroutine

    subroutine gpu_set_stream_c(handle, stream) bind(C, name='gpu_set_stream')
      import
      type(c_ptr) :: handle, stream
    end subroutine

    subroutine gpu_synchronize() bind(C)
      import
    end subroutine

    subroutine gpu_blas_create_c(handle) bind(C, name='gpu_blas_create')
      import
      type(c_ptr) :: handle
    end subroutine

    subroutine gpu_blas_destroy_c(handle) bind(C, name='gpu_blas_destroy')
      import
      type(c_ptr) :: handle
    end subroutine

    subroutine gpu_ddot_c(handle, n, dx, incx, dy, incy, res) bind(C, name='gpu_ddot')
      import
      type(c_ptr), intent(in)        :: handle
      integer(c_int64_t), value      :: n, incx, incy
      type(c_ptr), intent(in), value :: dx, dy
      real(c_double), intent(out)    :: res
    end subroutine

    subroutine gpu_sdot_c(handle, n, dx, incx, dy, incy, res) bind(C, name='gpu_sdot')
      import
      type(c_ptr), intent(in)        :: handle
      integer(c_int64_t), value      :: n, incx, incy
      type(c_ptr), intent(in), value :: dx, dy
      real(c_float), intent(out)  :: res
    end subroutine

    subroutine gpu_dgeam_c(handle, transa, transb, m, n, alpha, a, lda, beta, &
      b, ldb, c, ldc) bind(C, name='gpu_dgeam')
      import
      type(c_ptr), intent(in)        :: handle
      character(c_char), intent(in), value :: transa, transb
      integer(c_int64_t), intent(in), value :: m, n, lda, ldb, ldc
      real(c_double), intent(in), value :: alpha, beta
      type(c_ptr), value :: a, b, c
    end subroutine

    subroutine gpu_sgeam_c(handle, transa, transb, m, n, alpha, a, lda, beta, &
      b, ldb, c, ldc) bind(C, name='gpu_sgeam')
      import
      type(c_ptr), intent(in)        :: handle
      character(c_char), intent(in), value :: transa, transb
      integer(c_int64_t), intent(in), value :: m, n, lda, ldb, ldc
      real(c_float), intent(in), value :: alpha, beta
      type(c_ptr), value :: a, b, c
    end subroutine

  end interface


! Polymorphic interfaces
! ----------------------

  interface gpu_allocate
    procedure gpu_allocate_double1  &
             ,gpu_allocate_double2  &
             ,gpu_allocate_double3  &
             ,gpu_allocate_double4  &
             ,gpu_allocate_double5  &
             ,gpu_allocate_double6
  end interface gpu_allocate

  interface gpu_deallocate
    procedure gpu_deallocate_double1  &
             ,gpu_deallocate_double2  &
             ,gpu_deallocate_double3  &
             ,gpu_deallocate_double4  &
             ,gpu_deallocate_double5  &
             ,gpu_deallocate_double6
  end interface gpu_deallocate

  interface gpu_upload
    procedure gpu_upload_double1  &
             ,gpu_upload_double2  &
             ,gpu_upload_double3  &
             ,gpu_upload_double4  &
             ,gpu_upload_double5  &
             ,gpu_upload_double6
  end interface gpu_upload

  interface gpu_download
    procedure gpu_download_double1  &
             ,gpu_download_double2  &
             ,gpu_download_double3  &
             ,gpu_download_double4  &
             ,gpu_download_double5  &
             ,gpu_download_double6
  end interface gpu_download

  interface gpu_copy
    procedure gpu_copy_double1  &
             ,gpu_copy_double2  &
             ,gpu_copy_double3  &
             ,gpu_copy_double4  &
             ,gpu_copy_double5  &
             ,gpu_copy_double6
  end interface gpu_copy


  contains


! gpu_allocate
! ------------

    subroutine gpu_allocate_double1(ptr, s)
      implicit none
      type(gpu_double1), intent(inout) :: ptr
      integer, intent(in) :: s

      call gpu_allocate_c(ptr%c, s*8_8)
      call c_f_pointer(ptr%c, ptr%f, (/ s /))
    end subroutine

    subroutine gpu_allocate_double2(ptr, s1, s2)
      implicit none
      type(gpu_double2), intent(inout) :: ptr
      integer, intent(in) :: s1, s2

      call gpu_allocate_c(ptr%c, s1*s2*8_8)
      call c_f_pointer(ptr%c, ptr%f, (/ s1, s2 /))
    end subroutine

    subroutine gpu_allocate_double3(ptr, s1, s2, s3)
      implicit none
      type(gpu_double3), intent(inout) :: ptr
      integer, intent(in) :: s1, s2, s3

      call gpu_allocate_c(ptr%c, s1*s2*s3*8_8)
      call c_f_pointer(ptr%c, ptr%f, (/ s1, s2, s3 /))
    end subroutine

    subroutine gpu_allocate_double4(ptr, s1, s2, s3, s4)
      implicit none
      type(gpu_double4), intent(inout) :: ptr
      integer, intent(in) :: s1, s2, s3, s4

      call gpu_allocate_c(ptr%c, s1*s2*s3*s4*8_8)
      call c_f_pointer(ptr%c, ptr%f, (/ s1, s2, s3, s4 /))
    end subroutine

    subroutine gpu_allocate_double5(ptr, s1, s2, s3, s4, s5)
      implicit none
      type(gpu_double5), intent(inout) :: ptr
      integer, intent(in) :: s1, s2, s3, s4, s5

      call gpu_allocate_c(ptr%c, s1*s2*s3*s4*s5*8_8)
      call c_f_pointer(ptr%c, ptr%f, (/ s1, s2, s3, s4, s5 /))
    end subroutine

    subroutine gpu_allocate_double6(ptr, s1, s2, s3, s4, s5, s6)
      implicit none
      type(gpu_double6), intent(inout) :: ptr
      integer, intent(in) :: s1, s2, s3, s4, s5, s6

      call gpu_allocate_c(ptr%c, s1*s2*s3*s4*s5*s6*8_8)
      call c_f_pointer(ptr%c, ptr%f, (/ s1, s2, s3, s4, s5, s6 /))
    end subroutine


! gpu_deallocate
! --------------

    subroutine gpu_deallocate_double1(ptr)
      implicit none
      type(gpu_double1), intent(inout) :: ptr
      call gpu_deallocate_c(ptr%c)
      NULLIFY(ptr%f)
    end subroutine

    subroutine gpu_deallocate_double2(ptr)
      implicit none
      type(gpu_double2), intent(inout) :: ptr
      call gpu_deallocate_c(ptr%c)
      NULLIFY(ptr%f)
    end subroutine

    subroutine gpu_deallocate_double3(ptr)
      implicit none
      type(gpu_double3), intent(inout) :: ptr
      call gpu_deallocate_c(ptr%c)
      NULLIFY(ptr%f)
    end subroutine

    subroutine gpu_deallocate_double4(ptr)
      implicit none
      type(gpu_double4), intent(inout) :: ptr
      call gpu_deallocate_c(ptr%c)
      NULLIFY(ptr%f)
    end subroutine

    subroutine gpu_deallocate_double5(ptr)
      implicit none
      type(gpu_double5), intent(inout) :: ptr
      call gpu_deallocate_c(ptr%c)
      NULLIFY(ptr%f)
    end subroutine

    subroutine gpu_deallocate_double6(ptr)
      implicit none
      type(gpu_double6), intent(inout) :: ptr
      call gpu_deallocate_c(ptr%c)
      NULLIFY(ptr%f)
    end subroutine


! gpu_upload
! ----------

    subroutine gpu_upload_double1(cpu_ptr, gpu_ptr)
      implicit none
      double precision, intent(in)     :: cpu_ptr(:)
      type(gpu_double1), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, 8_8*size(gpu_ptr%f))
    end subroutine

    subroutine gpu_upload_double2(cpu_ptr, gpu_ptr)
      implicit none
      double precision, intent(in)     :: cpu_ptr(:,:)
      type(gpu_double2), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, product(shape(gpu_ptr%f)*1_8)*8_8)
    end subroutine

    subroutine gpu_upload_double3(cpu_ptr, gpu_ptr)
      implicit none
      double precision, intent(in)     :: cpu_ptr(:,:,:)
      type(gpu_double3), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, product(shape(gpu_ptr%f)*1_8)*8_8)
    end subroutine

    subroutine gpu_upload_double4(cpu_ptr, gpu_ptr)
      implicit none
      double precision, intent(in)     :: cpu_ptr(:,:,:,:)
      type(gpu_double4), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, product(shape(gpu_ptr%f)*1_8)*8_8)
    end subroutine

    subroutine gpu_upload_double5(cpu_ptr, gpu_ptr)
      implicit none
      double precision, intent(in)     :: cpu_ptr(:,:,:,:,:)
      type(gpu_double5), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, product(shape(gpu_ptr%f)*1_8)*8_8)
    end subroutine

    subroutine gpu_upload_double6(cpu_ptr, gpu_ptr)
      implicit none
      double precision, intent(in)     :: cpu_ptr(:,:,:,:,:,:)
      type(gpu_double6), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, product(shape(gpu_ptr%f)*1_8)*8_8)
    end subroutine


! gpu_download
! ------------

    subroutine gpu_download_double1(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double1), intent(in)  :: gpu_ptr
      double precision, intent(in)   :: cpu_ptr(:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*size(gpu_ptr%f))
    end subroutine

    subroutine gpu_download_double2(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double2), intent(in)  :: gpu_ptr
      double precision, intent(in)   :: cpu_ptr(:,:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*product(shape(gpu_ptr%f)*1_8))
    end subroutine

    subroutine gpu_download_double3(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double3), intent(in)  :: gpu_ptr
      double precision, intent(in)   :: cpu_ptr(:,:,:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*product(shape(gpu_ptr%f)*1_8))
    end subroutine

    subroutine gpu_download_double4(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double4), intent(in)  :: gpu_ptr
      double precision, intent(in)   :: cpu_ptr(:,:,:,:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*product(shape(gpu_ptr%f)*1_8))
    end subroutine

    subroutine gpu_download_double5(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double5), intent(in)  :: gpu_ptr
      double precision, intent(in)   :: cpu_ptr(:,:,:,:,:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*product(shape(gpu_ptr%f)*1_8))
    end subroutine

    subroutine gpu_download_double6(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double6), intent(in)  :: gpu_ptr
      double precision, intent(in)   :: cpu_ptr(:,:,:,:,:,:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*product(shape(gpu_ptr%f)*1_8))
    end subroutine

! gpu_copy
! --------

    subroutine gpu_copy_double1(gpu_ptr_src, gpu_ptr_dest)
      implicit none
      type(gpu_double1), intent(in)        :: gpu_ptr_src
      type(gpu_double1), intent(in)        :: gpu_ptr_dest
      call gpu_copy_c(gpu_ptr_src%c, gpu_ptr_dest%c, 8_8*size(gpu_ptr_dest%f))
    end subroutine

    subroutine gpu_copy_double2(gpu_ptr_src, gpu_ptr_dest)
      implicit none
      type(gpu_double2), intent(in)        :: gpu_ptr_src
      type(gpu_double2), intent(in)        :: gpu_ptr_dest
      call gpu_copy_c(gpu_ptr_src%c, gpu_ptr_dest%c, 8_8*product(shape(gpu_ptr_dest%f)*1_8))
    end subroutine

    subroutine gpu_copy_double3(gpu_ptr_src, gpu_ptr_dest)
      implicit none
      type(gpu_double3), intent(in)        :: gpu_ptr_src
      type(gpu_double3), intent(in)        :: gpu_ptr_dest
      call gpu_copy_c(gpu_ptr_src%c, gpu_ptr_dest%c, 8_8*product(shape(gpu_ptr_dest%f)*1_8))
    end subroutine

    subroutine gpu_copy_double4(gpu_ptr_src, gpu_ptr_dest)
      implicit none
      type(gpu_double4), intent(in)        :: gpu_ptr_src
      type(gpu_double4), intent(in)        :: gpu_ptr_dest
      call gpu_copy_c(gpu_ptr_src%c, gpu_ptr_dest%c, 8_8*product(shape(gpu_ptr_dest%f)*1_8))
    end subroutine

    subroutine gpu_copy_double5(gpu_ptr_src, gpu_ptr_dest)
      implicit none
      type(gpu_double5), intent(in)        :: gpu_ptr_src
      type(gpu_double5), intent(in)        :: gpu_ptr_dest
      call gpu_copy_c(gpu_ptr_src%c, gpu_ptr_dest%c, 8_8*product(shape(gpu_ptr_dest%f)*1_8))
    end subroutine

    subroutine gpu_copy_double6(gpu_ptr_src, gpu_ptr_dest)
      implicit none
      type(gpu_double6), intent(in)        :: gpu_ptr_src
      type(gpu_double6), intent(in)        :: gpu_ptr_dest
      call gpu_copy_c(gpu_ptr_src%c, gpu_ptr_dest%c, 8_8*product(shape(gpu_ptr_dest%f)*1_8))
    end subroutine


! gpu_stream
! ----------

    subroutine gpu_stream_create(stream)
      import
      type(gpu_stream) :: stream
      call gpu_stream_create_c(stream%c)
    end subroutine

    subroutine gpu_stream_destroy(stream)
      import
      type(gpu_stream) :: stream
      call gpu_stream_destroy_c(stream%c)
    end subroutine

    subroutine gpu_set_stream(handle, stream)
      import
      type(gpu_blas)   :: handle
      type(gpu_stream) :: stream
      call gpu_set_stream_c(handle%c, stream%c)
    end subroutine


! gpu_blas
! --------

    subroutine gpu_blas_create(handle)
      import
      type(gpu_blas) :: handle
      call gpu_blas_create_c(handle%c)
    end subroutine

    subroutine gpu_blas_destroy(handle)
      import
      type(gpu_blas) :: handle
      call gpu_blas_destroy_c(handle%c)
    end subroutine


end module



! dot
! ---

subroutine gpu_ddot(handle, n, dx, incx, dy, incy, res)
  use gpu
  type(gpu_blas), intent(in)     :: handle
  integer*8                      :: n, incx, incy
  double precision, intent(in)   :: dx(*), dy(*)
  double precision, intent(out)    :: res
  call gpu_ddot_c(handle%c, n, c_loc(dx), incx, c_loc(dy), incy, res)
end subroutine

subroutine gpu_sdot(handle, n, dx, incx, dy, incy, res)
  use gpu
  type(gpu_blas), intent(in)     :: handle
  integer*8                      :: n, incx, incy
  real, intent(in)               :: dx(*), dy(*)
  real, intent(out)              :: res
  call gpu_sdot_c(handle%c, n, c_loc(dx), incx, c_loc(dy), incy, res)
end subroutine


! geam
! ----

subroutine gpu_dgeam(handle, transa, transb, m, n, alpha, a, lda, beta, &
  b, ldb, c, ldc)
  use gpu
  type(gpu_blas), intent(in)   :: handle
  character, intent(in)        :: transa, transb
  integer*8, intent(in)        :: m, n, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  double precision             :: a(lda,*), b(ldb,*), c(ldc,*)
  call gpu_dgeam_c(handle%c, transa, transb, m, n, alpha, c_loc(a), lda, beta, &
        c_loc(b), ldb, c_loc(c), ldc)
end subroutine

subroutine gpu_sgeam(handle, transa, transb, m, n, alpha, a, lda, beta, &
  b, ldb, c, ldc)
 use gpu 
  type(gpu_blas), intent(in)   :: handle
  character, intent(in)        :: transa, transb
  integer*8, intent(in)        :: m, n, lda, ldb, ldc
  real, intent(in)             :: alpha, beta
  real                         :: a(lda,*), b(ldb,*), c(ldc,*)
  call gpu_sgeam_c(handle%c, transa, transb, m, n, alpha, c_loc(a), lda, beta, &
        c_loc(b), ldb, c_loc(c), ldc)
end subroutine

