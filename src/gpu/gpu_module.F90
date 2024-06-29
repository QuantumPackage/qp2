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
    logical(c_bool) function no_gpu() bind(C)
      import
    end function

    integer function gpu_ndevices() bind(C)
      import
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
      type(c_ptr), value :: handle, stream
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
      type(c_ptr), value, intent(in) :: handle
      integer(c_int64_t), value      :: n, incx, incy
      type(c_ptr), value             :: dx, dy
      real(c_double), intent(out)    :: res
    end subroutine

    subroutine gpu_sdot_c(handle, n, dx, incx, dy, incy, res) bind(C, name='gpu_sdot')
      import
      type(c_ptr), value, intent(in) :: handle
      integer(c_int64_t), value      :: n, incx, incy
      type(c_ptr), intent(in), value :: dx, dy
      real(c_float), intent(out)  :: res
    end subroutine

    subroutine gpu_dgeam_c(handle, transa, transb, m, n, alpha, a, lda, beta, &
      b, ldb, c, ldc) bind(C, name='gpu_dgeam')
      import
      type(c_ptr), value, intent(in)        :: handle
      character(c_char), intent(in), value  :: transa, transb
      integer(c_int64_t), intent(in), value :: m, n, lda, ldb, ldc
      real(c_double), intent(in), value :: alpha, beta
      type(c_ptr), value :: a, b, c
    end subroutine

    subroutine gpu_sgeam_c(handle, transa, transb, m, n, alpha, a, lda, beta, &
      b, ldb, c, ldc) bind(C, name='gpu_sgeam')
      import
      type(c_ptr), value, intent(in)        :: handle
      character(c_char), intent(in), value  :: transa, transb
      integer(c_int64_t), intent(in), value :: m, n, lda, ldb, ldc
      real(c_float), intent(in), value :: alpha, beta
      type(c_ptr), value :: a, b, c
    end subroutine

    subroutine gpu_dgemm_c(handle, transa, transb, m, n, k, alpha, a, lda, &
      b, ldb, beta, c, ldc) bind(C, name='gpu_dgemm')
      import
      type(c_ptr), value, intent(in)        :: handle
      character(c_char), intent(in), value  :: transa, transb
      integer(c_int64_t), intent(in), value :: m, n, k, lda, ldb, ldc
      real(c_double), intent(in), value :: alpha, beta
      type(c_ptr), value :: a, b, c
    end subroutine

    subroutine gpu_sgemm_c(handle, transa, transb, m, n, k, alpha, a, lda, &
      b, ldb, beta, c, ldc) bind(C, name='gpu_sgemm')
      import
      type(c_ptr), value, intent(in)        :: handle
      character(c_char), intent(in), value  :: transa, transb
      integer(c_int64_t), intent(in), value :: m, n, k, lda, ldb, ldc
      real(c_float), intent(in), value :: alpha, beta
      type(c_ptr), value :: a, b, c
    end subroutine

  end interface


! Polymorphic interfaces
! ----------------------

  interface gpu_allocate
    procedure gpu_allocate_double1     &
             ,gpu_allocate_double2     &
             ,gpu_allocate_double3     &
             ,gpu_allocate_double4     &
             ,gpu_allocate_double5     &
             ,gpu_allocate_double6     &
             ,gpu_allocate_double1_64  &
             ,gpu_allocate_double2_64  &
             ,gpu_allocate_double3_64  &
             ,gpu_allocate_double4_64  &
             ,gpu_allocate_double5_64  &
             ,gpu_allocate_double6_64
  end interface gpu_allocate

  interface gpu_deallocate
    procedure gpu_deallocate_double1     &
             ,gpu_deallocate_double2     &
             ,gpu_deallocate_double3     &
             ,gpu_deallocate_double4     &
             ,gpu_deallocate_double5     &
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


    subroutine gpu_allocate_double1_64(ptr, s)
      implicit none
      type(gpu_double1), intent(inout) :: ptr
      integer*8, intent(in) :: s

      call gpu_allocate_c(ptr%c, s)
      call c_f_pointer(ptr%c, ptr%f, (/ s /))
    end subroutine

    subroutine gpu_allocate_double2_64(ptr, s1, s2)
      implicit none
      type(gpu_double2), intent(inout) :: ptr
      integer*8, intent(in) :: s1, s2

      call gpu_allocate_c(ptr%c, s1*s2*8_8)
      call c_f_pointer(ptr%c, ptr%f, (/ s1, s2 /))
    end subroutine

    subroutine gpu_allocate_double3_64(ptr, s1, s2, s3)
      implicit none
      type(gpu_double3), intent(inout) :: ptr
      integer*8, intent(in) :: s1, s2, s3

      call gpu_allocate_c(ptr%c, s1*s2*s3*8_8)
      call c_f_pointer(ptr%c, ptr%f, (/ s1, s2, s3 /))
    end subroutine

    subroutine gpu_allocate_double4_64(ptr, s1, s2, s3, s4)
      implicit none
      type(gpu_double4), intent(inout) :: ptr
      integer*8, intent(in) :: s1, s2, s3, s4

      call gpu_allocate_c(ptr%c, s1*s2*s3*s4*8_8)
      call c_f_pointer(ptr%c, ptr%f, (/ s1, s2, s3, s4 /))
    end subroutine

    subroutine gpu_allocate_double5_64(ptr, s1, s2, s3, s4, s5)
      implicit none
      type(gpu_double5), intent(inout) :: ptr
      integer*8, intent(in) :: s1, s2, s3, s4, s5

      call gpu_allocate_c(ptr%c, s1*s2*s3*s4*s5*8_8)
      call c_f_pointer(ptr%c, ptr%f, (/ s1, s2, s3, s4, s5 /))
    end subroutine

    subroutine gpu_allocate_double6_64(ptr, s1, s2, s3, s4, s5, s6)
      implicit none
      type(gpu_double6), intent(inout) :: ptr
      integer*8, intent(in) :: s1, s2, s3, s4, s5, s6

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
      double precision, target, intent(in)     :: cpu_ptr(*)
      type(gpu_double1), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, 8_8*size(gpu_ptr%f))
    end subroutine

    subroutine gpu_upload_double2(cpu_ptr, gpu_ptr)
      implicit none
      double precision, target, intent(in)     :: cpu_ptr(:,:)
      type(gpu_double2), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, product(shape(gpu_ptr%f)*1_8)*8_8)
    end subroutine

    subroutine gpu_upload_double3(cpu_ptr, gpu_ptr)
      implicit none
      double precision, target, intent(in)     :: cpu_ptr(:,:,:)
      type(gpu_double3), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, product(shape(gpu_ptr%f)*1_8)*8_8)
    end subroutine

    subroutine gpu_upload_double4(cpu_ptr, gpu_ptr)
      implicit none
      double precision, target, intent(in)     :: cpu_ptr(:,:,:,:)
      type(gpu_double4), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, product(shape(gpu_ptr%f)*1_8)*8_8)
    end subroutine

    subroutine gpu_upload_double5(cpu_ptr, gpu_ptr)
      implicit none
      double precision, target, intent(in)     :: cpu_ptr(:,:,:,:,:)
      type(gpu_double5), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, product(shape(gpu_ptr%f)*1_8)*8_8)
    end subroutine

    subroutine gpu_upload_double6(cpu_ptr, gpu_ptr)
      implicit none
      double precision, target, intent(in)     :: cpu_ptr(:,:,:,:,:,:)
      type(gpu_double6), intent(in)    :: gpu_ptr
      call gpu_upload_c(c_loc(cpu_ptr), gpu_ptr%c, product(shape(gpu_ptr%f)*1_8)*8_8)
    end subroutine


! gpu_download
! ------------

    subroutine gpu_download_double1(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double1), intent(in)  :: gpu_ptr
      double precision, target, intent(in)   :: cpu_ptr(:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*size(gpu_ptr%f))
    end subroutine

    subroutine gpu_download_double2(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double2), intent(in)  :: gpu_ptr
      double precision, target, intent(in)   :: cpu_ptr(:,:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*product(shape(gpu_ptr%f)*1_8))
    end subroutine

    subroutine gpu_download_double3(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double3), intent(in)  :: gpu_ptr
      double precision, target, intent(in)   :: cpu_ptr(:,:,:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*product(shape(gpu_ptr%f)*1_8))
    end subroutine

    subroutine gpu_download_double4(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double4), intent(in)  :: gpu_ptr
      double precision, target, intent(in)   :: cpu_ptr(:,:,:,:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*product(shape(gpu_ptr%f)*1_8))
    end subroutine

    subroutine gpu_download_double5(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double5), intent(in)  :: gpu_ptr
      double precision, target, intent(in)   :: cpu_ptr(:,:,:,:,:)
      call gpu_download_c(gpu_ptr%c, c_loc(cpu_ptr), 8_8*product(shape(gpu_ptr%f)*1_8))
    end subroutine

    subroutine gpu_download_double6(gpu_ptr, cpu_ptr)
      implicit none
      type(gpu_double6), intent(in)  :: gpu_ptr
      double precision, target, intent(in)   :: cpu_ptr(:,:,:,:,:,:)
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
      type(gpu_stream) :: stream
      call gpu_stream_create_c(stream%c)
    end subroutine

    subroutine gpu_stream_destroy(stream)
      type(gpu_stream) :: stream
      call gpu_stream_destroy_c(stream%c)
    end subroutine

    subroutine gpu_set_stream(handle, stream)
      type(gpu_blas)   :: handle
      type(gpu_stream) :: stream
      call gpu_set_stream_c(handle%c, stream%c)
    end subroutine


! gpu_blas
! --------

    subroutine gpu_blas_create(handle)
      type(gpu_blas) :: handle
      call gpu_blas_create_c(handle%c)
    end subroutine

    subroutine gpu_blas_destroy(handle)
      type(gpu_blas) :: handle
      call gpu_blas_destroy_c(handle%c)
    end subroutine


end module



! dot
! ---

subroutine gpu_ddot(handle, n, dx, incx, dy, incy, res)
  use gpu
  type(gpu_blas), intent(in)     :: handle
  integer*4                      :: n, incx, incy
  type(gpu_double1), intent(in)  :: dx, dy
  double precision, intent(out)  :: res
  call gpu_ddot_c(handle%c, int(n,c_int64_t), dx%c, int(incx,c_int64_t), dy%c, int(incy,c_int64_t), res)
end subroutine

subroutine gpu_ddot_f(handle, n, dx, incx, dy, incy, res)
  use gpu
  type(gpu_blas), intent(in)     :: handle
  integer*4                      :: n, incx, incy
  double precision, target       :: dx(*), dy(*)
  double precision, intent(out)  :: res
  call gpu_ddot_c(handle%c, int(n,c_int64_t), c_loc(dx), int(incx,c_int64_t), c_loc(dy), int(incy,c_int64_t), res)
end subroutine


subroutine gpu_ddot_64(handle, n, dx, incx, dy, incy, res)
  use gpu
  type(gpu_blas), intent(in)     :: handle
  integer*8                      :: n, incx, incy
  type(gpu_double1), intent(in)  :: dx, dy
  double precision, intent(out)  :: res
  call gpu_ddot_c(handle%c, n, dx%c, incx, dy%c, incy, res)
end subroutine

subroutine gpu_ddot_f_64(handle, n, dx, incx, dy, incy, res)
  use gpu
  type(gpu_blas), intent(in)     :: handle
  integer*8                      :: n, incx, incy
  double precision, target       :: dx(*), dy(*)
  double precision, intent(out)  :: res
  call gpu_ddot_c(handle%c, n, c_loc(dx), incx, c_loc(dy), incy, res)
end subroutine


! geam
! ----

subroutine gpu_dgeam(handle, transa, transb, m, n, alpha, a, lda, beta, &
  b, ldb, c, ldc)
  use gpu
  type(gpu_blas), intent(in)   :: handle
  character, intent(in)        :: transa, transb
  integer*4, intent(in)        :: m, n, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  type(gpu_double2)            :: a, b, c
  call gpu_dgeam_c(handle%c, transa, transb, int(m,c_int64_t), int(n,c_int64_t), alpha, a%c, int(lda,c_int64_t), beta, &
        b%c, int(ldb,c_int64_t), c%c, int(ldc,c_int64_t))
end subroutine


subroutine gpu_dgeam_f(handle, transa, transb, m, n, alpha, a, lda, beta, &
  b, ldb, c, ldc)
  use gpu
  type(gpu_blas), intent(in)   :: handle
  character, intent(in)        :: transa, transb
  integer*4, intent(in)        :: m, n, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  double precision, target     :: a(*), b(*), c(*)
  call gpu_dgeam_c(handle%c, transa, transb, int(m,c_int64_t), int(n,c_int64_t), alpha, c_loc(a), int(lda,c_int64_t), beta, &
        c_loc(b), int(ldb,c_int64_t), c_loc(c), int(ldc,c_int64_t))
end subroutine


subroutine gpu_dgeam_64(handle, transa, transb, m, n, alpha, a, lda, beta, &
  b, ldb, c, ldc)
  use gpu
  type(gpu_blas), intent(in)   :: handle
  character, intent(in)        :: transa, transb
  integer*8, intent(in)        :: m, n, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  type(gpu_double2)            :: a, b, c
  call gpu_dgeam_c(handle%c, transa, transb, int(m,c_int64_t), int(n,c_int64_t), alpha, a%c, int(lda,c_int64_t), beta, &
        b%c, int(ldb,c_int64_t), c%c, int(ldc,c_int64_t))
end subroutine


subroutine gpu_dgeam_f_64(handle, transa, transb, m, n, alpha, a, lda, beta, &
  b, ldb, c, ldc)
  use gpu
  type(gpu_blas), intent(in)   :: handle
  character, intent(in)        :: transa, transb
  integer*8, intent(in)        :: m, n, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  double precision, target     :: a(*), b(*), c(*)
  call gpu_dgeam_c(handle%c, transa, transb, int(m,c_int64_t), int(n,c_int64_t), alpha, c_loc(a), int(lda,c_int64_t), beta, &
        c_loc(b), int(ldb,c_int64_t), c_loc(c), int(ldc,c_int64_t))
end subroutine


! gemm
! ----

subroutine gpu_dgemm(handle, transa, transb, m, n, k, alpha, a, lda, &
  b, ldb, beta, c, ldc)
  use gpu
  type(gpu_blas), intent(in)   :: handle
  character, intent(in)        :: transa, transb
  integer*4, intent(in)        :: m, n, k, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  type(gpu_double2)            :: a, b, c
  call gpu_dgemm_c(handle%c, transa, transb, int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), &
        alpha, a%c, int(lda,c_int64_t), &
        b%c, int(ldb,c_int64_t), beta, c%c, int(ldc,c_int64_t))
end subroutine

subroutine gpu_dgemm_64(handle, transa, transb, m, n, k, alpha, a, lda, &
  b, ldb, beta, c, ldc)
  use gpu
  type(gpu_blas), intent(in)   :: handle
  character, intent(in)        :: transa, transb
  integer*8, intent(in)        :: m, n, k, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  type(gpu_double2)            :: a, b, c
  call gpu_dgemm_c(handle%c, transa, transb, m, n, k, &
        alpha, a%c, lda, b%c, ldb, beta, c%c, ldc)
end subroutine

subroutine gpu_dgemm_f(handle, transa, transb, m, n, k, alpha, a, lda, &
  b, ldb, beta, c, ldc)
  use gpu
  type(gpu_blas), intent(in)   :: handle
  character, intent(in)        :: transa, transb
  integer*4, intent(in)        :: m, n, k, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  double precision, target     :: a(*), b(*), c(*)
  call gpu_dgemm_c(handle%c, transa, transb, int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), &
        alpha, c_loc(a), int(lda,c_int64_t), &
        c_loc(b), int(ldb,c_int64_t), beta, c_loc(c), int(ldc,c_int64_t))
end subroutine

subroutine gpu_dgemm_f_64(handle, transa, transb, m, n, k, alpha, a, lda, &
  b, ldb, beta, c, ldc)
  use gpu
  type(gpu_blas), intent(in)   :: handle
  character, intent(in)        :: transa, transb
  integer*8, intent(in)        :: m, n, k, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  double precision, target     :: a(*), b(*), c(*)
  call gpu_dgemm_c(handle%c, transa, transb, m, n, k, &
        alpha, c_loc(a), lda, c_loc(b), ldb, beta, c_loc(c), ldc)
end subroutine

