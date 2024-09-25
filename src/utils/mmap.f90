module mmap_module

  use iso_c_binding

  type mmap_type
    type(c_ptr)     :: ptr          ! Pointer to the data
    character*(128) :: filename     ! Name of the file
    integer*8       :: length       ! Size of the array in bytes
    integer         :: fd           ! File descriptor

    ! Pointers to data
    integer, pointer, dimension (:)       :: i1
    integer, pointer, dimension (:,:)     :: i2
    integer, pointer, dimension (:,:,:)   :: i3
    integer, pointer, dimension (:,:,:,:) :: i4

    ! Pointers to data
    double precision, pointer, dimension (:)       :: d1
    double precision, pointer, dimension (:,:)     :: d2
    double precision, pointer, dimension (:,:,:)   :: d3
    double precision, pointer, dimension (:,:,:,:) :: d4
  end type mmap_type

  interface

    ! File descriptors
    ! ----------------

    type(c_ptr) function c_mmap_fortran(filename, length, fd, read_only, single_node) bind(c,name='mmap_fortran')
      use iso_c_binding
      character(c_char), intent(in)  :: filename(*)
      integer(c_size_t), intent(in), value :: length
      integer(c_int), intent(out)    :: fd
      integer(c_int), intent(in), value    :: read_only
      integer(c_int), intent(in), value    :: single_node
    end function

    subroutine c_munmap_fortran(length, fd, map) bind(c,name='munmap_fortran')
      use iso_c_binding
      integer(c_size_t), intent(in), value :: length
      integer(c_int), intent(in), value :: fd
      type(c_ptr), intent(in), value    :: map
    end subroutine

    subroutine c_msync_fortran(length, fd, map) bind(c,name='msync_fortran')
      use iso_c_binding
      integer(c_size_t), intent(in), value :: length
      integer(c_int), intent(in), value :: fd
      type(c_ptr), intent(in), value    :: map
    end subroutine

  end interface

  contains

  subroutine mmap(filename, shape, bytes, fd, read_only, single_node, map)
      use iso_c_binding
      implicit none
      character*(*), intent(in)      :: filename   ! Name of the mapped file
      integer*8, intent(in)          :: shape(:)   ! Shape of the array to map
      integer, intent(in)            :: bytes      ! Number of bytes per element
      logical, intent(in)            :: read_only  ! If true, mmap is read-only
      logical, intent(in)            :: single_node! If true, mmap is on a single node
      integer, intent(out)           :: fd         ! File descriptor
      type(c_ptr), intent(out)       :: map        ! C Pointer

      integer(c_size_t)              :: length
      integer(c_int)                 :: fd_, read_only_, single_node_

      integer :: i

      read_only_ = 0
      single_node_ = 0
      if (read_only) read_only_ = 1
      if (single_node) single_node_ = 1

      length = int(bytes,8)
      do i=1,size(shape)
        length = length * shape(i)
      enddo

      map = c_mmap_fortran( trim(filename)//char(0), length, fd_, read_only_, single_node_)
      fd = fd_
  end subroutine

  subroutine munmap(shape, bytes, fd, map)
      use iso_c_binding
      implicit none
      integer*8, intent(in)          :: shape(:)  ! Shape of the array to map
      integer, intent(in)            :: bytes     ! Number of bytes per element
      integer, intent(in)            :: fd        ! File descriptor
      type(c_ptr), intent(in)        :: map       ! C pointer

      integer(c_size_t)              :: length
      integer(c_int)                 :: fd_

      integer :: i

      length = int(bytes,8)
      do i=1,size(shape)
        length = length * shape(i)
      enddo
      fd_ = fd
      call c_munmap_fortran(length, fd_, map)
  end subroutine

  subroutine msync(shape, bytes, fd, map)
      use iso_c_binding
      implicit none
      integer*8, intent(in)          :: shape(:)  ! Shape of the array to map
      integer, intent(in)            :: bytes     ! Number of bytes per element
      integer, intent(in)            :: fd        ! File descriptor
      type(c_ptr), intent(in)        :: map       ! C pointer

      integer(c_size_t)              :: length
      integer(c_int)                 :: fd_

      integer :: i

      length = int(bytes,8)
      do i=1,size(shape)
        length = length * shape(i)
      enddo
      fd_ = fd
      call c_msync_fortran( length, fd_, map)
  end subroutine


  ! Functions for the mmap_type

  subroutine mmap_create(filename, shape, bytes, read_only, single_node, map)
      implicit none
      character*(*), intent(in)      :: filename   ! Name of the mapped file
      integer*8, intent(in)          :: shape(:)   ! Shape of the array to map
      integer, intent(in)            :: bytes      ! Number of bytes per element
      logical, intent(in)            :: read_only  ! If true, mmap is read-only
      logical, intent(in)            :: single_node! If true, mmap is on a single node
      type(mmap_type), intent(out)   :: map        ! mmap

      integer :: i

      map%filename = filename
      map%length = int(bytes,8)
      do i=1,size(shape)
        map%length = map%length * shape(i)
      enddo
      call mmap(map%filename, &
                shape,        &
                bytes,        &
                map%fd,       &
                read_only,    &
                single_node,  &
                map%ptr)

      map%d1 => NULL()
      map%d2 => NULL()
      map%d3 => NULL()
      map%d4 => NULL()
      map%i1 => NULL()
      map%i2 => NULL()
      map%i3 => NULL()
      map%i4 => NULL()

  end

  subroutine mmap_create_d(filename, shape, bytes, read_only, single_node, map)
      implicit none
      character*(*), intent(in)      :: filename   ! Name of the mapped file
      integer*8, intent(in)          :: shape(:)   ! Shape of the array to map
      integer, intent(in)            :: bytes      ! Number of bytes per element
      logical, intent(in)            :: read_only  ! If true, mmap is read-only
      logical, intent(in)            :: single_node! If true, mmap is on a single node
      type(mmap_type), intent(out)   :: map        ! mmap

      call mmap_create(filename, shape, bytes, read_only, single_node, map)

      select case (size(shape))
        case (1)
          call c_f_pointer(map%ptr, map%d1, shape)
        case (2)
          call c_f_pointer(map%ptr, map%d2, shape)
        case (3)
          call c_f_pointer(map%ptr, map%d3, shape)
        case (4)
          call c_f_pointer(map%ptr, map%d4, shape)
        case default
          stop 'mmap: dimension not implemented'
      end select
  end subroutine

  subroutine mmap_create_i(filename, shape, bytes, read_only, single_node, map)
      implicit none
      character*(*), intent(in)      :: filename   ! Name of the mapped file
      integer*8, intent(in)          :: shape(:)   ! Shape of the array to map
      integer, intent(in)            :: bytes      ! Number of bytes per element
      logical, intent(in)            :: read_only  ! If true, mmap is read-only
      logical, intent(in)            :: single_node! If true, mmap is on a single node
      type(mmap_type), intent(out)   :: map        ! mmap

      call mmap_create(filename, shape, bytes, read_only, single_node, map)

      select case (size(shape))
        case (1)
          call c_f_pointer(map%ptr, map%i1, shape)
        case (2)
          call c_f_pointer(map%ptr, map%i2, shape)
        case (3)
          call c_f_pointer(map%ptr, map%i3, shape)
        case (4)
          call c_f_pointer(map%ptr, map%i4, shape)
        case default
          stop 'mmap: dimension not implemented'
      end select
  end subroutine

  subroutine mmap_destroy(map)
      implicit none
      type(mmap_type), intent(inout)    :: map

      call c_munmap_fortran(map%length, map%fd, map%ptr)

      map%ptr      = C_NULL_PTR
      map%filename = ''
      map%length   = 0
      map%fd       = 0
      map%d1 => NULL()
      map%d2 => NULL()
      map%d3 => NULL()
      map%d4 => NULL()
      map%i1 => NULL()
      map%i2 => NULL()
      map%i3 => NULL()
      map%i4 => NULL()
  end subroutine


  subroutine mmap_sync(map)
      implicit none
      type(mmap_type), intent(inout)    :: map

      call c_msync_fortran(map%length, map%fd, map%ptr)
  end subroutine

end module mmap_module


