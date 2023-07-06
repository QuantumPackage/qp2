module mmap_module

  use iso_c_binding

  interface

    ! File descriptors
    ! ----------------

    type(c_ptr) function c_mmap_fortran(filename, length, fd, read_only) bind(c,name='mmap_fortran')
      use iso_c_binding
      character(c_char), intent(in)  :: filename(*)
      integer(c_size_t), intent(in), value :: length
      integer(c_int), intent(out)    :: fd
      integer(c_int), intent(in), value    :: read_only
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

  subroutine mmap(filename, shape, bytes, fd, read_only, map)
      use iso_c_binding
      implicit none
      character*(*), intent(in)      :: filename   ! Name of the mapped file
      integer*8, intent(in)          :: shape(:)   ! Shape of the array to map
      integer, intent(in)            :: bytes      ! Number of bytes per element
      logical, intent(in)            :: read_only  ! If true, mmap is read-only
      integer, intent(out)           :: fd         ! File descriptor
      type(c_ptr), intent(out)       :: map        ! C Pointer

      integer(c_size_t)              :: length
      integer(c_int)                 :: fd_

      integer :: i

      length = int(bytes,8)
      do i=1,size(shape)
        length = length * shape(i)
      enddo

      if (read_only) then
          map = c_mmap_fortran( trim(filename)//char(0), length, fd_, 1)
      else
          map = c_mmap_fortran( trim(filename)//char(0), length, fd_, 0)
      endif
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
      call c_munmap_fortran( length, fd_, map)
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

end module mmap_module


