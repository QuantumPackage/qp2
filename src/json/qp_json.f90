!> Plain Fortran interface of the JSON result document (qp3).
!!
!! Call sites write one line:
!!     call json_add('/results/scf/energy', SCF_energy)
!! and, once per program that produces results,
!!     call json_init('scf')            ! src/json/json.irp.f
!!
!! No JSON state is visible here: every specific converts its argument
!! to the interoperable type of the C shim (`src/json/qp_json_shim.c`)
!! and the shim owns the document, the file and the throttle.  An
!! unsupported actual argument (single-precision real, integer(8), a
!! rank-2 array) is a compile error by design: no silent narrowing and
!! no silent flattening.
!!
!! Array specifics take assumed-shape, possibly strided dummies because
!! call sites pass sections such as extrapolated_energy(N_iter, 1:n);
!! they copy into a contiguous buffer before the bind(C) call.
module qp_json
  use iso_c_binding
  implicit none
  private

  public :: json_add
  public :: qp_json_begin

  interface
    subroutine qp_json_begin_c(kind, filename, ezfio_file, master) &
        bind(C, name='qp_json_begin')
      use iso_c_binding
      implicit none
      character(kind=c_char), intent(in) :: kind(*)
      character(kind=c_char), intent(in) :: filename(*)
      character(kind=c_char), intent(in) :: ezfio_file(*)
      integer(c_int), value, intent(in) :: master
    end subroutine qp_json_begin_c

    subroutine qp_json_add_int_c(path, v) bind(C, name='qp_json_add_int')
      use iso_c_binding
      implicit none
      character(kind=c_char), intent(in) :: path(*)
      integer(c_long_long), value, intent(in) :: v
    end subroutine qp_json_add_int_c

    subroutine qp_json_add_real_c(path, v) bind(C, name='qp_json_add_real')
      use iso_c_binding
      implicit none
      character(kind=c_char), intent(in) :: path(*)
      real(c_double), value, intent(in) :: v
    end subroutine qp_json_add_real_c

    subroutine qp_json_add_bool_c(path, v) bind(C, name='qp_json_add_bool')
      use iso_c_binding
      implicit none
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value, intent(in) :: v
    end subroutine qp_json_add_bool_c

    subroutine qp_json_add_str_c(path, v) bind(C, name='qp_json_add_str')
      use iso_c_binding
      implicit none
      character(kind=c_char), intent(in) :: path(*)
      character(kind=c_char), intent(in) :: v(*)
    end subroutine qp_json_add_str_c

    subroutine qp_json_add_arr_int_c(path, v, n) &
        bind(C, name='qp_json_add_arr_int')
      use iso_c_binding
      implicit none
      character(kind=c_char), intent(in) :: path(*)
      integer(c_long_long), intent(in) :: v(*)
      integer(c_int), value, intent(in) :: n
    end subroutine qp_json_add_arr_int_c

    subroutine qp_json_add_arr_real_c(path, v, n) &
        bind(C, name='qp_json_add_arr_real')
      use iso_c_binding
      implicit none
      character(kind=c_char), intent(in) :: path(*)
      real(c_double), intent(in) :: v(*)
      integer(c_int), value, intent(in) :: n
    end subroutine qp_json_add_arr_real_c
  end interface

  interface json_add
    module procedure json_add_int
    module procedure json_add_int64
    module procedure json_add_real
    module procedure json_add_bool
    module procedure json_add_str
    module procedure json_add_arr_int
    module procedure json_add_arr_int64
    module procedure json_add_arr_real
  end interface

contains

  !> Start the JSON logging of this run and write the first version of
  !! the file.  Called once per program, before any json_add.
  !! `master` is 1 on the MPI master, 0 elsewhere; non-master ranks do
  !! no work and no I/O afterwards.
  subroutine qp_json_begin(kind, filename, ezfio_file, master)
    character(len=*), intent(in) :: kind
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: ezfio_file
    integer(c_int), intent(in) :: master

    call qp_json_begin_c(trim(kind)//c_null_char, &
                         trim(filename)//c_null_char, &
                         trim(ezfio_file)//c_null_char, master)
  end subroutine qp_json_begin

  subroutine json_add_int(path, value)
    character(len=*), intent(in) :: path
    integer, intent(in) :: value

    call qp_json_add_int_c(trim(path)//c_null_char, int(value, c_long_long))
  end subroutine json_add_int

  subroutine json_add_int64(path, value)
    character(len=*), intent(in) :: path
    integer(c_int64_t), intent(in) :: value

    call qp_json_add_int_c(trim(path)//c_null_char, int(value, c_long_long))
  end subroutine json_add_int64

  subroutine json_add_real(path, value)
    character(len=*), intent(in) :: path
    real(c_double), intent(in) :: value

    call qp_json_add_real_c(trim(path)//c_null_char, value)
  end subroutine json_add_real

  subroutine json_add_bool(path, value)
    character(len=*), intent(in) :: path
    logical, intent(in) :: value

    call qp_json_add_bool_c(trim(path)//c_null_char, &
                            merge(1_c_int, 0_c_int, value))
  end subroutine json_add_bool

  subroutine json_add_str(path, value)
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: value

    call qp_json_add_str_c(trim(path)//c_null_char, &
                           trim(value)//c_null_char)
  end subroutine json_add_str

  subroutine json_add_arr_int(path, values)
    use iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: path
    integer, intent(in), dimension(:) :: values
    integer(c_long_long), allocatable :: buf(:)
    integer :: i, info

    allocate(buf(size(values)), stat=info)
    if (info /= 0) then
      write(error_unit, '(A)') &
        'qp_json: out of memory; an integer array was not logged'
      return
    end if
    do i = 1, size(values)
      buf(i) = int(values(i), c_long_long)
    end do
    call qp_json_add_arr_int_c(trim(path)//c_null_char, buf, &
                               int(size(buf), c_int))
    deallocate(buf)
  end subroutine json_add_arr_int

  subroutine json_add_arr_int64(path, values)
    character(len=*), intent(in) :: path
    integer(c_int64_t), intent(in), dimension(:) :: values

    call qp_json_add_arr_int_c(trim(path)//c_null_char, values, &
                               int(size(values), c_int))
  end subroutine json_add_arr_int64

  subroutine json_add_arr_real(path, values)
    use iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: path
    real(c_double), intent(in), dimension(:) :: values
    real(c_double), allocatable :: buf(:)
    integer :: info

    allocate(buf(size(values)), stat=info)
    if (info /= 0) then
      write(error_unit, '(A)') &
        'qp_json: out of memory; a real array was not logged'
      return
    end if
    buf = values
    call qp_json_add_arr_real_c(trim(path)//c_null_char, buf, &
                                int(size(buf), c_int))
    deallocate(buf)
  end subroutine json_add_arr_real

end module qp_json
