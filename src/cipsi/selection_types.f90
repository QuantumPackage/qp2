module selection_types
  type selection_buffer
    integer :: N, cur
    integer(8)      , pointer :: det(:,:,:)
    double precision, pointer :: val(:)
    double precision          :: mini
  endtype

  type pt2_type
    double precision, allocatable :: pt2(:)
    double precision, allocatable :: rpt2(:)
    double precision, allocatable :: variance(:)
    double precision, allocatable :: overlap(:,:)
    double precision, allocatable :: overlap_imag(:,:)
  endtype

  contains

  integer function pt2_type_size(N,has_imag)
    implicit none
    integer, intent(in) :: N
    logical, intent(in), optional :: has_imag
    logical :: has_imag_tmp
    if(present(has_imag)) then
      has_imag_tmp = has_imag
    else
      has_imag_tmp = .False.
    endif

    if (has_imag_tmp) then
      pt2_type_size = (3*n + 2*n*n)
    else
      pt2_type_size = (3*n + n*n)
    endif

  end function

end module

