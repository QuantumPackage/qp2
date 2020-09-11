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
  endtype

  contains

  integer function pt2_type_size(N)
    implicit none
    integer, intent(in) :: N
    if (is_complex) then
      pt2_type_size = (3*n + 2*n*n)
    else
      pt2_type_size = (3*n + n*n)
    endif

  end function

end module

