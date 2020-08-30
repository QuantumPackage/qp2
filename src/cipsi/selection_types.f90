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
    double precision, allocatable :: norm2(:)
    double precision, allocatable :: overlap(:,:)
  endtype

  contains

  integer function pt2_type_size(N)
    implicit none
    integer, intent(in) :: N
    pt2_type_size = (4*n + n*n)
  end function

end module

