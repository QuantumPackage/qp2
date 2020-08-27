module selection_types
  type selection_buffer
    integer :: N, cur
    integer(8)      , pointer :: det(:,:,:)
    double precision, pointer :: val(:)
    double precision          :: mini
  endtype

  type pt2_type
    double precision, allocatable :: pt2(:)
    double precision, allocatable :: variance(:)
    double precision, allocatable :: norm2(:)
  endtype
end module

