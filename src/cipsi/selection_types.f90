module selection_types
  type selection_buffer
    integer :: N, cur
    integer(8)      , pointer :: det(:,:,:)
    double precision, pointer :: val(:)
    double precision          :: mini
  endtype

  type pt2_type
    double precision, allocatable :: pt2(:)
    double precision, allocatable :: pt2_err(:)
    double precision, allocatable :: rpt2(:)
    double precision, allocatable :: rpt2_err(:)
    double precision, allocatable :: variance(:)
    double precision, allocatable :: variance_err(:)
    double precision, allocatable :: norm2(:)
    double precision, allocatable :: norm2_err(:)
    double precision, allocatable :: overlap(:,:)
    double precision, allocatable :: overlap_err(:,:)
  endtype


end module

