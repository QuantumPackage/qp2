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
    ! Per-call logging counters (JSON history); zeroed by pt2_alloc.
    integer          :: windows      ! windows walked, 0 if not run
    integer          :: tasks        ! (generator, subset) tasks, 0 if not run
    double precision :: samples      ! estimator draws, summed over states
    double precision :: time         ! seconds spent in this PT2 call
    double precision :: epoch_start  ! qp_epoch_c() at entry
    double precision :: epoch_end    ! qp_epoch_c() at exit
  endtype

  contains

  integer function pt2_type_size(N)
    implicit none
    integer, intent(in) :: N
    pt2_type_size = (3*n + n*n)
  end function

end module

