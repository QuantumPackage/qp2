subroutine pt2_alloc(pt2_data,N)
  implicit none
  use selection_types
  type(pt2_type), intent(inout) :: pt2_data
  integer, intent(in)           :: N
  integer :: k

  allocate(pt2_data % pt2(N)           &
          ,pt2_data % pt2_err(N)       &
          ,pt2_data % variance(N)      &
          ,pt2_data % variance_err(N)  &
          ,pt2_data % norm2(N)         &
          ,pt2_data % norm2_err(N)     &
          ,pt2_data % rpt2(N)          &
          ,pt2_data % rpt2_err(N)      &
          ,pt2_data % overlap(N,N)     &
          ,pt2_data % overlap_err(N,N) &
          )

  pt2_data % pt2(:)           = 0.d0
  pt2_data % pt2_err(:)       = 0.d0
  pt2_data % variance(:)      = 0.d0
  pt2_data % variance_err(:)  = 0.d0
  pt2_data % norm2(:)         = 0.d0
  pt2_data % norm2_err(:)     = 0.d0
  pt2_data % rpt2(:)          = 0.d0
  pt2_data % rpt2_err(:)      = 0.d0
  pt2_data % overlap(:,:)     = 0.d0
  pt2_data % overlap_err(:,:) = 0.d0

  do k=1,N
    pt2_data % overlap(k,k) = 1.d0
  enddo
end subroutine

subroutine pt2_dealloc(pt2_data)
  implicit none
  use selection_types
  type(pt2_type), intent(inout) :: pt2_data
  deallocate(pt2_data % pt2         &
            ,pt2_data % pt2_err     &
            ,pt2_data % variance    &
            ,pt2_data % variance_err&
            ,pt2_data % norm2       &
            ,pt2_data % norm2_err   &
            ,pt2_data % rpt2        &
            ,pt2_data % rpt2_err    &
            ,pt2_data % overlap     &
            ,pt2_data % overlap_err &
            )
end subroutine

subroutine pt2_add(p1, p2)
  implicit none
  use selection_types
  BEGIN_DOC
! p1 += p2
  END_DOC
  type(pt2_type), intent(inout) :: p1
  type(pt2_type), intent(in)    :: p2

  p1 % pt2(:)            += p2 % pt2(:)
  p1 % pt2_err(:)        += p2 % pt2_err(:)
  p1 % rpt2(:)           += p2 % rpt2(:)
  p1 % rpt2_err(:)       += p2 % rpt2_err(:)
  p1 % variance(:)       += p2 % variance(:)
  p1 % variance_err(:)   += p2 % variance_err(:)
  p1 % norm2(:)          += p2 % norm2(:)
  p1 % norm2_err(:)      += p2 % norm2_err(:)
  p1 % overlap(:,:)      += p2 % overlap(:,:)
  p1 % overlap_err(:,:)  += p2 % overlap_err(:,:)
end subroutine


