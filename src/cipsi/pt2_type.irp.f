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

subroutine pt2_add(p1, w, p2)
  implicit none
  use selection_types
  BEGIN_DOC
! p1 += w * p2
  END_DOC
  type(pt2_type), intent(inout) :: p1
  double precision, intent(in)  :: w
  type(pt2_type), intent(in)    :: p2

  if (w == 1.d0) then

    p1 % pt2(:)            = p1 % pt2(:)           + p2 % pt2(:)
    p1 % pt2_err(:)        = p1 % pt2_err(:)       + p2 % pt2_err(:)
    p1 % rpt2(:)           = p1 % rpt2(:)          + p2 % rpt2(:)
    p1 % rpt2_err(:)       = p1 % rpt2_err(:)      + p2 % rpt2_err(:)
    p1 % variance(:)       = p1 % variance(:)      + p2 % variance(:)
    p1 % variance_err(:)   = p1 % variance_err(:)  + p2 % variance_err(:)
    p1 % norm2(:)          = p1 % norm2(:)         + p2 % norm2(:)
    p1 % norm2_err(:)      = p1 % norm2_err(:)     + p2 % norm2_err(:)
    p1 % overlap(:,:)      = p1 % overlap(:,:)     + p2 % overlap(:,:)
    p1 % overlap_err(:,:)  = p1 % overlap_err(:,:) + p2 % overlap_err(:,:)

  else

    p1 % pt2(:)            = p1 % pt2(:)           + w * p2 % pt2(:)
    p1 % pt2_err(:)        = p1 % pt2_err(:)       + w * p2 % pt2_err(:)
    p1 % rpt2(:)           = p1 % rpt2(:)          + w * p2 % rpt2(:)
    p1 % rpt2_err(:)       = p1 % rpt2_err(:)      + w * p2 % rpt2_err(:)
    p1 % variance(:)       = p1 % variance(:)      + w * p2 % variance(:)
    p1 % variance_err(:)   = p1 % variance_err(:)  + w * p2 % variance_err(:)
    p1 % norm2(:)          = p1 % norm2(:)         + w * p2 % norm2(:)
    p1 % norm2_err(:)      = p1 % norm2_err(:)     + w * p2 % norm2_err(:)
    p1 % overlap(:,:)      = p1 % overlap(:,:)     + w * p2 % overlap(:,:)
    p1 % overlap_err(:,:)  = p1 % overlap_err(:,:) + w * p2 % overlap_err(:,:)

  endif

end subroutine


subroutine pt2_add2(p1, w, p2)
  implicit none
  use selection_types
  BEGIN_DOC
! p1 += w * p2**2
  END_DOC
  type(pt2_type), intent(inout) :: p1
  double precision, intent(in)  :: w
  type(pt2_type), intent(in)    :: p2

  if (w == 1.d0) then

    p1 % pt2(:)           = p1 % pt2(:)           + p2 % pt2(:)           * p2 % pt2(:)
    p1 % pt2_err(:)       = p1 % pt2_err(:)       + p2 % pt2_err(:)       * p2 % pt2_err(:)
    p1 % rpt2(:)          = p1 % rpt2(:)          + p2 % rpt2(:)          * p2 % rpt2(:)
    p1 % rpt2_err(:)      = p1 % rpt2_err(:)      + p2 % rpt2_err(:)      * p2 % rpt2_err(:)
    p1 % variance(:)      = p1 % variance(:)      + p2 % variance(:)      * p2 % variance(:)
    p1 % variance_err(:)  = p1 % variance_err(:)  + p2 % variance_err(:)  * p2 % variance_err(:)
    p1 % norm2(:)         = p1 % norm2(:)         + p2 % norm2(:)         * p2 % norm2(:)
    p1 % norm2_err(:)     = p1 % norm2_err(:)     + p2 % norm2_err(:)     * p2 % norm2_err(:)
    p1 % overlap(:,:)     = p1 % overlap(:,:)     + p2 % overlap(:,:)     * p2 % overlap(:,:)
    p1 % overlap_err(:,:) = p1 % overlap_err(:,:) + p2 % overlap_err(:,:) * p2 % overlap_err(:,:)

  else

    p1 % pt2(:)           = p1 % pt2(:)           + w * p2 % pt2(:)           * p2 % pt2(:)
    p1 % pt2_err(:)       = p1 % pt2_err(:)       + w * p2 % pt2_err(:)       * p2 % pt2_err(:)
    p1 % rpt2(:)          = p1 % rpt2(:)          + w * p2 % rpt2(:)          * p2 % rpt2(:)
    p1 % rpt2_err(:)      = p1 % rpt2_err(:)      + w * p2 % rpt2_err(:)      * p2 % rpt2_err(:)
    p1 % variance(:)      = p1 % variance(:)      + w * p2 % variance(:)      * p2 % variance(:)
    p1 % variance_err(:)  = p1 % variance_err(:)  + w * p2 % variance_err(:)  * p2 % variance_err(:)
    p1 % norm2(:)         = p1 % norm2(:)         + w * p2 % norm2(:)         * p2 % norm2(:)
    p1 % norm2_err(:)     = p1 % norm2_err(:)     + w * p2 % norm2_err(:)     * p2 % norm2_err(:)
    p1 % overlap(:,:)     = p1 % overlap(:,:)     + w * p2 % overlap(:,:)     * p2 % overlap(:,:)
    p1 % overlap_err(:,:) = p1 % overlap_err(:,:) + w * p2 % overlap_err(:,:) * p2 % overlap_err(:,:)

  endif

end subroutine


subroutine pt2_serialize(pt2_data, n, x)
  implicit none
  use selection_types
  type(pt2_type), intent(in)    :: pt2_data
  integer, intent(in)           :: n
  double precision, intent(out) :: x(*)

  integer :: i,k,n2

  n2 = n*n
  x(1:n)        =  pt2_data % pt2(1:n)
  x(n+1:2*n)    =  pt2_data % pt2_err(1:n)
  x(2*n+1:3*n)  =  pt2_data % rpt2(1:n)
  x(3*n+1:4*n)  =  pt2_data % rpt2_err(1:n)
  x(4*n+1:5*n)  =  pt2_data % variance(1:n)
  x(5*n+1:6*n)  =  pt2_data % variance_err(1:n)
  x(6*n+1:7*n)  =  pt2_data % norm2(1:n)
  x(7*n+1:8*n)  =  pt2_data % norm2_err(1:n)
  k=8*n
  x(k+1:k+n2)   =  reshape(pt2_data % overlap(1:n,1:n), (/ n2 /))
  k=8*n+n2
  x(k+1:k+n2)   =  reshape(pt2_data % overlap_err(1:n,1:n), (/ n2 /))

end

subroutine pt2_deserialize(pt2_data, n, x)
  implicit none
  use selection_types
  type(pt2_type), intent(inout) :: pt2_data
  integer, intent(in)           :: n
  double precision, intent(in)  :: x(*)

  integer :: i,k,n2

  n2 = n*n
  pt2_data % pt2(1:n)           =   x(1:n)
  pt2_data % pt2_err(1:n)       =   x(n+1:2*n)
  pt2_data % rpt2(1:n)          =   x(2*n+1:3*n)
  pt2_data % rpt2_err(1:n)      =   x(3*n+1:4*n)
  pt2_data % variance(1:n)      =   x(4*n+1:5*n)
  pt2_data % variance_err(1:n)  =   x(5*n+1:6*n)
  pt2_data % norm2(1:n)         =   x(6*n+1:7*n)
  pt2_data % norm2_err(1:n)     =   x(7*n+1:8*n)
  k=8*n
  pt2_data % overlap(1:n,1:n) = reshape(x(k+1:k+n2), (/ n, n /))
  k=8*n+n2
  pt2_data % overlap_err(1:n,1:n) = reshape(x(k+1:k+n2), (/ n, n /))

end
