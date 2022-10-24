subroutine pt2_alloc(pt2_data,N)
  implicit none
  use selection_types
  type(pt2_type), intent(inout) :: pt2_data
  integer, intent(in)           :: N
  integer :: k

  allocate(pt2_data % pt2(N)           &
          ,pt2_data % variance(N)      &
          ,pt2_data % rpt2(N)          &
          ,pt2_data % overlap(N,N)     &
          )

  pt2_data % pt2(:)           = 0.d0
  pt2_data % variance(:)      = 0.d0
  pt2_data % rpt2(:)          = 0.d0
  pt2_data % overlap(:,:)     = 0.d0

end subroutine

subroutine pt2_dealloc(pt2_data)
  implicit none
  use selection_types
  type(pt2_type), intent(inout) :: pt2_data
  deallocate(pt2_data % pt2         &
            ,pt2_data % variance    &
            ,pt2_data % rpt2        &
            ,pt2_data % overlap     &
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
    p1 % rpt2(:)           = p1 % rpt2(:)          + p2 % rpt2(:)
    p1 % variance(:)       = p1 % variance(:)      + p2 % variance(:)
    p1 % overlap(:,:)      = p1 % overlap(:,:)     + p2 % overlap(:,:)

  else

    p1 % pt2(:)            = p1 % pt2(:)           + w * p2 % pt2(:)
    p1 % rpt2(:)           = p1 % rpt2(:)          + w * p2 % rpt2(:)
    p1 % variance(:)       = p1 % variance(:)      + w * p2 % variance(:)
    p1 % overlap(:,:)      = p1 % overlap(:,:)     + w * p2 % overlap(:,:)

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
    p1 % rpt2(:)          = p1 % rpt2(:)          + p2 % rpt2(:)          * p2 % rpt2(:)
    p1 % variance(:)      = p1 % variance(:)      + p2 % variance(:)      * p2 % variance(:)
    p1 % overlap(:,:)     = p1 % overlap(:,:)     + p2 % overlap(:,:)     * p2 % overlap(:,:)

  else

    p1 % pt2(:)           = p1 % pt2(:)           + w * p2 % pt2(:)           * p2 % pt2(:)
    p1 % rpt2(:)          = p1 % rpt2(:)          + w * p2 % rpt2(:)          * p2 % rpt2(:)
    p1 % variance(:)      = p1 % variance(:)      + w * p2 % variance(:)      * p2 % variance(:)
    p1 % overlap(:,:)     = p1 % overlap(:,:)     + w * p2 % overlap(:,:)     * p2 % overlap(:,:)

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
  x(1:n)           =  pt2_data % pt2(1:n)
  k=n
  x(k+1:k+n)     =  pt2_data % rpt2(1:n)
  k=k+n
  x(k+1:k+n)     =  pt2_data % variance(1:n)
  k=k+n
  x(k+1:k+n2)  =  reshape(pt2_data % overlap(1:n,1:n), (/ n2 /))

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
  k=n
  pt2_data % rpt2(1:n)          =   x(k+1:k+n)
  k=k+n
  pt2_data % variance(1:n)      =   x(k+1:k+n)
  k=k+n
  pt2_data % overlap(1:n,1:n) = reshape(x(k+1:k+n2), (/ n, n /))

end
