BEGIN_TEMPLATE

 subroutine sorted_$Xnumber(x,isize,n)
  implicit none
  BEGIN_DOC
! Returns the number of sorted elements
  END_DOC
  integer, intent(in)            :: isize
  $type, intent(in)           :: x(isize)
  integer, intent(out)           :: n
  integer :: i
  n=1

  if (isize < 2) then
    return
  endif

  do i=2,isize
    if (x(i-1) <= x(i)) then
      n=n+1
    endif
  enddo

 end

SUBST [ X, type ]
   ; real ;;
 d ; double precision ;;
 i ; integer ;;
 i8 ; integer*8 ;;
 i2 ; integer*2 ;;
END_TEMPLATE



BEGIN_TEMPLATE
 subroutine $Xset_order(x,iorder,isize)
  implicit none
  BEGIN_DOC
  ! array A has already been sorted, and iorder has contains the new order of
  ! elements of A. This subroutine changes the order of x to match the new order of A.
  END_DOC
  integer                        :: isize
  $type                          :: x(*)
  $type,allocatable              :: xtmp(:)
  integer                        :: iorder(*)
  integer                        :: i

  allocate(xtmp(isize))
  do i=1,isize
    xtmp(i) = x(iorder(i))
  enddo

  do i=1,isize
    x(i) = xtmp(i)
  enddo
  deallocate(xtmp)
 end

 subroutine $Xset_order_big(x,iorder,isize)
  implicit none
  BEGIN_DOC
  ! array A has already been sorted, and iorder has contains the new order of
  ! elements of A. This subroutine changes the order of x to match the new order of A.
  ! This is a version for very large arrays where the indices need
  ! to be in integer*8 format
  END_DOC
  integer*8                      :: isize
  $type                          :: x(*)
  $type, allocatable             :: xtmp(:)
  integer*8                      :: iorder(*)
  integer*8                      :: i
  allocate(xtmp(isize))
  do i=1_8,isize
    xtmp(i) = x(iorder(i))
  enddo

  do i=1_8,isize
    x(i) = xtmp(i)
  enddo
  deallocate(xtmp)
 end

SUBST [ X, type ]
   ; real ;;
 d ; double precision ;;
 i ; integer ;;
 i8; integer*8 ;;
 i2; integer*2 ;;
END_TEMPLATE


