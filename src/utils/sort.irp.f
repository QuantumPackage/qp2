BEGIN_TEMPLATE
 subroutine insertion_$Xsort (x,iorder,isize)
  implicit none
  BEGIN_DOC
  ! Sort array x(isize) using the insertion sort algorithm.
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  END_DOC
  integer,intent(in)             :: isize
  $type,intent(inout)            :: x(isize)
  integer,intent(inout)          :: iorder(isize)
  $type                          :: xtmp
  integer                        :: i, i0, j, jmax

  do i=2,isize
    xtmp = x(i)
    i0 = iorder(i)
    j=i-1
    do while (j>0)
      if ((x(j) <= xtmp)) exit
      x(j+1) = x(j)
      iorder(j+1) = iorder(j)
      j=j-1
    enddo
    x(j+1) = xtmp
    iorder(j+1) = i0
  enddo
 end subroutine insertion_$Xsort

 subroutine quick_$Xsort(x, iorder, isize)
  implicit none
  BEGIN_DOC
  ! Sort array x(isize) using the quicksort algorithm.
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  END_DOC
  integer,intent(in)             :: isize
  $type,intent(inout)            :: x(isize)
  integer,intent(inout)          :: iorder(isize)
  integer, external              :: omp_get_num_threads
  call rec_$X_quicksort(x,iorder,isize,1,isize,nproc)
 end

 recursive subroutine rec_$X_quicksort(x, iorder, isize, first, last, level)
  implicit none
  integer, intent(in)            :: isize, first, last, level
  integer,intent(inout)          :: iorder(isize)
  $type, intent(inout)           :: x(isize)
  $type                          :: c, tmp
  integer                        :: itmp
  integer                        :: i, j

  if(isize<2)return

  c = x( shiftr(first+last,1) )
  i = first
  j = last
  do
    do while (x(i) < c)
      i=i+1
    end do
    do while (c < x(j))
      j=j-1
    end do
    if (i >= j) exit
    tmp  = x(i)
    x(i) = x(j)
    x(j) = tmp
    itmp      = iorder(i)
    iorder(i) = iorder(j)
    iorder(j) = itmp
    i=i+1
    j=j-1
  enddo
  if ( ((i-first <= 10000).and.(last-j <= 10000)).or.(level<=0) ) then
    if (first < i-1) then
      call rec_$X_quicksort(x, iorder, isize, first, i-1,level/2)
    endif
    if (j+1 < last) then
      call rec_$X_quicksort(x, iorder, isize, j+1, last,level/2)
    endif
  else
    if (first < i-1) then
      call rec_$X_quicksort(x, iorder, isize, first, i-1,level/2)
    endif
    if (j+1 < last) then
      call rec_$X_quicksort(x, iorder, isize, j+1, last,level/2)
    endif
  endif
 end

 subroutine heap_$Xsort(x,iorder,isize)
  implicit none
  BEGIN_DOC
  ! Sort array x(isize) using the heap sort algorithm.
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  END_DOC
  integer,intent(in)             :: isize
  $type,intent(inout)            :: x(isize)
  integer,intent(inout)          :: iorder(isize)

  integer                        :: i, k, j, l, i0
  $type                          :: xtemp

  l = isize/2+1
  k = isize
  do while (.True.)
    if (l>1) then
      l=l-1
      xtemp = x(l)
      i0 = iorder(l)
    else
      xtemp = x(k)
      i0 = iorder(k)
      x(k) = x(1)
      iorder(k) = iorder(1)
      k = k-1
      if (k == 1) then
        x(1) = xtemp
        iorder(1) = i0
        exit
      endif
    endif
    i=l
    j = shiftl(l,1)
    do while (j<k)
      if ( x(j) < x(j+1) ) then
        j=j+1
      endif
      if (xtemp < x(j)) then
        x(i) = x(j)
        iorder(i) = iorder(j)
        i = j
        j = shiftl(j,1)
      else
        j = k+1
      endif
    enddo
    if (j==k) then
      if (xtemp < x(j)) then
        x(i) = x(j)
        iorder(i) = iorder(j)
        i = j
        j = shiftl(j,1)
      else
        j = k+1
      endif
    endif
    x(i) = xtemp
    iorder(i) = i0
  enddo
 end subroutine heap_$Xsort

 subroutine heap_$Xsort_big(x,iorder,isize)
  implicit none
  BEGIN_DOC
  ! Sort array x(isize) using the heap sort algorithm.
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  ! This is a version for very large arrays where the indices need
  ! to be in integer*8 format
  END_DOC
  integer*8,intent(in)           :: isize
  $type,intent(inout)            :: x(isize)
  integer*8,intent(inout)        :: iorder(isize)

  integer*8                      :: i, k, j, l, i0
  $type                          :: xtemp

  l = isize/2+1
  k = isize
  do while (.True.)
    if (l>1) then
      l=l-1
      xtemp = x(l)
      i0 = iorder(l)
    else
      xtemp = x(k)
      i0 = iorder(k)
      x(k) = x(1)
      iorder(k) = iorder(1)
      k = k-1
      if (k == 1) then
        x(1) = xtemp
        iorder(1) = i0
        exit
      endif
    endif
    i=l
    j = shiftl(l,1)
    do while (j<k)
      if ( x(j) < x(j+1) ) then
        j=j+1
      endif
      if (xtemp < x(j)) then
        x(i) = x(j)
        iorder(i) = iorder(j)
        i = j
        j = shiftl(j,1)
      else
        j = k+1
      endif
    enddo
    if (j==k) then
      if (xtemp < x(j)) then
        x(i) = x(j)
        iorder(i) = iorder(j)
        i = j
        j = shiftl(j,1)
      else
        j = k+1
      endif
    endif
    x(i) = xtemp
    iorder(i) = i0
  enddo

 end subroutine heap_$Xsort_big

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


!---------------------- INTEL
IRP_IF INTEL

BEGIN_TEMPLATE
 subroutine $Xsort(x,iorder,isize)
  use intel
  implicit none
  BEGIN_DOC
  ! Sort array x(isize).
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  END_DOC
  integer,intent(in)             :: isize
  $type,intent(inout)            :: x(isize)
  integer,intent(inout)          :: iorder(isize)
  integer                        :: n
  character, allocatable         :: tmp(:)
  if (isize < 2) return
  call ippsSortRadixIndexGetBufferSize(isize, $ippsz, n)
  allocate(tmp(n))
  call ippsSortRadixIndexAscend_$ityp(x, $n, iorder, isize, tmp)
  deallocate(tmp)
  iorder(1:isize) = iorder(1:isize)+1
  call $Xset_order(x,iorder,isize)
 end

 subroutine $Xsort_noidx(x,isize)
  use intel
  implicit none
  BEGIN_DOC
  ! Sort array x(isize).
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  END_DOC
  integer,intent(in)             :: isize
  $type,intent(inout)            :: x(isize)
  integer                        :: n
  character, allocatable         :: tmp(:)
  if (isize < 2) return
  call ippsSortRadixIndexGetBufferSize(isize, $ippsz, n)
  allocate(tmp(n))
  call ippsSortRadixAscend_$ityp_I(x, isize, tmp)
  deallocate(tmp)
 end

SUBST [ X, type, ityp, n, ippsz ]
   ; real ; 32f ; 4 ; 13 ;;
 i ; integer ; 32s ; 4 ; 11 ;;
 i2 ; integer*2 ; 16s ; 2 ; 7 ;;
END_TEMPLATE

BEGIN_TEMPLATE

 subroutine $Xsort(x,iorder,isize)
  implicit none
  BEGIN_DOC
  ! Sort array x(isize).
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  END_DOC
  integer,intent(in)             :: isize
  $type,intent(inout)            :: x(isize)
  integer,intent(inout)          :: iorder(isize)
  integer                        :: n
  if (isize < 2) then
    return
  endif
!  call sorted_$Xnumber(x,isize,n)
!  if (isize == n) then
!    return
!  endif
  if ( isize < 32) then
    call insertion_$Xsort(x,iorder,isize)
  else
!    call heap_$Xsort(x,iorder,isize)
    call quick_$Xsort(x,iorder,isize)
  endif
 end subroutine $Xsort

SUBST [ X, type ]
 d ; double precision ;;
END_TEMPLATE

BEGIN_TEMPLATE

 subroutine $Xsort(x,iorder,isize)
  implicit none
  BEGIN_DOC
  ! Sort array x(isize).
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  END_DOC
  integer,intent(in)             :: isize
  $type,intent(inout)            :: x(isize)
  integer,intent(inout)          :: iorder(isize)
  integer                        :: n
  if (isize < 2) then
    return
  endif
  call sorted_$Xnumber(x,isize,n)
  if (isize == n) then
    return
  endif
  if ( isize < 32) then
    call insertion_$Xsort(x,iorder,isize)
  else
!    call $Xradix_sort(x,iorder,isize,-1)
    call quick_$Xsort(x,iorder,isize)
  endif
 end subroutine $Xsort

SUBST [ X, type ]
 i8 ; integer*8 ;;
END_TEMPLATE

!---------------------- END INTEL
IRP_ELSE
!---------------------- NON-INTEL
BEGIN_TEMPLATE

 subroutine $Xsort_noidx(x,isize)
  implicit none
  BEGIN_DOC
  ! Sort array x(isize).
  END_DOC
  integer,intent(in)             :: isize
  $type,intent(inout)            :: x(isize)
  integer, allocatable           :: iorder(:)
  integer                        :: i
  allocate(iorder(isize))
  do i=1,isize
   iorder(i)=i
  enddo
  call $Xsort(x,iorder,isize)
  deallocate(iorder)
 end subroutine $Xsort_noidx

SUBST [ X, type ]
   ; real ;;
 d ; double precision ;;
 i ; integer ;;
 i8 ; integer*8 ;;
 i2 ; integer*2 ;;
END_TEMPLATE

BEGIN_TEMPLATE

 subroutine $Xsort(x,iorder,isize)
  implicit none
  BEGIN_DOC
  ! Sort array x(isize).
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  END_DOC
  integer,intent(in)             :: isize
  $type,intent(inout)            :: x(isize)
  integer,intent(inout)          :: iorder(isize)
  integer                        :: n
  if (isize < 2) then
    return
  endif
!  call sorted_$Xnumber(x,isize,n)
!  if (isize == n) then
!    return
!  endif
  if ( isize < 32) then
    call insertion_$Xsort(x,iorder,isize)
  else
!    call heap_$Xsort(x,iorder,isize)
    call quick_$Xsort(x,iorder,isize)
  endif
 end subroutine $Xsort

SUBST [ X, type ]
   ; real ;;
 d ; double precision ;;
END_TEMPLATE

BEGIN_TEMPLATE

 subroutine $Xsort(x,iorder,isize)
  implicit none
  BEGIN_DOC
  ! Sort array x(isize).
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  END_DOC
  integer,intent(in)             :: isize
  $type,intent(inout)            :: x(isize)
  integer,intent(inout)          :: iorder(isize)
  integer                        :: n
  if (isize < 2) then
    return
  endif
  call sorted_$Xnumber(x,isize,n)
  if (isize == n) then
    return
  endif
  if ( isize < 32) then
    call insertion_$Xsort(x,iorder,isize)
  else
!    call $Xradix_sort(x,iorder,isize,-1)
    call quick_$Xsort(x,iorder,isize)
  endif
 end subroutine $Xsort

SUBST [ X, type ]
 i ; integer ;;
 i8 ; integer*8 ;;
 i2 ; integer*2 ;;
END_TEMPLATE

IRP_ENDIF
!---------------------- END NON-INTEL



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

SUBST [ X, type ]
   ; real ;;
 d ; double precision ;;
 i ; integer ;;
 i8; integer*8 ;;
 i2; integer*2 ;;
END_TEMPLATE


BEGIN_TEMPLATE
 subroutine insertion_$Xsort_big (x,iorder,isize)
  implicit none
  BEGIN_DOC
  ! Sort array x(isize) using the insertion sort algorithm.
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  ! This is a version for very large arrays where the indices need
  ! to be in integer*8 format
  END_DOC
  integer*8,intent(in)           :: isize
  $type,intent(inout)            :: x(isize)
  integer*8,intent(inout)        :: iorder(isize)
  $type                          :: xtmp
  integer*8                      :: i, i0, j, jmax

  do i=2_8,isize
    xtmp = x(i)
    i0 = iorder(i)
    j = i-1_8
    do while (j>0_8)
      if (x(j)<=xtmp) exit
      x(j+1_8) = x(j)
      iorder(j+1_8) = iorder(j)
      j = j-1_8
    enddo
    x(j+1_8) = xtmp
    iorder(j+1_8) = i0
  enddo

 end subroutine insertion_$Xsort_big

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


BEGIN_TEMPLATE

recursive subroutine $Xradix_sort$big(x,iorder,isize,iradix)
  implicit none

  BEGIN_DOC
  ! Sort integer array x(isize) using the radix sort algorithm.
  ! iorder in input should be (1,2,3,...,isize), and in output
  ! contains the new order of the elements.
  ! iradix should be -1 in input.
  END_DOC
  integer*$int_type, intent(in)  :: isize
  integer*$int_type, intent(inout) :: iorder(isize)
  integer*$type, intent(inout)   :: x(isize)
  integer, intent(in)            :: iradix
  integer                        :: iradix_new
  integer*$type, allocatable     :: x2(:), x1(:)
  integer*$type                  :: i4               ! data type
  integer*$int_type, allocatable :: iorder1(:),iorder2(:)
  integer*$int_type              :: i0, i1, i2, i3, i ! index type
  integer*$type                  :: mask
  integer                        :: err
  !DIR$ ATTRIBUTES ALIGN : 128   :: iorder1,iorder2, x2, x1

  if (isize < 2) then
    return
  endif

  if (iradix == -1) then ! Sort Positive and negative

    allocate(x1(isize),iorder1(isize), x2(isize),iorder2(isize),stat=err)
    if (err /= 0) then
      print *,  irp_here, ': Unable to allocate arrays'
      stop
    endif

    i1=1_$int_type
    i2=1_$int_type
    do i=1_$int_type,isize
      if (x(i) < 0_$type) then
        iorder1(i1) = iorder(i)
        x1(i1) = -x(i)
        i1 = i1+1_$int_type
      else
        iorder2(i2) = iorder(i)
        x2(i2) = x(i)
        i2 = i2+1_$int_type
      endif
    enddo
    i1=i1-1_$int_type
    i2=i2-1_$int_type

    do i=1_$int_type,i2
      iorder(i1+i) = iorder2(i)
      x(i1+i) = x2(i)
    enddo
    deallocate(x2,iorder2,stat=err)
    if (err /= 0) then
      print *,  irp_here, ': Unable to deallocate arrays x2, iorder2'
      stop
    endif


    if (i1 > 1_$int_type) then
      call $Xradix_sort$big(x1,iorder1,i1,-2)
      do i=1_$int_type,i1
        x(i) = -x1(1_$int_type+i1-i)
        iorder(i) = iorder1(1_$int_type+i1-i)
      enddo
    endif

    if (i2>1_$int_type) then
      call $Xradix_sort$big(x(i1+1_$int_type),iorder(i1+1_$int_type),i2,-2)
    endif

    deallocate(x1,iorder1,stat=err)
    if (err /= 0) then
      print *,  irp_here, ': Unable to deallocate arrays x1, iorder1'
      stop
    endif
    return

  else if (iradix == -2) then ! Positive

    ! Find most significant bit

    i0 = 0_$int_type
    i4 = maxval(x)

    iradix_new = max($integer_size-1-leadz(i4),1)
    mask = ibset(0_$type,iradix_new)

    allocate(x1(isize),iorder1(isize), x2(isize),iorder2(isize),stat=err)
    if (err /= 0) then
      print *,  irp_here, ': Unable to allocate arrays'
      stop
    endif

    i1=1_$int_type
    i2=1_$int_type

    do i=1_$int_type,isize
      if (iand(mask,x(i)) == 0_$type) then
        iorder1(i1) = iorder(i)
        x1(i1) = x(i)
        i1 = i1+1_$int_type
      else
        iorder2(i2) = iorder(i)
        x2(i2) = x(i)
        i2 = i2+1_$int_type
      endif
    enddo
    i1=i1-1_$int_type
    i2=i2-1_$int_type

    do i=1_$int_type,i1
      iorder(i0+i) = iorder1(i)
      x(i0+i) = x1(i)
    enddo
    i0 = i0+i1
    i3 = i0
    deallocate(x1,iorder1,stat=err)
    if (err /= 0) then
      print *,  irp_here, ': Unable to deallocate arrays x1, iorder1'
      stop
    endif


    do i=1_$int_type,i2
      iorder(i0+i) = iorder2(i)
      x(i0+i) = x2(i)
    enddo
    i0 = i0+i2
    deallocate(x2,iorder2,stat=err)
    if (err /= 0) then
      print *,  irp_here, ': Unable to deallocate arrays x2, iorder2'
      stop
    endif


    if (i3>1_$int_type) then
      call $Xradix_sort$big(x,iorder,i3,iradix_new-1)
    endif

    if (isize-i3>1_$int_type) then
      call $Xradix_sort$big(x(i3+1_$int_type),iorder(i3+1_$int_type),isize-i3,iradix_new-1)
    endif

    return
  endif

  ASSERT (iradix >= 0)

  if (isize < 48) then
    call insertion_$Xsort$big(x,iorder,isize)
    return
  endif


  allocate(x2(isize),iorder2(isize),stat=err)
  if (err /= 0) then
    print *,  irp_here, ': Unable to allocate arrays x1, iorder1'
    stop
  endif


  mask = ibset(0_$type,iradix)
  i0=1_$int_type
  i1=1_$int_type

  do i=1_$int_type,isize
    if (iand(mask,x(i)) == 0_$type) then
      iorder(i0) = iorder(i)
      x(i0) = x(i)
      i0 = i0+1_$int_type
    else
      iorder2(i1) = iorder(i)
      x2(i1) = x(i)
      i1 = i1+1_$int_type
    endif
  enddo
  i0=i0-1_$int_type
  i1=i1-1_$int_type

  do i=1_$int_type,i1
    iorder(i0+i) = iorder2(i)
    x(i0+i) = x2(i)
  enddo

  deallocate(x2,iorder2,stat=err)
  if (err /= 0) then
    print *,  irp_here, ': Unable to allocate arrays x2, iorder2'
    stop
  endif


  if (iradix == 0) then
    return
  endif


  if (i1>1_$int_type) then
    call $Xradix_sort$big(x(i0+1_$int_type),iorder(i0+1_$int_type),i1,iradix-1)
  endif
  if (i0>1) then
    call $Xradix_sort$big(x,iorder,i0,iradix-1)
  endif

 end

SUBST [ X, type, integer_size, is_big, big, int_type ]
 i  ; 4 ; 32 ; .False. ;      ; 4 ;;
 i8 ; 8 ; 64 ; .False. ;      ; 4 ;;
 i2 ; 2 ; 16 ; .False. ;      ; 4 ;;
 i  ; 4 ; 32 ; .True.  ; _big ; 8 ;;
 i8 ; 8 ; 64 ; .True.  ; _big ; 8 ;;
END_TEMPLATE



