
subroutine create_selection_buffer(N, size_in, res)
  use selection_types
  implicit none
  BEGIN_DOC
! Allocates the memory for a selection buffer.
! The arrays have dimension size_in and the maximum number of elements is N
  END_DOC

  integer, intent(in) :: N, size_in
  type(selection_buffer), intent(out) :: res

  integer :: siz
  siz = max(size_in,1)

  double precision :: rss
  double precision, external :: memory_of_double
  rss = memory_of_double(siz)*(N_int*2+1)
  call check_mem(rss,irp_here)

  allocate(res%det(N_int, 2, siz), res%val(siz))

  res%val(:) = 0d0
  res%det(:,:,:) = 0_8
  res%N = N
  res%mini = 0d0
  res%cur = 0
end subroutine

subroutine delete_selection_buffer(b)
  use selection_types
  implicit none
  type(selection_buffer), intent(inout) :: b
  if (associated(b%det)) then
    deallocate(b%det)
  endif
  if (associated(b%val)) then
    deallocate(b%val)
  endif
  NULLIFY(b%det)
  NULLIFY(b%val)
  b%cur = 0
  b%mini = 0.d0
  b%N = 0
end


subroutine add_to_selection_buffer(b, det, val)
  use selection_types
  implicit none

  type(selection_buffer), intent(inout) :: b
  integer(bit_kind), intent(in) :: det(N_int, 2)
  double precision, intent(in) :: val
  integer :: i

  if(b%N > 0 .and. val <= b%mini) then
    b%cur += 1
    b%det(1:N_int,1:2,b%cur) = det(1:N_int,1:2)
    b%val(b%cur) = val
    if(b%cur == size(b%val)) then
      call sort_selection_buffer(b)
    end if
  end if
end subroutine

subroutine merge_selection_buffers(b1, b2)
  use selection_types
  implicit none
  BEGIN_DOC
! Merges the selection buffers b1 and b2 into b2
  END_DOC
  type(selection_buffer), intent(inout) :: b1
  type(selection_buffer), intent(inout) :: b2
  integer(bit_kind), pointer     :: detmp(:,:,:)
  double precision, pointer      :: val(:)
  integer                        :: i, i1, i2, k, nmwen, sze
  if (b1%cur == 0) return
  do while (b1%val(b1%cur) > b2%mini)
    b1%cur = b1%cur-1
    if (b1%cur == 0) then
      return
    endif
  enddo
  nmwen = min(b1%N, b1%cur+b2%cur)
  double precision :: rss
  double precision, external :: memory_of_double
  sze = max(size(b1%val), size(b2%val))
  rss = memory_of_double(sze) + 2*N_int*memory_of_double(sze)
  call check_mem(rss,irp_here)
  allocate(val(sze), detmp(N_int, 2, sze))
  i1=1
  i2=1
  do i=1,nmwen
    if ( (i1 > b1%cur).and.(i2 > b2%cur) ) then
      exit
    else if (i1 > b1%cur) then
        val(i) = b2%val(i2)
        detmp(1:N_int,1,i) = b2%det(1:N_int,1,i2)
        detmp(1:N_int,2,i) = b2%det(1:N_int,2,i2)
        i2=i2+1
    else if (i2 > b2%cur) then
        val(i) = b1%val(i1)
        detmp(1:N_int,1,i) = b1%det(1:N_int,1,i1)
        detmp(1:N_int,2,i) = b1%det(1:N_int,2,i1)
        i1=i1+1
    else
      if (b1%val(i1) <= b2%val(i2)) then
        val(i) = b1%val(i1)
        detmp(1:N_int,1,i) = b1%det(1:N_int,1,i1)
        detmp(1:N_int,2,i) = b1%det(1:N_int,2,i1)
        i1=i1+1
      else
        val(i) = b2%val(i2)
        detmp(1:N_int,1,i) = b2%det(1:N_int,1,i2)
        detmp(1:N_int,2,i) = b2%det(1:N_int,2,i2)
        i2=i2+1
      endif
    endif
  enddo
  deallocate(b2%det, b2%val)
  do i=nmwen+1,b2%N
    val(i) = 0.d0
    detmp(1:N_int,1:2,i) = 0_bit_kind
  enddo
  b2%det => detmp
  b2%val => val
  b2%mini = min(b2%mini,b2%val(b2%N))
  b2%cur = nmwen
end


subroutine sort_selection_buffer(b)
  use selection_types
  implicit none

  type(selection_buffer), intent(inout) :: b
  integer, allocatable :: iorder(:)
  integer(bit_kind), pointer :: detmp(:,:,:)
  integer :: i, nmwen
  logical, external :: detEq
  if (b%N == 0 .or. b%cur == 0) return
  nmwen = min(b%N, b%cur)

  double precision :: rss
  double precision, external :: memory_of_double, memory_of_int
  rss = memory_of_int(b%cur) + 2*N_int*memory_of_double(size(b%det,3))
  call check_mem(rss,irp_here)
  allocate(iorder(b%cur), detmp(N_int, 2, size(b%det,3)))
  do i=1,b%cur
    iorder(i) = i
  end do
  call dsort(b%val, iorder, b%cur)
  do i=1, nmwen
    detmp(1:N_int,1,i) = b%det(1:N_int,1,iorder(i))
    detmp(1:N_int,2,i) = b%det(1:N_int,2,iorder(i))
  end do
  deallocate(b%det,iorder)
  b%det => detmp
  b%mini = min(b%mini,b%val(b%N))
  b%cur = nmwen
end subroutine

subroutine make_selection_buffer_s2(b)
  use selection_types
  type(selection_buffer), intent(inout) :: b

  integer(bit_kind), allocatable :: o(:,:,:)
  double precision, allocatable  :: val(:)

  integer :: n_d
  integer :: i,k,sze,n_alpha,j,n
  logical                        :: dup

  ! Sort
  integer, allocatable           :: iorder(:)
  integer*8, allocatable         :: bit_tmp(:)
  integer*8, external            :: configuration_search_key
  integer(bit_kind), allocatable :: tmp_array(:,:,:)
  logical, allocatable           :: duplicate(:)

  n_d = b%cur
  double precision :: rss
  double precision, external :: memory_of_double
  rss = (4*N_int+4)*memory_of_double(n_d)
  call check_mem(rss,irp_here)
  allocate(o(N_int,2,n_d), iorder(n_d), duplicate(n_d), bit_tmp(n_d), &
           tmp_array(N_int,2,n_d), val(n_d) )

  do i=1,n_d
    do k=1,N_int
      o(k,1,i) = ieor(b%det(k,1,i), b%det(k,2,i))
      o(k,2,i) = iand(b%det(k,1,i), b%det(k,2,i))
    enddo
    iorder(i) = i
    bit_tmp(i) = configuration_search_key(o(1,1,i),N_int)
  enddo

  deallocate(b%det)

  call i8sort(bit_tmp,iorder,n_d)

  do i=1,n_d
    do k=1,N_int
      tmp_array(k,1,i) = o(k,1,iorder(i))
      tmp_array(k,2,i) = o(k,2,iorder(i))
    enddo
    val(i) = b%val(iorder(i))
    duplicate(i) = .False.
  enddo

  ! Find duplicates
  do i=1,n_d-1
    if (duplicate(i)) then
      cycle
    endif
    j = i+1
    do while (bit_tmp(j)==bit_tmp(i))
      if (duplicate(j)) then
        j+=1
        if (j>n_d) then
          exit
        endif
        cycle
      endif
      dup = .True.
      do k=1,N_int
        if ( (tmp_array(k,1,i) /= tmp_array(k,1,j))                   &
              .or. (tmp_array(k,2,i) /= tmp_array(k,2,j)) ) then
          dup = .False.
          exit
        endif
      enddo
      if (dup) then
        val(i) = max(val(i), val(j))
        duplicate(j) = .True.
      endif
      j+=1
      if (j>n_d) then
        exit
      endif
    enddo
  enddo

  deallocate (b%val)
  ! Copy filtered result
  integer :: n_p
  n_p=0
  do i=1,n_d
    if (duplicate(i)) then
      cycle
    endif
    n_p = n_p + 1
    do k=1,N_int
      o(k,1,n_p) = tmp_array(k,1,i)
      o(k,2,n_p) = tmp_array(k,2,i)
    enddo
    val(n_p) = val(i)
  enddo

  ! Sort by importance
  do i=1,n_p
    iorder(i) = i
  end do
  call dsort(val,iorder,n_p)
  do i=1,n_p
    do k=1,N_int
      tmp_array(k,1,i) = o(k,1,iorder(i))
      tmp_array(k,2,i) = o(k,2,iorder(i))
    enddo
  enddo
  do i=1,n_p
    do k=1,N_int
      o(k,1,i) = tmp_array(k,1,i)
      o(k,2,i) = tmp_array(k,2,i)
    enddo
  enddo

  ! Create determinants
  n_d = 0
  do i=1,n_p
    call configuration_to_dets_size(o(1,1,i),sze,elec_alpha_num,N_int)
    n_d = n_d + sze
    if (n_d > b%cur) then
!      if (n_d - b%cur > b%cur - n_d + sze) then
!        n_d = n_d - sze
!      endif
      exit
    endif
  enddo

  rss = (4*N_int+2)*memory_of_double(n_d)
  call check_mem(rss,irp_here)
  allocate(b%det(N_int,2,2*n_d), b%val(2*n_d))
  k=1
  do i=1,n_p
    n=n_d
    call configuration_to_dets_size(o(1,1,i),n,elec_alpha_num,N_int)
    call configuration_to_dets(o(1,1,i),b%det(1,1,k),n,elec_alpha_num,N_int)
    do j=k,k+n-1
      b%val(j) = val(i)
    enddo
    k = k+n
    if (k > n_d) exit
  enddo
  deallocate(o)
  b%cur = n_d
  b%N = n_d
end




subroutine remove_duplicates_in_selection_buffer(b)
  use selection_types
  type(selection_buffer), intent(inout) :: b

  integer(bit_kind), allocatable :: o(:,:,:)
  double precision, allocatable  :: val(:)

  integer :: n_d
  integer :: i,k,sze,n_alpha,j,n
  logical                        :: dup

  ! Sort
  integer, allocatable           :: iorder(:)
  integer*8, allocatable         :: bit_tmp(:)
  integer*8, external            :: det_search_key
  integer(bit_kind), allocatable :: tmp_array(:,:,:)
  logical, allocatable           :: duplicate(:)

  n_d = b%cur
  logical                        :: found_duplicates
  double precision               :: rss
  double precision, external     :: memory_of_double
  rss = (4*N_int+4)*memory_of_double(n_d)
  call check_mem(rss,irp_here)

  found_duplicates = .False.
  allocate(iorder(n_d), duplicate(n_d), bit_tmp(n_d), &
           tmp_array(N_int,2,n_d), val(n_d) )

  do i=1,n_d
    iorder(i) = i
    bit_tmp(i) = det_search_key(b%det(1,1,i),N_int)
  enddo

  call i8sort(bit_tmp,iorder,n_d)

  do i=1,n_d
    do k=1,N_int
      tmp_array(k,1,i) = b%det(k,1,iorder(i))
      tmp_array(k,2,i) = b%det(k,2,iorder(i))
    enddo
    val(i) = b%val(iorder(i))
    duplicate(i) = .False.
  enddo

  ! Find duplicates
  do i=1,n_d-1
    if (duplicate(i)) then
      cycle
    endif
    j = i+1
    do while (bit_tmp(j)==bit_tmp(i))
      if (duplicate(j)) then
        j+=1
        if (j>n_d) then
          exit
        endif
        cycle
      endif
      dup = .True.
      do k=1,N_int
        if ( (tmp_array(k,1,i) /= tmp_array(k,1,j))                   &
              .or. (tmp_array(k,2,i) /= tmp_array(k,2,j)) ) then
          dup = .False.
          exit
        endif
      enddo
      if (dup) then
        duplicate(j) = .True.
        found_duplicates = .True.
      endif
      j+=1
      if (j>n_d) then
        exit
      endif
    enddo
  enddo

  if (found_duplicates) then

    ! Copy filtered result
    integer :: n_p
    n_p=0
    do i=1,n_d
      if (duplicate(i)) then
        cycle
      endif
      n_p = n_p + 1
      do k=1,N_int
        b%det(k,1,n_p) = tmp_array(k,1,i)
        b%det(k,2,n_p) = tmp_array(k,2,i)
      enddo
      val(n_p) = val(i)
    enddo
    b%cur=n_p
    b%N=n_p

  endif

end



