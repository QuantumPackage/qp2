! phase

subroutine get_phase_general(det1,det2,phase,degree,Nint)
  implicit none

  integer, intent(in)           :: Nint
  integer(bit_kind), intent(in) :: det1(Nint,2), det2(Nint,2)
  double precision, intent(out) :: phase
  integer, intent(out)          :: degree
  integer :: n(2)
  integer, allocatable :: list_anni(:,:), list_crea(:,:)

  allocate(list_anni(N_int*bit_kind_size,2))
  allocate(list_crea(N_int*bit_kind_size,2))

  call get_excitation_general(det1,det2,degree,n,list_anni,list_crea,phase,Nint)
end

! Get excitation general

subroutine get_excitation_general(det1,det2,degree,n,list_anni,list_crea,phase,Nint)

  use bitmasks

  implicit none

  integer, intent(in) :: Nint
  integer(bit_kind), intent(in)  :: det1(Nint,2), det2(Nint,2)
  double precision, intent(out)  :: phase
  integer, intent(out)           :: list_crea(Nint*bit_kind_size,2)
  integer, intent(out)           :: list_anni(Nint*bit_kind_size,2)
  integer, intent(out)           :: degree, n(2)
  
  integer, allocatable           :: l1(:,:), l2(:,:) 
  integer(bit_kind), allocatable :: det_crea(:,:), det_anni(:,:)
  integer, allocatable           :: pos_anni(:,:), pos_crea(:,:)

  integer :: n1(2),n2(2),n_crea(2),n_anni(2),i,j,k,d

  allocate(l1(Nint*bit_kind_size,2))
  allocate(l2(Nint*bit_kind_size,2))
  allocate(det_crea(Nint,2),det_anni(Nint,2))

  ! 1      111010
  ! 2      110101
  !
  !not 1-> 000101
  !    2   110101
  !and     000101 -> crea
  !
  !    1   111010
  !not 2-> 001010
  !        001010 -> anni

  do j = 1, 2
    do i = 1, Nint
      det_crea(i,j) = iand(not(det1(i,j)),det2(i,j))
    enddo
  enddo
  
  do j = 1, 2
    do i = 1, Nint
      det_anni(i,j) = iand(det1(i,j),not(det2(i,j)))
    enddo
  enddo

  call bitstring_to_list_ab(det1,l1,n1,Nint)
  call bitstring_to_list_ab(det2,l2,n2,Nint)
  call bitstring_to_list_ab(det_crea,list_crea,n_crea,Nint)
  call bitstring_to_list_ab(det_anni,list_anni,n_anni,Nint)

  do i = 1, 2
    if (n_crea(i) /= n_anni(i)) then
      print*,'Well, it seems we have a problem here...'
      call abort
    endif
  enddo

  !1    11110011001  1 2 3 4 7 8  11
  !pos               1 2 3 4 5 6  7 
  !2    11100101011  1 2 3 6 8 10 11
  !anni 00010010000  4 7
  !pos               4 5
  !crea 00000100010  6 10
  !pos               4 6
  !4 -> 6  pos(4 -> 4)
  !7 -> 10 pos(5 -> 6)

  n = n_anni
  degree = n_anni(1) + n_anni(2)

  allocate(pos_anni(max(n(1),n(2)),2))
  allocate(pos_crea(max(n(1),n(2)),2))
  
  ! Search pos anni
  do j = 1, 2
    k = 1
    do i = 1, n1(j)
       if (k > n_anni(j)) exit
       if (l1(i,j) /= list_anni(k,j)) cycle
       pos_anni(k,j) = i
       k = k + 1
    enddo
  enddo

  ! Search pos crea
  do j = 1, 2
    k = 1
    do i = 1, n2(j)
       if (k > n_crea(j)) exit
       if (l2(i,j) /= list_crea(k,j)) cycle
       pos_crea(k,j) = i
       k = k + 1
    enddo
  enddo

  ! Distance between the ith anni and the ith crea op
  ! By doing so there is no crossing between the different pairs of anni/crea
  ! and the phase is determined by the sum of the distances
  ! -> (-1)^{sum of the distances}
  d = 0
  do j = 1, 2
    do i = 1, n(j)
      d = d + abs(pos_anni(i,j) - pos_crea(i,j))
    enddo
  enddo
  
  phase = dble((-1)**d)

  ! Debug
  !print*,l2(1:n2(1),1)
  !print*,l2(1:n2(2),2)
  !!call print_det(det1,Nint)
  !!call print_det(det2,Nint)
  !print*,phase
  !print*,''
end
