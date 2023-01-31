  use bitmasks

subroutine filter_connected_array(key1,key2,ld,Nint,sze,idx)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Filters out the determinants that are not connected by H
  !
  ! returns the array idx which contains the index of the
  !
  ! determinants in the array key1 that interact
  !
  ! via the H operator with key2.
  !
  ! idx(0) is the number of determinants that interact with key1
  END_DOC
  integer, intent(in)            :: Nint, ld,sze
  integer(bit_kind), intent(in)  :: key1(Nint,2,ld)
  integer(bit_kind), intent(in)  :: key2(Nint,2)
  integer, intent(out)           :: idx(0:sze)

  integer                        :: i,j,l
  integer                        :: degree_x2

  ASSERT (Nint > 0)
  ASSERT (sze >= 0)

  l=1

  if (Nint==1) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      degree_x2 = popcnt(    xor( key1(1,1,i), key2(1,1))) &
                + popcnt(    xor( key1(1,2,i), key2(1,2)))
!      print*,degree_x2
      if (degree_x2 > 4) then
        cycle
      else
        idx(l) = i
        l = l+1
      endif
    enddo

  else if (Nint==2) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      degree_x2 =  popcnt(xor( key1(1,1,i), key2(1,1))) +            &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2)))
      if (degree_x2 > 4) then
        cycle
      else
        idx(l) = i
        l = l+1
      endif
    enddo

  else if (Nint==3) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      degree_x2 = popcnt(xor( key1(1,1,i), key2(1,1))) +             &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2))) +                     &
          popcnt(xor( key1(3,1,i), key2(3,1))) +                     &
          popcnt(xor( key1(3,2,i), key2(3,2)))
      if (degree_x2 > 4) then
        cycle
      else
        idx(l) = i
        l = l+1
      endif
    enddo

  else

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      degree_x2 = 0
      !DIR$ LOOP COUNT MIN(4)
      do j=1,Nint
        degree_x2 = degree_x2+ popcnt(xor( key1(j,1,i), key2(j,1))) +&
            popcnt(xor( key1(j,2,i), key2(j,2)))
        if (degree_x2 > 4) then
          exit
        endif
      enddo
      if (degree_x2 <= 5) then
        idx(l) = i
        l = l+1
      endif
    enddo

  endif
  idx(0) = l-1
!  print*,'idx(0) = ',idx(0)
end


 BEGIN_PROVIDER [ integer, n_sparse_mat]
&BEGIN_PROVIDER [ integer, n_connected_per_det, (N_det)]
&BEGIN_PROVIDER [ integer, n_max_connected_per_det]
 implicit none
 BEGIN_DOC
! n_sparse_mat = total number of connections in the CI matrix 
!
! n_connected_per_det(i) = number of connected determinants to the determinant psi_det(1,1,i)
!
! n_max_connected_per_det = maximum number of connected determinants 
 END_DOC
 integer, allocatable :: idx(:)
 allocate(idx(0:N_det))
 integer :: i
 n_sparse_mat = 0
 do i = 1, N_det
  call filter_connected_array(psi_det_sorted,psi_det_sorted(1,1,i),psi_det_size,N_int,N_det,idx)
  n_connected_per_det(i) = idx(0)
  n_sparse_mat += idx(0)
 enddo
 n_max_connected_per_det = maxval(n_connected_per_det)
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), connected_det_per_det, (N_int,2,n_max_connected_per_det,N_det)]
&BEGIN_PROVIDER [ integer(bit_kind), list_connected_det_per_det, (n_max_connected_per_det,N_det)]
 implicit none
 BEGIN_DOC
! connected_det_per_det(:,:,j,i) = jth connected determinant to the determinant psi_det(:,:,i)
!
! list_connected_det_per_det(j,i) = index of jth determinant in psi_det which is connected to psi_det(:,:,i)
 END_DOC
 integer, allocatable :: idx(:)
 allocate(idx(0:N_det))
 integer :: i,j
 do i = 1, N_det
  call filter_connected_array(psi_det_sorted,psi_det_sorted(1,1,i),psi_det_size,N_int,N_det,idx)
  do j = 1, idx(0)
   connected_det_per_det(1:N_int,1:2,j,i) = psi_det_sorted(1:N_int,1:2,idx(j))
   list_connected_det_per_det(j,i) = idx(j)
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, sparse_h_mat, (n_max_connected_per_det, N_det)]
 implicit none
 BEGIN_DOC
! sparse matrix format 
!
! sparse_h_mat(j,i) = matrix element between the jth connected determinant and psi_det(:,:,i)
 END_DOC
 integer :: i,j 
 double precision :: hij
 do i = 1, N_det
  do j = 1, n_connected_per_det(i)
   call i_H_j(psi_det(1,1,i),connected_det_per_det(1,1,j,i),N_int,hij)
   sparse_h_mat(j,i) = hij
  enddo
 enddo

END_PROVIDER 

