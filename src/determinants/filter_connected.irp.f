
subroutine filter_not_connected(key1,key2,Nint,sze,idx)
  use bitmasks
 implicit none
  BEGIN_DOC
  ! Returns the array idx which contains the index of the
  !
  ! determinants in the array key1 that DO NOT interact
  !
  ! via the H operator with key2.
  !
  ! idx(0) is the number of determinants that DO NOT interact with key1
  END_DOC
  integer, intent(in)            :: Nint, sze
  integer(bit_kind), intent(in)  :: key1(Nint,2,sze)
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
      if (degree_x2 > 4) then
        idx(l) = i
        l = l+1
      else
        cycle
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
        idx(l) = i
        l = l+1
      else
        cycle
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
        idx(l) = i
        l = l+1
      else
        cycle
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
        idx(l) = i
        l = l+1
        endif
      enddo
      if (degree_x2 <= 5) then
          exit
      endif
    enddo

  endif
  idx(0) = l-1
end

subroutine filter_connected(key1,key2,Nint,sze,idx)
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
  integer, intent(in)            :: Nint, sze
  integer(bit_kind), intent(in)  :: key1(Nint,2,sze)
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
end


subroutine getMobiles(key,key_mask, mobiles,Nint)
  use bitmasks
  implicit none
  integer(bit_kind),intent(in) :: key(Nint,2), key_mask(Nint,2)
  integer,intent(out) :: mobiles(2)
  integer,intent(in) :: Nint

  integer(bit_kind) :: mobileMask(Nint,2)
  integer :: list(Nint*bit_kind_size), nel,j

  do j=1,Nint
    mobileMask(j,1) = xor(key(j,1), key_mask(j,1))
    mobileMask(j,2) = xor(key(j,2), key_mask(j,2))
  end do

  call bitstring_to_list(mobileMask(1,1), list, nel, Nint)
  if(nel == 2) then
    mobiles(1) = list(1)
    mobiles(2) = list(2)
  else if(nel == 1) then
    mobiles(1) = list(1)
    call bitstring_to_list(mobileMask(1,2), list, nel, Nint)
    mobiles(2) = list(1) + mo_num
  else
    call bitstring_to_list(mobileMask(1,2), list, nel, Nint)
    mobiles(1) = list(1) + mo_num
    mobiles(2) = list(2) + mo_num
  end if
end subroutine


subroutine create_microlist(minilist, N_minilist, key_mask, microlist, idx_microlist, N_microlist, ptr_microlist, Nint)
  use bitmasks
  implicit none
  integer, intent(in) :: Nint, N_minilist
  integer(bit_kind), intent(in) :: minilist(Nint,2,N_minilist), key_mask(Nint,2)

  integer, intent(out) :: N_microlist(0:mo_num*2), ptr_microlist(0:mo_num*2+1), idx_microlist(N_minilist*4)
  integer(bit_kind), intent(out) :: microlist(Nint,2,N_minilist*4)

  integer :: i,j,k,nt,n_element(2)
  integer :: list(Nint*bit_kind_size,2)
  integer, allocatable :: cur_microlist(:)
  allocate (cur_microlist(0:mo_num*2+1))
  integer(bit_kind) :: key_mask_neg(Nint,2), mobileMask(Nint,2)
  integer :: mo_num_2
  mo_num_2 = mo_num+mo_num


  do i=1,Nint
    key_mask_neg(i,1) = not(key_mask(i,1))
    key_mask_neg(i,2) = not(key_mask(i,2))
  end do

  do i=0,mo_num_2
    N_microlist(i) = 0
  enddo

  do i=1, N_minilist
    do j=1,Nint
      mobileMask(j,1) = iand(key_mask_neg(j,1), minilist(j,1,i))
      mobileMask(j,2) = iand(key_mask_neg(j,2), minilist(j,2,i))
    end do

    call bitstring_to_list(mobileMask(1,1), list(1,1), n_element(1), Nint)
    call bitstring_to_list(mobileMask(1,2), list(1,2), n_element(2), Nint)

    if(n_element(1) + n_element(2) /= 4) then
      N_microlist(0) = N_microlist(0) + 1
    else
      do j=1,n_element(1)
        nt = list(j,1)
        N_microlist(nt) = N_microlist(nt) + 1
      end do

      do j=1,n_element(2)
        nt = list(j,2) + mo_num
        N_microlist(nt) = N_microlist(nt) + 1
      end do
    end if
  end do

  ptr_microlist(0) = 1
  do i=1,mo_num_2+1
    ptr_microlist(i) = ptr_microlist(i-1) + N_microlist(i-1)
  end do

  do i=0,mo_num_2+1
    cur_microlist(i) = ptr_microlist(i)
  end do


  do i=1, N_minilist
    do j=1,Nint
      mobileMask(j,1) = iand(key_mask_neg(j,1), minilist(j,1,i))
      mobileMask(j,2) = iand(key_mask_neg(j,2), minilist(j,2,i))
    end do

    call bitstring_to_list(mobileMask(1,1), list(1,1), n_element(1), Nint)
    call bitstring_to_list(mobileMask(1,2), list(1,2), n_element(2), Nint)


    if(n_element(1) + n_element(2) /= 4) then
      idx_microlist(cur_microlist(0)) = i
      do k=1,Nint
        microlist(k,1,cur_microlist(0)) = minilist(k,1,i)
        microlist(k,2,cur_microlist(0)) = minilist(k,2,i)
      enddo
      cur_microlist(0) = cur_microlist(0) + 1
    else
      do j=1,n_element(1)
        nt = list(j,1)
        idx_microlist(cur_microlist(nt)) = i
        ! TODO : Page faults
        do k=1,Nint
          microlist(k,1,cur_microlist(nt)) = minilist(k,1,i)
          microlist(k,2,cur_microlist(nt)) = minilist(k,2,i)
        enddo
        cur_microlist(nt) = cur_microlist(nt) + 1
      end do

      do j=1,n_element(2)
        nt = list(j,2) + mo_num
        idx_microlist(cur_microlist(nt)) = i
        do k=1,Nint
          microlist(k,1,cur_microlist(nt)) = minilist(k,1,i)
          microlist(k,2,cur_microlist(nt)) = minilist(k,2,i)
        enddo
        cur_microlist(nt) = cur_microlist(nt) + 1
      end do
    end if
  end do
  deallocate(cur_microlist)
end subroutine


subroutine filter_connected_i_H_psi0(key1,key2,Nint,sze,idx)
  use bitmasks
  BEGIN_DOC
  ! Returns the array idx which contains the index of the
  !
  ! determinants in the array key1 that interact
  !
  ! via the H operator with key2.
  !
  ! idx(0) is the number of determinants that interact with key1
  END_DOC
  implicit none
  integer, intent(in)            :: Nint, sze
  integer(bit_kind), intent(in)  :: key1(Nint,2,sze)
  integer(bit_kind), intent(in)  :: key2(Nint,2)
  integer, intent(out)           :: idx(0:sze)

  integer                        :: i,l,m
  integer                        :: degree_x2

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sze > 0)

  l=1

  if (Nint==1) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      degree_x2 = popcnt(xor( key1(1,1,i), key2(1,1))) +             &
          popcnt(xor( key1(1,2,i), key2(1,2)))
      if (degree_x2 <= 4) then
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
      if (degree_x2 <= 4) then
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
      if (degree_x2 <= 4) then
        idx(l) = i
        l = l+1
      endif
    enddo

  else


    !DIR$ LOOP COUNT (1000)
    outer: do i=1,sze
      degree_x2 = 0
      !DIR$ LOOP COUNT MIN(4)
      do m=1,Nint
        if ( key1(m,1,i) /=  key2(m,1)) then
          degree_x2 = degree_x2+ popcnt(xor( key1(m,1,i), key2(m,1)))
        endif
        if ( key1(m,2,i) /=  key2(m,2)) then
          degree_x2 = degree_x2+ popcnt(xor( key1(m,2,i), key2(m,2)))
        endif
        if (degree_x2 > 4) then
          cycle outer
        endif
      enddo
      idx(l) = i
      l = l+1
    enddo outer

  endif
  idx(0) = l-1
end


