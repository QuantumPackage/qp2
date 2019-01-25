logical function det_inf(key1, key2, Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
! Ordering function for determinants.
  END_DOC
  integer,intent(in)                 :: Nint
  integer(bit_kind),intent(in)       :: key1(Nint, 2), key2(Nint, 2)
  integer                            :: i,j

  det_inf = .false.

  do i=1,2
    do j=Nint,1,-1
      if(key1(j,i) < key2(j,i)) then
        det_inf = .true.
        return
      else if(key1(j,i) > key2(j,i)) then
        return
      end if
    end do
  end do
end function


subroutine tamiser(key, idx, no, n, Nint, N_key)
  use bitmasks
  implicit none
  integer,intent(in)                    :: no, n, Nint, N_key
  integer(bit_kind),intent(inout)       :: key(Nint, 2, N_key)
  integer,intent(inout)                 :: idx(N_key)
  integer                               :: k,j,tmpidx
  integer(bit_kind)                     :: tmp(Nint, 2)
  logical                               :: det_inf
  integer                               :: ni

  k = no
  j = 2*k
  do while(j <= n)
    if(j < n) then
      if (det_inf(key(1,1,j), key(1,1,j+1), Nint)) then
        j = j+1
      endif
    endif
    if(det_inf(key(1,1,k), key(1,1,j), Nint)) then
      do ni=1,Nint
        tmp(ni,1)   = key(ni,1,k)
        tmp(ni,2)   = key(ni,2,k)
        key(ni,1,k) = key(ni,1,j)
        key(ni,2,k) = key(ni,2,j)
        key(ni,1,j) = tmp(ni,1)
        key(ni,2,j) = tmp(ni,2)
      enddo
      tmpidx = idx(k)
      idx(k) = idx(j)
      idx(j) = tmpidx
      k = j
      j = k+k
    else
      return
    endif
  enddo
end subroutine


subroutine sort_dets_ba_v(key_in, key_out, idx, shortcut, version, N_key, Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
! Deprecated routine
  END_DOC
  integer, intent(in)            :: Nint, N_key
  integer(bit_kind),intent(in)   :: key_in(Nint,2,N_key)
  integer(bit_kind),intent(out)  :: key_out(Nint,N_key)
  integer,intent(out)            :: idx(N_key)
  integer,intent(out)            :: shortcut(0:N_key+1)
  integer(bit_kind),intent(out)  :: version(Nint,N_key+1)
  integer(bit_kind), allocatable :: key(:,:,:)
  integer                        :: i,ni

  allocate ( key(Nint,2,N_key) )
  do i=1,N_key
    do ni=1,Nint
      key(ni,1,i) = key_in(ni,2,i)
      key(ni,2,i) = key_in(ni,1,i)
    enddo
  enddo

  call sort_dets_ab_v(key, key_out, idx, shortcut, version, N_key, Nint)
  deallocate ( key )
end subroutine



subroutine sort_dets_ab_v(key_in, key_out, idx, shortcut, version, N_key, Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
! Deprecated routine
  END_DOC
  integer, intent(in)                   :: Nint, N_key
  integer(bit_kind),intent(in)          :: key_in(Nint,2,N_key)
  integer(bit_kind),intent(out)         :: key_out(Nint,N_key)
  integer,intent(out)                   :: idx(N_key)
  integer,intent(out)                   :: shortcut(0:N_key+1)
  integer(bit_kind),intent(out)         :: version(Nint,N_key+1)
  integer(bit_kind), allocatable        :: key(:,:,:)
  integer(bit_kind)                     :: tmp(Nint, 2)
  integer                               :: tmpidx,i,ni

  allocate (key(Nint,2,N_key))
  do i=1,N_key
    do ni=1,Nint
      key(ni,1,i) = key_in(ni,1,i)
      key(ni,2,i) = key_in(ni,2,i)
    enddo
    idx(i) = i
  end do

  do i=N_key/2,1,-1
    call tamiser(key, idx, i, N_key, Nint, N_key)
  end do

  do i=N_key,2,-1
    do ni=1,Nint
      tmp(ni,1) = key(ni,1,i)
      tmp(ni,2) = key(ni,2,i)
      key(ni,1,i) = key(ni,1,1)
      key(ni,2,i) = key(ni,2,1)
      key(ni,1,1) = tmp(ni,1)
      key(ni,2,1) = tmp(ni,2)
    enddo
    tmpidx = idx(i)
    idx(i) = idx(1)
    idx(1) = tmpidx
    call tamiser(key, idx, 1, i-1, Nint, N_key)
  end do

  shortcut(0) = 1
  shortcut(1) = 1
  do ni=1,Nint
    version(ni,1) = key(ni,1,1)
  enddo
  do i=2,N_key
    do ni=1,nint
      if(key(ni,1,i) /= key(ni,1,i-1)) then
        shortcut(0) = shortcut(0) + 1
        shortcut(shortcut(0)) = i
        version(:,shortcut(0)) = key(:,1,i)
        exit
      end if
    end do
  end do
  shortcut(shortcut(0)+1) = N_key+1
  do i=1,N_key
    do ni=1,Nint
      key_out(ni,i) = key(ni,2,i)
    enddo
  enddo
  deallocate (key)
end subroutine


subroutine sort_dets_ab(key, idx, shortcut, N_key, Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
! Deprecated routine
  END_DOC
  integer, intent(in)                   :: Nint, N_key
  integer(bit_kind),intent(inout)       :: key(Nint,2,N_key)
  integer,intent(inout)                   :: idx(N_key)
  integer,intent(inout)                   :: shortcut(0:N_key+1)
  integer(bit_kind)                     :: tmp(Nint, 2)
  integer                               :: tmpidx,i,ni

  do i=1,N_key
    idx(i) = i
  end do

  do i=N_key/2,1,-1
    call tamiser(key, idx, i, N_key, Nint, N_key)
  end do

  do i=N_key,2,-1
    do ni=1,Nint
      tmp(ni,1) = key(ni,1,i)
      tmp(ni,2) = key(ni,2,i)
      key(ni,1,i) = key(ni,1,1)
      key(ni,2,i) = key(ni,2,1)
      key(ni,1,1) = tmp(ni,1)
      key(ni,2,1) = tmp(ni,2)
    enddo

    tmpidx = idx(i)
    idx(i) = idx(1)
    idx(1) = tmpidx
    call tamiser(key, idx, 1, i-1, Nint, N_key)
  end do

  shortcut(0) = 1
  shortcut(1) = 1
  do i=2,N_key
    do ni=1,nint
      if(key(ni,1,i) /= key(ni,1,i-1)) then
        shortcut(0) = shortcut(0) + 1
        shortcut(shortcut(0)) = i
        exit
      end if
    end do
  end do
  shortcut(shortcut(0)+1) = N_key+1
end subroutine


