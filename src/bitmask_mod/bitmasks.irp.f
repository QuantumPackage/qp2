use bitmasks

BEGIN_PROVIDER [ integer, N_int ]
  implicit none
  include 'utils/constants.include.F'
  BEGIN_DOC
  ! Number of 64-bit integers needed to represent determinants as binary strings
  END_DOC
  N_int = (ao_num-1)/bit_kind_size + 1
  call write_int(6,N_int, 'N_int')
  if (N_int > N_int_max) then
    stop 'N_int > N_int_max'
  endif
  
END_PROVIDER


BEGIN_PROVIDER [ integer(bit_kind), full_ijkl_bitmask, (N_int) ]
  implicit none
  BEGIN_DOC
  ! Bitmask to include all possible AOs
  END_DOC
  
  integer                        :: i,j,k
  k=0
  do j=1,N_int
    full_ijkl_bitmask(j) = 0_bit_kind
    do i=0,bit_kind_size-1
      k=k+1
      full_ijkl_bitmask(j) = ibset(full_ijkl_bitmask(j),i)
      if (k == ao_num) exit
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), full_ijkl_bitmask_4, (N_int,4) ]
  implicit none
  integer                        :: i
  do i=1,N_int
    full_ijkl_bitmask_4(i,1) = full_ijkl_bitmask(i)
    full_ijkl_bitmask_4(i,2) = full_ijkl_bitmask(i)
    full_ijkl_bitmask_4(i,3) = full_ijkl_bitmask(i)
    full_ijkl_bitmask_4(i,4) = full_ijkl_bitmask(i)
  enddo
END_PROVIDER

