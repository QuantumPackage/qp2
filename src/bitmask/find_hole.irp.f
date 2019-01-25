logical function is_the_hole_in_det(key_in,ispin,i_hole)
  use bitmasks
  ! returns true if the electron ispin is absent from i_hole
 implicit none
 integer, intent(in) :: i_hole,ispin
 integer(bit_kind), intent(in) :: key_in(N_int,2)
 integer(bit_kind) :: key_tmp(N_int)
 integer(bit_kind) :: itest(N_int)
 integer :: i,j,k
 do i = 1, N_int
  itest(i) = 0_bit_kind
 enddo
 k = shiftr(i_hole-1,bit_kind_shift)+1
 j = i_hole-shiftl(k-1,bit_kind_shift)-1
 itest(k) = ibset(itest(k),j)
 j = 0
 do i = 1, N_int
  key_tmp(i) = iand(itest(i),key_in(i,ispin))
  j += popcnt(key_tmp(i))
 enddo
 if(j==0)then
  is_the_hole_in_det = .True.
 else
  is_the_hole_in_det = .False.
 endif

end

logical function is_the_particl_in_det(key_in,ispin,i_particl)
  use bitmasks
  ! returns true if the electron ispin is absent from i_particl
 implicit none
 integer, intent(in) :: i_particl,ispin
 integer(bit_kind), intent(in) :: key_in(N_int,2)
 integer(bit_kind) :: key_tmp(N_int)
 integer(bit_kind) :: itest(N_int)
 integer :: i,j,k
 do i = 1, N_int
  itest(i) = 0_bit_kind
 enddo
 k = shiftr(i_particl-1,bit_kind_shift)+1
 j = i_particl-shiftl(k-1,bit_kind_shift)-1
 itest(k) = ibset(itest(k),j)
 j = 0
 do i = 1, N_int
  key_tmp(i) = iand(itest(i),key_in(i,ispin))
  j += popcnt(key_tmp(i))
 enddo
 if(j==0)then
  is_the_particl_in_det = .False.
 else
  is_the_particl_in_det = .True.
 endif

end
