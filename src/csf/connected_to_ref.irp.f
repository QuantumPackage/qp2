integer*8 function configuration_search_key(cfg,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns an integer*8 corresponding to a determinant index for searching.
  ! The left-most 8 bits contain the number of open shells+1. This ensures that the CSF
  ! are packed with the same seniority.
  END_DOC
  integer, intent(in) :: Nint
  integer(bit_kind), intent(in) :: cfg(Nint,2)
  integer :: i, n_open_shells
  integer*8 :: mask

  i = shiftr(elec_alpha_num, bit_kind_shift)+1
  configuration_search_key = int(shiftr(ior(cfg(i,1),cfg(i,2)),1)+sum(cfg),8)

  mask = int(Z'00FFFFFFFFFFFFFF',8)
  configuration_search_key = iand(mask,configuration_search_key)

  n_open_shells = 1
  do i=1,Nint
    if (cfg(i,1) == 0_bit_kind) cycle
    n_open_shells = n_open_shells + popcnt(cfg(i,1))
  enddo
  mask = n_open_shells
  mask = shiftl(mask,56)
  configuration_search_key = ior (mask,configuration_search_key)

end
