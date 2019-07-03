subroutine only_act_bitmask
 implicit none
  integer :: i,j,k
  do k = 1, N_generators_bitmask
   do j = 1, 6
    do i = 1, N_int
    generators_bitmask(i,1,j,k) = act_bitmask(i,1)
    generators_bitmask(i,2,j,k) = act_bitmask(i,2)
    enddo
   enddo
  enddo
  touch generators_bitmask
end

