program print_bitmasks
  call routine
end

subroutine routine
 implicit none
 integer :: i
 do i = 1, mo_num
  print*,i ,'',mo_class(i)
 enddo
 print*, 'core'
 do i = 1, n_core_orb
  print*, list_core(i)
 enddo
 call debug_det(core_bitmask, N_int)
 print*, 'inact'
 do i = 1, n_inact_orb
  print*, list_inact(i)
 enddo
 call debug_det(inact_bitmask, N_int)
 print*, 'act'
 do i = 1, n_act_orb
  print*, list_act(i)
 enddo
 call debug_det(act_bitmask, N_int)
 print*, 'virt'
 do i = 1, n_virt_orb
  print*, list_virt(i)
 enddo
 call debug_det(virt_bitmask, N_int)

end

