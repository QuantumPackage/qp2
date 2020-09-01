program print_mo_2e_integrals
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer ::i,j,k,l

  provide mo_two_e_integrals_in_map
  complex*16 :: get_two_e_integral_complex, tmp_cmplx
  do i=1,mo_num 
    do j=1,mo_num 
      do k=1,mo_num 
        do l=1,mo_num 
          tmp_cmplx = get_two_e_integral_complex(i,j,k,l,mo_integrals_map,mo_integrals_map_2)
          if (cdabs(tmp_cmplx).gt. 1d-12) then
            print'(4(I4),2(E23.15))',i,j,k,l,tmp_cmplx
          endif
        enddo
      enddo
    enddo
  enddo
!  print*,'map1'
!  do i=0,mo_integrals_map%map_size
!    print*,i,mo_integrals_map%map(i)%value(:)
!    print*,i,mo_integrals_map%map(i)%key(:)
!  enddo
!  print*,'map2'
!  do i=0,mo_integrals_map_2%map_size
!    print*,i,mo_integrals_map_2%map(i)%value(:)
!    print*,i,mo_integrals_map_2%map(i)%key(:)
!  enddo
end
