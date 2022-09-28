program print_ao_2e_integrals
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer ::i,j,k,l

  provide ao_two_e_integrals_in_map
  complex*16 :: get_ao_two_e_integral_complex, tmp_cmplx
  do i=1,ao_num 
    do j=1,ao_num 
      do k=1,ao_num 
        do l=1,ao_num 
          tmp_cmplx = get_ao_two_e_integral_complex(i,j,k,l,ao_integrals_map,ao_integrals_map_2)
          if (cdabs(tmp_cmplx) .gt. 1E-10) then
            print'(4(I4),2(E23.15))',i,k,j,l,tmp_cmplx
          endif
        enddo
      enddo
    enddo
  enddo
  !print*,'map1'
  !do i=0,ao_integrals_map%map_size
  !  print*,i,ao_integrals_map%map(i)%value(:)
  !  print*,i,ao_integrals_map%map(i)%key(:)
  !enddo
  !print*,'map2'
  !do i=0,ao_integrals_map_2%map_size
  !  print*,i,ao_integrals_map_2%map(i)%value(:)
  !  print*,i,ao_integrals_map_2%map(i)%key(:)
  !enddo
end
