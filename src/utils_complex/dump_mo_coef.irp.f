program print_mo_coef
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer ::i,j,k,l

  provide mo_coef_complex
  complex*16 :: tmp_cmplx
!  complex*16 :: get_two_e_integral_complex, tmp_cmplx
  do i=1,ao_num
    do j=1,mo_num
      tmp_cmplx = mo_coef_complex(i,j)
      if (cdabs(tmp_cmplx).gt.1.d-10) then
        print'(2(I4),2(E23.15))',i,j,tmp_cmplx
      endif
    enddo
  enddo
end
