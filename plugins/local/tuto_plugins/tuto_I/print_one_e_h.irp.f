program my_program_to_print_stuffs
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
 integer :: i,j
 print*,'AO integrals '
 do i = 1, ao_num
  do j = 1, ao_num
   print*,j,i,ao_one_e_integrals(j,i)
  enddo
 enddo

 print*,'MO integrals '
 do i = 1, mo_num
  do j = 1, mo_num
   print*,j,i,mo_one_e_integrals(j,i)
  enddo
 enddo
end
