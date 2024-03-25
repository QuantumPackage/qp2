program my_program_to_print_stuffs
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
 integer :: i,j,k,l
 double precision :: integral 
 double precision :: get_ao_two_e_integral, get_two_e_integral ! declaration of the functions 
 print*,'AO integrals, physicist notations : <i j|k l>'
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     integral = get_ao_two_e_integral(i, j, k, l, ao_integrals_map)
     print*,i,j,k,l,integral 
    enddo
   enddo
  enddo
 enddo

 print*,'MO integrals, physicist notations : <i j|k l>'
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     integral = get_two_e_integral(i, j, k, l, mo_integrals_map)
     print*,i,j,k,l,integral 
    enddo
   enddo
  enddo
 enddo
end
