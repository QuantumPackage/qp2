program my_program_to_print_stuffs
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
 integer :: i,j,k,l,m
 double precision :: integral, accu, accu_tot, integral_cholesky
 double precision :: get_ao_two_e_integral, get_two_e_integral ! declaration of the functions 
 print*,'AO integrals, physicist notations : <i j|k l>'
 accu_tot = 0.D0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     integral = get_ao_two_e_integral(i, j, k, l, ao_integrals_map)
     integral_cholesky = 0.D0
     do m = 1, cholesky_ao_num
      integral_cholesky += cholesky_ao_transp(m,i,k) * cholesky_ao_transp(m,j,l)
     enddo
     accu = dabs(integral_cholesky-integral)
     accu_tot += accu
     if(accu.gt.1.d-10)then
      print*,i,j,k,l
      print*,accu, integral, integral_cholesky
     endif
    enddo
   enddo
  enddo
 enddo
 print*,'accu_tot',accu_tot

 print*,'MO integrals, physicist notations : <i j|k l>'
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     integral = get_two_e_integral(i, j, k, l, mo_integrals_map)
     accu = 0.D0
     integral_cholesky = 0.D0
     do m = 1, cholesky_mo_num
      integral_cholesky += cholesky_mo_transp(m,i,k) * cholesky_mo_transp(m,j,l)
     enddo
     accu = dabs(integral_cholesky-integral)
     accu_tot += accu
     if(accu.gt.1.d-10)then
      print*,i,j,k,l
      print*,accu, integral, integral_cholesky
     endif
    enddo
   enddo
  enddo
 enddo
end
