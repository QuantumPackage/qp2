 BEGIN_PROVIDER [double precision, big_array_coulomb_integrals, (mo_num,mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, big_array_exchange_integrals,(mo_num,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! big_array_coulomb_integrals(j,i,k)  = <ij|kj> = (ik|jj)
 !
 ! big_array_exchange_integrals(i,j,k) = <ij|jk> = (ij|kj)
 END_DOC
 integer :: i,j,k,l
 double precision :: get_two_e_integral
 double precision :: integral

 do k = 1, mo_num
  do i = 1, mo_num
   do j = 1, mo_num
     l = j
     integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
     big_array_coulomb_integrals(j,i,k) = integral
     l = j
     integral = get_two_e_integral(i,j,l,k,mo_integrals_map)
     big_array_exchange_integrals(j,i,k) = integral
   enddo
  enddo
 enddo

END_PROVIDER

