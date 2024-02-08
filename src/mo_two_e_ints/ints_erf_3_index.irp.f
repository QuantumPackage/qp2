 BEGIN_PROVIDER [double precision, int_erf_3_index, (mo_num,mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, int_erf_3_index_exc,(mo_num,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! int_erf_3_index(i,j)     = <ij|ij> = (ii|jj) with the erf interaction
 !
 ! int_erf_3_index_exc(i,j) = <ij|ji> = (ij|ij) with the erf interaction
 END_DOC
 integer :: i,j,k,l
 double precision :: get_mo_two_e_integral_erf
 double precision :: integral

 do k = 1, mo_num
  do i = 1, mo_num
   do j = 1, mo_num
     l = j
     integral = get_mo_two_e_integral_erf(i,j,k,l,mo_integrals_erf_map)
     int_erf_3_index(j,i,k) = integral
     l = j
     integral = get_mo_two_e_integral_erf(i,j,l,k,mo_integrals_erf_map)
     int_erf_3_index_exc(j,i,k) = integral
   enddo
  enddo
 enddo


END_PROVIDER

