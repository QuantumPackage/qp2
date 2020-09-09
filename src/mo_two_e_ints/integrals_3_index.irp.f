 BEGIN_PROVIDER [double precision, big_array_coulomb_integrals, (mo_num,mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, big_array_exchange_integrals,(mo_num,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! big_array_coulomb_integrals(i,j)  = <ij|ij> = (ii|jj)
 !
 ! big_array_exchange_integrals(i,j) = <ij|ji> = (ij|ij)
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

 BEGIN_PROVIDER [complex*16, big_array_coulomb_integrals_complex, (mo_num,mo_num, mo_num)]
&BEGIN_PROVIDER [complex*16, big_array_exchange_integrals_complex,(mo_num,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! big_array_coulomb_integrals(j,i,k)  = <ij|kj> = (ik|jj)
 ! big_array_exchange_integrals(j,i,k) = <ij|jk> = (ij|jk)
 !    for both of these, i and k must be from same kpt for integral to be nonzero
 ! TODO: only loop over half, and assign two elements:
 !        b_a_coul_int(j,i,k) = b_a_coul_int(j,k,i)*
 !        b_a_exch_int(j,i,k) = b_a_exch_int(j,k,i)*
 END_DOC
 integer :: i,j,k,l
 complex*16 :: get_two_e_integral_complex
 complex*16 :: integral

 do k = 1, mo_num
  do i = 1, mo_num
   do j = 1, mo_num
     l = j
     integral = get_two_e_integral_complex(i,j,k,l,mo_integrals_map,mo_integrals_map_2)
     big_array_coulomb_integrals_complex(j,i,k) = integral
     l = j
     integral = get_two_e_integral_complex(i,j,l,k,mo_integrals_map,mo_integrals_map_2)
     big_array_exchange_integrals_complex(j,i,k) = integral
   enddo
  enddo
 enddo

END_PROVIDER

 BEGIN_PROVIDER [complex*16, big_array_coulomb_integrals_kpts, (mo_num_per_kpt,kpt_num,mo_num_per_kpt, mo_num_per_kpt,kpt_num)]
&BEGIN_PROVIDER [complex*16, big_array_exchange_integrals_kpts,(mo_num_per_kpt,kpt_num,mo_num_per_kpt, mo_num_per_kpt,kpt_num)]
 implicit none
 BEGIN_DOC
 ! big_array_coulomb_integrals(j,kj,i,k,ki)  = <ij|kj> = (ik|jj)
 ! big_array_exchange_integrals(j,kj,i,k,ki) = <ij|jk> = (ij|jk)
 !    for both of these, i and k must be from same kpt for integral to be nonzero
 ! TODO: only loop over half, and assign two elements:
 !        b_a_coul_int(j,i,k) = b_a_coul_int(j,k,i)*
 !        b_a_exch_int(j,i,k) = b_a_exch_int(j,k,i)*
 END_DOC
 integer :: i,j,k,l
 integer :: ki,kj,kk,kl
 complex*16 :: get_two_e_integral_kpts
 complex*16 :: integral
 
 do ki = 1,kpt_num
  kk=ki
  do k = 1, mo_num_per_kpt
   do i = 1, mo_num_per_kpt
    do kj = 1,kpt_num
     kl=kj
     do j = 1, mo_num_per_kpt
       l = j
       integral = get_two_e_integral_kpts(i,j,k,l,ki,kj,kk,kl,mo_integrals_map,mo_integrals_map_2)
       big_array_coulomb_integrals_kpts(j,kj,i,k,ki) = integral
       l = j
       integral = get_two_e_integral_kpts(i,j,l,k,ki,kj,kl,kk,mo_integrals_map,mo_integrals_map_2)
       big_array_exchange_integrals_kpts(j,kj,i,k,ki) = integral
     enddo
    enddo
   enddo
  enddo
 enddo

END_PROVIDER

