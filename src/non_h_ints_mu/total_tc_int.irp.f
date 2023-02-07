
! ---

BEGIN_PROVIDER [double precision, ao_tc_int_chemist, (ao_num, ao_num, ao_num, ao_num)]

  implicit none
  integer          :: i, j, k, l
  double precision :: wall1, wall0

  print *, ' providing ao_tc_int_chemist ...'
  call wall_time(wall0)
  
  if(test_cycle_tc)then
   ao_tc_int_chemist = ao_tc_int_chemist_test
  else
   do j = 1, ao_num
     do l = 1, ao_num
       do i = 1, ao_num
         do k = 1, ao_num
           ao_tc_int_chemist(k,i,l,j) = tc_grad_square_ao(k,i,l,j) + tc_grad_and_lapl_ao(k,i,l,j) + ao_two_e_coul(k,i,l,j)
         enddo
       enddo
     enddo
   enddo
  endif

  call wall_time(wall1)
  print *, ' wall time for ao_tc_int_chemist ', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, ao_tc_int_chemist_test, (ao_num, ao_num, ao_num, ao_num)]

  implicit none
  integer          :: i, j, k, l
  double precision :: wall1, wall0

  print *, ' providing ao_tc_int_chemist_test ...'
  call wall_time(wall0)

   do j = 1, ao_num
     do l = 1, ao_num
       do i = 1, ao_num
         do k = 1, ao_num
           ao_tc_int_chemist_test(k,i,l,j) = tc_grad_square_ao_test(k,i,l,j) + tc_grad_and_lapl_ao_test(k,i,l,j) + ao_two_e_coul(k,i,l,j)
         enddo
       enddo
     enddo
   enddo

  call wall_time(wall1)
  print *, ' wall time for ao_tc_int_chemist_test ', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, ao_two_e_coul, (ao_num, ao_num, ao_num, ao_num) ]

  BEGIN_DOC
  !
  ! ao_two_e_coul(k,i,l,j) = ( k i | 1/r12 | l j ) = < l k | 1/r12 | j i > 
  !
  END_DOC

  integer                    :: i, j, k, l
  double precision           :: integral
  double precision, external :: get_ao_two_e_integral

  PROVIDE ao_integrals_map

  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num

          !  < 1:k, 2:l | 1:i, 2:j > 
          integral = get_ao_two_e_integral(i, j, k, l, ao_integrals_map)

          ao_two_e_coul(k,i,l,j) = integral
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

