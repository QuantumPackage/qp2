
! TODO
! remove ao_two_e_coul and use map directly

! ---

BEGIN_PROVIDER [double precision, ao_vartc_int_chemist, (ao_num, ao_num, ao_num, ao_num)]

  implicit none
  integer          :: i, j, k, l
  double precision :: wall1, wall0

  print *, ' providing ao_vartc_int_chemist ...'
  call wall_time(wall0)
  
  if(test_cycle_tc) then

    PROVIDE j1b_type
    if(j1b_type .ne. 3) then
      print*, ' TC integrals with cycle can not be used for j1b_type =', j1b_type
      stop
    endif

    do j = 1, ao_num
      do l = 1, ao_num
        do i = 1, ao_num
          do k = 1, ao_num
            ao_vartc_int_chemist(k,i,l,j) = tc_grad_square_ao_test(k,i,l,j) + ao_two_e_coul(k,i,l,j)
          enddo
        enddo
      enddo
    enddo

  else

    do j = 1, ao_num
      do l = 1, ao_num
        do i = 1, ao_num
          do k = 1, ao_num
            ao_vartc_int_chemist(k,i,l,j) = tc_grad_square_ao(k,i,l,j) + ao_two_e_coul(k,i,l,j)
          enddo
        enddo
      enddo
    enddo

  endif

  call wall_time(wall1)
  print *, ' wall time for ao_vartc_int_chemist ', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, ao_tc_int_chemist, (ao_num, ao_num, ao_num, ao_num)]

  implicit none
  integer          :: i, j, k, l
  double precision :: wall1, wall0

  print *, ' providing ao_tc_int_chemist ...'
  call wall_time(wall0)
  
  if(test_cycle_tc) then

    PROVIDE j1b_type
    if(j1b_type .ne. 3) then
      print*, ' TC integrals with cycle can not be used for j1b_type =', j1b_type
      stop
    endif

    ao_tc_int_chemist = ao_tc_int_chemist_test

  else

    PROVIDE tc_grad_square_ao tc_grad_and_lapl_ao ao_two_e_coul

    do j = 1, ao_num
      do l = 1, ao_num
        do i = 1, ao_num
          do k = 1, ao_num
            ao_tc_int_chemist(k,i,l,j) = tc_grad_square_ao(k,i,l,j) + tc_grad_and_lapl_ao(k,i,l,j) + ao_two_e_coul(k,i,l,j)
!            ao_tc_int_chemist(k,i,l,j) = ao_two_e_coul(k,i,l,j)
          enddo
        enddo
      enddo
    enddo
  endif

  FREE tc_grad_square_ao tc_grad_and_lapl_ao ao_two_e_coul

  call wall_time(wall1)
  print *, ' wall time for ao_tc_int_chemist ', wall1 - wall0
  call print_memory_usage()

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, ao_tc_int_chemist_no_cycle, (ao_num, ao_num, ao_num, ao_num)]

  implicit none
  integer          :: i, j, k, l
  double precision :: wall1, wall0

  print *, ' providing ao_tc_int_chemist_no_cycle ...'
  call wall_time(wall0)

  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          ao_tc_int_chemist_no_cycle(k,i,l,j) = tc_grad_square_ao(k,i,l,j) + tc_grad_and_lapl_ao(k,i,l,j) + ao_two_e_coul(k,i,l,j)
          !ao_tc_int_chemist(k,i,l,j) = ao_two_e_coul(k,i,l,j)
        enddo
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print *, ' wall time for ao_tc_int_chemist_no_cycle ', wall1 - wall0

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
!           ao_tc_int_chemist_test(k,i,l,j) =  ao_two_e_coul(k,i,l,j)
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

  integer                       :: i, j, k, l
  double precision              :: integral
  double precision, allocatable :: tmp(:)
  double precision, external    :: get_ao_two_e_integral

  PROVIDE ao_integrals_map

  !$OMP PARALLEL DEFAULT(NONE)                          &
  !$OMP SHARED(ao_num, ao_two_e_coul, ao_integrals_map) &
  !$OMP PRIVATE(i, j, k, l)
  !$OMP DO
  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          !  < 1:k, 2:l | 1:i, 2:j > 
          ao_two_e_coul(k,i,l,j) = get_ao_two_e_integral(i, j, k, l, ao_integrals_map)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


! TODO
!  allocate(tmp(ao_num))
!
!  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,l,j,k,tmp)
!  do j = 1, ao_num
!    do l = 1, ao_num
!      do i = 1, ao_num
!        call get_ao_two_e_integrals(i, l, l, ao_num, tmp(1))
!        do k = 1, ao_num
!          ao_two_e_coul(k,i,l,j) = tmp(k)
!        enddo
!      enddo
!    enddo
!  enddo
!  !$OMP END PARALLEL DO
!
!  deallocate(tmp)

END_PROVIDER 

! ---

