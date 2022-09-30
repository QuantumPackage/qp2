
! ---

program test_cosgtos

  implicit none
  integer :: iao, jao, kao, lao

  call init_expo()

!  call test_coef()
  call test_2e()

  iao = 1
  jao = 1
  kao = 1
  lao = 21
!  call test_2e_cpx (iao, jao, kao, lao)
!  call test_2e_real(iao, jao, kao, lao)

end

! ---

subroutine init_expo()

  implicit none

  integer                       :: i, j
  double precision, allocatable :: expo_im(:,:)

  allocate(expo_im(ao_num, ao_prim_num_max))

  do j = 1, ao_prim_num_max
    do i = 1, ao_num
      ao_expoim_cosgtos(i,j) = 0.d0
    enddo
  enddo

  call ezfio_set_cosgtos_ao_int_ao_expoim_cosgtos(expo_im)

  deallocate(expo_im)

end subroutine init_expo

! ---

subroutine test_coef()

  implicit none

  integer          :: i, j
  double precision :: coef, coef_gtos, coef_cosgtos 
  double precision :: delta, accu_abs

  print*, ' check coefs' 

  accu_abs = 0.d0
  accu_abs = 0.d0
  do i = 1, ao_num
    do j = 1, ao_prim_num(i)

      coef         = ao_coef(i,j)
      coef_gtos    = 1.d0 * ao_coef_normalized_ordered_transp(j,i)
      coef_cosgtos = 2.d0 * ao_coef_norm_ord_transp_cosgtos  (j,i)

      delta = dabs(coef_gtos - coef_cosgtos)
      accu_abs += delta

      if(delta .gt. 1.d-10) then
        print*, ' problem on: '
        print*, i, j
        print*, coef_gtos, coef_cosgtos, delta
        print*, coef 
        stop
      endif

    enddo
  enddo

  print*, 'accu_abs = ', accu_abs 

end subroutine test_coef


! ---

subroutine test_2e()

  implicit none

  integer          :: iao, jao, kao, lao
  double precision :: integral_gtos, integral_cosgtos
  double precision :: delta, accu_abs

  double precision :: ao_two_e_integral, ao_two_e_integral_cosgtos

  print*, ' check integrals' 

  accu_abs = 0.d0
  accu_abs = 0.d0

 ! iao = 1
 ! jao = 1
 ! kao = 1
 ! lao = 24

  do iao = 1, ao_num ! r1
    do jao = 1, ao_num ! r2
      do kao = 1, ao_num ! r1
        do lao = 1, ao_num ! r2

          integral_gtos    = ao_two_e_integral        (iao, kao, jao, lao)
          integral_cosgtos = ao_two_e_integral_cosgtos(iao, kao, jao, lao)

          delta = dabs(integral_gtos - integral_cosgtos)
          accu_abs += delta

          if(delta .gt. 1.d-7) then
            print*, ' problem on: '
            print*, iao, jao, kao, lao
            print*, integral_gtos, integral_cosgtos, delta
            !stop
          endif

        enddo
      enddo
    enddo
  enddo

 print*,'accu_abs = ', accu_abs

end subroutine test_2e

! ---

subroutine test_2e_cpx(iao, jao, kao, lao)

  implicit none
  integer, intent(in) :: iao, jao, kao, lao
  double precision    :: integral
  double precision    :: ao_two_e_integral_cosgtos

  integral = ao_two_e_integral_cosgtos(iao, kao, jao, lao)
  print *, ' cosgtos: ', integral

end subroutine test_2e_cpx

! ---

subroutine test_2e_real(iao, jao, kao, lao)

  implicit none
  integer, intent(in) :: iao, jao, kao, lao
  double precision    :: integral
  double precision    :: ao_two_e_integral

  integral = ao_two_e_integral(iao, kao, jao, lao)
  print *, ' gtos: ', integral

end subroutine test_2e_real

! ---


