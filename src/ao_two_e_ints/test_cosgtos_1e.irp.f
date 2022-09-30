
! ---

program test_cosgtos

  implicit none
  integer :: i, j

  call init_expo()

!  call test_coef()
  call test_1e_kin()
  call test_1e_coul()

  i = 1
  j = 1
!  call test_1e_coul_real(i, j)
!  call test_1e_coul_cpx (i, j)

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

subroutine test_1e_kin()

  implicit none

  integer          :: i, j
  double precision :: integral_gtos, integral_cosgtos
  double precision :: delta, accu_abs

  print*, ' check kin 1e integrals' 

  accu_abs = 0.d0
  accu_abs = 0.d0

  do j = 1, ao_num
    do i = 1, ao_num

      integral_gtos    = ao_kinetic_integrals        (i,j) 
      integral_cosgtos = ao_kinetic_integrals_cosgtos(i,j) 


      delta = dabs(integral_gtos - integral_cosgtos)
      accu_abs += delta

      if(delta .gt. 1.d-7) then
        print*, ' problem on: '
        print*, i, j
        print*, integral_gtos, integral_cosgtos, delta
        !stop
      endif

    enddo
  enddo

 print*,'accu_abs = ', accu_abs

end subroutine test_1e_kin

! ---

subroutine test_1e_coul()

  implicit none

  integer          :: i, j
  double precision :: integral_gtos, integral_cosgtos
  double precision :: delta, accu_abs

  print*, ' check Coulomb 1e integrals' 

  accu_abs = 0.d0
  accu_abs = 0.d0

  do j = 1, ao_num
    do i = 1, ao_num

      integral_gtos    = ao_integrals_n_e        (i,j) 
      integral_cosgtos = ao_integrals_n_e_cosgtos(i,j) 

      delta = dabs(integral_gtos - integral_cosgtos)
      accu_abs += delta

      if(delta .gt. 1.d-7) then
        print*, ' problem on: '
        print*, i, j
        print*, integral_gtos, integral_cosgtos, delta
        !stop
      endif

    enddo
  enddo

 print*,'accu_abs = ', accu_abs

end subroutine test_1e_coul

! ---

subroutine test_1e_coul_cpx(i, j)

  implicit none

  integer, intent(in) :: i, j
  double precision    :: integral

  integral = ao_integrals_n_e_cosgtos(i,j) 

  print*, ' cpx Coulomb 1e integrals', integral

end subroutine test_1e_coul_cpx

! ---

subroutine test_1e_coul_real(i, j)

  implicit none

  integer, intent(in) :: i, j
  double precision    :: integral

  integral = ao_integrals_n_e(i,j) 

  print*, ' real Coulomb 1e integrals', integral

end subroutine test_1e_coul_real

! ---
