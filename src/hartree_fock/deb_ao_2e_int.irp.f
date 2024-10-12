
program deb_ao_2e_int

  implicit none

  call check_ao_two_e_integral_cosgtos()
  !call check_crint1()
  !call check_crint2()
  !call check_crint3()

end

! ---

subroutine check_ao_two_e_integral_cosgtos()

  implicit none

  integer                    :: i, j, k, l
  double precision           :: acc, nrm, dif
  double precision           :: tmp1, tmp2

  double precision, external :: ao_two_e_integral
  double precision, external :: ao_two_e_integral_cosgtos

  acc = 0.d0
  nrm = 0.d0

  i = 11
  j = 100
  k = 74
  l = 104
 ! do i = 1, ao_num
 !   do k = 1, ao_num
 !     do j = 1, ao_num
 !       do l = 1, ao_num

          tmp1 = ao_two_e_integral        (i, j, k, l)
          tmp2 = ao_two_e_integral_cosgtos(i, j, k, l)

          dif = abs(tmp1 - tmp2)
          !if(dif .gt. 1d-10) then
            print*, ' error on:', i, j, k, l
            print*, tmp1, tmp2, dif
            !stop
          !endif

          acc += dif
          nrm += abs(tmp1)
 !       enddo
 !     enddo
 !   enddo
 ! enddo

  print *, ' acc (%) = ', dif * 100.d0 / nrm

end

! ---

subroutine check_crint1()

  implicit none
  integer                    :: i, n, i_rho
  double precision           :: dif_thr
  double precision           :: dif_re, dif_im, acc_re, nrm_re, acc_im, nrm_im
  complex*16                 :: rho_test(1:10) = (/ (1d-12, 0.d0),  &
                                                    (+1d-9, +1d-6), &
                                                    (-1d-6, -1d-5), &
                                                    (+1d-3, -1d-2), &
                                                    (-1d-1, +1d-1), &
                                                    (+1d-0, +1d-1), &
                                                    (-1d+1, +1d+1), &
                                                    (+1d+2, +1d+1), &
                                                    (-1d+3, +1d+2), &
                                                    (+1d+4, +1d+4) /)
  complex*16                 :: rho
  complex*16                 :: int_an, int_nm

  double precision, external :: rint
  complex*16, external       :: crint_1, crint_2

  n = 10
  dif_thr = 1d-7

  do i_rho = 8, 10
  !do i_rho = 7, 7
  
    !rho = (-10.d0, 0.1d0)
    !rho = (+10.d0, 0.1d0)
    rho = rho_test(i_rho)
    print*, "rho = ", real(rho), aimag(rho)

    acc_re = 0.d0
    nrm_re = 0.d0
    acc_im = 0.d0
    nrm_im = 0.d0
    do i = 0, n
      !int_an = crint_1(i, rho)
      int_an = crint_2(i, rho)
      call crint_quad_1(i, rho, 100000000, int_nm)
      
      dif_re = dabs(real(int_an) - real(int_nm))
      dif_im = dabs(aimag(int_an) - aimag(int_nm))

      if((dif_re .gt. dif_thr) .or. (dif_im .gt. dif_thr)) then
        print*, ' error on i =', i
        print*, real(int_an), real(int_nm), dif_re
        print*, aimag(int_an), aimag(int_nm), dif_im
        !print*, rint(i, real(rho))
        print*, crint_1(i, rho)
        !print*, crint_2(i, rho)
        stop
      endif
      acc_re += dif_re
      nrm_re += dabs(real(int_nm))
      acc_im += dif_im
      nrm_im += dabs(aimag(int_nm))
    enddo
  
    print*, "accuracy on real part (%):", 100.d0 * acc_re / (nrm_re+1d-15)
    print*, "accuracy on imag part (%):", 100.d0 * acc_im / (nrm_im+1d-15)
  enddo

end

! ---

subroutine check_crint2()

  implicit none

  integer              :: i, n, i_rho
  double precision     :: dif_thr
  double precision     :: dif_re, dif_im, acc_re, nrm_re, acc_im, nrm_im
  complex*16           :: rho_test(1:10) = (/ (1d-12, 0.d0),  &
                                              (+1d-9, +1d-6), &
                                              (-1d-6, -1d-5), &
                                              (+1d-3, -1d-2), &
                                              (-1d-1, +1d-1), &
                                              (+1d-0, +1d-1), &
                                              (-1d+1, +1d+1), &
                                              (+1d+2, +1d+1), &
                                              (-1d+3, +1d+2), &
                                              (+1d+4, +1d+4) /)
  complex*16           :: rho
  complex*16           :: int_an, int_nm
  complex*16, external :: crint_1, crint_2

  n = 30
  dif_thr = 1d-12

  do i_rho = 1, 10
    rho = rho_test(i_rho)
    print*, "rho = ", real(rho), aimag(rho)

    acc_re = 0.d0
    nrm_re = 0.d0
    acc_im = 0.d0
    nrm_im = 0.d0
    do i = 0, n
      int_an = crint_1(i, rho)
      int_nm = crint_2(i, rho)
      
      dif_re = dabs(real(int_an) - real(int_nm))
      !if(dif_re .gt. dif_thr) then
      !  print*, ' error in real part:', i
      !  print*, real(int_an), real(int_nm), dif_re
      !  stop
      !endif
      acc_re += dif_re
      nrm_re += dabs(real(int_nm))
  
      dif_im = dabs(aimag(int_an) - aimag(int_nm))
      !if(dif_im .gt. dif_thr) then
      !  print*, ' error in imag part:', i
      !  print*, aimag(int_an), aimag(int_nm), dif_im
      !  stop
      !endif
      acc_im += dif_im
      nrm_im += dabs(aimag(int_nm))
    enddo
  
    print*, "accuracy on real part (%):", 100.d0 * acc_re / (nrm_re+1d-15)
    print*, "accuracy on imag part (%):", 100.d0 * acc_im / (nrm_im+1d-15)
  enddo

end

! ---

subroutine check_crint3()

  implicit none

  integer                    :: i_test, n_test
  integer                    :: nx, ny, n, n_quad
  integer                    :: i, seed_size, clock_time
  double precision           :: xr(1:4), x
  double precision           :: yr(1:4), y
  double precision           :: dif_re, dif_im, acc_re, nrm_re, acc_im, nrm_im
  double precision           :: delta_ref
  double precision           :: t1, t2, t_int1, t_int2
  complex*16                 :: rho
  complex*16                 :: int1_old, int1_ref, int2_old, int2_ref
  integer, allocatable       :: seed(:)

  complex*16, external       :: crint_2

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock_time)                  
  seed = clock_time + 37 * (/ (i, i=0, seed_size-1) /) 
  !seed = 123456789
  call random_seed(put=seed)


  t_int1 = 0.d0
  t_int2 = 0.d0

  n_test = 5

  acc_re = 0.d0
  nrm_re = 0.d0
  acc_im = 0.d0
  nrm_im = 0.d0
  do i_test = 1, n_test

    ! Re(rho)
    call random_number(xr)
    x = xr(1)
    if(xr(2) .gt. 0.5d0) x = -x
    nx = int(15.d0 * xr(3))
    if(xr(4) .gt. 0.5d0) nx = -nx
    x = x * 10.d0**nx

    ! Im(rho)
    call random_number(yr)
    y = yr(1)
    if(yr(2) .gt. 0.5d0) y = -y
    ny = int(5.d0 * yr(3))
    if(yr(4) .gt. 0.5d0) ny = -ny
    y = y * 10.d0**ny

    rho = x + (0.d0, 1.d0) * y

    call random_number(x)
    x = 31.d0 * x
    n = int(x)
    !if(n.eq.0) cycle

    print*, " n = ", n
    print*, " rho = ", real(rho), aimag(rho)

    call wall_time(t1)
    int1_old = crint_2(n, rho)
    !n_quad = 10000000
    !call crint_quad_1(n, rho, n_quad, int1_old)
    !!delta_ref = 1.d0
    !!do while(delta_ref .gt. 1d-12)
    !!  n_quad = n_quad * 10
    !!  !print*, " delta = ", delta_ref
    !!  !print*, " increasing n_quad to:", n_quad
    !!  call crint_quad_1(n, rho, n_quad, int1_ref)
    !!  delta_ref = abs(int1_ref - int1_old)
    !!  int1_old = int1_ref
    !!  if(n_quad .ge. 1000000000) then
    !!    print*, ' convergence was not reached for crint_quad_1'
    !!    print*, " delta = ", delta_ref
    !!    exit
    !!  endif
    !!enddo
    call wall_time(t2)
    t_int1 = t_int1 + t2 - t1
    !print*, " n_quad for crint_quad_1:", n_quad

    call wall_time(t1)
    n_quad = 10000000
    call crint_quad_12(n, rho, n_quad, int2_old)
    !delta_ref = 1.d0
    !do while(delta_ref .gt. 1d-12)
    !  n_quad = n_quad * 10
    !  !print*, " delta = ", delta_ref
    !  !print*, " increasing n_quad to:", n_quad
    !  call crint_quad_12(n, rho, n_quad, int2_ref)
    !  delta_ref = abs(int2_ref - int2_old)
    !  int2_old = int2_ref
    !  if(n_quad .ge. 1000000000) then
    !    print*, ' convergence was not reached for crint_quad_2'
    !    print*, " delta = ", delta_ref
    !    exit
    !  endif
    !enddo
    call wall_time(t2)
    t_int2 = t_int2 + t2 - t1
    !print*, " n_quad for crint_quad_2:", n_quad

    dif_re = dabs(real(int1_old) - real(int2_old))
    dif_im = dabs(aimag(int1_old) - aimag(int2_old))
    if((dif_re .gt. 1d-10) .or. (dif_im .gt. 1d-10)) then
      print*, ' important error found: '
      print*, " n = ", n
      print*, " rho = ", real(rho), aimag(rho)
      print*, real(int1_old), real(int2_old), dif_re
      print*, aimag(int1_old), aimag(int2_old), dif_im
      !stop
    endif

    if((real(int1_old) /= real(int1_old)) .or. (aimag(int1_old) /= aimag(int1_old)) .or. &
       (real(int2_old) /= real(int2_old)) .or. (aimag(int2_old) /= aimag(int2_old)) ) then
      cycle
    else
      acc_re += dif_re
      acc_im += dif_im
      nrm_re += dabs(real(int1_old))
      nrm_im += dabs(aimag(int1_old))
    endif
  enddo

  print*, "accuracy on real part (%):", 100.d0 * acc_re / (nrm_re + 1d-15)
  print*, "accuracy on imag part (%):", 100.d0 * acc_im / (nrm_im + 1d-15)

  print*, "crint_quad_1 wall time (sec) = ", t_int1
  print*, "crint_quad_2 wall time (sec) = ", t_int2


  deallocate(seed)

end

! ---

