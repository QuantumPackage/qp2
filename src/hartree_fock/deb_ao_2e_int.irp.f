
program deb_ao_2e_int

  implicit none

  !call main()
  call check_ao_one_e_integral_cgtos()
  !call check_ao_two_e_integral_cgtos()
  !call check_crint1()
  !call check_crint2()
  !call check_crint3()
  !call check_crint4()
  !call check_crint5()
  !call check_crint6()

end

! ---

subroutine main()

  implicit none

  integer :: i, j


  PROVIDE ao_overlap
  PROVIDE ao_kinetic_integrals
  PROVIDE ao_integrals_n_e

  print*, "ao_overlap:"
  do i = 1, ao_num
    print*, (ao_overlap(i,j), j=1, ao_num)
  enddo

  print*, "ao_kinetic_integrals:"
  do i = 1, ao_num
    print*, (ao_kinetic_integrals(i,j), j=1, ao_num)
  enddo

  print*, "ao_integrals_n_e:"
  do i = 1, ao_num
    print*, (ao_integrals_n_e(i,j), j=1, ao_num)
  enddo

  return
end

! ---

subroutine check_ao_one_e_integral_cgtos()

  implicit none

  integer          :: i, j
  double precision :: acc, nrm, dif
  double precision :: tmp1, tmp2
  double precision :: t1, t2, tt

  PROVIDE ao_overlap ao_overlap_cgtos
  PROVIDE ao_integrals_n_e ao_integrals_n_e_cgtos
  PROVIDE ao_kinetic_integrals ao_kinetic_integrals_cgtos

  ! ---

!  print *, "overlap:"
!  acc = 0.d0
!  nrm = 0.d0
!  do i = 1, ao_num
!    do j = 1, ao_num
!      tmp1 = ao_overlap      (i,j)
!      tmp2 = ao_overlap_cgtos(i,j)
!      dif = abs(tmp1 - tmp2)
!      if(dif .gt. 1d-10) then
!        print*, ' error on:', i, j
!        print*, tmp1, tmp2, dif
!        !stop
!      endif
!      acc += dif
!      nrm += abs(tmp1)
!    enddo
!  enddo
!  print *, ' acc (%) = ', 100.d0 * acc / nrm
!
!  ! ---
!
!  print *, "kinetic:"
!  acc = 0.d0
!  nrm = 0.d0
!  do i = 1, ao_num
!    do j = 1, ao_num
!      tmp1 = ao_kinetic_integrals      (i,j)
!      tmp2 = ao_kinetic_integrals_cgtos(i,j)
!      dif = abs(tmp1 - tmp2)
!      if(dif .gt. 1d-10) then
!        print*, ' error on:', i, j
!        print*, tmp1, tmp2, dif
!        !stop
!      endif
!      acc += dif
!      nrm += abs(tmp1)
!    enddo
!  enddo
!  print *, ' acc (%) = ', 100.d0 * acc / nrm

  ! ---

  print *, "NAI:"
  acc = 0.d0
  nrm = 0.d0
  do i = 1, ao_num
  !do i = 9, 9
    do j = 1, ao_num
    !do j = 16, 16
      tmp1 = ao_integrals_n_e      (i,j)
      tmp2 = ao_integrals_n_e_cgtos(i,j)
      dif = dabs(tmp1 - tmp2)
      if(dif .gt. 1d-10) then
        print*, ' error on:', i, j
        print*, tmp1, tmp2, dif
        stop
      endif
      acc += dif
      nrm += dabs(tmp1)
    enddo
  enddo
  print *, ' acc (%) = ', 100.d0 * acc / nrm

end

! ---


subroutine check_ao_two_e_integral_cgtos()

  implicit none

  integer                    :: i, j, k, l
  double precision           :: acc, nrm, dif
  double precision           :: tmp1, tmp2
  double precision           :: t1, t2, tt

  double precision, external :: ao_two_e_integral
  double precision, external :: ao_two_e_integral_cgtos

  acc = 0.d0
  nrm = 0.d0

  tt = 0.d0
  do i = 1, ao_num
  !do i = 1, 1
    call wall_time(t1)
    do j = 1, ao_num
    !do j = 1, 1
      do k = 1, ao_num
      !do k = 1, 1
        do l = 1, ao_num
        !do l = 21, 21

          !call deb_ao_2eint_cgtos(i, j, k, l)

          tmp1 = ao_two_e_integral      (i, j, k, l)
          tmp2 = ao_two_e_integral_cgtos(i, j, k, l)

          dif = abs(tmp1 - tmp2)
          if(dif .gt. 1d-10) then
            print*, ' error on:', i, j, k, l
            print*, tmp1, tmp2, dif
            !stop
          endif
          acc += dif
          nrm += abs(tmp1)
        enddo
      enddo
    enddo
    call wall_time(t2)
    tt += t2 - t1
    print*, " % done = ", 100.d0 * dble(i) / ao_num
    print*, ' ellapsed time (sec) =', tt
  enddo

  !print *, ' acc (%) = ', 100.d0 * acc / nrm

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

  n_test = 1

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

    n = 0
    !rho = (-6.83897018210218d0, -7.24479852507338d0)
    rho = (-9.83206247355480d0, 0.445269582329036d0)

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

subroutine check_crint4()

  implicit none

  integer                    :: i_test, n_test
  integer                    :: i, seed_size, clock_time
  double precision           :: xr(1), x, shift
  double precision           :: yr(1), y
  double precision           :: dif_re, dif_im, acc_re, nrm_re, acc_im, nrm_im
  double precision           :: t1, t2, t_int1, t_int2
  complex*16                 :: rho
  complex*16                 :: int1, int2, int3
  integer, allocatable       :: seed(:)



  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock_time)                  
  seed = clock_time + 37 * (/ (i, i=0, seed_size-1) /) 
  !seed = 123456789
  call random_seed(put=seed)


  t_int1 = 0.d0
  t_int2 = 0.d0

  n_test = 100
  shift = 15.d0

  acc_re = 0.d0
  nrm_re = 0.d0
  acc_im = 0.d0
  nrm_im = 0.d0
  do i_test = 1, n_test

    call random_number(xr)
    call random_number(yr)

    x = 1.d0 * (2.d0 * shift * xr(1) - shift)
    y = 1.d0 * (2.d0 * shift * yr(1) - shift)

    rho = x + (0.d0, 1.d0) * y

    call wall_time(t1)
    call zboysfun00_1(rho, int1)
    call wall_time(t2)
    t_int1 = t_int1 + t2 - t1

    call wall_time(t1)
    call zboysfun00_2(rho, int2)
    call wall_time(t2)
    t_int2 = t_int2 + t2 - t1

    dif_re = dabs(real(int1) - real(int2))
    dif_im = dabs(aimag(int1) - aimag(int2))
    if((dif_re .gt. 1d-10) .or. (dif_im .gt. 1d-10)) then
      print*, ' important error found: '
      print*, " rho = ", x, y
      print*, real(int1), real(int2), dif_re
      print*, aimag(int1), aimag(int2), dif_im
      call crint_quad_12(0, rho, 10000000, int3)
      if(zabs(int1 - int3) .lt. zabs(int2 - int3)) then
        print*, ' implementation 2 seems to be wrong'
      else
        print*, ' implementation 1 seems to be wrong'
        print*, ' quad 10000000:', real(int3), aimag(int3)
        call crint_quad_12(0, rho, 100000000, int3)
        print*, ' quad 100000000:', real(int3), aimag(int3)
      endif
      !print*, ' quad:', real(int3), aimag(int3)
      !stop
    endif

    if((real(int1) /= real(int1)) .or. (aimag(int1) /= aimag(int1)) .or. &
       (real(int2) /= real(int2)) .or. (aimag(int2) /= aimag(int2)) ) then
      cycle
    else
      acc_re += dif_re
      acc_im += dif_im
      nrm_re += dabs(real(int1))
      nrm_im += dabs(aimag(int1))
    endif
  enddo

  print*, "accuracy on real part (%):", 100.d0 * acc_re / (nrm_re + 1d-15)
  print*, "accuracy on imag part (%):", 100.d0 * acc_im / (nrm_im + 1d-15)

  print*, "zerf_1 wall time (sec) = ", t_int1
  print*, "zerf_2 wall time (sec) = ", t_int2


  deallocate(seed)

end

! ---

subroutine check_crint5()

  implicit none

  integer                    :: i_test, n_test
  integer                    :: i, seed_size, clock_time
  integer                    :: n
  double precision           :: xr(1), yr(1), nr(1), x, shift, y
  double precision           :: dif1_re, dif1_im, acc1_re, acc1_im
  double precision           :: dif2_re, dif2_im, acc2_re, acc2_im
  double precision           :: nrm_re, nrm_im
  double precision           :: t1, t2, t_int1, t_int2
  complex*16                 :: rho
  complex*16                 :: int1, int2, int_ref
  integer, allocatable       :: seed(:)

  complex*16, external       :: crint_1, crint_2



  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock_time)                  
  seed = clock_time + 37 * (/ (i, i=0, seed_size-1) /) 
  !seed = 123456789
  call random_seed(put=seed)


  t_int1 = 0.d0
  t_int2 = 0.d0

  n_test = 100

  acc1_re = 0.d0
  acc1_im = 0.d0
  acc2_re = 0.d0
  acc2_im = 0.d0
  nrm_re = 0.d0
  nrm_im = 0.d0
  do i_test = 1, n_test

    call random_number(xr)
    call random_number(yr)
    call random_number(nr)

    x = 1.d+1 * (30.d0 * xr(1) - 15.d0)
    y = 1.d+1 * (30.d0 * yr(1) - 15.d0)
    n = int(16.d0 * nr(1))

    rho = x + (0.d0, 1.d0) * y

    call wall_time(t1)
    int1 = crint_1(n, rho)
    call wall_time(t2)
    t_int1 = t_int1 + t2 - t1

    call wall_time(t1)
    int2 = crint_2(n, rho)
    call wall_time(t2)
    t_int2 = t_int2 + t2 - t1

    call crint_quad_12(n, rho, 10000000, int_ref)

    dif1_re = dabs(real(int1) - real(int_ref))
    dif1_im = dabs(aimag(int1) - aimag(int_ref))

    dif2_re = dabs(real(int2) - real(int_ref))
    dif2_im = dabs(aimag(int2) - aimag(int_ref))

    if((dif2_re .gt. 1d-7) .or. (dif2_im .gt. 1d-7)) then
      print*, ' important error found: '
      print*, " n, rho = ", n, x, y
      print*, real(int1), real(int2), real(int_ref)
      print*, aimag(int1), aimag(int2), aimag(int_ref)
      !stop
    endif

    acc1_re += dif1_re
    acc1_im += dif1_im

    acc2_re += dif2_re
    acc2_im += dif2_im

    nrm_re += dabs(real(int_ref))
    nrm_im += dabs(aimag(int_ref))
  enddo

  print*, "accuracy on boys_1 (%):", 100.d0 * acc1_re / (nrm_re + 1d-15), 100.d0 * acc1_im / (nrm_im + 1d-15)
  print*, "accuracy on boys_2 (%):", 100.d0 * acc1_re / (nrm_re + 1d-15), 100.d0 * acc2_im / (nrm_im + 1d-15)

  print*, "boys_1 wall time (sec) = ", t_int1
  print*, "boys_2 wall time (sec) = ", t_int2


  deallocate(seed)

end

! ---

subroutine check_crint6()

  implicit none

  integer                    :: i_test, n_test
  integer                    :: i, seed_size, clock_time
  integer                    :: n
  double precision           :: xr(1), yr(1), nr(1), x, shift, y
  double precision           :: dif_re, dif_im, acc_re, acc_im
  double precision           :: nrm_re, nrm_im
  double precision           :: t1, t2, t_int1, t_int2
  complex*16                 :: rho
  complex*16                 :: int1, int2, int3
  integer, allocatable       :: seed(:)

  complex*16, external       :: crint_1, crint_2



  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock_time)                  
  seed = clock_time + 37 * (/ (i, i=0, seed_size-1) /) 
  !seed = 123456789
  call random_seed(put=seed)


  t_int1 = 0.d0
  t_int2 = 0.d0

  n_test = 100

  acc_re = 0.d0
  acc_im = 0.d0
  nrm_re = 0.d0
  nrm_im = 0.d0
  do i_test = 1, n_test

    call random_number(xr)
    call random_number(yr)
    call random_number(nr)

    x = 1.d0 * (30.d0 * xr(1) - 15.d0)
    y = 1.d0 * (30.d0 * yr(1) - 15.d0)
    n = int(16.d0 * nr(1))

    rho = x + (0.d0, 1.d0) * y

    call wall_time(t1)
    int1 = crint_1(n, rho)
    call wall_time(t2)
    t_int1 = t_int1 + t2 - t1

    call wall_time(t1)
    int2 = crint_2(n, rho)
    call wall_time(t2)
    t_int2 = t_int2 + t2 - t1

    dif_re = dabs(real(int1) - real(int2))
    dif_im = dabs(aimag(int1) - aimag(int2))

    if((dif_re .gt. 1d-10) .or. (dif_im .gt. 1d-10)) then
      print*, ' important error found: '
      print*, " n, rho = ", n, x, y
      print*, real(int1), real(int2), dif_re
      print*, aimag(int1), aimag(int2), dif_im
      call crint_quad_12(n, rho, 100000000, int3)
      print*, ' quad 100000000:', real(int3), aimag(int3)
      !print*, ' quad 100000000:', dabs(real(int1) - real(int3)), dabs(aimag(int1) - aimag(int3))
      !stop
    endif

    acc_re += dif_re
    acc_im += dif_im
    nrm_re += dabs(real(int1))
    nrm_im += dabs(aimag(int1))
  enddo

  print*, "diff (%):", 100.d0 * acc_re / (nrm_re + 1d-15), 100.d0 * acc_im / (nrm_im + 1d-15)

  print*, "boys_1 wall time (sec) = ", t_int1
  print*, "boys_2 wall time (sec) = ", t_int2


  deallocate(seed)

end

! ---


! ---

subroutine deb_ao_2eint_cgtos(i, j, k, l)

  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in)        :: i, j, k, l
                             
  integer                    :: p, q, r, s
  integer                    :: num_i, num_j, num_k, num_l, dim1, I_power(3), J_power(3), K_power(3), L_power(3)
  integer                    :: iorder_p1(3), iorder_p2(3), iorder_q1(3), iorder_q2(3)
  complex*16                 :: I_center(3), J_center(3), K_center(3), L_center(3)
  complex*16                 :: expo1, expo2, expo3, expo4
  complex*16                 :: P1_center(3), pp1
  complex*16                 :: P2_center(3), pp2
  complex*16                 :: Q1_center(3), qq1
  complex*16                 :: Q2_center(3), qq2



  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)

  if(num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k) then

    !print*, ao_prim_num(i), ao_prim_num(j), ao_prim_num(k), ao_prim_num(l)

    do p = 1, 3
      I_power(p)  = ao_power(i,p)
      J_power(p)  = ao_power(j,p)
      K_power(p)  = ao_power(k,p)
      L_power(p)  = ao_power(l,p)
      I_center(p) = nucl_coord(num_i,p) * (1.d0, 0.d0) 
      J_center(p) = nucl_coord(num_j,p) * (1.d0, 0.d0)
      K_center(p) = nucl_coord(num_k,p) * (1.d0, 0.d0)
      L_center(p) = nucl_coord(num_l,p) * (1.d0, 0.d0)
    enddo

    do p = 1, ao_prim_num(i)
      expo1 = ao_expo_cgtos_ord_transp(p,i) 
      !print*, "expo1 = ", expo1
      !print*, "center1 = ", I_center

      do q = 1, ao_prim_num(j)
        expo2 = ao_expo_cgtos_ord_transp(q,j) 
        !print*, "expo2 = ", expo2
        !print*, "center2 = ", J_center

        pp1 = expo1 + expo2
        P1_center(1:3) = (expo1 * I_center(1:3) + expo2 * J_center(1:3)) / pp1
        iorder_p1(1:3) = I_power(1:3) + J_power(1:3)

        pp2 = conjg(expo1) + expo2
        P2_center(1:3) = (conjg(expo1) * I_center(1:3) + expo2 * J_center(1:3)) / pp2
        iorder_p2(1:3) = I_power(1:3) + J_power(1:3)

        do r = 1, ao_prim_num(k)
          expo3 = ao_expo_cgtos_ord_transp(r,k) 
          !print*, "expo3 = ", expo3
          !print*, "center3 = ", K_center

          do s = 1, ao_prim_num(l)
            expo4 = ao_expo_cgtos_ord_transp(s,l) 
            !print*, "expo4 = ", expo4
            !print*, "center4 = ", L_center

            qq1 = expo3 + expo4
            Q1_center(1:3) = (expo3 * K_center(1:3) + expo4 * L_center(1:3)) / qq1
            iorder_q1(1:3) = K_power(1:3) + L_power(1:3)

            qq2 = conjg(expo3) + expo4
            Q2_center(1:3) = (conjg(expo3) * K_center(1:3) + expo4 * L_center(1:3)) / qq2
            iorder_q2(1:3) = K_power(1:3) + L_power(1:3)

            call deb_cboys(P1_center, pp1, iorder_p1, Q1_center, qq1, iorder_q1)
            call deb_cboys(P1_center, pp1, iorder_p1, Q2_center, qq2, iorder_q2)
            call deb_cboys(P2_center, pp2, iorder_p2, Q1_center, qq1, iorder_q1)
            call deb_cboys(P2_center, pp2, iorder_p2, Q2_center, qq2, iorder_q2)
            call deb_cboys(conjg(P2_center), conjg(pp2), iorder_p2, Q1_center, qq1, iorder_q1)
            call deb_cboys(conjg(P2_center), conjg(pp2), iorder_p2, Q2_center, qq2, iorder_q2)
            call deb_cboys(conjg(P1_center), conjg(pp1), iorder_p1, Q1_center, qq1, iorder_q1)
            call deb_cboys(conjg(P1_center), conjg(pp1), iorder_p1, Q2_center, qq2, iorder_q2)
          enddo ! s
        enddo ! r
      enddo ! q
    enddo ! p

  endif ! same centers

  return
end

! ---

subroutine deb_cboys(P_center, p, iorder_p, Q_center, q, iorder_q)


  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in) :: iorder_p(3), iorder_q(3)
  complex*16, intent(in) :: P_center(3), p
  complex*16, intent(in) :: Q_center(3), q

  integer                :: iorder, n
  complex*16             :: dist, rho
  complex*16             :: int1, int2

  complex*16, external   :: crint_2


  dist = (P_center(1) - Q_center(1)) * (P_center(1) - Q_center(1)) &
       + (P_center(2) - Q_center(2)) * (P_center(2) - Q_center(2)) &
       + (P_center(3) - Q_center(3)) * (P_center(3) - Q_center(3))
  rho = dist * p * q / (p + q)

  if(real(rho) .lt. -5.d0) then
    print*, 'warning ! impotant negative rho: ', rho
  endif

  !if(abs(rho) .lt. 1d-15) return

  iorder = 2*iorder_p(1)+2*iorder_q(1) + 2*iorder_p(2)+2*iorder_q(2) + 2*iorder_p(3)+2*iorder_q(3)
  n = shiftr(iorder, 1)
  
  !write(33,*) n, real(rho), aimag(rho)
  !print*, n, real(rho), aimag(rho)

  int1 = crint_2(n, rho)
  call crint_quad_12(n, rho, 1000000, int2)

  if(abs(int1 - int2) .gt. 1d-5) then
    print*, ' important error found: '
    print*, p!, P_center
    print*, q!, Q_center
    print*, dist
    print*, " n, tho = ", n, real(rho), aimag(rho)
    print*, real(int1), real(int2), dabs(real(int1-int2))
    print*, aimag(int1), aimag(int2), dabs(aimag(int1-int2))
    stop
  endif

end

! ---

