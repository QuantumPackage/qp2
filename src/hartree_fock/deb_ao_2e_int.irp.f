
program deb_ao_2e_int

  !call check_ao_two_e_integral_cosgtos()
  call check_crint1()
  !call check_crint2()

end

! ---

subroutine check_ao_two_e_integral_cosgtos()

  implicit none

  integer                    :: i, j, k, l
  double precision           :: tmp1, tmp2
  double precision           :: acc, nrm, dif

  double precision, external :: ao_two_e_integral
  double precision, external :: ao_two_e_integral_cosgtos

  acc = 0.d0
  nrm = 0.d0

  i = 1
  j = 6
  k = 1
  l = 16
!  do i = 1, ao_num
!    do k = 1, ao_num
!      do j = 1, ao_num
!        do l = 1, ao_num

          tmp1 = ao_two_e_integral        (i, j, k, l)
          tmp2 = ao_two_e_integral_cosgtos(i, j, k, l)

          dif  = dabs(tmp1 - tmp2)
          if(dif .gt. 1d-12) then
            print*, ' error on:', i, j, k, l
            print*, tmp1, tmp2, dif
            stop
          endif

!          acc += dif
!          nrm += dabs(tmp1)
!        enddo
!      enddo
!    enddo
!  enddo

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
  complex*16, external       :: crint_1, crint_2, crint_quad

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
      !int_an = crint_1   (i, rho)
      int_an = crint_2   (i, rho)
      int_nm = crint_quad(i, rho)
      
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


