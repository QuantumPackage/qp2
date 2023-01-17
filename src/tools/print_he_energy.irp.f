
! ---

program print_he_energy

  implicit none

  call print_overlap()

  !call print_energy1()
  !call print_energy2()

end 

! ---

subroutine print_overlap()

  implicit none
  integer          :: i, j, k, l
  double precision :: S_ij

  print *, ' ao_overlap:'
  do i = 1, ao_num
    do j = 1, ao_num
      print *, j, i, ao_overlap(j,i)
    enddo
  enddo


  print *, ' mo_overlap:'
  do i = 1, mo_num
    do j = 1, mo_num

      S_ij = 0.d0
      do k = 1, ao_num
        do l = 1, ao_num
          S_ij += mo_coef(k,i) * ao_overlap(k,l) * mo_coef(l,j)
        enddo
      enddo
      
      print *, i, j, S_ij
    enddo
  enddo

end subroutine print_overlap

! ---

subroutine print_energy1()

  implicit none
  integer                    :: i, j, k, l
  double precision           :: e, n, e_tmp, n_tmp
  double precision, external :: ao_two_e_integral

  e = 0.d0
  n = 0.d0

  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  ! < phi_1 phi_1 | h1 | phi_1 phi_1 > 

  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      e_tmp += mo_coef(i,1) * ao_one_e_integrals(i,j) * mo_coef(j,1)
      n_tmp += mo_coef(i,1) * ao_overlap(i,j)         * mo_coef(j,1)
    enddo
  enddo

  e += e_tmp * n_tmp

  ! ---

  ! < phi_1 phi_1 | h2 | phi_1 phi_1 > 

  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      n_tmp += mo_coef(i,1) * ao_overlap(i,j)         * mo_coef(j,1)
      e_tmp += mo_coef(i,1) * ao_one_e_integrals(i,j) * mo_coef(j,1)
    enddo
  enddo

  e += e_tmp * n_tmp

  ! ---

  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  ! ---

  do i = 1, ao_num
    do j = 1, ao_num
      do k = 1, ao_num
        do l = 1, ao_num

          ! ( phi_1 phi_1 | phi_1 phi_1 )
          e += mo_coef(i,1) * mo_coef(j,1) * ao_two_e_integral(i,j,k,l) * mo_coef(k,1) * mo_coef(l,1)

        enddo
      enddo
    enddo
  enddo

  ! ---

  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  ! ---

  ! < phi_1 phi_1 | phi_1 phi_1 >
  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      e_tmp += mo_coef(i,1) * ao_overlap(i,j) * mo_coef(j,1)
      n_tmp += mo_coef(i,1) * ao_overlap(i,j) * mo_coef(j,1)
    enddo
  enddo

  n += e_tmp * n_tmp

  ! ---

  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  e = e / n 
  print *, ' energy = ', e

end subroutine print_energy1

! ---

subroutine print_energy2()

  implicit none
  integer                    :: i, j, k, l
  double precision           :: e, n, e_tmp, n_tmp
  double precision, external :: ao_two_e_integral

  e = 0.d0
  n = 0.d0

  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  ! < phi_1 phi_2 | h1 | phi_1 phi_2 > 

  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      e_tmp += mo_coef(i,1) * ao_one_e_integrals(i,j) * mo_coef(j,1)
      n_tmp += mo_coef(i,2) * ao_overlap(i,j)         * mo_coef(j,2)
    enddo
  enddo

  e += e_tmp * n_tmp

  ! ---

  ! < phi_1 phi_2 | h1 | phi_2 phi_1 > 

  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      e_tmp += mo_coef(i,1) * ao_one_e_integrals(i,j) * mo_coef(j,2)
      n_tmp += mo_coef(i,2) * ao_overlap(i,j)         * mo_coef(j,1)
    enddo
  enddo

  e -= e_tmp * n_tmp

  ! ---

  ! < phi_2 phi_1 | h1 | phi_1 phi_2 > 

  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      e_tmp += mo_coef(i,2) * ao_one_e_integrals(i,j) * mo_coef(j,1)
      n_tmp += mo_coef(i,1) * ao_overlap(i,j)         * mo_coef(j,2)
    enddo
  enddo

  e -= e_tmp * n_tmp

  ! ---

  ! < phi_2 phi_1 | h1 | phi_2 phi_1 >

  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      e_tmp += mo_coef(i,2) * ao_one_e_integrals(i,j) * mo_coef(j,2)
      n_tmp += mo_coef(i,1) * ao_overlap(i,j)         * mo_coef(j,1)
    enddo
  enddo

  e += e_tmp * n_tmp

  ! ---

  ! < phi_1 phi_2 | h2 | phi_1 phi_2 > 

  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      n_tmp += mo_coef(i,1) * ao_overlap(i,j)         * mo_coef(j,1)
      e_tmp += mo_coef(i,2) * ao_one_e_integrals(i,j) * mo_coef(j,2)
    enddo
  enddo

  e += e_tmp * n_tmp

  ! ---

  ! < phi_1 phi_2 | h2 | phi_2 phi_1 >

  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      n_tmp += mo_coef(i,1) * ao_overlap(i,j)         * mo_coef(j,2)
      e_tmp += mo_coef(i,2) * ao_one_e_integrals(i,j) * mo_coef(j,1)
    enddo
  enddo

  e -= e_tmp * n_tmp

  ! ---

  ! < phi_2 phi_1 | h2 | phi_1 phi_2 >

  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      n_tmp += mo_coef(i,2) * ao_overlap(i,j)         * mo_coef(j,1)
      e_tmp += mo_coef(i,1) * ao_one_e_integrals(i,j) * mo_coef(j,2)
    enddo
  enddo

  e -= e_tmp * n_tmp

  ! ---

  ! < phi_2 phi_1 | h2 | phi_2 phi_1 > 

  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      n_tmp += mo_coef(i,2) * ao_overlap(i,j)         * mo_coef(j,2)
      e_tmp += mo_coef(i,1) * ao_one_e_integrals(i,j) * mo_coef(j,1)
    enddo
  enddo

  e += e_tmp * n_tmp

  ! ---

  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  ! ---

  do i = 1, ao_num
    do j = 1, ao_num
      do k = 1, ao_num
        do l = 1, ao_num

          ! ( phi_1 phi_1 | phi_2 phi_2 )
          e += mo_coef(i,1) * mo_coef(j,1) * ao_two_e_integral(i,j,k,l) * mo_coef(k,2) * mo_coef(l,2)

          ! ( phi_1 phi_2 | phi_2 phi_1 )
          e -= mo_coef(i,1) * mo_coef(j,2) * ao_two_e_integral(i,j,k,l) * mo_coef(k,2) * mo_coef(l,1)

          ! ( phi_2 phi_1 | phi_1 phi_2 )
          e -= mo_coef(i,2) * mo_coef(j,1) * ao_two_e_integral(i,j,k,l) * mo_coef(k,1) * mo_coef(l,2)

          ! ( phi_2 phi_2 | phi_1 phi_1 )
          e += mo_coef(i,2) * mo_coef(j,2) * ao_two_e_integral(i,j,k,l) * mo_coef(k,1) * mo_coef(l,1)

        enddo
      enddo
    enddo
  enddo

  ! ---

  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  ! ---

  ! < phi_1 phi_2 | phi_1 phi_2 >
  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      e_tmp += mo_coef(i,1) * ao_overlap(i,j) * mo_coef(j,1)
      n_tmp += mo_coef(i,2) * ao_overlap(i,j) * mo_coef(j,2)
    enddo
  enddo

  n += e_tmp * n_tmp

  ! ---

  ! < phi_1 phi_2 | phi_2 phi_1 >
  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      e_tmp += mo_coef(i,1) * ao_overlap(i,j) * mo_coef(j,2)
      n_tmp += mo_coef(i,2) * ao_overlap(i,j) * mo_coef(j,1)
    enddo
  enddo

  n -= e_tmp * n_tmp

  ! ---

  ! < phi_2 phi_1 | phi_1 phi_2 >
  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      e_tmp += mo_coef(i,2) * ao_overlap(i,j) * mo_coef(j,1)
      n_tmp += mo_coef(i,1) * ao_overlap(i,j) * mo_coef(j,2)
    enddo
  enddo

  n -= e_tmp * n_tmp

  ! ---

  ! < phi_2 phi_1 | phi_2 phi_1 >
  e_tmp = 0.d0
  n_tmp = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      e_tmp += mo_coef(i,2) * ao_overlap(i,j) * mo_coef(j,2)
      n_tmp += mo_coef(i,1) * ao_overlap(i,j) * mo_coef(j,1)
    enddo
  enddo

  n += e_tmp * n_tmp

  ! ---

  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  e = e / n 
  print *, ' energy = ', e

end subroutine print_energy2

! ---
