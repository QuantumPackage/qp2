
program print_scf_int

  call main()

end

subroutine main()

  implicit none
  integer ::  i, j, k, l

  print *, " Hcore:"
  do j = 1, ao_num
    do i = 1, ao_num
      print *, i, j, ao_one_e_integrals(i,j)
    enddo
  enddo

  print *, " P:"
  do j = 1, ao_num
    do i = 1, ao_num
      print *, i, j, SCF_density_matrix_ao_alpha(i,j)
    enddo
  enddo


  double precision :: integ, density_a, density_b, density
  double precision :: J_scf(ao_num, ao_num)
  double precision :: K_scf(ao_num, ao_num)


  double precision, external :: get_ao_two_e_integral
  PROVIDE ao_integrals_map

  print *, " J:"
  !do j = 1, ao_num
  !  do l = 1, ao_num
  !    do i = 1, ao_num
  !      do k = 1, ao_num
  !        !  < 1:k, 2:l | 1:i, 2:j > 
  !        print *, '< k l | i j >', k, l, i, j
  !        print *, get_ao_two_e_integral(i, j, k, l, ao_integrals_map)
  !      enddo
  !    enddo
  !  enddo
  !enddo

  !do k = 1, ao_num
  !  do i = 1, ao_num
  !    do j = 1, ao_num
  !      do l = 1, ao_num
  !        !  ( 1:k, 1:i | 2:l, 2:j ) 
  !        print *, '(k i | l j)', k, i, l, j
  !        print *, get_ao_two_e_integral(l, j, k, i, ao_integrals_map)
  !      enddo
  !    enddo
  !    print *, ''
  !  enddo
  !enddo

  J_scf = 0.d0
  K_scf = 0.d0
  do i = 1, ao_num
    do k = 1, ao_num
      do j = 1, ao_num
        do l = 1, ao_num

          density_a = SCF_density_matrix_ao_alpha(l,j)
          density_b = SCF_density_matrix_ao_beta (l,j)
          density   = density_a + density_b

          integ = get_ao_two_e_integral(l, j, k, i, ao_integrals_map)
          J_scf(k,i) += density   * integ
          integ = get_ao_two_e_integral(l, i, k, j, ao_integrals_map)
          K_scf(k,i) -= density_a * integ
        enddo
      enddo
    enddo
  enddo

  print *, 'J x P'
  do i = 1, ao_num
    do k = 1, ao_num
      print *, k, i, J_scf(k,i)
    enddo
  enddo

  print *, ''
  print *, 'K x P'
  do i = 1, ao_num
    do k = 1, ao_num
      print *, k, i, K_scf(k,i)
    enddo
  enddo

  print *, ''
  print *, 'F in AO'
  do i = 1, ao_num
    do k = 1, ao_num
      print *, k, i, Fock_matrix_ao(k,i)
    enddo
  enddo

  print *, ''
  print *, 'F in MO'
  do i = 1, ao_num
    do k = 1, ao_num
      print *, k, i, 2.d0 * Fock_matrix_mo_alpha(k,i)
    enddo
  enddo

end

