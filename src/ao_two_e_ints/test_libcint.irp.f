program test_libcint

  implicit none

  BEGIN_DOC
  ! Print every |AO| integral produced by the active engine, with tags, so
  ! that two runs differing only by `ao_basis.use_libcint` can be diffed.
  !
  ! All integral lines are of the form
  !
  !     TAG i j value        (one-electron)
  !     TAG i j k l value    (two-electron)
  !
  ! and the accompanying `INFO` lines are not compared.
  END_DOC

  integer                    :: i, j, k, l
  double precision           :: dmax, s
  double precision, external :: ao_two_e_integral
  double precision, external :: ao_two_e_integral_erf

  PROVIDE ao_num shell_num nucl_num ao_l_max use_libcint
  PROVIDE use_libcint_effective libcint_available libcint_nuc_available
  PROVIDE ao_overlap ao_kinetic_integrals ao_integrals_n_e
  PROVIDE mu_erf

  write(6,'(A,1X,I0)') 'INFO ao_num          ', ao_num
  write(6,'(A,1X,I0)') 'INFO shell_num       ', shell_num
  write(6,'(A,1X,I0)') 'INFO nucl_num        ', nucl_num
  write(6,'(A,1X,I0)') 'INFO ao_l_max        ', ao_l_max
  write(6,'(A,1X,L1)') 'INFO use_libcint     ', use_libcint
  write(6,'(A,1X,L1)') 'INFO libcint_avail   ', libcint_available
  write(6,'(A,1X,L1)') 'INFO libcint_nuc_avail', libcint_nuc_available
  write(6,'(A,1X,ES24.16)') 'INFO mu_erf          ', mu_erf

  do j = 1, ao_num
    do i = 1, ao_num
      write(6,'(A,1X,2I5,1X,ES24.16)') 'OVLP', i, j, ao_overlap(i,j)
    enddo
  enddo

  do j = 1, ao_num
    do i = 1, ao_num
      write(6,'(A,1X,2I5,1X,ES24.16)') 'KIN', i, j, ao_kinetic_integrals(i,j)
    enddo
  enddo

  do j = 1, ao_num
    do i = 1, ao_num
      write(6,'(A,1X,2I5,1X,ES24.16)') 'NAI', i, j, ao_integrals_n_e(i,j)
    enddo
  enddo

  do l = 1, ao_num
    do k = 1, l
      do j = 1, ao_num
        do i = 1, j
          write(6,'(A,1X,4I5,1X,ES24.16)') 'ERI', i, j, k, l, &
              ao_two_e_integral(i,j,k,l)
        enddo
      enddo
    enddo
  enddo

  ! Long-range operator erf(mu_erf*r12)/r12, same index ordering
  do l = 1, ao_num
    do k = 1, l
      do j = 1, ao_num
        do i = 1, j
          write(6,'(A,1X,4I5,1X,ES24.16)') 'ERF', i, j, k, l, &
              ao_two_e_integral_erf(i,j,k,l)
        enddo
      enddo
    enddo
  enddo

  if (use_libcint_effective) then
    PROVIDE ao_overlap_qp2
    dmax = 0.d0
    do j = 1, ao_num
      do i = 1, ao_num
        dmax = max(dmax, abs(ao_overlap(i,j) - ao_overlap_qp2(i,j)))
      enddo
    enddo
    write(6,'(A,1X,ES24.16)') 'INFO max|ovlp-libcint - ovlp-qp2|', dmax
  endif

  s = 0.d0
  do j = 1, ao_num
    do i = 1, ao_num
      s = s + ao_two_e_integral(i,i,j,j)
    enddo
  enddo
  write(6,'(A,1X,ES24.16)') 'INFO sum(iijj)', s

end program test_libcint
