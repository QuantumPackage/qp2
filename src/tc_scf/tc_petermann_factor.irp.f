
! ---

program tc_petermann_factor

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call main()

end

! ---

subroutine main()

  implicit none
  integer                       :: i, j
  double precision              :: Pf_diag_av
  double precision, allocatable :: Sl(:,:), Sr(:,:), Pf(:,:)

  allocate(Sl(mo_num,mo_num), Sr(mo_num,mo_num), Pf(mo_num,mo_num))


  call LTxSxR(ao_num, mo_num, mo_l_coef, ao_overlap, mo_r_coef, Sl)
  !call dgemm( "T", "N", mo_num, mo_num, ao_num, 1.d0                       &
  !          , mo_l_coef, size(mo_l_coef, 1), mo_l_coef, size(mo_l_coef, 1) &
  !          , 0.d0, Sl, size(Sl, 1) )

  print *, ''
  print *, ' left-right orthog matrix:'
  do i = 1, mo_num
    write(*,'(100(F8.4,X))') Sl(:,i)
  enddo

  call LTxSxR(ao_num, mo_num, mo_l_coef, ao_overlap, mo_l_coef, Sl)
  !call dgemm( "T", "N", mo_num, mo_num, ao_num, 1.d0                       &
  !          , mo_l_coef, size(mo_l_coef, 1), mo_l_coef, size(mo_l_coef, 1) &
  !          , 0.d0, Sl, size(Sl, 1) )

  print *, ''
  print *, ' left-orthog matrix:'
  do i = 1, mo_num
    write(*,'(100(F8.4,X))') Sl(:,i)
  enddo

  call LTxSxR(ao_num, mo_num, mo_r_coef, ao_overlap, mo_r_coef, Sr)
!  call dgemm( "T", "N", mo_num, mo_num, ao_num, 1.d0                       &
!            , mo_r_coef, size(mo_r_coef, 1), mo_r_coef, size(mo_r_coef, 1) &
!            , 0.d0, Sr, size(Sr, 1) )

  print *, ''
  print *, ' right-orthog matrix:'
  do i = 1, mo_num
    write(*,'(100(F8.4,X))') Sr(:,i)
  enddo

  print *, ''
  print *, ' Petermann matrix:'
  do i = 1, mo_num
    do j = 1, mo_num
      Pf(j,i) = Sl(j,i) * Sr(j,i)
    enddo
    write(*,'(100(F8.4,X))') Pf(:,i)
  enddo

  Pf_diag_av = 0.d0
  do i = 1, mo_num
    Pf_diag_av = Pf_diag_av + Pf(i,i)
  enddo
  Pf_diag_av = Pf_diag_av / dble(mo_num)

  print *, ''
  print *, ' mean of the diagonal Petermann factor = ', Pf_diag_av

  deallocate(Sl, Sr, Pf)

  return
end subroutine

! ---

