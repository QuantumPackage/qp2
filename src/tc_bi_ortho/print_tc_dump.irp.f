program tc_bi_ortho

  BEGIN_DOC
  ! TODO
  END_DOC
  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call ERI_dump()
  call KMat_tilde_dump()
  call LMat_tilde_dump()

end

! ---

subroutine KMat_tilde_dump()

  implicit none
  integer :: i, j, k, l

  PROVIDE mo_bi_ortho_tc_two_e_chemist

  print *, ' Kmat_tilde in chem notation'

  open(33, file='Kmat_tilde.dat', action='write')  
    do l = 1, mo_num
      do k = 1, mo_num
        do j = 1, mo_num
          do i = 1, mo_num
            write(33, '(4(I4, 2X), 4X, E15.7)') i, j, k, l, mo_bi_ortho_tc_two_e_chemist(i,j,k,l)
            ! TCHint convention
            !write(33, '(4(I4, 2X), 4X, E15.7)') i, j, k, l, mo_bi_ortho_tc_two_e_chemist(j,i,l,k)
          enddo
        enddo
      enddo
    enddo
  close(33)

  return
end subroutine KMat_tilde_dump

! ---

subroutine ERI_dump()

  implicit none
  integer                       :: i, j, k, l
  double precision, allocatable :: a1(:,:,:,:), a2(:,:,:,:)

  PROVIDE mo_r_coef mo_l_coef

  allocate(a2(ao_num,ao_num,ao_num,mo_num))

  call dgemm( 'T', 'N', ao_num*ao_num*ao_num, mo_num, ao_num, 1.d0   &
            , ao_two_e_coul(1,1,1,1), ao_num, mo_l_coef(1,1), ao_num &
            , 0.d0 , a2(1,1,1,1), ao_num*ao_num*ao_num)

  allocate(a1(ao_num,ao_num,mo_num,mo_num))

  call dgemm( 'T', 'N', ao_num*ao_num*mo_num, mo_num, ao_num, 1.d0 &
            , a2(1,1,1,1), ao_num, mo_r_coef(1,1), ao_num          &
            , 0.d0, a1(1,1,1,1), ao_num*ao_num*mo_num)

  deallocate(a2)
  allocate(a2(ao_num,mo_num,mo_num,mo_num))

  call dgemm( 'T', 'N', ao_num*mo_num*mo_num, mo_num, ao_num, 1.d0 &
            , a1(1,1,1,1), ao_num, mo_l_coef(1,1), ao_num          &
            , 0.d0, a2(1,1,1,1), ao_num*mo_num*mo_num)

  deallocate(a1)
  allocate(a1(mo_num,mo_num,mo_num,mo_num))

  call dgemm( 'T', 'N', mo_num*mo_num*mo_num, mo_num, ao_num, 1.d0 &
            , a2(1,1,1,1), ao_num, mo_r_coef(1,1), ao_num          &
            , 0.d0, a1(1,1,1,1), mo_num*mo_num*mo_num)

  deallocate(a2)

  open(33, file='ERI.dat', action='write')  
    do l = 1, mo_num
      do k = 1, mo_num
        do j = 1, mo_num
          do i = 1, mo_num
            write(33, '(4(I4, 2X), 4X, E15.7)') i, j, k, l, a1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo
  close(33)

  deallocate(a1)

  return
end subroutine ERI_dump

! ---

subroutine LMat_tilde_dump()

  implicit none
  integer          :: i, j, k, l, m, n
  double precision :: integral

  PROVIDE mo_bi_ortho_tc_two_e_chemist

  print *, ' Lmat_tilde in phys notation'

  open(33, file='Lmat_tilde.dat', action='write')  
    do n = 1, mo_num
      do m = 1, mo_num
        do l = 1, mo_num
          do k = 1, mo_num
            do j = 1, mo_num
              do i = 1, mo_num
                ! < i j k | -L | l m n > with a BI-ORTHONORMAL MOLECULAR ORBITALS
                call give_integrals_3_body_bi_ort(i, j, k, l, m, n, integral)
                write(33, '(6(I4, 2X), 4X, E15.7)') i, j, k, l, m, n, integral
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  close(33)

  return
end subroutine LMat_tilde_dump

! ---
