program tc_bi_ortho

  BEGIN_DOC
  ! TODO
  END_DOC
  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call KMat_tilde_dump()
end

! ---

subroutine KMat_tilde_dump()

  implicit none
  integer                       :: i, j, k, l

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
