program print_wf_qp_edit

  implicit none

  BEGIN_DOC
  ! Print the psi_det wave function up to n_det_qp_edit
  END_DOC

  read_wf = .true.
  touch read_wf

  call print_wf_qp_edit_subroutine()

end

subroutine print_wf_qp_edit_subroutine()

  implicit none
   BEGIN_DOC
   ! Prints the first n_det_qp_edit determinants of the wave
   ! function stored in psi_coef and psi_det: for each determinant
   ! its index, its coefficients psi_coef(i,:) and the determinant
   ! itself, via print_det.
   END_DOC

  integer :: i

  do i = 1, n_det_qp_edit
    print*,i
    write(*,'(100(1pE12.4))') psi_coef(i,:)
    call print_det(psi_det(1,1,i),N_int)
    print*,''
  enddo

end
