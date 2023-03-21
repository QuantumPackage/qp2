program run

  implicit none

  read_wf = .true.
  touch read_wf

  call print_wf_qp_edit()

end

subroutine print_wf_qp_edit()
  
  implicit none

  BEGIN_DOC
  ! Print the psi_det wave function up to n_det_qp_edit
  END_DOC
  
  integer :: i

  do i = 1, n_det_qp_edit
    print*,i
    write(*,'(100(1pE12.4))') psi_coef(i,:)
    call print_det(psi_det(1,1,i),N_int)
    print*,''
  enddo

end
