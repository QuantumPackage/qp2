program print_sorted_wf_coef
 implicit none
 BEGIN_DOC
 ! Writes to the file <ezfio_filename>.wf_sorted the
 ! coefficients of the first state of the wave function stored
 ! in the |EZFIO| directory, in the psi_coef_sorted order
 ! (decreasing state-averaged norm contribution).
 END_DOC
 read_wf =.True.
 call routine

end

subroutine routine
 implicit none
 BEGIN_DOC
 ! Writes one line per determinant to
 ! <ezfio_filename>.wf_sorted: the index i, the magnitude of
 ! psi_coef_sorted(i,1) and its ratio to psi_coef_sorted(1,1).
 END_DOC
 integer :: i
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.wf_sorted'
 i_unit_output = getUnitAndOpen(output,'w')
 do i= 1, N_det
  write(i_unit_output,*)i,dabs(psi_coef_sorted(i,1)),dabs(psi_coef_sorted(i,1)/psi_coef_sorted(1,1))
 enddo

end
