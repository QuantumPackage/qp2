program sort_wf
 implicit none
 BEGIN_DOC
 ! Writes to the file <ezfio_filename>.wf_sorted the amplitudes
 ! |c_i/c_1| of the wave function sorted by contribution to the
 ! norm, read from the |EZFIO| directory, one line per determinant.
 END_DOC
 read_wf = .true.
 touch read_wf 
 call routine

end

subroutine routine
 implicit none
 BEGIN_DOC
 ! Writes <ezfio_filename>.wf_sorted: one line per determinant, with
 ! its index i and |psi_coef_sorted(i,1)|/|psi_coef_sorted(1,1)|.
 END_DOC
 integer :: i
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.wf_sorted'
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, N_det
  write(i_unit_output, *)i,dabs(psi_coef_sorted(i,1))/dabs(psi_coef_sorted(1,1))
 enddo

end
