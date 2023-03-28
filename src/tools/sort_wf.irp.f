program sort_wf
 implicit none
 read_wf = .true.
 touch read_wf 
 call routine

end

subroutine routine
 implicit none
 integer :: i
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.wf_sorted'
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, N_det
  write(i_unit_output, *)i,dabs(psi_coef_sorted(i,1))/dabs(psi_coef_sorted(1,1))
 enddo

end
