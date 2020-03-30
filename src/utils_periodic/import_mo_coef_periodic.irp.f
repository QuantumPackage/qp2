program import_mo_coef_periodic

  PROVIDE ezfio_filename
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j
  double precision :: int_re, int_im

  
  iunit = getunitandopen('C.qp','r')
  do 
    read (iunit,*,end=10) i,j, mo_coef(i,j), mo_coef_imag(i,j)
  enddo
  10 continue
  close(iunit)
  mo_label = "None"
  call save_mos

end
