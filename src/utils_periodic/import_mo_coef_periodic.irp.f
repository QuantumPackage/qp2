program import_mo_coef_periodic
  BEGIN_DOC
  ! Imports the MO coefficients (real and imaginary parts)
  ! from the text file C.qp in the current directory into
  ! the |EZFIO| database named by ezfio_filename.
  ! TODO(doc): mo_num, mo_occ, mo_class and ao_md5 are not
  ! set in this program; are they read from the database by
  ! save_mos?
  END_DOC

  PROVIDE ezfio_filename
  call run
end

subroutine run
  use map_module
  implicit none
  BEGIN_DOC
  ! Reads C.qp line by line: a line "i j re im" is stored
  ! in mo_coef(i,j) and mo_coef_imag(i,j) (real and
  ! imaginary parts). Sets mo_label to "None" and writes
  ! the MO basis to the |EZFIO| database via save_mos.
  ! TODO(doc): is mo_coef_imag stored in the database?
  ! save_mos writes mo_coef but not mo_coef_imag.
  END_DOC
  
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
