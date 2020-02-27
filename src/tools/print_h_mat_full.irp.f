program hdump
  read_wf=.True.
  touch read_wf
  call hmatdump
end
subroutine hmatdump
  implicit none
  BEGIN_DOC
!
  END_DOC
  character*(128) :: output
  integer :: iunit,getUnitAndOpen

  integer :: i,j,k,l
  integer :: i1,j1,k1,l1
  integer :: i2,j2,k2,l2
  integer*8 :: m
  character*(2), allocatable :: A(:)
  PROVIDE h_matrix_all_dets

  if(N_det.ge.10000)then
    print*,'Warning !!!'
    print*,'Number of determinants is ',N_det
    print*,'It means that the H matrix will be enormous !'
    print*,'stoppping ..'
    stop
  endif
  output=trim(ezfio_filename)//'.HDUMP'
  iunit = getUnitAndOpen(output,'w')

 do i = 1, N_det
  write(iunit,'(10000(F25.16))') H_matrix_all_dets(i,:)
 enddo
  close(iunit)
end
