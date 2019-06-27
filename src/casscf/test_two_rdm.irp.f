program print_two_rdm
  implicit none
  integer                        :: i,j,k,l
  read_wf = .True.
  TOUCH read_wf

  double precision, parameter    :: thr = 1.d-15

  double precision :: accu,twodm
  accu = 0.d0
  do i=1,n_act_orb
    do j=1,n_act_orb
      do k=1,n_act_orb
        do l=1,n_act_orb
         twodm = coussin_peter_two_rdm_mo(list_act(i),list_act(j),list_act(k),list_act(l),1)
            if(dabs(twodm - P0tuvx(i,j,k,l)).gt.thr)then
             print*,''
             print*,'sum'
             write(*,'(3X,4(I2,X),3(F16.13,X))'),  i, j, k, l, twodm,P0tuvx(i,j,k,l),dabs(twodm - P0tuvx(i,j,k,l))
             print*,''
            endif
            accu += dabs(twodm - P0tuvx(i,j,k,l))
        enddo
      enddo
    enddo
  enddo
  print*,'accu = ',accu
  print*,'<accu> ',accu / dble(mo_num**4)

end
