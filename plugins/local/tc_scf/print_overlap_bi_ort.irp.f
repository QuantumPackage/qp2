program print_overlap_bi_ort
 implicit none
 integer :: i,j
 print*,'Overlap between right orbitals '
 do i = 1, mo_num
  write(*,'(40(F10.4,X))')overlap_mo_read_r(i,:)
 enddo
 print*,'Overlap between left orbitals '
 do i = 1, mo_num
  write(*,'(40(F10.4,X))')overlap_mo_read_l(i,:)
 enddo
 
end
