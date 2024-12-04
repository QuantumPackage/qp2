program extra_basis
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
 integer :: i
 print*,'extra_nucl_num = ',extra_nucl_num
 do i = 1, extra_nucl_num
  print*,'i = ',i
  print*,'extra_nucl_label  = ',extra_nucl_label(i)
  print*,'extra_nucl_charge = ',extra_nucl_charge(i)
  print*,'extra_nucl_coord  = '
  print*,extra_nucl_coord(i,1:3)
 enddo
end
