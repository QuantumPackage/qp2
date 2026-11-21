program check_orthonormality
 implicit none
 BEGIN_DOC
 ! Checks the orthonormality of the MOs by printing the deviation
 ! of the average diagonal element of the MO overlap matrix from
 ! 1 and the average of the absolute values of the off-diagonal
 ! elements.
 END_DOC
 integer :: i,j

 call do_print
 print *,  '--'
! call orthonormalize_mos
! call do_print
end

subroutine do_print
 implicit none
 BEGIN_DOC
 ! Accumulates the absolute values of the elements of the MO
 ! overlap matrix mo_overlap and prints the deviation of the
 ! average diagonal element from 1, and the average of the
 ! absolute values of the off-diagonal elements.
 END_DOC
 integer :: i,j
 double precision :: off_diag, diag

 off_diag = 0.d0
 diag = 0.d0
 do j=1,mo_num
   do i=1,mo_num
     off_diag += abs(mo_overlap(i,j))
   enddo
   diag += abs(mo_overlap(j,j))
   off_diag -= abs(mo_overlap(j,j))
 enddo
 print *,  'Diag      = ', abs(1.d0-diag/dble(mo_num))
 print *,  'Off-Diag  = ', off_diag/(dble(mo_num)**2-dble(mo_num))
 print *,  '--'

end


