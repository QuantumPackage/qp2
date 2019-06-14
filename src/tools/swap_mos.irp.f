program swap_mos
  implicit none
  BEGIN_DOC
  ! Swaps the indices of two molecular orbitals
  END_DOC
  integer                        :: i,j, i1, i2
  double precision               :: x
  print *,  'MOs to swap?'
  read(*,*) i1, i2
  do i=1,ao_num
    x = mo_coef(i,i1)
    mo_coef(i,i1) = mo_coef(i,i2)
    mo_coef(i,i2) = x
  enddo
  call save_mos
  
end
