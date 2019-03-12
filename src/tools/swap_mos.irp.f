program swap_mos
  implicit none
  integer :: i,j, i1, i2
  double precision :: x
  print *,  'MOs to swap?'
  read(*,*) i1, i2
  do i=1,ao_num
    x = mo_coef(i,i1)
    mo_coef(i,i1) = mo_coef(i,i2)
    mo_coef(i,i2) = x
  enddo
  call save_mos

end
