program swap_mos
  implicit none
  BEGIN_DOC
  ! Swaps the indices of two molecular orbitals
  END_DOC
  integer                        :: i,j, i1, i2
  double precision               :: x
  print *,  'MOs to swap?'
  read(*,*) i1, i2
  if (is_complex) then
    complex*16 :: xc
    do i=1,ao_num
      xc = mo_coef_complex(i,i1)
      mo_coef_complex(i,i1) = mo_coef_complex(i,i2)
      mo_coef_complex(i,i2) = xc
    enddo
  else
    do i=1,ao_num
      x = mo_coef(i,i1)
      mo_coef(i,i1) = mo_coef(i,i2)
      mo_coef(i,i2) = x
    enddo
  endif
  call save_mos
  
end
