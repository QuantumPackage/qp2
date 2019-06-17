program rotate_mos
  implicit none
  BEGIN_DOC
  ! Rotates molecular orbitals i and j by combining them as
  ! $1/\sqrt{2} ( \phi_i + \phi_j )$ and
  ! $1/\sqrt{2} ( \phi_i - \phi_j )$.
  END_DOC
  integer                        :: iorb,jorb
  integer                        :: i,j
  double precision               :: dsqrt2_inv
  double precision, allocatable  :: mo_coef_tmp(:,:)

  read(5,*)iorb,jorb

  allocate(mo_coef_tmp(ao_num,mo_num))
  mo_coef_tmp = mo_coef

  dsqrt2_inv = 1.d0/dsqrt(2.d0)
  do i = 1, ao_num
    mo_coef(i,iorb) = dsqrt2_inv * ( mo_coef_tmp(i,iorb) + mo_coef_tmp(i,jorb) )
    mo_coef(i,jorb) = dsqrt2_inv * ( mo_coef_tmp(i,iorb) - mo_coef_tmp(i,jorb) )
  enddo

  touch mo_coef
  call save_mos
  
end
