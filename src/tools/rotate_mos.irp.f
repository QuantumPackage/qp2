program rotate_mos
 implicit none
 integer :: iorb,jorb
 read(5,*)iorb,jorb
 double precision, allocatable :: mo_coef_tmp(:,:)
 allocate(mo_coef_tmp(ao_num,mo_num))
 mo_coef_tmp = mo_coef
 integer :: i,j
 double precision :: dsqrt2_inv
 dsqrt2_inv = 1.d0/dsqrt(2.d0)
 do i = 1, ao_num
  mo_coef(i,iorb) = dsqrt2_inv * ( mo_coef_tmp(i,iorb) + mo_coef_tmp(i,jorb) )
  mo_coef(i,jorb) = dsqrt2_inv * ( mo_coef_tmp(i,iorb) - mo_coef_tmp(i,jorb) )
 enddo
 touch mo_coef 
 call save_mos

end
