program print_mos
 implicit none
 integer :: i,nx
 double precision :: r(3), xmax, dx, accu
 double precision, allocatable :: mos_array(:)
 double precision:: alpha,envelop,dm_a,dm_b
 allocate(mos_array(mo_num))
 xmax = 5.d0
 nx = 1000
 dx=xmax/dble(nx)
 r = 0.d0
 alpha = 0.5d0
 do i = 1, nx
  call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
  call give_all_mos_at_r(r,mos_array)
  accu = mos_array(3)**2+mos_array(4)**2+mos_array(5)**2
  accu = dsqrt(accu)
  envelop = (1.d0 - dexp(-alpha * r(3)**2))
  write(33,'(100(F16.10,X))')r(3), mos_array(1), mos_array(2), accu, dm_a+dm_b, envelop
  r(3) += dx
 enddo

end

