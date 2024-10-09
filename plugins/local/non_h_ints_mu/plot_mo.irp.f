program plot_mo
 implicit none
 integer :: i,npt
 double precision :: xmin,xmax,dx,r(3)
 double precision,allocatable :: mos_array(:)
 allocate(mos_array(mo_num))
 npt = 10000
 xmin =0.d0
 xmax =10.d0
 dx=(xmax-xmin)/dble(npt)
 r=0.d0
 r(1) = xmin
 do i = 1, npt
  call give_all_mos_at_r(r,mos_array)
  write(33,'(100(F16.10,X))')r(1),mos_array(1),mos_array(2),mos_array(3)
  r(1) += dx
 enddo

end
