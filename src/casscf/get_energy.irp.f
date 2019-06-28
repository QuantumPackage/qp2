program print_2rdm
 implicit none
 read_wf = .True.
 touch read_wf
 integer :: i,j,k,l
 double precision :: accu(4),twodm,thr,act_twodm2,integral,get_two_e_integral
 thr = 1.d-10

 accu = 0.d0
 do l = 1, mo_num
  do k = 1, mo_num
   do j = 1, mo_num
    do i = 1, mo_num
     integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
     accu(1) += act_two_rdm_spin_trace_mo(i,j,k,l) * integral
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu(1)
end
