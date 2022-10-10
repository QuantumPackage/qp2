program test_non_h
 implicit none
  my_grid_becke  = .True.
  my_n_pt_r_grid = 50
  my_n_pt_a_grid = 74
!  my_n_pt_r_grid = 10 ! small grid for quick debug
!  my_n_pt_a_grid = 26 ! small grid for quick debug
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
!call routine_grad_squared
 call routine_fit
end

subroutine routine_lapl_grad
 implicit none
 integer :: i,j,k,l
 double precision :: grad_lapl, get_ao_tc_sym_two_e_pot,new,accu,contrib
 double precision :: ao_two_e_integral_erf,get_ao_two_e_integral,count_n,accu_relat
! !!!!!!!!!!!!!!!!!!!!! WARNING
! THIS ROUTINE MAKES SENSE ONLY IF HAND MODIFIED coef_gauss_eff_pot(1:n_max_fit_slat) = 0. to cancel (1-erf(mu*r12))^2
 accu = 0.d0
 accu_relat = 0.d0
 count_n = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     grad_lapl  = get_ao_tc_sym_two_e_pot(i,j,k,l,ao_tc_sym_two_e_pot_map) ! pure gaussian part : comes from Lapl
     grad_lapl += ao_two_e_integral_erf(i, k, j, l)                        ! erf(mu r12)/r12    : comes from Lapl
     grad_lapl += ao_non_hermit_term_chemist(k,i,l,j)                      ! \grad u(r12) . grad
     new        = tc_grad_and_lapl_ao(k,i,l,j)
     new       += get_ao_two_e_integral(i,j,k,l,ao_integrals_map)
     contrib    = dabs(new - grad_lapl)
     if(dabs(grad_lapl).gt.1.d-12)then
      count_n += 1.d0
      accu_relat += 2.0d0 * contrib/dabs(grad_lapl+new)
     endif
     if(contrib.gt.1.d-10)then
      print*,i,j,k,l
      print*,grad_lapl,new,contrib
      print*,2.0d0*contrib/dabs(grad_lapl+new+1.d-12)
     endif 
     accu += contrib
    enddo
   enddo
  enddo
 enddo
 print*,'accu      = ',accu/count_n
 print*,'accu/rel  = ',accu_relat/count_n

end

subroutine routine_grad_squared
 implicit none
 integer :: i,j,k,l
 double precision :: grad_squared, get_ao_tc_sym_two_e_pot,new,accu,contrib
 double precision :: count_n,accu_relat
! !!!!!!!!!!!!!!!!!!!!! WARNING
! THIS ROUTINE MAKES SENSE ONLY IF HAND MODIFIED coef_gauss_eff_pot(n_max_fit_slat:n_max_fit_slat+1) = 0. to cancel exp(-'mu*r12)^2)
 accu = 0.d0
 accu_relat = 0.d0
 count_n = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     grad_squared  = get_ao_tc_sym_two_e_pot(i,j,k,l,ao_tc_sym_two_e_pot_map) ! pure gaussian part : comes from Lapl
     new        = tc_grad_square_ao(k,i,l,j)
     contrib    = dabs(new - grad_squared)
     if(dabs(grad_squared).gt.1.d-12)then
      count_n += 1.d0
      accu_relat += 2.0d0 * contrib/dabs(grad_squared+new)
     endif
     if(contrib.gt.1.d-10)then
      print*,i,j,k,l
      print*,grad_squared,new,contrib
      print*,2.0d0*contrib/dabs(grad_squared+new+1.d-12)
     endif 
     accu += contrib
    enddo
   enddo
  enddo
 enddo
 print*,'accu      = ',accu/count_n
 print*,'accu/rel  = ',accu_relat/count_n

end

subroutine routine_fit
 implicit none
 integer :: i,nx
 double precision :: dx,xmax,x,j_mu,j_mu_F_x_j,j_mu_fit_gauss
 nx = 500
 xmax = 5.d0
 dx = xmax/dble(nx)
 x = 0.d0
 print*,'coucou',mu_erf
 do i = 1, nx
  write(33,'(100(F16.10,X))') x,j_mu(x),j_mu_F_x_j(x),j_mu_fit_gauss(x)
  x += dx
 enddo

end
