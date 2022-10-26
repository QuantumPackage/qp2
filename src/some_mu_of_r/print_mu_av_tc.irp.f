program print_mu_av_tc
 implicit none
 read_wf = .True.
 touch read_wf
   print*,'average_mu_lda       = ',average_mu_lda
!   print*,'average_mu_rs        = ',average_mu_rs 
   print*,'average_mu_rs_c      = ',average_mu_rs_c
   print*,'average_mu_rs_c_lda  = ',average_mu_rs_c_lda
   call plot_mu_tc
end



subroutine plot_mu_tc
 implicit none
 integer :: i,nx,m
 double precision :: xmin,xmax,dx
 double precision :: r(3)
 double precision :: sqpi
 double precision :: rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n,grad_n_a(3), grad_n_b(3)
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs,mu_rs, mu_rs_c, mu_lda, mu_grad_n 
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 if (no_core_density)then
  output=trim(ezfio_filename)//'.tc_mu_fc'
 else
  output=trim(ezfio_filename)//'.tc_mu'
 endif
 i_unit_output = getUnitAndOpen(output,'w')

 nx = 5000
 xmin = -10.d0
 xmax = 10.d0
 dx   = (xmax - xmin)/dble(nx)
 r = 0.d0
 r(3) = xmin

 sqpi = dsqrt(dacos(-1.d0))
 cst_rs   = (4.d0 * dacos(-1.d0)/3.d0)**(-1.d0/3.d0)
 alpha_rs = 2.d0 * dsqrt((9.d0 * dacos(-1.d0)/4.d0)**(-1.d0/3.d0)) / sqpi

  write(i_unit_output,*)'#       r(3)              rho_hf          mu_rs           mu_rs_c           mu_lda         mu_grad_n'
 do i = 1, nx
  rho_a_hf = 0.d0
  grad_n   = 0.d0
  call density_and_grad_alpha_beta(r,rho_a_hf,rho_b_hf, grad_n_a, grad_n_b)
  do m = 1, 3
   grad_n += grad_n_a(m)*grad_n_a(m) + grad_n_b(m)*grad_n_b(m) + 2.d0 * grad_n_a(m) * grad_n_b(m)
  enddo
  rho_hf = rho_a_hf + rho_b_hf
  grad_n = dsqrt(grad_n)
  grad_n = grad_n/(4.d0 * rho_hf)
  rs = cst_rs * rho_hf**(-1.d0/3.d0)
  g0 = g0_UEG_mu_inf(rho_a_hf,rho_b_hf)

  mu_rs     = 1.d0/rs 
  mu_rs_c   = alpha_rs/dsqrt(rs) 
  mu_lda    =  - 1.d0 / (dlog(2.d0 * g0) * sqpi)  
  mu_grad_n = grad_n 
  write(i_unit_output,'(100(F16.10,X))')r(3),rho_hf,mu_rs,mu_rs_c,mu_lda,mu_grad_n
  r(3) += dx
 enddo


end
