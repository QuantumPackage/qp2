double precision function mu_rs_c(rho)
 implicit none
 double precision, intent(in) :: rho
 include 'constants.include.F'
 double precision :: cst_rs,alpha_rs,rs
 cst_rs   = (4.d0 * dacos(-1.d0)/3.d0)**(-1.d0/3.d0)
 alpha_rs = 2.d0 * dsqrt((9.d0 * dacos(-1.d0)/4.d0)**(-1.d0/3.d0)) / sqpi
  
 rs = cst_rs * rho**(-1.d0/3.d0)
 mu_rs_c =  alpha_rs/dsqrt(rs)

end

double precision function mu_grad_rho_func(r)
 implicit none
 double precision , intent(in) :: r(3)
 integer :: m
 double precision :: rho, dm_a, dm_b, grad_dm_a(3), grad_dm_b(3)
 double precision :: eta, grad_rho(3), grad_sqr
 eta = 0.135d0
  call density_and_grad_alpha_beta(r,dm_a,dm_b, grad_dm_a, grad_dm_b)  
  rho = dm_a + dm_b
  do m = 1,3
   grad_rho(m) = grad_dm_a(m) + grad_dm_b(m)
  enddo
  grad_sqr=0.d0
  do m = 1,3
   grad_sqr=grad_sqr+grad_rho(m)*grad_rho(m)
  enddo
  grad_sqr = dsqrt(grad_sqr)
  if (rho<1.d-12) then
   mu_grad_rho_func = 1.d-10
  else
   mu_grad_rho_func = eta * grad_sqr / rho
  endif
  
end
