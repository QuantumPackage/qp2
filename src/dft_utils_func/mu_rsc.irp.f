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

