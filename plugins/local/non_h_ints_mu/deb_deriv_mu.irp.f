program test_j_mu_of_r
 implicit none
 double precision :: x,mu_min,dmu,mu_max, mu, mu_p, mu_m
 double precision :: j_simple,j_p, j_m,numeric_d_mu,d_dx_mu
 double precision :: accu
 integer :: npt,i
 npt = 1000
 mu_min = 0.3d0
 mu_max = 10.d0
 dmu = (mu_max - mu_min)/dble(npt)
 x = 0.7d0
 mu = mu_min
 do i = 1, npt
  call get_deriv_mu_j12(x,mu,d_dx_mu)
  mu_p = mu + dmu
  mu_m = mu - dmu
  j_p = j_simple(x,mu_p)
  j_m = j_simple(x,mu_m)
  numeric_d_mu = 0.5d0 * (j_p - j_m)/dmu
  print*,mu
  print*,numeric_d_mu,d_dx_mu,dabs(d_dx_mu-numeric_d_mu)
  accu += dabs(d_dx_mu-numeric_d_mu)
  mu += dmu
 enddo
 accu *= dmu
 print*,'accu = ',accu
end

