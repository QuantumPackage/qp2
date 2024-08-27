subroutine get_j_sum_mu_of_r(r1,r2,jast)
 implicit none
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: jast
 double precision :: mu_r1, dm_r1, grad_mu_r1(3), grad_dm_r1(3), j_mu_r1
 double precision :: mu_r2, dm_r2, grad_mu_r2(3), grad_dm_r2(3), j_mu_r2
 double precision :: j12_mu_input
 jast = 0.d0
 call grad_mu_of_r_mean_field(r1,mu_r1, dm_r1, grad_mu_r1, grad_dm_r1) 
 call grad_mu_of_r_mean_field(r2,mu_r2, dm_r2, grad_mu_r2, grad_dm_r2) 
 j_mu_r1 = j12_mu_input(r1, r2, mu_r1)
 j_mu_r2 = j12_mu_input(r1, r2, mu_r2)
 if(dm_r1 + dm_r2.lt.1.d-7)return
 jast = (dm_r1 * j_mu_r1 + dm_r2 * j_mu_r2) / (dm_r1 + dm_r2)

end

subroutine grad_j_sum_mu_of_r(r1,r2,jast,grad_jast)
 implicit none
 BEGIN_DOC
 END_DOC
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: jast, grad_jast(3)
 jast = 0.d0
 grad_jast = 0.d0
 double precision :: num, denom, grad_num(3), grad_denom(3)
 double precision :: j_r1, grad_j_r1(3),j_r2, grad_j_r2(3)
 double precision :: dm_r1, grad_dm_r1(3), grad_jmu_r2(3)
 double precision :: dm_r2, grad_dm_r2(3),mu_r2, grad_mu_r2(3)
 double precision :: j12_mu_input

 if(murho_type==0)then
! J(r1,r2) = [rho(r1) * j(mu(r1),r12) + rho(r2) * j(mu(r2),r12)] / [rho(r1) + rho(r2)]
! 
!          = num(r1,r2) / denom(r1,r2)
!
! d/dx1 J(r1,r2) = [denom(r1,r2) X d/dx1 num(r1,r2) - num(r1,r2) X d/dx1 denom(r1,r2) ] / denom(r1,r2)^2
!
! d/dx1 num(r1,r2) =  j(mu(r1),r12)*d/dx1 rho(r1) + rho(r1) * d/dx1 j(mu(r1),r12) 
!                   + rho(r2) d/dx1 j(mu(r2),r12)
! d/dx1 denom(r1,r2) = d/dx1 rho(r1)
  call grad_j_mu_of_r_1(r1,r2,j_r1, grad_j_r1,dm_r1, grad_dm_r1)
  call grad_mu_of_r_mean_field(r2,mu_r2, dm_r2, grad_mu_r2, grad_dm_r2) 
  j_r2 = j12_mu_input(r1, r2, mu_r2) ! j(mu(r2),r1,r2)
  num = dm_r1 * j_r1 + dm_r2 * j_r2
  denom = dm_r1 + dm_r2
  if(denom.lt.1.d-7)return
  jast = num / denom
 
  grad_denom = grad_dm_r1
 
  call grad_j12_mu_input(r1, r2, mu_r2, grad_jmu_r2)
  grad_num =  j_r1 * grad_dm_r1 + dm_r1 * grad_j_r1 + dm_r2 * grad_jmu_r2
  
  grad_jast = (grad_num * denom - num * grad_denom)/(denom*denom)
 else if(murho_type==-1)then
! J(r1,r2) = 0.5 * [j(mu(r1),r12) + j(mu(r2),r12)] 
!
! d/dx1 J(r1,r2)  = 0.5 * (d/dx1 j(mu(r1),r12) + d/dx1 j(mu(r2),r12))
  call grad_j_mu_of_r_1(r1,r2,j_r1, grad_j_r1,dm_r1, grad_dm_r1)
  call grad_mu_of_r_mean_field(r2,mu_r2, dm_r2, grad_mu_r2, grad_dm_r2) 
  j_r2 = j12_mu_input(r1, r2, mu_r2) ! j(mu(r2),r1,r2)
  call grad_j12_mu_input(r1, r2, mu_r2, grad_jmu_r2)
  jast = 0.5d0 * (j_r1 + j_r2)
  grad_jast = 0.5d0 * (grad_j_r1 + grad_jmu_r2)
  
 endif

end

subroutine grad_j_mu_of_r_1(r1,r2,jast, grad_jast, dm_r1, grad_dm_r1)
 implicit none
  include 'constants.include.F'
 BEGIN_DOC
! grad_r1 of j(mu(r1),r12)
  !
  ! d/dx1 j(mu(r1),r12) = exp(-(mu(r1)*r12)**2) /(2 *sqrt(pi) * mu(r1)**2 ) d/dx1 mu(r1) 
  !                     + d/dx1 j(mu(r1),r12) 
  !
  ! with 
  !
  !           j(mu,r12) = 1/2 r12 (1 - erf(mu r12)) - 1/2 (sqrt(pi) * mu) e^{-(mu*r12)^2}
  !
  ! and d/dx1 j(mu,r12) = 0.5 * (1 - erf(mu *r12))/r12 * (x1 - x2)
  !
  ! here mu(r1) is obtained by MU MEAN FIELD 
 END_DOC
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: jast, grad_jast(3),dm_r1, grad_dm_r1(3)
 double precision :: dx, dy, dz,  r12, mu_der(3)
 double precision :: mu_tmp, tmp, grad(3), mu_val
 jast = 0.d0
 grad = 0.d0

 dx  = r1(1) - r2(1)
 dy  = r1(2) - r2(2)
 dz  = r1(3) - r2(3)
 r12 = dsqrt(dx * dx + dy * dy + dz * dz)
 ! get mu(r1) == mu_val and its gradient d/dx1 mu(r1) == mu_der 
 call grad_mu_of_r_mean_field(r1,mu_val, dm_r1, mu_der, grad_dm_r1) 
 mu_tmp  = mu_val * r12
 ! evalulation of the jastrow j(mu(r1),r12)
 jast = 0.5d0 * r12 * (1.d0 - derf(mu_tmp)) - inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / mu_val

 ! tmp = exp(-(mu(r1)*r12)**2) /(2 *sqrt(pi) * mu(r1)**2 )
 tmp     = inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / (mu_val * mu_val)
 ! grad = 
 grad(1) = tmp * mu_der(1)
 grad(2) = tmp * mu_der(2)
 grad(3) = tmp * mu_der(3)

 if(r12 .lt. 1d-10) return
 tmp     = 0.5d0 * (1.d0 - derf(mu_tmp)) / r12 ! d/dx1 j(mu(r1),r12) 
 grad(1) = grad(1) + tmp * dx
 grad(2) = grad(2) + tmp * dy
 grad(3) = grad(3) + tmp * dz

 grad_jast = grad 
end

! ---

double precision function j12_mu_input(r1, r2, mu)

  BEGIN_DOC
  ! j(mu,r12) = 1/2 r12 (1 - erf(mu r12)) - 1/2 (sqrt(pi) * mu) e^{-(mu*r12)^2}
  END_DOC
  include 'constants.include.F'

  implicit none
  double precision, intent(in) :: r1(3), r2(3), mu
  double precision             :: mu_tmp, r12

    r12 = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
               + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
               + (r1(3) - r2(3)) * (r1(3) - r2(3)) )
    mu_tmp = mu * r12

    j12_mu_input = 0.5d0 * r12 * (1.d0 - derf(mu_tmp)) - inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / mu

end

subroutine grad_j12_mu_input(r1, r2, mu, grad_jmu)
 implicit none
 BEGIN_DOC
 ! grad_jmu = 
 END_DOC
 double precision, intent(in) :: r1(3), r2(3), mu
 double precision, intent(out):: grad_jmu(3)
 double precision             :: mu_tmp, r12, dx, dy, dz, grad(3), tmp
 grad_jmu = 0.d0
 dx  = r1(1) - r2(1)
 dy  = r1(2) - r2(2)
 dz  = r1(3) - r2(3)
 r12 = dsqrt(dx * dx + dy * dy + dz * dz)
 if(r12 .lt. 1d-10) return
 mu_tmp  = mu * r12
 tmp     = 0.5d0 * (1.d0 - derf(mu_tmp)) / r12 ! d/dx1 j(mu(r1),r12) 
 grad(1) = tmp * dx
 grad(2) = tmp * dy
 grad(3) = tmp * dz

 grad_jmu = grad 
end
