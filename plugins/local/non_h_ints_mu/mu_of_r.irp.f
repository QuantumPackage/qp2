
subroutine grad_mu_of_r_mean_field(r,mu_mf, dm, grad_mu_mf, grad_dm)
 implicit none
 BEGIN_DOC
 ! returns the value and gradients of the mu(r) mean field, together with the HF density and its gradients. 
 END_DOC
 include 'constants.include.F'
 double precision, intent(in) :: r(3)
 double precision, intent(out):: grad_mu_mf(3), grad_dm(3)
 double precision, intent(out):: mu_mf, dm
 double precision :: grad_f_mf_ab(3), grad_two_bod_dens(3),grad_dm_a(3), grad_dm_b(3)
 double precision :: f_mf_ab,two_bod_dens, dm_a, dm_b
 
 double precision :: dist
 call get_grad_f_mf_ab(r,grad_f_mf_ab, grad_two_bod_dens,f_mf_ab,two_bod_dens, dm_a, dm_b,grad_dm_a, grad_dm_b)
 
 dm = dm_a + dm_b
 grad_dm(1:3) = grad_dm_a(1:3) + grad_dm_b(1:3)

 if(dabs(two_bod_dens).lt.1.d-10)then
  mu_mf = 1.d+10
  grad_mu_mf = 0.d0
 else
  if(mu_of_r_tc=="Erfmu")then
   mu_mf           = 0.3333333333d0 * sqpi * (f_mf_ab/two_bod_dens + 0.25d0)
   grad_mu_mf(1:3) = 0.3333333333d0 * sqpi * (grad_f_mf_ab(1:3) * two_bod_dens - f_mf_ab * grad_two_bod_dens(1:3))& 
                                    /(two_bod_dens*two_bod_dens)
  else if(mu_of_r_tc=="Standard")then
   mu_mf = 0.5d0 * sqpi * f_mf_ab/two_bod_dens
   grad_mu_mf(1:3) = 0.5d0 * sqpi * (grad_f_mf_ab(1:3) * two_bod_dens - f_mf_ab * grad_two_bod_dens(1:3))& 
                                    /(two_bod_dens*two_bod_dens)
  else if(mu_of_r_tc=="Erfmugauss")then
   mu_mf = (f_mf_ab/two_bod_dens + 0.25d0)/c_mu_gauss_tot 
   grad_mu_mf(1:3) = 1.d0/c_mu_gauss_tot* (grad_f_mf_ab(1:3) * two_bod_dens - f_mf_ab * grad_two_bod_dens(1:3))& 
                                    /(two_bod_dens*two_bod_dens)
  else 
   print*,'Wrong value for mu_of_r_tc !'
   stop
  endif
 endif 

end

