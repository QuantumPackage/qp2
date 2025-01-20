double precision function coul_full_ao_pq_r_1s(p,q,R,R_p,R_q)
 implicit none
 include 'constants.include.F'
 BEGIN_DOC
 ! coul_full_pq_r_1s(p,q,r) = \int d^3r phi_p(r) phi_q(r) 1/(r-R)
 !
 ! where phi_q and phi_p are centered in R_q and R_p. 
 ! 
 ! WARNING :: works only for purely 1s extra basis  !!
 END_DOC
 double precision, intent(in) :: R(3),R_p(3),R_q(3)
 integer, intent(in) :: p,q
 double precision :: coef,dist,P_pq(3),coefaos
 coefaos= ao_extra_coef_normalized(p,1) * ao_extra_coef_normalized(q,1) 
 coef = inv_pi_gamma_pq_3_2_ao_extra(p,q) * E_pq_ao_extra(p,q)
 P_pq = ao_extra_expo(p,1) * R_p + ao_extra_expo(q,1) * R_q
 P_pq = P_pq * inv_gamma_pq_ao_extra(q,p)
 dist = (P_pq(1)-R(1)) * (P_pq(1)-R(1))
 dist+= (P_pq(2)-R(2)) * (P_pq(2)-R(2))
 dist+= (P_pq(3)-R(3)) * (P_pq(3)-R(3))
 dist = dsqrt(dist)
 if(dist.gt.1.d-10)then
  coul_full_ao_pq_r_1s = coefaos * coef * derf(sqrt_gamma_pq_ao_extra(q,p) * dist)/dist
 else
  coul_full_ao_pq_r_1s = coefaos * coef * 2.d0 * sqrt_gamma_pq_ao_extra(q,p) * inv_sq_pi
 endif

end

double precision function coul_pq_r_1s(p,q,R,R_p,R_q)
 implicit none
 include 'constants.include.F'
 BEGIN_DOC
 ! coul_full_pq_r_1s(p,q,r) = \int d^3r exp(-alpha_p (r-R_p)^2) exp(-alpha_q (r-R_q)^2) 1/|r-R|
 !
 ! where alpha_p and alpha_q are the 1s extra basis
 ! 
 ! WARNING :: works only for purely 1s extra basis  !!
 END_DOC
 double precision, intent(in) :: R(3),R_p(3),R_q(3)
 integer, intent(in) :: p,q
 double precision :: dist,P_pq(3)
 P_pq = ao_extra_expo(p,1) * R_p + ao_extra_expo(q,1) * R_q
 P_pq = P_pq * inv_gamma_pq_ao_extra(q,p)
 dist = (P_pq(1)-R(1)) * (P_pq(1)-R(1))
 dist+= (P_pq(2)-R(2)) * (P_pq(2)-R(2))
 dist+= (P_pq(3)-R(3)) * (P_pq(3)-R(3))
 dist = dsqrt(dist)
 if(dist.gt.1.d-10)then
  coul_pq_r_1s = derf(sqrt_gamma_pq_ao_extra(q,p) * dist)/dist
 else
  coul_pq_r_1s = 2.d0 * sqrt_gamma_pq_ao_extra(q,p) * inv_sq_pi
 endif

end
