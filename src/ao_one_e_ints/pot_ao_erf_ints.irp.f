subroutine give_all_erf_kl_ao(integrals_ao,mu_in,C_center)
  implicit none
  BEGIN_DOC
  ! Subroutine that returns all integrals over $r$ of type
  ! $\frac{ \erf(\mu * |r-R_C|) }{ |r-R_C| }$
  END_DOC
  double precision, intent(in)   :: mu_in,C_center(3)
  double precision, intent(out)  :: integrals_ao(ao_num,ao_num)
  double precision               :: NAI_pol_mult_erf_ao
  integer                        :: i,j,l,k,m
  do k = 1, ao_num
    do m = 1, ao_num
      integrals_ao(m,k) = NAI_pol_mult_erf_ao(m,k,mu_in,C_center)
    enddo
  enddo
end


double precision function NAI_pol_mult_erf_ao(i_ao,j_ao,mu_in,C_center)
  implicit none
  BEGIN_DOC
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) \frac{\erf(\mu |r-R_C|)}{|r-R_C|}$.
  END_DOC
  integer, intent(in)            :: i_ao,j_ao
  double precision, intent(in)   :: mu_in, C_center(3)
  integer                        :: i,j,num_A,num_B, power_A(3), power_B(3), n_pt_in
  double precision               :: A_center(3), B_center(3),integral, alpha,beta
  double precision               :: NAI_pol_mult_erf
  num_A = ao_nucl(i_ao)
  power_A(1:3)= ao_power(i_ao,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)
  num_B = ao_nucl(j_ao)
  power_B(1:3)= ao_power(j_ao,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)
  n_pt_in = n_pt_max_integrals
  NAI_pol_mult_erf_ao = 0.d0
  do i = 1, ao_prim_num(i_ao)
    alpha = ao_expo_ordered_transp(i,i_ao)
    do j = 1, ao_prim_num(j_ao)
      beta = ao_expo_ordered_transp(j,j_ao)
      integral =  NAI_pol_mult_erf(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in,mu_in)
      NAI_pol_mult_erf_ao += integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
    enddo
  enddo

end



double precision function NAI_pol_mult_erf(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in,mu_in)
  BEGIN_DOC
  ! Computes the following integral :
  ! $\int dr (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !  \frac{\erf(\mu |r-R_C|)}{|r-R_c|}$.
  END_DOC

  implicit none
  integer, intent(in)            :: n_pt_in
  double precision,intent(in)    :: C_center(3),A_center(3),B_center(3),alpha,beta,mu_in
  integer, intent(in)            :: power_A(3),power_B(3)
  integer                        :: i,j,k,l,n_pt
  double precision               :: P_center(3)

  double precision               :: d(0:n_pt_in),pouet,coeff,dist,const,pouet_2,factor
  double precision               :: I_n_special_exact,integrate_bourrin,I_n_bibi
  double precision               :: V_e_n,const_factor,dist_integral,tmp
  double precision               :: accu,rint,p_inv,p,rho,p_inv_2
  integer                        :: n_pt_out,lmax
  include 'utils/constants.include.F'
  p = alpha + beta
  p_inv = 1.d0/p
  p_inv_2 = 0.5d0 * p_inv
  rho = alpha * beta * p_inv

  dist = 0.d0
  dist_integral = 0.d0
  do i = 1, 3
    P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
    dist += (A_center(i) - B_center(i))*(A_center(i) - B_center(i))
    dist_integral += (P_center(i) - C_center(i))*(P_center(i) - C_center(i))
  enddo
  const_factor = dist*rho
  if(const_factor > 80.d0)then
    NAI_pol_mult_erf = 0.d0
    return
  endif
  double precision               :: p_new
  p_new = mu_in/dsqrt(p+ mu_in * mu_in)
  factor = dexp(-const_factor)
  coeff = dtwo_pi * factor * p_inv * p_new
  lmax = 20

  !  print*, "b"
  do i = 0, n_pt_in
    d(i) = 0.d0
  enddo
  n_pt =  2 * ( (power_A(1) + power_B(1)) +(power_A(2) + power_B(2)) +(power_A(3) + power_B(3)) )
  const = p * dist_integral * p_new * p_new
  if (n_pt == 0) then
    pouet = rint(0,const)
    NAI_pol_mult_erf = coeff * pouet
    return
  endif

  ! call give_polynomial_mult_center_one_e_erf(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out,mu_in)
  p_new = p_new * p_new
  call give_polynomial_mult_center_one_e_erf_opt(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out,mu_in,p,p_inv,p_inv_2,p_new,P_center)


  if(n_pt_out<0)then
    NAI_pol_mult_erf = 0.d0
    return
  endif
  accu = 0.d0

  ! sum of integrals of type : int {t,[0,1]}  exp-(rho.(P-Q)^2 * t^2) * t^i
  do i =0 ,n_pt_out,2
    accu +=  d(i) * rint(i/2,const)
  enddo
  NAI_pol_mult_erf = accu * coeff

end


subroutine give_polynomial_mult_center_one_e_erf_opt(A_center,B_center,alpha,beta,&
      power_A,power_B,C_center,n_pt_in,d,n_pt_out,mu_in,p,p_inv,p_inv_2,p_new,P_center)
  BEGIN_DOC
  ! Returns the explicit polynomial in terms of the $t$ variable of the following polynomial:
  !
  ! $I_{x1}(a_x, d_x,p,q) \times I_{x1}(a_y, d_y,p,q) \times I_{x1}(a_z, d_z,p,q)$.
  END_DOC
  implicit none
  integer, intent(in)            :: n_pt_in
  integer,intent(out)            :: n_pt_out
  double precision, intent(in)   :: A_center(3), B_center(3),C_center(3),p,p_inv,p_inv_2,p_new,P_center(3)
  double precision, intent(in)   :: alpha,beta,mu_in
  integer, intent(in)            :: power_A(3), power_B(3)
  integer                        :: a_x,b_x,a_y,b_y,a_z,b_z
  double precision               :: d(0:n_pt_in)
  double precision               :: d1(0:n_pt_in)
  double precision               :: d2(0:n_pt_in)
  double precision               :: d3(0:n_pt_in)
  double precision               :: accu
  accu = 0.d0
  ASSERT (n_pt_in > 1)

  double precision               :: R1x(0:2), B01(0:2), R1xp(0:2),R2x(0:2)
  R1x(0)  = (P_center(1) - A_center(1))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(1) - C_center(1))* p_new
  ! R1x = (P_x - A_x) - (P_x - C_x) ( t * mu/sqrt(p+mu^2) )^2
  R1xp(0)  = (P_center(1) - B_center(1))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(1) - C_center(1))* p_new
  !R1xp = (P_x - B_x) - (P_x - C_x) ( t * mu/sqrt(p+mu^2) )^2
  R2x(0)  =  p_inv_2
  R2x(1)  = 0.d0
  R2x(2)  = -p_inv_2* p_new
  !R2x  = 0.5 / p - 0.5/p ( t * mu/sqrt(p+mu^2) )^2
  do i = 0,n_pt_in
    d(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d1(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d2(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d3(i) = 0.d0
  enddo
  integer                        :: n_pt1,n_pt2,n_pt3,dim,i
  n_pt1 = n_pt_in
  n_pt2 = n_pt_in
  n_pt3 = n_pt_in
  a_x = power_A(1)
  b_x = power_B(1)
  call I_x1_pol_mult_one_e(a_x,b_x,R1x,R1xp,R2x,d1,n_pt1,n_pt_in)
  if(n_pt1<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif

  R1x(0)  = (P_center(2) - A_center(2))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(2) - C_center(2))* p_new
  ! R1x = (P_x - A_x) - (P_x - C_x) ( t * mu/sqrt(p+mu^2) )^2
  R1xp(0)  = (P_center(2) - B_center(2))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(2) - C_center(2))* p_new
  !R1xp = (P_x - B_x) - (P_x - C_x) ( t * mu/sqrt(p+mu^2) )^2
  a_y = power_A(2)
  b_y = power_B(2)
  call I_x1_pol_mult_one_e(a_y,b_y,R1x,R1xp,R2x,d2,n_pt2,n_pt_in)
  if(n_pt2<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif


  R1x(0)  = (P_center(3) - A_center(3))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(3) - C_center(3))* p_new
  ! R1x = (P_x - A_x) - (P_x - C_x) ( t * mu/sqrt(p+mu^2) )^2
  R1xp(0)  = (P_center(3) - B_center(3))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(3) - C_center(3))* p_new
  !R2x  = 0.5 / p - 0.5/p ( t * mu/sqrt(p+mu^2) )^2
  a_z = power_A(3)
  b_z = power_B(3)

  call I_x1_pol_mult_one_e(a_z,b_z,R1x,R1xp,R2x,d3,n_pt3,n_pt_in)
  if(n_pt3<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif
  integer                        :: n_pt_tmp
  n_pt_tmp = 0
  call multiply_poly(d1,n_pt1,d2,n_pt2,d,n_pt_tmp)
  do i = 0,n_pt_tmp
    d1(i) = 0.d0
  enddo
  n_pt_out = 0
  call multiply_poly(d ,n_pt_tmp ,d3,n_pt3,d1,n_pt_out)
  do i = 0, n_pt_out
    d(i) = d1(i)
  enddo

end




subroutine give_polynomial_mult_center_one_e_erf(A_center,B_center,alpha,beta,&
      power_A,power_B,C_center,n_pt_in,d,n_pt_out,mu_in)
  BEGIN_DOC
  ! Returns the explicit polynomial in terms of the $t$ variable of the following polynomial:
  !
  ! $I_{x1}(a_x, d_x,p,q) \times I_{x1}(a_y, d_y,p,q) \times I_{x1}(a_z, d_z,p,q)$.
  END_DOC
  implicit none
  integer, intent(in)            :: n_pt_in
  integer,intent(out)            :: n_pt_out
  double precision, intent(in)   :: A_center(3), B_center(3),C_center(3)
  double precision, intent(in)   :: alpha,beta,mu_in
  integer, intent(in)            :: power_A(3), power_B(3)
  integer                        :: a_x,b_x,a_y,b_y,a_z,b_z
  double precision               :: d(0:n_pt_in)
  double precision               :: d1(0:n_pt_in)
  double precision               :: d2(0:n_pt_in)
  double precision               :: d3(0:n_pt_in)
  double precision               :: accu,  pq_inv, p10_1, p10_2, p01_1, p01_2
  double precision               :: p,P_center(3),rho,p_inv,p_inv_2
  accu = 0.d0
  !COMPTEUR irp_rdtsc1 = irp_rdtsc()
  ASSERT (n_pt_in > 1)
  p = alpha+beta
  p_inv = 1.d0/p
  p_inv_2 = 0.5d0/p
  do i =1, 3
    P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
  enddo

  double precision               :: R1x(0:2), B01(0:2), R1xp(0:2),R2x(0:2)
  R1x(0)  = (P_center(1) - A_center(1))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(1) - C_center(1))* mu_in**2 / (p+mu_in*mu_in)
  ! R1x = (P_x - A_x) - (P_x - C_x) ( t * mu/sqrt(p+mu^2) )^2
  R1xp(0)  = (P_center(1) - B_center(1))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(1) - C_center(1))* mu_in**2 / (p+mu_in*mu_in)
  !R1xp = (P_x - B_x) - (P_x - C_x) ( t * mu/sqrt(p+mu^2) )^2
  R2x(0)  =  p_inv_2
  R2x(1)  = 0.d0
  R2x(2)  = -p_inv_2* mu_in**2 / (p+mu_in*mu_in)
  !R2x  = 0.5 / p - 0.5/p ( t * mu/sqrt(p+mu^2) )^2
  do i = 0,n_pt_in
    d(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d1(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d2(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d3(i) = 0.d0
  enddo
  integer                        :: n_pt1,n_pt2,n_pt3,dim,i
  n_pt1 = n_pt_in
  n_pt2 = n_pt_in
  n_pt3 = n_pt_in
  a_x = power_A(1)
  b_x = power_B(1)
  call I_x1_pol_mult_one_e(a_x,b_x,R1x,R1xp,R2x,d1,n_pt1,n_pt_in)
  ! print*,'passed the first I_x1'
  if(n_pt1<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif

  R1x(0)  = (P_center(2) - A_center(2))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(2) - C_center(2))* mu_in**2 / (p+mu_in*mu_in)
  ! R1x = (P_x - A_x) - (P_x - C_x) ( t * mu/sqrt(p+mu^2) )^2
  R1xp(0)  = (P_center(2) - B_center(2))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(2) - C_center(2))* mu_in**2 / (p+mu_in*mu_in)
  !R1xp = (P_x - B_x) - (P_x - C_x) ( t * mu/sqrt(p+mu^2) )^2
  a_y = power_A(2)
  b_y = power_B(2)
  call I_x1_pol_mult_one_e(a_y,b_y,R1x,R1xp,R2x,d2,n_pt2,n_pt_in)
  ! print*,'passed the second I_x1'
  if(n_pt2<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif


  R1x(0)  = (P_center(3) - A_center(3))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(3) - C_center(3))* mu_in**2 / (p+mu_in*mu_in)
  ! R1x = (P_x - A_x) - (P_x - C_x) ( t * mu/sqrt(p+mu^2) )^2
  R1xp(0)  = (P_center(3) - B_center(3))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(3) - C_center(3))* mu_in**2 / (p+mu_in*mu_in)
  !R2x  = 0.5 / p - 0.5/p ( t * mu/sqrt(p+mu^2) )^2
  a_z = power_A(3)
  b_z = power_B(3)

  ! print*,'a_z = ',a_z
  ! print*,'b_z = ',b_z
  call I_x1_pol_mult_one_e(a_z,b_z,R1x,R1xp,R2x,d3,n_pt3,n_pt_in)
  ! print*,'passed the third I_x1'
  if(n_pt3<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif
  integer                        :: n_pt_tmp
  n_pt_tmp = 0
  call multiply_poly(d1,n_pt1,d2,n_pt2,d,n_pt_tmp)
  do i = 0,n_pt_tmp
    d1(i) = 0.d0
  enddo
  n_pt_out = 0
  call multiply_poly(d ,n_pt_tmp ,d3,n_pt3,d1,n_pt_out)
  do i = 0, n_pt_out
    d(i) = d1(i)
  enddo

end

