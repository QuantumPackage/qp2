subroutine give_all_stuffs_in_r_for_lyp_88(r,rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: rho_a(N_states),rho_b(N_states),grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_2(N_states),rho(N_states)
 double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states),grad_rho_a_b(N_states)
 double precision :: grad_aos_array(3,ao_num),aos_array(ao_num)

 call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
 integer :: i,istate
 rho = rho_a + rho_b
 grad_rho_a_2 = 0.d0
 grad_rho_b_2 = 0.d0
 grad_rho_a_b = 0.d0
 do istate = 1, N_states
  do i = 1, 3 
   grad_rho_a_2(istate) += grad_rho_a(i,istate) * grad_rho_a(i,istate)
   grad_rho_b_2(istate) += grad_rho_b(i,istate) * grad_rho_b(i,istate)
   grad_rho_a_b(istate) += grad_rho_a(i,istate) * grad_rho_b(i,istate)
  enddo
 enddo
 grad_rho_2 = grad_rho_a_2 + grad_rho_b_2 + 2.d0 * grad_rho_a_b

end


double precision function ec_lyp_88(rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)

 implicit none

 BEGIN_DOC
! LYP functional of the Lee, Yan, Parr, Phys. Rev B 1988, Vol 37, page 785.  
! The expression used is the one by Miehlich, Savin, Stoll, Preuss, CPL, 1989 which gets rid of the laplacian of the density
 END_DOC

 include 'constants.include.F'

! Input variables

 double precision, intent(in) :: rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2

! Local variables

 double precision :: a,b,c,d,c_f,omega,delta
 double precision :: rho_13,rho_inv_13,rho_83,rho_113,rho_inv_113,denom
 double precision :: thr,huge_num
 double precision :: cst_2_113,cst_8_3,rho_2,rho_a_2,rho_b_2
 double precision :: tmp1,tmp2,tmp3,tmp4
 double precision :: big1,big2,big3

! Output variables


! Constants of the LYP correlation functional

 a = 0.04918d0
 b = 0.132d0
 c = 0.2533d0
 d = 0.349d0

 thr = 1d-15
 huge_num = 1.d0/thr

 rho_13  = rho**(1d0/3d0)
 rho_113 = rho**(11d0/3d0)

 if(abs(rho_13) < thr) then
  rho_inv_13 = huge_num
 else
  rho_inv_13 = 1.d0/rho_13
 endif

 if (abs(rho_113) < thr) then
  rho_inv_113 = huge_num
 else
  rho_inv_113 = 1d0/rho_113
 endif

! Useful quantities to predefine

 denom = 1d0/(1d0 + d*rho_inv_13)
 omega = rho_inv_113*exp(-c*rho_inv_13)*denom
 delta = c*rho_inv_13 + d*rho_inv_13*denom
 c_f   = 0.3d0*(3d0*pi*pi)**(2d0/3d0)

 rho_2   = rho  *rho
 rho_a_2 = rho_a*rho_a
 rho_b_2 = rho_b*rho_b

 cst_2_113 = 2d0**(11d0/3d0)
 cst_8_3   = 8d0/3d0

 ! first term in the equation (2) of Preuss CPL, 1989

 big1 = 4d0*denom*rho_a*rho_b/rho

 tmp1 = cst_2_113*c_f*(rho_a**cst_8_3 + rho_b**cst_8_3)
 tmp2 = (47d0/18d0 - 7d0/18d0*delta)*grad_rho_2
 tmp3 = - (5d0/2d0 - 1.d0/18d0*delta)*(grad_rho_a_2 + grad_rho_b_2)
 tmp4 = - (delta - 11d0)/9d0*(rho_a/rho*grad_rho_a_2 + rho_b/rho*grad_rho_b_2)
 big2 = rho_a*rho_b*(tmp1 + tmp2 + tmp3 + tmp4)

 tmp1 = -2d0/3d0*rho_2*grad_rho_2
 tmp2 = grad_rho_b_2*(2d0/3d0*rho_2 - rho_a_2)
 tmp3 = grad_rho_a_2*(2d0/3d0*rho_2 - rho_b_2)
 big3 = tmp1 + tmp2 + tmp3


 ec_lyp_88 = -a*big1 -a*b*omega*big2 -a*b*omega*big3

end

