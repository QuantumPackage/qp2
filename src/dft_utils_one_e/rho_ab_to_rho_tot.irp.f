subroutine rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)
 implicit none
 double precision, intent(in)  :: rho_a,rho_b
 double precision, intent(out) :: rho_o,rho_c
 rho_c=rho_a+rho_b
 rho_o=rho_a-rho_b
end

subroutine rho_oc_to_rho_ab(rho_o,rho_c,rho_a,rho_b)
 implicit none
 double precision, intent(in)  :: rho_o,rho_c
 double precision, intent(out) :: rho_a,rho_b
 rho_a= 0.5d0*(rho_c+rho_o)
 rho_b= 0.5d0*(rho_c-rho_o)
end



subroutine grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_o_2,grad_rho_c_2,grad_rho_o_c)
 implicit none
 double precision, intent(in)  :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision, intent(out) :: grad_rho_o_2,grad_rho_c_2,grad_rho_o_c
 grad_rho_c_2 = grad_rho_a_2 + grad_rho_b_2 + 2d0*grad_rho_a_b
 grad_rho_o_2 = grad_rho_a_2 + grad_rho_b_2 - 2d0*grad_rho_a_b
 grad_rho_o_c = grad_rho_a_2 - grad_rho_b_2
end



subroutine v_rho_ab_to_v_rho_oc(v_rho_a,v_rho_b,v_rho_o,v_rho_c)
 implicit none
 double precision, intent(in)  :: v_rho_a,v_rho_b
 double precision, intent(out) :: v_rho_o,v_rho_c
 v_rho_c = 0.5d0*(v_rho_a + v_rho_b)
 v_rho_o = 0.5d0*(v_rho_a - v_rho_b)
end

subroutine v_rho_oc_to_v_rho_ab(v_rho_o,v_rho_c,v_rho_a,v_rho_b)
 implicit none
 double precision, intent(in)  :: v_rho_o,v_rho_c
 double precision, intent(out) :: v_rho_a,v_rho_b
 v_rho_a = v_rho_c + v_rho_o
 v_rho_b = v_rho_c - v_rho_o
end



subroutine v_grad_rho_oc_to_v_grad_rho_ab(v_grad_rho_o_2,v_grad_rho_c_2,v_grad_rho_o_c,v_grad_rho_a_2,v_grad_rho_b_2,v_grad_rho_a_b)
 implicit none
 double precision, intent(in)  :: v_grad_rho_o_2,v_grad_rho_c_2,v_grad_rho_o_c
 double precision, intent(out) :: v_grad_rho_a_2,v_grad_rho_b_2,v_grad_rho_a_b
 v_grad_rho_a_2 = v_grad_rho_o_2 + v_grad_rho_c_2 + v_grad_rho_o_c
 v_grad_rho_b_2 = v_grad_rho_o_2 + v_grad_rho_c_2 - v_grad_rho_o_c
 v_grad_rho_a_b = -2d0 * v_grad_rho_o_2 + 2d0 * v_grad_rho_c_2
end



















