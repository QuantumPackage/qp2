
subroutine GGA_sr_type_functionals(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
                                ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &
                                ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b )
 implicit none
 BEGIN_DOC
 ! routine that helps in building the x/c potentials on the AO basis for a GGA functional with a short-range interaction
 END_DOC
 double precision, intent(in)  :: mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision, intent(out) :: ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b
 double precision, intent(out) :: ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b
 integer          :: istate
 double precision :: r2(3),dr2(3), local_potential,r12,dx2

  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo

   ! exhange energy and potentials 
   call ex_pbe_sr(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b)

   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a,rho_b,rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,sigmaoo,sigmacc,sigmaco)

   ! correlation energy and potentials 
   call ec_pbe_sr(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)

   ! convertion from (closed, open) formalism to (alpha,beta) formalism 
   call v_rho_oc_to_v_rho_ab(vrhoo,vrhoc,vc_rho_a,vc_rho_b)
   call v_grad_rho_oc_to_v_grad_rho_ab(vsigmaoo,vsigmacc,vsigmaco,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b)
end


