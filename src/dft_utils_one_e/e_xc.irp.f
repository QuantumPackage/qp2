

 BEGIN_PROVIDER[double precision, old_energy_x_lda, (N_states) ]
&BEGIN_PROVIDER[double precision, old_energy_c_lda, (N_states) ]
 implicit none
 BEGIN_DOC
! exchange/correlation energy with the short range lda functional
 END_DOC
 integer :: istate,i,j
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: e_c,vc_a,vc_b,e_x,vx_a,vx_b
 double precision, allocatable :: rhoa(:),rhob(:)
 allocate(rhoa(N_states), rhob(N_states))
 old_energy_x_lda = 0.d0
 old_energy_c_lda = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)
   rhoa(istate) = one_e_dm_alpha_at_r(i,istate)
   rhob(istate) = one_e_dm_beta_at_r(i,istate)
   call ec_lda(rhoa(istate),rhob(istate),e_c,vc_a,vc_b)
   call ex_lda(rhoa(istate),rhob(istate),e_x,vx_a,vx_b)
   old_energy_x_lda(istate) += weight * e_x
   old_energy_c_lda(istate) += weight * e_c
  enddo
 enddo

 END_PROVIDER

 BEGIN_PROVIDER[double precision, old_energy_x_pbe, (N_states) ]
&BEGIN_PROVIDER[double precision, old_energy_c_pbe, (N_states) ]
 implicit none
 BEGIN_DOC
! exchange/correlation energy with the short range pbe functional
 END_DOC
 integer :: istate,i,j,m
 double precision :: r(3)
 double precision :: mu,weight
 double precision, allocatable :: ex(:), ec(:)
 double precision, allocatable :: rho_a(:),rho_b(:),grad_rho_a(:,:),grad_rho_b(:,:),grad_rho_a_2(:),grad_rho_b_2(:),grad_rho_a_b(:)
 double precision, allocatable :: contrib_grad_xa(:,:),contrib_grad_xb(:,:),contrib_grad_ca(:,:),contrib_grad_cb(:,:)
 double precision, allocatable :: vc_rho_a(:), vc_rho_b(:), vx_rho_a(:), vx_rho_b(:)
 double precision, allocatable :: vx_grad_rho_a_2(:), vx_grad_rho_b_2(:), vx_grad_rho_a_b(:), vc_grad_rho_a_2(:), vc_grad_rho_b_2(:), vc_grad_rho_a_b(:)
 allocate(vc_rho_a(N_states), vc_rho_b(N_states), vx_rho_a(N_states), vx_rho_b(N_states))
 allocate(vx_grad_rho_a_2(N_states), vx_grad_rho_b_2(N_states), vx_grad_rho_a_b(N_states), vc_grad_rho_a_2(N_states), vc_grad_rho_b_2(N_states), vc_grad_rho_a_b(N_states))


 allocate(rho_a(N_states), rho_b(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states))
 allocate(grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states), ex(N_states), ec(N_states))
 old_energy_x_pbe = 0.d0
 old_energy_c_pbe = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)
   rho_a(istate) =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b(istate) =  one_e_dm_and_grad_beta_in_r(4,i,istate)
   grad_rho_a(1:3,istate) =  one_e_dm_and_grad_alpha_in_r(1:3,i,istate)
   grad_rho_b(1:3,istate) =  one_e_dm_and_grad_beta_in_r(1:3,i,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate) * grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate) * grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate) * grad_rho_b(m,istate)
   enddo

                             ! inputs
   call GGA_type_functionals(r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,                 &  ! outputs exchange
                             ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  ! outputs correlation
                             ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b  )
   old_energy_x_pbe += ex * weight
   old_energy_c_pbe += ec * weight
  enddo
 enddo


END_PROVIDER

