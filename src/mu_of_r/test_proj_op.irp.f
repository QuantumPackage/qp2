program projected_operators
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
  ! You specify that you want to avoid any contribution from 
  ! orbitals coming from core 
  no_core_density = .True.
  touch no_core_density
  mu_of_r_potential = "cas_full"
  touch mu_of_r_potential 
  print*,'Using Valence Only functions'
!  call test_f_HF_valence_ab
!  call routine_full_mos
!   call test_f_ii_valence_ab
!   call test_f_ia_valence_ab
!  call test_f_ii_ia_aa_valence_ab
! call test
!  call test_f_mean_field
 call test_grad_f_mean_field
end


subroutine test
 implicit none
 integer :: i_point
 double precision :: ref, new, accu, weight
 accu = 0.d0
 do i_point = 1, n_points_final_grid
  ref = f_hf_cholesky_sparse(i_point)
  new = f_hf_cholesky_sparse_bis(i_point)
  weight = final_weight_at_r_vector(i_point)
  accu += dabs(ref - new) * weight
 enddo
 print*,'accu = ',accu

end

subroutine test_f_mean_field
 implicit none
 integer :: i_point
 double precision :: weight,r(3)
 double precision :: ref_f, new_f, accu_f
 double precision :: ref_two_dens, new_two_dens, accu_two_dens, dm_a, dm_b
 accu_f = 0.d0
 accu_two_dens = 0.d0
 do i_point = 1, n_points_final_grid
  r(1:3)   = final_grid_points(1:3,i_point)
  weight = final_weight_at_r_vector(i_point)
  call get_f_mf_ab(r,new_f,new_two_dens, dm_a, dm_b)
  call f_HF_valence_ab(r,r,ref_f,ref_two_dens)
  accu_f += weight * dabs(new_f- ref_f)
  accu_two_dens += weight * dabs(new_two_dens - ref_two_dens)
 enddo
 print*,'accu_f        = ',accu_f
 print*,'accu_two_dens = ',accu_two_dens

end

subroutine test_grad_f_mean_field
 implicit none
 integer :: i_point,k
 double precision :: weight,r(3)
 double precision :: grad_f_mf_ab(3), grad_two_bod_dens(3)
 double precision :: grad_dm_a(3), grad_dm_b(3)
 double precision :: f_mf_ab,two_bod_dens, dm_a, dm_b

 double precision :: num_grad_f_mf_ab(3), num_grad_two_bod_dens(3)
 double precision :: num_grad_dm_a(3), num_grad_dm_b(3)
 double precision :: f_mf_ab_p,f_mf_ab_m
 double precision :: two_bod_dens_p, two_bod_dens_m
 double precision :: dm_a_p, dm_a_m
 double precision :: dm_b_p, dm_b_m
 double precision :: rbis(3), dr
 double precision :: accu_grad_f_mf_ab(3),accu_grad_two_bod_dens(3)
 double precision :: accu_grad_dm_a(3),accu_grad_dm_b(3)
 double precision :: accu_f_mf_ab, accu_two_bod_dens, accu_dm_a, accu_dm_b
 dr = 0.00001d0
 accu_f_mf_ab = 0.d0 
 accu_two_bod_dens = 0.d0 
 accu_dm_a = 0.d0 
 accu_dm_b = 0.d0

 accu_grad_f_mf_ab = 0.d0
 accu_grad_two_bod_dens = 0.d0
 accu_grad_dm_a = 0.d0
 accu_grad_dm_b = 0.d0
 do i_point = 1, n_points_final_grid
  r(1:3)   = final_grid_points(1:3,i_point)
  weight = final_weight_at_r_vector(i_point)
  call get_grad_f_mf_ab(r,grad_f_mf_ab, grad_two_bod_dens,f_mf_ab,two_bod_dens, dm_a, dm_b,grad_dm_a, grad_dm_b)
  call get_f_mf_ab(r,f_mf_ab_p,two_bod_dens_p, dm_a_p, dm_b_p)
  accu_f_mf_ab += weight * dabs(f_mf_ab - f_mf_ab_p)
  accu_two_bod_dens += weight * dabs(two_bod_dens - two_bod_dens_p)
  accu_dm_a += weight*dabs(dm_a - dm_a_p)
  accu_dm_b += weight*dabs(dm_b - dm_b_p)
  do k = 1, 3
   rbis = r
   rbis(k) += dr
   call get_f_mf_ab(rbis,f_mf_ab_p,two_bod_dens_p, dm_a_p, dm_b_p)
   rbis = r
   rbis(k) -= dr
   call get_f_mf_ab(rbis,f_mf_ab_m,two_bod_dens_m, dm_a_m, dm_b_m)
   num_grad_f_mf_ab(k) = (f_mf_ab_p - f_mf_ab_m)/(2.d0*dr)
   num_grad_two_bod_dens(k) = (two_bod_dens_p - two_bod_dens_m)/(2.d0*dr)
   num_grad_dm_a(k) = (dm_a_p - dm_a_m)/(2.d0*dr)
   num_grad_dm_b(k) = (dm_b_p - dm_b_m)/(2.d0*dr)
  enddo
  do k = 1, 3
   accu_grad_f_mf_ab(k) += weight * dabs(grad_f_mf_ab(k) - num_grad_f_mf_ab(k))
   accu_grad_two_bod_dens(k) += weight * dabs(grad_two_bod_dens(k) - num_grad_two_bod_dens(k))
   accu_grad_dm_a(k) += weight * dabs(grad_dm_a(k) - num_grad_dm_a(k))
   accu_grad_dm_b(k) += weight * dabs(grad_dm_b(k) - num_grad_dm_b(k))
  enddo
 enddo
 print*,'accu_f_mf_ab = ',accu_f_mf_ab
 print*,'accu_two_bod_dens = ',accu_two_bod_dens
 print*,'accu_dm_a = ',accu_dm_a
 print*,'accu_dm_b = ',accu_dm_b
 print*,'accu_grad_f_mf_ab = '
 print*,accu_grad_f_mf_ab
 print*,'accu_grad_two_bod_dens = '
 print*,accu_grad_two_bod_dens
 print*,'accu_dm_a = '
 print*,accu_grad_dm_a
 print*,'accu_dm_b = '
 print*,accu_grad_dm_b

end
