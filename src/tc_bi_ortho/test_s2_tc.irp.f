program test_tc 
 implicit none
 read_wf = .True.
 my_grid_becke = .True.
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 read_wf = .True.
 touch read_wf
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 call routine_test_s2
 call routine_test_s2_davidson
end

subroutine routine_test_s2
 implicit none
 logical           :: do_right 
 integer :: sze ,i, N_st, j
 double precision :: sij, accu_e, accu_s, accu_e_0, accu_s_0
 double precision, allocatable :: v_0_ref(:,:),u_0(:,:),s_0_ref(:,:)
 double precision, allocatable :: v_0_new(:,:),s_0_new(:,:)
 sze = N_det
 N_st = 1
 allocate(v_0_ref(N_det,1),u_0(N_det,1),s_0_ref(N_det,1),s_0_new(N_det,1),v_0_new(N_det,1))
 print*,'Checking first the Left '
 do_right = .False.
 do i = 1, sze
  u_0(i,1) = psi_l_coef_bi_ortho(i,1)
 enddo
 call H_tc_u_0_nstates_openmp(v_0_ref,u_0,N_st,sze, do_right)
 s_0_ref = 0.d0
 do i = 1, sze
  do j = 1, sze
    call get_s2(psi_det(1,1,i),psi_det(1,1,j),N_int,sij)
    s_0_ref(i,1) += u_0(j,1) * sij
  enddo
 enddo
 call H_tc_s2_u_0_nstates_openmp(v_0_new,s_0_new,u_0,N_st,sze, do_right)
 accu_e = 0.d0
 accu_s = 0.d0
 accu_e_0 = 0.d0
 accu_s_0 = 0.d0
 do i = 1, sze
  accu_e_0 += v_0_ref(i,1) * psi_r_coef_bi_ortho(i,1)
  accu_s_0 += s_0_ref(i,1) * psi_r_coef_bi_ortho(i,1)
  accu_e += dabs(v_0_ref(i,1) - v_0_new(i,1))
  accu_s += dabs(s_0_ref(i,1) - s_0_new(i,1))
 enddo
 print*,'accu_e   = ',accu_e
 print*,'accu_s   = ',accu_s
 print*,'accu_e_0 = ',accu_e_0
 print*,'accu_s_0 = ',accu_s_0

 print*,'Checking then the right '
 do_right = .True.
 do i = 1, sze
  u_0(i,1) = psi_r_coef_bi_ortho(i,1)
 enddo
 call H_tc_u_0_nstates_openmp(v_0_ref,u_0,N_st,sze, do_right)
 s_0_ref = 0.d0
 do i = 1, sze
  do j = 1, sze
    call get_s2(psi_det(1,1,i),psi_det(1,1,j),N_int,sij)
    s_0_ref(i,1) += u_0(j,1) * sij
  enddo
 enddo
 call H_tc_s2_u_0_nstates_openmp(v_0_new,s_0_new,u_0,N_st,sze, do_right)
 accu_e = 0.d0
 accu_s = 0.d0
 accu_e_0 = 0.d0
 accu_s_0 = 0.d0
 do i = 1, sze
  accu_e_0 += v_0_ref(i,1) * psi_l_coef_bi_ortho(i,1)
  accu_s_0 += s_0_ref(i,1) * psi_l_coef_bi_ortho(i,1)
  accu_e += dabs(v_0_ref(i,1) - v_0_new(i,1))
  accu_s += dabs(s_0_ref(i,1) - s_0_new(i,1))
 enddo
 print*,'accu_e   = ',accu_e
 print*,'accu_s   = ',accu_s
 print*,'accu_e_0 = ',accu_e_0
 print*,'accu_s_0 = ',accu_s_0


end

subroutine routine_test_s2_davidson
 implicit none
 double precision, allocatable :: H_jj(:),vec_tmp(:,:), energies(:) , s2(:)
 integer :: i,istate
 logical :: converged 
 external H_tc_s2_dagger_u_0_opt
 external H_tc_s2_u_0_opt
 allocate(H_jj(N_det),vec_tmp(N_det,n_states_diag),energies(n_states_diag), s2(n_states_diag))
 do i = 1, N_det
   call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,i), psi_det(1,1,i), N_int, H_jj(i))
 enddo
 ! Preparing the left-eigenvector
 print*,'Computing the left-eigenvector '
 vec_tmp = 0.d0
 do istate = 1, N_states
  vec_tmp(1:N_det,istate) = psi_l_coef_bi_ortho(1:N_det,istate)
 enddo
 do istate = N_states+1, n_states_diag
  vec_tmp(istate,istate) = 1.d0
 enddo
 do istate = 1, N_states
  leigvec_tc_bi_orth(1:N_det,istate) = vec_tmp(1:N_det,istate)
 enddo
 integer :: n_it_max
 n_it_max = 1
 call davidson_hs2_nonsym_b1space(vec_tmp, H_jj, s2, energies, N_det, n_states, n_states_diag, n_it_max, converged, H_tc_s2_dagger_u_0_opt)
 double precision, allocatable :: v_0_new(:,:),s_0_new(:,:)
 integer :: sze,N_st
 logical           :: do_right 
 sze = N_det
 N_st = 1
 do_right = .False.
 allocate(s_0_new(N_det,1),v_0_new(N_det,1))
 call H_tc_s2_u_0_nstates_openmp(v_0_new,s_0_new,vec_tmp,N_st,sze, do_right)
 double precision :: accu_e_0, accu_s_0
 accu_e_0 = 0.d0
 accu_s_0 = 0.d0
 do i = 1, sze
  accu_e_0 += v_0_new(i,1) * vec_tmp(i,1)
  accu_s_0 += s_0_new(i,1) * vec_tmp(i,1)
 enddo
 print*,'energies = ',energies
 print*,'s2       = ',s2
 print*,'accu_e_0',accu_e_0
 print*,'accu_s_0',accu_s_0

 ! Preparing the right-eigenvector
 print*,'Computing the right-eigenvector '
 vec_tmp = 0.d0
 do istate = 1, N_states
  vec_tmp(1:N_det,istate) = psi_r_coef_bi_ortho(1:N_det,istate)
 enddo
 do istate = N_states+1, n_states_diag
  vec_tmp(istate,istate) = 1.d0
 enddo
 do istate = 1, N_states
  leigvec_tc_bi_orth(1:N_det,istate) = vec_tmp(1:N_det,istate)
 enddo
 n_it_max = 1
 call davidson_hs2_nonsym_b1space(vec_tmp, H_jj, s2, energies, N_det, n_states, n_states_diag, n_it_max, converged, H_tc_s2_u_0_opt)
 sze = N_det
 N_st = 1
 do_right = .True.
 v_0_new = 0.d0
 s_0_new = 0.d0
 call H_tc_s2_u_0_nstates_openmp(v_0_new,s_0_new,vec_tmp,N_st,sze, do_right)
 accu_e_0 = 0.d0
 accu_s_0 = 0.d0
 do i = 1, sze
  accu_e_0 += v_0_new(i,1) * vec_tmp(i,1)
  accu_s_0 += s_0_new(i,1) * vec_tmp(i,1)
 enddo
 print*,'energies = ',energies
 print*,'s2       = ',s2
 print*,'accu_e_0',accu_e_0
 print*,'accu_s_0',accu_s_0

end
