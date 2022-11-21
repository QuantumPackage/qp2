 BEGIN_PROVIDER [ double precision, reigvec_tc_cisd_sc2_bi_ortho, (N_det,N_states)]
&BEGIN_PROVIDER [ double precision, leigvec_tc_cisd_sc2_bi_ortho, (N_det,N_states)]
&BEGIN_PROVIDER [ double precision, eigval_tc_cisd_sc2_bi_ortho, (N_states)]
 implicit none
 integer :: it,n_real,degree,i,istate
 double precision :: e_before, e_current,thr, hmono,htwoe,hthree,accu
 double precision, allocatable :: e_corr_dets(:),h0j(:), h_sc2(:,:), dressing_dets(:)
 double precision, allocatable :: leigvec_tc_bi_orth_tmp(:,:),reigvec_tc_bi_orth_tmp(:,:),eigval_right_tmp(:)
 allocate(leigvec_tc_bi_orth_tmp(N_det,N_det),reigvec_tc_bi_orth_tmp(N_det,N_det),eigval_right_tmp(N_det))
 allocate(e_corr_dets(N_det),h0j(N_det),h_sc2(N_det,N_det),dressing_dets(N_det))
 allocate(H_jj(N_det),vec_tmp(N_det,n_states_diag),eigval_tmp(N_states))
 dressing_dets = 0.d0
 do i = 1, N_det
  call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,i), psi_det(1,1,i), N_int, H_jj(i))
  call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
  if(degree == 1 .or. degree == 2)then
   call htilde_mu_mat_bi_ortho(HF_bitmask,psi_det(1,1,i),N_int,hmono,htwoe,hthree,h0j(i))
  endif
 enddo
 reigvec_tc_bi_orth_tmp = 0.d0
 do i = 1, N_det 
  reigvec_tc_bi_orth_tmp(i,1) = psi_r_coef_bi_ortho(i,1) 
 enddo
 vec_tmp = 0.d0
 do istate = 1, N_states
  vec_tmp(:,istate) = reigvec_tc_bi_orth_tmp(:,istate)
 enddo
 do istate = N_states+1, n_states_diag
  vec_tmp(istate,istate) = 1.d0
 enddo
 print*,'Diagonalizing the TC CISD '
 call davidson_general_diag_dressed_ext_rout_nonsym_b1space(vec_tmp, H_jj, dressing_dets,eigval_tmp, N_det, n_states, n_states_diag, converged, htc_bi_ortho_calc_tdav)
 do i = 1, N_det 
  e_corr_dets(i) = reigvec_tc_bi_orth_tmp(i,1) * h0j(i)/reigvec_tc_bi_orth_tmp(1,1)
 enddo
 E_before = eigval_tmp(1)
 print*,'Starting from ',E_before

 e_current = 10.d0
 thr = 1.d-5
 it = 0
 dressing_dets = 0.d0
  double precision, allocatable :: H_jj(:),vec_tmp(:,:),eigval_tmp(:)
  external                         htc_bi_ortho_calc_tdav
  external                         htcdag_bi_ortho_calc_tdav
  logical                       :: converged
 do while (dabs(E_before-E_current).gt.thr)
  it += 1
  E_before = E_current
!  h_sc2 = htilde_matrix_elmt_bi_ortho
  call get_cisd_sc2_dressing(psi_det,e_corr_dets,N_det,dressing_dets)
  do i = 1, N_det
!   print*,'dressing_dets(i) = ',dressing_dets(i)
   h_sc2(i,i) += dressing_dets(i)
  enddo
  print*,'********************'
  print*,'iteration       ',it
!  call non_hrmt_real_diag(N_det,h_sc2,& 
!       leigvec_tc_bi_orth_tmp,reigvec_tc_bi_orth_tmp,& 
!       n_real,eigval_right_tmp)
!  print*,'eigval_right_tmp(1)',eigval_right_tmp(1)
  vec_tmp = 0.d0
  do istate = 1, N_states
   vec_tmp(:,istate) = reigvec_tc_bi_orth_tmp(:,istate)
  enddo
  do istate = N_states+1, n_states_diag
   vec_tmp(istate,istate) = 1.d0
  enddo
  call davidson_general_diag_dressed_ext_rout_nonsym_b1space(vec_tmp, H_jj, dressing_dets,eigval_tmp, N_det, n_states, n_states_diag, converged, htc_bi_ortho_calc_tdav)
  print*,'outside Davidson'
  print*,'eigval_tmp(1) = ',eigval_tmp(1)
  do i = 1, N_det 
   reigvec_tc_bi_orth_tmp(i,1) = vec_tmp(i,1)
   e_corr_dets(i) = reigvec_tc_bi_orth_tmp(i,1) * h0j(i)/reigvec_tc_bi_orth_tmp(1,1)
  enddo
!  E_current = eigval_right_tmp(1)
  E_current = eigval_tmp(1)
  print*,'it, E(SC)^2 = ',it,E_current
 enddo
 eigval_tc_cisd_sc2_bi_ortho(1:N_states) = eigval_right_tmp(1:N_states)
 reigvec_tc_cisd_sc2_bi_ortho(1:N_det,1:N_states) = reigvec_tc_bi_orth_tmp(1:N_det,1:N_states)
 leigvec_tc_cisd_sc2_bi_ortho(1:N_det,1:N_states) = leigvec_tc_bi_orth_tmp(1:N_det,1:N_states)
 
END_PROVIDER 

subroutine get_cisd_sc2_dressing(dets,e_corr_dets,ndet,dressing_dets)
 implicit none
  use bitmasks
 integer, intent(in) :: ndet
 integer(bit_kind), intent(in)  :: dets(N_int,2,ndet)
 double precision, intent(in)   :: e_corr_dets(ndet)
 double precision, intent(out) :: dressing_dets(ndet)
 integer, allocatable  :: degree(:),hole(:,:),part(:,:),spin(:,:)
 integer(bit_kind), allocatable :: hole_part(:,:,:)
 integer :: i,j,k, exc(0:2,2,2),h1,p1,h2,p2,s1,s2
 integer(bit_kind) :: xorvec(2,N_int)

 double precision :: phase
 dressing_dets = 0.d0
 allocate(degree(ndet),hole(2,ndet),part(2,ndet), spin(2,ndet),hole_part(N_int,2,ndet))
 do i = 2, ndet
  call get_excitation_degree(HF_bitmask,dets(1,1,i),degree(i),N_int)
  do j = 1, N_int
   hole_part(j,1,i) = xor( HF_bitmask(j,1), dets(j,1,i))
   hole_part(j,2,i) = xor( HF_bitmask(j,2), dets(j,2,i))
  enddo
  if(degree(i) == 1)then
   call get_single_excitation(HF_bitmask,psi_det(1,1,i),exc,phase,N_int)
  else if(degree(i) == 2)then
   call get_double_excitation(HF_bitmask,psi_det(1,1,i),exc,phase,N_int)
  endif
  call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
  hole(1,i) = h1
  hole(2,i) = h2
  part(1,i) = p1
  part(2,i) = p2
  spin(1,i) = s1
  spin(2,i) = s2
 enddo
 
 integer :: same
 if(elec_alpha_num+elec_beta_num<3)return
 do i = 2, ndet
  do j = i+1, ndet
   same = 0
   if(degree(i) == degree(j) .and. degree(i)==1)cycle
   do k = 1, N_int
    xorvec(k,1) = iand(hole_part(k,1,i),hole_part(k,1,j))
    xorvec(k,2) = iand(hole_part(k,2,i),hole_part(k,2,j))
    same += popcnt(xorvec(k,1)) + popcnt(xorvec(k,2)) 
   enddo
!   print*,'i,j',i,j
!   call debug_det(dets(1,1,i),N_int) 
!   call debug_det(hole_part(1,1,i),N_int) 
!   call debug_det(dets(1,1,j),N_int) 
!   call debug_det(hole_part(1,1,j),N_int) 
!   print*,'same = ',same
   if(same.eq.0)then
    dressing_dets(i) += e_corr_dets(j)  
    dressing_dets(j) += e_corr_dets(i)  
   endif
  enddo
 enddo
 
end
