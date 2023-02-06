
BEGIN_PROVIDER [ double precision, fock_3_w_kk_sum, (n_points_final_grid,3)]
 implicit none
 integer :: mm, ipoint,k
 double precision :: w_kk
 fock_3_w_kk_sum = 0.d0
 do k = 1, elec_beta_num
  do mm = 1, 3
   do ipoint = 1, n_points_final_grid
    w_kk   = x_W_ij_erf_rk(ipoint,mm,k,k) 
    fock_3_w_kk_sum(ipoint,mm) += w_kk
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_ki_mos_k, (n_points_final_grid,3,mo_num)]
 implicit none
 integer :: mm, ipoint,k,i
 double precision :: w_ki, mo_k
 fock_3_w_ki_mos_k = 0.d0
 do i = 1, mo_num
  do k = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i) 
     mo_k = mos_in_r_array(k,ipoint)
     fock_3_w_ki_mos_k(ipoint,mm,i) += w_ki * mo_k
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_kl_w_kl, (n_points_final_grid,3)]
 implicit none
 integer :: k,j,ipoint,mm
 double precision :: w_kj
 fock_3_w_kl_w_kl = 0.d0
 do j = 1, elec_beta_num
  do k = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     w_kj   = x_W_ij_erf_rk(ipoint,mm,k,j) 
     fock_3_w_kl_w_kl(ipoint,mm) += w_kj * w_kj
    enddo
   enddo
  enddo
 enddo


END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_rho_beta, (n_points_final_grid)]
 implicit none
 integer :: ipoint,k
 fock_3_rho_beta = 0.d0
 do ipoint = 1, n_points_final_grid
  do k = 1, elec_beta_num
   fock_3_rho_beta(ipoint) += mos_in_r_array(k,ipoint) * mos_in_r_array(k,ipoint)
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_kl_mo_k_mo_l, (n_points_final_grid,3)]
 implicit none
 integer :: ipoint,k,l,mm
 double precision :: mos_k, mos_l, w_kl
 fock_3_w_kl_mo_k_mo_l = 0.d0
 do k = 1, elec_beta_num
  do l = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     mos_k  = mos_in_r_array_transp(ipoint,k) 
     mos_l  = mos_in_r_array_transp(ipoint,l) 
     w_kl   = x_W_ij_erf_rk(ipoint,mm,l,k)
     fock_3_w_kl_mo_k_mo_l(ipoint,mm) += w_kl * mos_k * mos_l 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_ki_wk_a, (n_points_final_grid,3,mo_num, mo_num)]
 implicit none
 integer :: ipoint,i,a,k,mm
 double precision :: w_ki,w_ka
 fock_3_w_ki_wk_a = 0.d0
 do i = 1, mo_num
  do a = 1, mo_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     do k = 1, elec_beta_num
      w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i)
      w_ka   = x_W_ij_erf_rk(ipoint,mm,k,a)
      fock_3_w_ki_wk_a(ipoint,mm,a,i) += w_ki * w_ka
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_trace_w_tilde, (n_points_final_grid,3)]
 implicit none
 integer :: ipoint,k,mm
 fock_3_trace_w_tilde = 0.d0
 do k = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     fock_3_trace_w_tilde(ipoint,mm) += fock_3_w_ki_wk_a(ipoint,mm,k,k)
    enddo
   enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_kl_wla_phi_k, (n_points_final_grid,3,mo_num)]
 implicit none
 integer :: ipoint,a,k,mm,l
 double precision :: w_kl,w_la, mo_k
 fock_3_w_kl_wla_phi_k = 0.d0
 do a = 1, mo_num
  do k = 1, elec_beta_num 
   do l = 1, elec_beta_num
    do mm = 1, 3
     do ipoint = 1, n_points_final_grid
      w_kl   = x_W_ij_erf_rk(ipoint,mm,l,k)
      w_la   = x_W_ij_erf_rk(ipoint,mm,l,a)
      mo_k  = mos_in_r_array_transp(ipoint,k) 
      fock_3_w_kl_wla_phi_k(ipoint,mm,a) += w_kl * w_la * mo_k
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

