  use bitmasks ! you need to include the bitmasks_module.f90 features
 BEGIN_PROVIDER [ double precision, e_tilde_00]
 implicit none
 double precision :: hmono,htwoe,hthree,htot
 call htilde_mu_mat_bi_ortho(HF_bitmask,HF_bitmask,N_int,hmono,htwoe,hthree,htot)
 e_tilde_00 = htot
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, e_pt2_tc_bi_orth]
&BEGIN_PROVIDER [ double precision, e_pt2_tc_bi_orth_single]
&BEGIN_PROVIDER [ double precision, e_pt2_tc_bi_orth_double]
 implicit none 
 integer :: i,degree
 double precision :: hmono,htwoe,hthree,htilde_ij,coef_pt1,e_i0,delta_e
 e_pt2_tc_bi_orth = 0.d0
 e_pt2_tc_bi_orth_single = 0.d0
 e_pt2_tc_bi_orth_double = 0.d0
 do i = 1, N_det
  call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
  if(degree == 1 .or. degree == 2)then
   call htilde_mu_mat_bi_ortho(psi_det(1,1,i),HF_bitmask,N_int,hmono,htwoe,hthree,htilde_ij)
   call htilde_mu_mat_bi_ortho(psi_det(1,1,i),psi_det(1,1,i),N_int,hmono,htwoe,hthree,e_i0)
   delta_e = e_tilde_00 - e_i0
   coef_pt1 = htilde_ij / delta_e
   call htilde_mu_mat_bi_ortho(HF_bitmask,psi_det(1,1,i),N_int,hmono,htwoe,hthree,htilde_ij)
   e_pt2_tc_bi_orth += coef_pt1 * htilde_ij
   if(degree == 1)then
    e_pt2_tc_bi_orth_single += coef_pt1 * htilde_ij
   else 
!    print*,'coef_pt1, e_pt2',coef_pt1,coef_pt1 * htilde_ij
    e_pt2_tc_bi_orth_double += coef_pt1 * htilde_ij
   endif
  endif
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, e_tilde_bi_orth_00]
 implicit none
 double precision :: hmono,htwoe,hthree,htilde_ij
 call htilde_mu_mat_bi_ortho(HF_bitmask,HF_bitmask,N_int,hmono,htwoe,hthree,e_tilde_bi_orth_00)
 e_tilde_bi_orth_00 += nuclear_repulsion
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, e_corr_bi_orth ]
&BEGIN_PROVIDER [ double precision, e_corr_bi_orth_proj ]
&BEGIN_PROVIDER [ double precision, e_corr_single_bi_orth ]
&BEGIN_PROVIDER [ double precision, e_corr_double_bi_orth ]
 implicit none 
 integer :: i,degree
 double precision :: hmono,htwoe,hthree,htilde_ij
 
 e_corr_bi_orth = 0.d0
 e_corr_single_bi_orth = 0.d0
 e_corr_double_bi_orth = 0.d0
 do i = 1, N_det
  call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
  call htilde_mu_mat_bi_ortho(HF_bitmask,psi_det(1,1,i),N_int,hmono,htwoe,hthree,htilde_ij)
  if(degree == 1)then
   e_corr_single_bi_orth += reigvec_tc_bi_orth(i,1) * htilde_ij/reigvec_tc_bi_orth(1,1)
  else if(degree == 2)then
   e_corr_double_bi_orth += reigvec_tc_bi_orth(i,1) * htilde_ij/reigvec_tc_bi_orth(1,1)
!   print*,'coef_wf , e_cor',reigvec_tc_bi_orth(i,1)/reigvec_tc_bi_orth(1,1), reigvec_tc_bi_orth(i,1) * htilde_ij/reigvec_tc_bi_orth(1,1)
  endif
 enddo
 e_corr_bi_orth_proj = e_corr_single_bi_orth + e_corr_double_bi_orth
 e_corr_bi_orth = eigval_right_tc_bi_orth(1) - e_tilde_bi_orth_00
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, e_tc_left_right ]
 implicit none
 integer :: i,j
 double precision :: hmono,htwoe,hthree,htilde_ij,accu
 e_tc_left_right = 0.d0
 accu = 0.d0
 do i = 1, N_det
  accu += reigvec_tc_bi_orth(i,1) * leigvec_tc_bi_orth(i,1)
  do j = 1, N_det
   call htilde_mu_mat_bi_ortho(psi_det(1,1,j),psi_det(1,1,i),N_int,hmono,htwoe,hthree,htilde_ij)
   e_tc_left_right += htilde_ij * reigvec_tc_bi_orth(i,1) * leigvec_tc_bi_orth(j,1)
  enddo
 enddo
 e_tc_left_right *= 1.d0/accu 
 e_tc_left_right += nuclear_repulsion

 END_PROVIDER 


BEGIN_PROVIDER [ double precision, coef_pt1_bi_ortho, (N_det)]
 implicit none
 integer :: i,degree
 double precision :: hmono,htwoe,hthree,htilde_ij,coef_pt1,e_i0,delta_e
 do i = 1, N_det
  call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
  if(degree==0)then
   coef_pt1_bi_ortho(i) = 1.d0
  else
   call htilde_mu_mat_bi_ortho(psi_det(1,1,i),HF_bitmask,N_int,hmono,htwoe,hthree,htilde_ij)
   call htilde_mu_mat_bi_ortho(psi_det(1,1,i),psi_det(1,1,i),N_int,hmono,htwoe,hthree,e_i0)
   delta_e = e_tilde_00 - e_i0
   coef_pt1 = htilde_ij / delta_e
   coef_pt1_bi_ortho(i)= coef_pt1
  endif
 enddo
END_PROVIDER
