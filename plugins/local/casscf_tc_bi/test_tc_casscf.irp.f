program tc_bi_ortho

  BEGIN_DOC
  !
  ! TODO : Reads psi_det in the EZFIO folder and prints out the left- and right-eigenvectors together 
  !        with the energy. Saves the left-right wave functions at the end. 
  !
  END_DOC

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  print*, ' nb of states = ', N_states
  print*, ' nb of det    = ', N_det

!  call routine_i_h_psi
!  call routine_grad
 call routine_grad_num_dm_one_body
end

subroutine routine_i_h_psi
 implicit none
 integer :: i,j
 double precision :: i_H_chi_array(0:3,N_states),i_H_phi_array(0:3,N_states)
 double precision               :: hmono, htwoe, hthree, htot
 double precision               :: accu_l_hmono, accu_l_htwoe, accu_l_hthree, accu_l_htot
 double precision               :: accu_r_hmono, accu_r_htwoe, accu_r_hthree, accu_r_htot
 double precision               :: test_l_hmono, test_l_htwoe, test_l_hthree, test_l_htot
 double precision               :: test_r_hmono, test_r_htwoe, test_r_hthree, test_r_htot

 test_l_hmono = 0.d0
 test_l_htwoe = 0.d0
 test_l_hthree= 0.d0
 test_l_htot  = 0.d0
 test_r_hmono = 0.d0
 test_r_htwoe = 0.d0
 test_r_hthree= 0.d0
 test_r_htot  = 0.d0

 do i = 1, N_det
   call i_H_tc_psi_phi(psi_det(1,1,i),psi_det,psi_l_coef_bi_ortho,psi_r_coef_bi_ortho,&
                       N_int,N_det,N_det,N_states,i_H_chi_array,i_H_phi_array)
   accu_l_hmono = 0.d0
   accu_l_htwoe = 0.d0
   accu_l_hthree= 0.d0
   accu_l_htot  = 0.d0
   accu_r_hmono = 0.d0
   accu_r_htwoe = 0.d0
   accu_r_hthree= 0.d0
   accu_r_htot  = 0.d0
   do j = 1, N_det
    call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
    accu_l_hmono  += psi_l_coef_bi_ortho(j,1) * hmono
    accu_l_htwoe  += psi_l_coef_bi_ortho(j,1) * htwoe
    accu_l_hthree += psi_l_coef_bi_ortho(j,1) * hthree
    accu_l_htot   += psi_l_coef_bi_ortho(j,1) * htot
    call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,i), psi_det(1,1,j), N_int, hmono, htwoe, hthree, htot)
    accu_r_hmono  += psi_r_coef_bi_ortho(j,1) * hmono
    accu_r_htwoe  += psi_r_coef_bi_ortho(j,1) * htwoe
    accu_r_hthree += psi_r_coef_bi_ortho(j,1) * hthree
    accu_r_htot   += psi_r_coef_bi_ortho(j,1) * htot
   enddo
   test_l_htot   += dabs(i_H_chi_array(0,1)-accu_l_htot)
   test_l_hmono  += dabs(i_H_chi_array(1,1)-accu_l_hmono)
   test_l_htwoe  += dabs(i_H_chi_array(2,1)-accu_l_htwoe)
   test_l_hthree += dabs(i_H_chi_array(3,1)-accu_l_hthree)

   test_r_htot   += dabs(i_H_phi_array(0,1)-accu_r_htot)
   test_r_hmono  += dabs(i_H_phi_array(1,1)-accu_r_hmono)
   test_r_htwoe  += dabs(i_H_phi_array(2,1)-accu_r_htwoe)
   test_r_hthree += dabs(i_H_phi_array(3,1)-accu_r_hthree)

 enddo
 
 test_l_htot   *= 1.D0/dble(N_det)
 test_l_hmono  *= 1.D0/dble(N_det)
 test_l_htwoe  *= 1.D0/dble(N_det)
 test_l_hthree *= 1.D0/dble(N_det)

 test_r_htot   *= 1.D0/dble(N_det)
 test_r_hmono  *= 1.D0/dble(N_det)
 test_r_htwoe  *= 1.D0/dble(N_det)
 test_r_hthree *= 1.D0/dble(N_det)
 
 print*,'**************************'
 print*,'test_l_htot    = ',test_l_htot
 print*,'test_l_hmono   = ',test_l_hmono
 print*,'test_l_htwoe   = ',test_l_htwoe
 print*,'test_l_hthree  = ',test_l_hthree
 print*,'**************************'
 print*,'test_r_htot    = ',test_r_htot
 print*,'test_r_hmono   = ',test_r_hmono
 print*,'test_r_htwoe   = ',test_r_htwoe
 print*,'test_r_hthree  = ',test_r_hthree

end

subroutine routine_grad_num
 implicit none
 integer :: indx,ihole,ipart
 integer :: p,q
 double precision :: accu_l, accu_r
 double precision :: contrib_l, contrib_r

 accu_l = 0.d0
 accu_r = 0.d0
  do indx=1,nMonoEx
   q = excit(1,indx)
   p = excit(2,indx)
   contrib_l  =  dabs(dabs(gradvec_detail_left_old(0,indx)) - 2.D0 * dabs( Fock_matrix_tc_mo_tot(q,p)))
   contrib_r  =  dabs(dabs(gradvec_detail_right_old(0,indx)) -2.D0 * dabs( Fock_matrix_tc_mo_tot(p,q)))
   if(contrib_l.gt.1.d-10.or.contrib_r.gt.1.d-10)then
    print*,indx,q,p
    print*,gradvec_detail_left_old(0,indx),gradvec_detail_right_old(0,indx) 
    print*,2.D0* Fock_matrix_tc_mo_tot(q,p), 2.d0* Fock_matrix_tc_mo_tot(p,q)
   endif
   accu_l += contrib_l
   accu_r += contrib_r
  enddo
 print*,'accu_l,accu_r'
 print*,accu_l,accu_r
 
! do i = 1, nMonoEx
!  print*,i,gradvec_old(i)
! enddo

end

subroutine routine_grad_num_dm_one_body
 implicit none
 integer :: indx,ii,i,a,aa,tt,t,ibody
 double precision :: accu_l, accu_r,ref_r, new_r, ref_l, new_l
 double precision :: contrib_l, contrib_r
 double precision :: res_l(0:3),res_r(0:3)

 ibody = 2 ! check only the two-body term
 provide gradvec_detail_left_old gradvec_tc_l 
 if(.True.)then
  print*,'**************************'
  print*,'**************************'
  print*,'testing inactive-->virtual'
  accu_l = 0.d0
  accu_r = 0.d0
  do ii = 1, n_core_inact_orb
   do aa = 1, n_virt_orb
    indx = mat_idx_c_v(ii,aa) 
    ref_l = gradvec_detail_left_old(ibody,indx)
    new_l = gradvec_tc_l(ibody,indx) 
    contrib_l  =  dabs(dabs(ref_l) - dabs(new_l))
    ref_r = gradvec_detail_right_old(ibody,indx)
    new_r = gradvec_tc_r(ibody,indx) 
    contrib_r  =  dabs(dabs(ref_r) - dabs(new_r))
    i = list_core_inact(ii)
    a = list_virt(aa)
!    if(i==1.and.a==9)then
!     print*,i,a,ref_r, new_r
!     stop
!    endif
    if(contrib_l.gt.1.d-10.or.contrib_r.gt.1.d-10)then
     print*,'---------'
     print*,'warning !'
     print*,indx,i,a,ii,aa
     print*,ref_l, new_l, contrib_l
     print*,ref_r, new_r, contrib_r
     print*,gradvec_detail_left_old(0,indx),gradvec_tc_l(0,indx)
     print*,gradvec_detail_right_old(0,indx),gradvec_tc_r(0,indx)
     print*,'---------'
    endif
    accu_l += contrib_l
    accu_r += contrib_r
   enddo
  enddo
  print*,'accu_l,accu_r'
  print*,accu_l,accu_r
  print*,'**************************'
  print*,'**************************'
 endif
 
 ibody = 2 ! check only the two-body term
 if(.True.)then
  print*,'**************************'
  print*,'**************************'
  print*,'testing inactive-->active'
  accu_l = 0.d0
  accu_r = 0.d0
  do ii = 1, n_core_inact_orb
   do tt = 1, n_act_orb
    indx = mat_idx_c_a(ii,tt) 
    ref_l = gradvec_detail_left_old(ibody,indx)
    new_l = gradvec_tc_l(ibody,indx) 
    contrib_l  =  dabs(dabs(ref_l) - dabs(new_l))
    ref_r = gradvec_detail_right_old(ibody,indx)
    new_r = gradvec_tc_r(ibody,indx) 
    contrib_r  =  dabs(dabs(ref_r) - dabs(new_r))
    if(contrib_l.gt.1.d-10.or.contrib_r.gt.1.d-10)then
     print*,'---------'
     print*,'warning !'
     i = list_core_inact(ii)
     t = list_act(tt)
     print*,indx,i,t
     print*,ref_l, new_l, contrib_l
     print*,ref_r, new_r, contrib_r
     print*,'---------'
    endif
    accu_l += contrib_l
    accu_r += contrib_r
   enddo
  enddo
  print*,'accu_l,accu_r'
  print*,accu_l,accu_r
 endif

 if(.True.)then
  print*,'**************************'
  print*,'**************************'
  print*,'testing active-->virtual '
  accu_l = 0.d0
  accu_r = 0.d0
  do tt = 1, n_act_orb
   do aa = 1, n_virt_orb
    indx = mat_idx_a_v(tt,aa) 
    ref_l = gradvec_detail_left_old(ibody,indx)
    new_l = gradvec_tc_l(ibody,indx) 
    contrib_l  =  dabs(dabs(ref_l) - dabs(new_l))
    ref_r = gradvec_detail_right_old(ibody,indx)
    new_r = gradvec_tc_r(ibody,indx) 
    contrib_r  =  dabs(dabs(ref_r) - dabs(new_r))
    if(contrib_l.gt.1.d-10.or.contrib_r.gt.1.d-10)then
     print*,'---------'
     print*,'warning !'
     a = list_virt(aa)
     t = list_act(tt)
     print*,indx,t,a
     print*,ref_l, new_l, contrib_l
     print*,ref_r, new_r, contrib_r
!     print*,gradvec_detail_right_old(0,indx),gradvec_tc_r(0,indx)
     print*,'---------'
    endif
    accu_l += contrib_l
    accu_r += contrib_r
   enddo
  enddo
  print*,'accu_l,accu_r'
  print*,accu_l,accu_r
 endif
 

end
