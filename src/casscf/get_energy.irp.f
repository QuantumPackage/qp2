program print_2rdm
 implicit none
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 integer :: i,j,k,l
 integer :: ii,jj,kk,ll
 double precision :: accu(4),twodm,thr,act_twodm2,integral,get_two_e_integral
 thr = 1.d-10


 accu = 0.d0
 do ll = 1, n_act_orb
  l = list_act(ll)
  do kk = 1, n_act_orb
   k = list_act(kk)
   do jj = 1, n_act_orb
    j = list_act(jj)
    do ii = 1, n_act_orb
     i = list_act(ii)
     integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
     accu(1) += act_two_rdm_spin_trace_mo(ii,jj,kk,ll) * integral
    !if(dabs(act_two_rdm_spin_trace_mo(ii,jj,kk,ll)).gt.thr)then
    !print*,'',ii,jj,kk,ll,act_two_rdm_spin_trace_mo(ii,jj,kk,ll)*integral
    !print*,'accu',accu(1)
    !endif
    enddo
   enddo
  enddo
 enddo
 print*,'accu             = ',accu(1)
 print*,'psi_energy_two_e = ',psi_energy_two_e
!double precision :: hij
!call i_H_j_double_alpha_beta(psi_det(1,1,1),psi_det(1,1,2),N_int,hij)
!print*,'hij * 2',hij * psi_coef(1,1) * psi_coef(2,1) * 2.d0
!print*,'psi diag         = ',psi_energy_two_e - hij * psi_coef(1,1) * psi_coef(2,1) * 2.d0
end
