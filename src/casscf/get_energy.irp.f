program print_2rdm
 implicit none
 BEGIN_DOC
 ! get the active part of the bielectronic energy on a given wave function.
 !
 ! useful to test the active part of the spin trace 2 rdms
 END_DOC
 no_vvvv_integrals = .True.
 read_wf = .True.
 touch read_wf no_vvvv_integrals
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
     accu(1) += state_av_act_two_rdm_spin_trace_mo(ii,jj,kk,ll) * integral
    enddo
   enddo
  enddo
 enddo
 print*,'accu             = ',accu(1)

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
     accu(1) += state_av_act_two_rdm_openmp_spin_trace_mo(ii,jj,kk,ll) * integral
    enddo
   enddo
  enddo
 enddo
 print*,'accu             = ',accu(1)
 print*,'psi_energy_two_e = ',psi_energy_two_e

 print *,  psi_energy_with_nucl_rep
end
