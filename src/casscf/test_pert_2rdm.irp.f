program test_pert_2rdm
 implicit none
 read_wf = .True.
 touch read_wf
!call get_pert_2rdm
 integer :: i,j,k,l,ii,jj,kk,ll
 double precision :: accu , get_two_e_integral, integral
 accu = 0.d0
 print*,'n_orb_pert_rdm = ',n_orb_pert_rdm
 do ii = 1, n_orb_pert_rdm
  i = list_orb_pert_rdm(ii) 
  do jj = 1, n_orb_pert_rdm
   j = list_orb_pert_rdm(jj) 
   do kk = 1, n_orb_pert_rdm
    k= list_orb_pert_rdm(kk) 
    do ll = 1, n_orb_pert_rdm
     l = list_orb_pert_rdm(ll) 
     integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
!    if(dabs(pert_2rdm_provider(ii,jj,kk,ll) * integral).gt.1.d-12)then
!     print*,i,j,k,l
!     print*,pert_2rdm_provider(ii,jj,kk,ll) * integral,pert_2rdm_provider(ii,jj,kk,ll), pert_2rdm_provider(ii,jj,kk,ll), integral
!    endif
     accu    += pert_2rdm_provider(ii,jj,kk,ll) * integral
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu
end
