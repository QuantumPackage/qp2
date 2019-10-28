program test_pert_2rdm
 implicit none
 read_wf = .True.
 touch read_wf
 pert_2rdm = .True.
 touch pert_2rdm

!provide is_pert_2rdm_provided
 call test
 end

 subroutine test
 implicit none
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
     if(dabs(pert_2rdm_provider(ii,jj,kk,ll) * integral).gt.1.d-12)then
      print*,i,j,k,l
      print*,pert_2rdm_provider(ii,jj,kk,ll) , integral, pert_2rdm_provider(ii,jj,kk,ll)* integral
     endif
     accu    += pert_2rdm_provider(ii,jj,kk,ll) * integral
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu
 do ii = 1, n_orb_pert_rdm
  i = list_orb_pert_rdm(ii) 
  do kk = 1, n_orb_pert_rdm
   k= list_orb_pert_rdm(kk) 
   do jj = 1, n_core_inact_orb
    j = list_core_inact(jj)
    l = j
    integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
    double precision :: exchange 
    exchange = get_two_e_integral(i,j,l,k,mo_integrals_map)
    accu    +=  pert_1rdm_provider(kk,ii) * ( 2.d0 * integral - exchange)
   enddo
  enddo
 enddo
 print*,'accu = ',accu
 double precision :: accu_1
 accu_1 = 0.d0
 do ii = 1, n_orb_pert_rdm
  i = list_orb_pert_rdm(ii) 
  do jj = 1, n_orb_pert_rdm
   j = list_orb_pert_rdm(jj) 
   if(dabs(pert_1rdm_provider(jj,ii) - pert_1rdm_provider_bis(jj,ii)).gt.1.d-10)then
    print*,ii,jj,pert_1rdm_provider(jj,ii),pert_1rdm_provider_bis(jj,ii)
   endif
   if(dabs(pert_1rdm_provider(jj,ii) * mo_one_e_integrals(j,i)).gt.1.d-10)then
    print*,j,i,pert_1rdm_provider(jj,ii) , mo_one_e_integrals(j,i),pert_1rdm_provider(jj,ii) * mo_one_e_integrals(j,i)
   endif
   accu_1 +=  pert_1rdm_provider(jj,ii) * mo_one_e_integrals(j,i)
  enddo
 enddo
 print*,'accu_1 = ',accu_1
 print*,'*********'
 print*,'*********'
 print*,'accu = ',accu + accu_1

 print*,'pt2  = ',pt2_pert_2rdm
end
