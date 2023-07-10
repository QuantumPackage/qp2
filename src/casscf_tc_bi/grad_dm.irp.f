 BEGIN_PROVIDER [real*8, gradvec_tc_r, (0:3,nMonoEx)]
&BEGIN_PROVIDER [real*8, gradvec_tc_l,  (0:3,nMonoEx)]
 implicit none
 integer                        :: ii,tt,aa,indx
 integer                        :: i,t,a,fff
 double precision :: res_l(0:3), res_r(0:3)
 gradvec_tc_l = 0.d0
 gradvec_tc_r = 0.d0
  do i=1,n_core_inact_orb
    ii=list_core_inact(i)
    do t=1,n_act_orb
      tt=list_act(t)
      indx = mat_idx_c_a(i,t) 
      call gradvec_tc_it(ii,tt,res_l,res_r)
      do fff = 0,3
       gradvec_tc_l(fff,indx)=res_l(fff)
       gradvec_tc_r(fff,indx)=res_r(fff)
      enddo
    end do
  end do
  
  do i=1,n_core_inact_orb
    ii=list_core_inact(i)
    do a=1,n_virt_orb
      indx = mat_idx_c_v(i,a) 
      aa=list_virt(a)
      call gradvec_tc_ia(ii,aa,res_l)
      call gradvec_tc_ia(aa,ii,res_r)
      do fff = 0,3
       gradvec_tc_l(fff,indx)=res_l(fff)
       gradvec_tc_r(fff,indx)=res_r(fff)
      enddo
    end do
  end do
  
  do t=1,n_act_orb
    do a=1,n_virt_orb
      indx = mat_idx_a_v(i,a) 
!      gradvec_tc_l(indx)=gradvec_ta(t,a)
    end do
  end do
END_PROVIDER 

subroutine gradvec_tc_ia(i,a,res)
 implicit none
 BEGIN_DOC
! doubly occupied --> virtual TC gradient 
!
! Corresponds to <X0|H E_i^a|Phi_0>
 END_DOC
 integer, intent(in) :: i,a
 double precision, intent(out) :: res(0:3)
 res = 0.d0
 res(1) = -2 * mo_bi_ortho_tc_one_e(i,a)
 
end

subroutine gradvec_tc_it(i,t,res_l, res_r)
 implicit none
 BEGIN_DOC
! doubly occupied --> active TC gradient 
!
! Corresponds to res_r = <X0|H E_i^t|Phi_0>
!
!                res_l = <X0|E_i^t H |Phi_0>
 END_DOC
 integer, intent(in) :: i,t
 double precision, intent(out) :: res_l(0:3),res_r(0:3)
 integer :: rr,r,ss,s,m
 double precision :: dm
 res_r = 0.d0
 do m = 1, mo_num
  res_r(1) += mo_bi_ortho_tc_one_e(i,m) * tc_transition_matrix_mo(t,m,1,1) & 
             -mo_bi_ortho_tc_one_e(m,t) * tc_transition_matrix_mo(m,i,1,1)
 enddo
 res_l = 0.d0
 do m = 1, mo_num
  res_l(1) += mo_bi_ortho_tc_one_e(t,m) * tc_transition_matrix_mo(i,m,1,1) & 
             -mo_bi_ortho_tc_one_e(m,i) * tc_transition_matrix_mo(m,t,1,1)
 enddo
 
end
