 BEGIN_PROVIDER [real*8, gradvec_tc_r, (0:3,nMonoEx)]
&BEGIN_PROVIDER [real*8, gradvec_tc_l,  (0:3,nMonoEx)]
 BEGIN_DOC
! gradvec_tc_r(0:3,i) = <Chi_0| H E_q^p(i) |Phi_0> 
! 
! gradvec_tc_l(0:3,i) = <Chi_0| E_p^q(i) H |Phi_0> 
!
! where the indices "i" corresponds to E_q^p(i) 
! 
! i = mat_idx_c_a(q,p) 
!
! and gradvec_tc_r/l(0) = full matrix element
!
!     gradvec_tc_r/l(1) = one-body part 

!     gradvec_tc_r/l(2) = two-body part 

!     gradvec_tc_r/l(3) = three-body part 
 END_DOC
 implicit none
 integer                        :: ii,tt,aa,indx
 integer                        :: i,t,a,fff
 double precision :: res_l(0:3), res_r(0:3)
 gradvec_tc_l = 0.d0
 gradvec_tc_r = 0.d0
 ! computing the core/inactive --> virtual orbitals gradients
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
      call gradvec_tc_ia(ii,aa,res_l,res_r)
      do fff = 0,3
       gradvec_tc_l(fff,indx)=res_l(fff)
       gradvec_tc_r(fff,indx)=res_r(fff)
      enddo
    end do
  end do
  
  do t=1,n_act_orb
    tt=list_act(t)
    do a=1,n_virt_orb
      aa=list_virt(a)
      indx = mat_idx_a_v(t,a) 
      call gradvec_tc_ta(tt,aa,res_l, res_r)
      do fff = 0,3
       gradvec_tc_l(fff,indx)=res_l(fff)
       gradvec_tc_r(fff,indx)=res_r(fff)
      enddo
    end do
  end do
END_PROVIDER 

subroutine gradvec_tc_ia(i,a,res_l, res_r)
 implicit none
 BEGIN_DOC
! doubly occupied --> virtual TC gradient 
!
! Corresponds to res_r = <X0|H E_i^a|Phi_0>, 
!
!                res_l = <X0|E_a^i H|Phi_0>
 END_DOC
 integer, intent(in) :: i,a
 double precision, intent(out) :: res_l(0:3), res_r(0:3)
 res_l = 0.d0
 res_r = 0.d0
 res_l(1) = -2 * mo_bi_ortho_tc_one_e(a,i)
 res_r(1) = -2 * mo_bi_ortho_tc_one_e(i,a)
 integer :: j,t,r,jj,tt,rr
 do jj = 1, n_core_inact_orb
  j = list_core_inact(jj)
  res_r(2) += -2.d0 * ( 2.d0 * mo_bi_ortho_tc_two_e(j,i,j,a) - mo_bi_ortho_tc_two_e(i,j,j,a))
  res_l(2) += -2.d0 * ( 2.d0 * mo_bi_ortho_tc_two_e(j,a,j,i) - mo_bi_ortho_tc_two_e(j,a,i,j))
 enddo
 do tt = 1, n_act_orb
  t = list_act(tt)
  do rr = 1, n_act_orb
   r = list_act(rr)
   res_r(2) += -0.5d0 * (                                                                                    &
   tc_transition_matrix_mo(r,t,1,1) *(2.d0 * mo_bi_ortho_tc_two_e(r,i,t,a) - mo_bi_ortho_tc_two_e(i,r,t,a))  & 
  +tc_transition_matrix_mo(t,r,1,1) *(2.d0 * mo_bi_ortho_tc_two_e(t,i,r,a) - mo_bi_ortho_tc_two_e(i,t,r,a))  &
   )
   res_l(2) += -0.5d0 * (                                                                                    &
   tc_transition_matrix_mo(t,r,1,1) *(2.d0 * mo_bi_ortho_tc_two_e(t,a,r,i) - mo_bi_ortho_tc_two_e(t,a,i,r))  & 
  +tc_transition_matrix_mo(t,r,1,1) *(2.d0 * mo_bi_ortho_tc_two_e(r,a,t,i) - mo_bi_ortho_tc_two_e(r,a,i,t))  &
   )
  enddo
 enddo
end

subroutine gradvec_tc_it(i,t,res_l, res_r)
 implicit none
 BEGIN_DOC
! doubly occupied --> active TC gradient 
!
! Corresponds to res_r = <X0|H E_i^t|Phi_0>
!
!                res_l = <X0|E_t^i H |Phi_0>
 END_DOC
 integer, intent(in) :: i,t
 double precision, intent(out) :: res_l(0:3),res_r(0:3)
 integer :: rr,r,ss,s,m,mm
 double precision :: dm
 res_r = 0.d0
 res_l = 0.d0
 res_r(1)  += -2.d0 * mo_bi_ortho_tc_one_e(i,t) 
 res_l(1)  +=  2.D0 * mo_bi_ortho_tc_one_e(t,i) 

 do rr = 1, n_act_orb
  r = list_act(rr)
  res_r(1) +=  mo_bi_ortho_tc_one_e(i,r) * tc_transition_matrix_mo(t,r,1,1)   
  res_l(1) += -mo_bi_ortho_tc_one_e(r,i) * tc_transition_matrix_mo(r,t,1,1)
 enddo
 
end

subroutine gradvec_tc_ta(t,a,res_l, res_r)
 implicit none
 BEGIN_DOC
! active --> virtual TC gradient 
!
! Corresponds to res_r = <X0|H E_t^a|Phi_0>
!
!                res_l = <X0|E_a^t H |Phi_0>
 END_DOC
 integer, intent(in) :: t,a
 double precision, intent(out) :: res_l(0:3),res_r(0:3)
 integer :: rr,r,m
 double precision :: dm
 res_r = 0.d0
 res_l = 0.d0
 do rr = 1, n_act_orb
  r = list_act(rr) 
  res_l(1) +=  mo_bi_ortho_tc_one_e(a,r) * tc_transition_matrix_mo(t,r,1,1)
  res_r(1) += -mo_bi_ortho_tc_one_e(r,a) * tc_transition_matrix_mo(r,t,1,1)
 enddo
 
end
